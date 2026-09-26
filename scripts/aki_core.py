"""Pure numerical core for the predeclared serum-creatinine label experiment.

This is a study-specific operational algorithm, not a complete clinical AKI
diagnosis or a claim to reproduce INSPIRE's unpublished implementation.
There are no file readers, patient loaders, bootstrap or cohort selection here.

INTERFACE
---------
``prepare_patient([(seconds, mg_dl), ...], aneend=..., discharge=...)``
    Accepts finite decimal strings, Decimal, integers or finite floats. Float
    inputs mean exactly their ``str`` representation: binary arithmetic already
    performed upstream cannot be undone. Prefer decimal strings from the source.
    Seconds must have one common origin; neither endpoint need be zero.
    Removes exact (time,value) duplicates, then excludes ENTIRE times with >1
    distinct value. Values must already satisfy the caller's quality/cohort
    rules; nonpositive/nonfinite values raise instead of becoming missing.
    Prior windows are [t-48h,t) and [t-7d,t). Postoperative times are
    (aneend,min(aneend+7d,discharge)]. The returned immutable plan holds fixed
    measurement and prior-window indices, including both E48 and E7. E48 is the
    primary observable set; E7 is the secondary ratio-only observable set.
    The caller implements preoperative coverage, eligibility and exclusions.

``fit_transforms(reference_values)``
    Takes the predeclared pool of eligible, deduplicated observations, preserving
    repeated numeric values from different observations. Returns four callables:
    identity, tail_only, tail_plus_bins, and lower_bounded_round_0p2.
    All reference and transformed values must be positive. Type-1 quantiles are
    empirical inverse-CDF quantiles. Midranks use the ORIGINAL pool. Strict tail
    inequalities clip; either exact tail endpoint is kept at that endpoint.
    Interior values choose nearest .05,.10,...,.95 midrank node, lower on ties.
    The .2 control is explicitly lower-bounded half-up rounding:
    .2 * max(1, floor(x/.2 + .5)); it is NOT ordinary unbounded rounding.

``evaluate_patient(plan, transform=None)``
    Reuses the plan's exact measurement/prior indices for every numeric version.
    At E48 times, the composite is current-min48 >= .3 OR current >= 1.5*min7.
    At E7 times, the secondary ratio-only rule is current >= 1.5*min7.
    Comparisons use exact rational arithmetic, not floating tolerance. Patient
    labels are any eligible positive, False only if an eligible set is nonempty,
    and None if that set is empty. Point outputs preserve both criteria.
"""

from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from typing import Callable, Iterable

Number = str | int | float | Decimal
H48 = Decimal(172800)
D7 = Decimal(604800)
Transform = Callable[[Decimal], Decimal]


def decimal(value: Number) -> Decimal:
    """Parse finite numbers without accepting bools or silently rounding."""
    if isinstance(value, bool):
        raise ValueError("Boolean is not a numeric measurement")
    try:
        result = value if isinstance(value, Decimal) else Decimal(str(value))
    except (InvalidOperation, ValueError, TypeError) as exc:
        raise ValueError(f"Invalid decimal: {value!r}") from exc
    if not result.is_finite():
        raise ValueError("Nonfinite decimal")
    return result


def positive_decimal(value: Number) -> Decimal:
    result = decimal(value)
    if result <= 0:
        raise ValueError("Creatinine values must be positive before/after transform")
    return result


@dataclass(frozen=True)
class Measurement:
    time: Decimal
    value: Decimal


@dataclass(frozen=True)
class CleanSeries:
    measurements: tuple[Measurement, ...]
    input_count: int
    exact_duplicates_removed: int
    conflicting_times: tuple[Decimal, ...]
    conflicting_distinct_records_removed: int


@dataclass(frozen=True)
class ObservablePoint:
    index: int
    prior48: tuple[int, ...]
    prior7: tuple[int, ...]


@dataclass(frozen=True)
class PatientPlan:
    cleaned: CleanSeries
    aneend: Decimal
    observation_end: Fraction
    postoperative_indices: tuple[int, ...]
    primary_points: tuple[ObservablePoint, ...]
    secondary_points: tuple[ObservablePoint, ...]


@dataclass(frozen=True)
class PointResult:
    time: Decimal
    value: Decimal
    min48: Decimal | None
    min7: Decimal
    absolute_positive: bool | None
    ratio_positive: bool
    composite_positive: bool | None


@dataclass(frozen=True)
class PatientResult:
    primary_label: bool | None
    secondary_ratio_label: bool | None
    primary_points: tuple[PointResult, ...]
    secondary_points: tuple[PointResult, ...]


def clean_series(records: Iterable[tuple[Number, Number]]) -> CleanSeries:
    by_time: dict[Decimal, set[Decimal]] = {}
    input_count = 0
    for time, value in records:
        input_count += 1
        by_time.setdefault(decimal(time), set()).add(positive_decimal(value))
    unique_count = sum(len(values) for values in by_time.values())
    conflicting = tuple(sorted(t for t, values in by_time.items() if len(values) > 1))
    return CleanSeries(
        measurements=tuple(
            Measurement(t, next(iter(by_time[t])))
            for t in sorted(by_time)
            if len(by_time[t]) == 1
        ),
        input_count=input_count,
        exact_duplicates_removed=input_count - unique_count,
        conflicting_times=conflicting,
        conflicting_distinct_records_removed=sum(len(by_time[t]) for t in conflicting),
    )


def prepare_patient(
    records: Iterable[tuple[Number, Number]], *, aneend: Number, discharge: Number
) -> PatientPlan:
    cleaned = clean_series(records)
    end = decimal(aneend)
    dis = decimal(discharge)
    if dis < end:
        raise ValueError("Discharge precedes anesthesia end")
    # Fractions retain subsecond precision even for huge epoch-like seconds.
    upper = min(Fraction(end) + Fraction(D7), Fraction(dis))
    times = [Fraction(m.time) for m in cleaned.measurements]
    postop = tuple(i for i, t in enumerate(times) if Fraction(end) < t <= upper)
    primary: list[ObservablePoint] = []
    secondary: list[ObservablePoint] = []
    for i in postop:
        t = times[i]
        prior48 = tuple(range(bisect_left(times, t - Fraction(H48)), i))
        prior7 = tuple(range(bisect_left(times, t - Fraction(D7)), i))
        point = ObservablePoint(i, prior48, prior7)
        if prior48:
            primary.append(point)
        if prior7:
            secondary.append(point)
    return PatientPlan(cleaned, end, upper, postop, tuple(primary), tuple(secondary))


def evaluate_patient(plan: PatientPlan, transform: Transform | None = None) -> PatientResult:
    transform = transform or (lambda x: x)
    values = tuple(positive_decimal(transform(m.value)) for m in plan.cleaned.measurements)

    def evaluate_point(point: ObservablePoint) -> PointResult:
        current = values[point.index]
        min7 = min(values[i] for i in point.prior7)
        min48 = min((values[i] for i in point.prior48), default=None)
        ratio = 2 * Fraction(current) >= 3 * Fraction(min7)
        absolute = (
            Fraction(current) - Fraction(min48) >= Fraction(3, 10)
            if min48 is not None else None
        )
        composite = (absolute or ratio) if absolute is not None else None
        return PointResult(
            plan.cleaned.measurements[point.index].time,
            current, min48, min7, absolute, ratio, composite,
        )

    primary = tuple(evaluate_point(p) for p in plan.primary_points)
    secondary = tuple(evaluate_point(p) for p in plan.secondary_points)
    return PatientResult(
        any(p.composite_positive for p in primary) if primary else None,
        any(p.ratio_positive for p in secondary) if secondary else None,
        primary, secondary,
    )


def type1_quantile(sorted_values: tuple[Decimal, ...], probability: Fraction) -> Decimal:
    """Exact empirical inverse CDF; input must already be sorted and nonempty."""
    if not sorted_values:
        raise ValueError("Empty quantile reference")
    if not Fraction(0) <= probability <= Fraction(1):
        raise ValueError("Probability outside [0,1]")
    rank = probability * len(sorted_values)
    # ceil numerator/denominator, bounded below so p=0 returns the minimum.
    k = max(1, (rank.numerator + rank.denominator - 1) // rank.denominator)
    return sorted_values[k - 1]


def _fifths_decimal(steps: int) -> Decimal:
    """Construct exact .2*steps even under a short Decimal precision context."""
    coefficient = str(2 * steps)
    return Decimal((0, tuple(int(d) for d in coefficient), -1))


def lower_bounded_round_0p2(value: Number) -> Decimal:
    value = positive_decimal(value)
    shifted = Fraction(value) * 5 + Fraction(1, 2)
    steps = max(1, shifted.numerator // shifted.denominator)
    return _fifths_decimal(steps)


@dataclass(frozen=True)
class NumericTransforms:
    reference: tuple[Decimal, ...]
    q025: Decimal
    q975: Decimal
    nodes: tuple[tuple[Fraction, Decimal], ...]

    def identity(self, value: Number) -> Decimal:
        return positive_decimal(value)

    def tail_only(self, value: Number) -> Decimal:
        value = positive_decimal(value)
        return min(max(value, self.q025), self.q975)

    def tail_plus_bins(self, value: Number) -> Decimal:
        value = positive_decimal(value)
        if value < self.q025:
            return self.q025
        if value > self.q975:
            return self.q975
        if value == self.q025 or value == self.q975:
            return value
        less = bisect_left(self.reference, value)
        less_equal = bisect_right(self.reference, value)
        midrank = Fraction(less + less_equal, 2 * len(self.reference))
        _, representative = min(self.nodes, key=lambda node: (abs(midrank - node[0]), node[0]))
        return representative

    def versions(self) -> dict[str, Transform]:
        return {
            "identity": self.identity,
            "tail_only": self.tail_only,
            "tail_plus_bins": self.tail_plus_bins,
            "lower_bounded_round_0p2": lower_bounded_round_0p2,
        }


def fit_transforms(reference_values: Iterable[Number]) -> NumericTransforms:
    reference = tuple(sorted(positive_decimal(x) for x in reference_values))
    if not reference:
        raise ValueError("Empty transform reference")
    return NumericTransforms(
        reference,
        type1_quantile(reference, Fraction(1, 40)),
        type1_quantile(reference, Fraction(39, 40)),
        tuple((Fraction(k, 20), type1_quantile(reference, Fraction(k, 20))) for k in range(1, 20)),
    )
