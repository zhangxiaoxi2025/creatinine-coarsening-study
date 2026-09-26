"""Exact-cent paired summaries and patient-cluster bootstrap; array API only.

No participant file readers, real-data CLI, cohort construction or timestamp
logic are included. The caller supplies an already frozen cohort and E48/E7
points, with reference labs from the same cohort's full eligible time window.

Raw values must be EXACT positive integer cents (0.01 mg/dL). Do not silently
round higher-precision source values to use this fast path: use the general
aki_core implementation or revise the common unit before analysis. Cent values
are never filtered after transformation; e.g. 9999 -> 10000 is retained.

Each reference observation is repeated according to its patient's bootstrap
multiplicity. Midrank and type-1 quantile decisions use integer arithmetic.
Since all four transforms are nondecreasing, transforming a previously
computed original minimum equals the minimum of the transformed window.

The bootstrap intervals describe empirical patient-resampling stability.
Ordinary percentile coverage for this discrete, fitted-threshold statistic is
not established; a degenerate zero interval is not evidence of zero population
risk. No simple binomial interval is silently substituted.
"""
from __future__ import annotations

from dataclasses import dataclass
from copy import deepcopy
import hashlib
from typing import Mapping

import numpy as np

VERSIONS = ("identity", "tail_only", "tail_plus_bins", "lower_bounded_round_0p2")
CENTER_ORDER = ("MOVER", "VitalDB")
STUDY_CENTER_ORDER = ("mover", "vitaldb")
ANALYSIS_ORDER = ("main", "parallel7", "S1_no_diagnosis", "S2_target_after5min",
                  "S3_exclude_extremes", "S4_exclude_over24h",
                  "S5_exclude_conflicts", "S6_all_postop_E48")
RATE_NAMES = ("discordance", "upward", "downward", "net_change",
              "original_positive_rate", "transformed_positive_rate")
# np.bincount's weighted implementation uses binary64. Integer additions are
# exact up to this bound; checking it makes conversion back to int64 exact.
MAX_EXACT_COUNT = 2**53 - 1
MAX_SAFE_CENTS = np.iinfo(np.int64).max // 4


def _integer_vector(value, name: str, *, positive: bool = False) -> np.ndarray:
    raw = np.asarray(value)
    if raw.ndim != 1 or raw.dtype.kind not in "iu":
        raise ValueError(f"{name} must be a one-dimensional integer array; no rounding")
    if raw.size and (np.any(raw < (1 if positive else 0)) or np.any(raw > np.iinfo(np.int64).max)):
        raise ValueError(f"{name} has out-of-range values")
    result = raw.astype(np.int64, copy=True)
    result.setflags(write=False)
    return result


def _readonly(value: np.ndarray) -> np.ndarray:
    value.setflags(write=False)
    return value


def _metric_rates(counts: np.ndarray) -> dict[str, np.ndarray | float]:
    # Last axis is always n00,n01,n10,n11. Inputs are exact integer counts.
    n = counts.sum(axis=-1)
    return {
        "discordance": (counts[..., 1] + counts[..., 2]) / n,
        "upward": counts[..., 1] / n,
        "downward": counts[..., 2] / n,
        "net_change": (counts[..., 1] - counts[..., 2]) / n,
        "original_positive_rate": (counts[..., 2] + counts[..., 3]) / n,
        "transformed_positive_rate": (counts[..., 1] + counts[..., 3]) / n,
    }


def center_rngs(seed: int = 20260921) -> dict[str, np.random.Generator]:
    """Generic two-center streams for standalone synthetic checks only.

    This does NOT implement the frozen study hierarchy. Study analyses must
    use study_rngs()[center][analysis] instead.
    """
    children = np.random.SeedSequence(seed).spawn(len(CENTER_ORDER))
    return {center: np.random.Generator(np.random.PCG64(child))
            for center, child in zip(CENTER_ORDER, children)}


def study_rngs(seed: int = 20260921) -> dict[str, dict[str, np.random.Generator]]:
    """Frozen two-level PCG64 streams, independent of analysis execution order.

    Centers and analyses exactly match analysis-spec-v1.0.json. Create this
    collection once per run, then select its named stream for every analysis.
    Returning all children together avoids call-order-dependent spawn behavior.
    """
    centers = np.random.SeedSequence(seed).spawn(len(STUDY_CENTER_ORDER))
    return {
        center: {
            analysis: np.random.Generator(np.random.PCG64(child))
            for analysis, child in zip(ANALYSIS_ORDER, center_seed.spawn(len(ANALYSIS_ORDER)))
        }
        for center, center_seed in zip(STUDY_CENTER_ORDER, centers)
    }


def draw_patient_multiplicities(n: int, rng: np.random.Generator) -> np.ndarray:
    if not isinstance(n, (int, np.integer)) or isinstance(n, (bool, np.bool_)) or n < 1:
        raise ValueError("n must be a positive integer")
    if not isinstance(rng, np.random.Generator):
        raise ValueError("an explicit NumPy Generator is required")
    return rng.multinomial(int(n), np.full(int(n), 1.0 / n))


@dataclass(frozen=True)
class MappingFit:
    reference_counts: np.ndarray
    q025: int
    q975: int
    node_values: np.ndarray
    mappings: Mapping[str, np.ndarray]


@dataclass(frozen=True)
class Evaluation:
    counts: np.ndarray  # versions x (n00,n01,n10,n11)
    labels: Mapping[str, np.ndarray]
    mapping_fit: MappingFit

    def tables(self) -> dict[str, dict[str, int | float]]:
        result = {}
        for k, name in enumerate(VERSIONS):
            c = self.counts[k]
            row = {key: int(v) for key, v in zip(("n00", "n01", "n10", "n11"), c)}
            row.update(n=int(c.sum()), discordance_count=int(c[1] + c[2]))
            row.update({key: float(v) for key, v in _metric_rates(c).items()})
            result[name] = row
        return result


@dataclass(frozen=True)
class BootstrapResult:
    point_estimate: Evaluation
    replicate_counts: np.ndarray  # resamples x versions x four exact cells
    intervals: Mapping[str, Mapping[str, tuple[float, float]]]
    diagnostics: Mapping[str, object]
    initial_rng_state: Mapping[str, object]
    final_rng_state: Mapping[str, object]


class PairedDesign:
    """Frozen reference labs and already observable points for one center/mode.

    Patient IDs are dense indices 0..n_patients-1. Every patient must have at
    least one reference observation and one fixed observable point. All point
    values/minima must occur in the reference pool. Input arrays are copied and
    marked read-only; reference multiplicity, not patient outcome, fits maps.

    mode='composite' requires min48_cents and min7_cents. mode='ratio' uses only
    min7_cents and expects min48_cents=None. Under the frozen labels, original
    min7 <= min48 must hold for every composite point.
    """

    def __init__(self, *, n_patients: int, reference_patient, reference_cents,
                 point_patient, current_cents, min7_cents, min48_cents=None,
                 mode: str = "composite"):
        if not isinstance(n_patients, (int, np.integer)) or isinstance(n_patients, (bool, np.bool_)) or n_patients < 1:
            raise ValueError("n_patients must be a positive integer")
        if mode not in ("composite", "ratio"):
            raise ValueError("mode must be composite or ratio")
        self.n_patients = int(n_patients)
        self.mode = mode
        rp = _integer_vector(reference_patient, "reference_patient")
        rv = _integer_vector(reference_cents, "reference_cents", positive=True)
        pp = _integer_vector(point_patient, "point_patient")
        cv = _integer_vector(current_cents, "current_cents", positive=True)
        m7 = _integer_vector(min7_cents, "min7_cents", positive=True)
        if len(rp) != len(rv) or not len(rp):
            raise ValueError("nonempty reference arrays must have matching lengths")
        if len(pp) != len(cv) or len(pp) != len(m7) or not len(pp):
            raise ValueError("nonempty point arrays must have matching lengths")
        if np.any(rp >= self.n_patients) or np.any(pp >= self.n_patients):
            raise ValueError("patient index outside frozen cohort")
        if np.any(rv > MAX_SAFE_CENTS) or np.any(cv > MAX_SAFE_CENTS) or np.any(m7 > MAX_SAFE_CENTS):
            raise ValueError("cent values exceed safe exact int64 arithmetic range")
        self._reference_per_patient = np.bincount(rp, minlength=self.n_patients).astype(np.int64)
        points_per_patient = np.bincount(pp, minlength=self.n_patients)
        if np.any(self._reference_per_patient == 0) or np.any(points_per_patient == 0):
            raise ValueError("every frozen patient must have reference labs and an observable point")
        if mode == "composite":
            if min48_cents is None:
                raise ValueError("composite mode requires min48_cents")
            m48 = _integer_vector(min48_cents, "min48_cents", positive=True)
            if len(m48) != len(pp) or np.any(m48 > MAX_SAFE_CENTS):
                raise ValueError("invalid min48_cents length/range")
            if np.any(m7 > m48):
                raise ValueError("7-day minimum cannot exceed nested 48-hour minimum")
        else:
            if min48_cents is not None:
                raise ValueError("ratio mode must not carry a composite 48-hour minimum")
            m48 = None
        # Stable sort permits fast complete patient-level OR without resampling
        # individual points or changing their within-patient membership.
        order = np.argsort(pp, kind="stable")
        self.point_patient = _readonly(pp[order])
        self.reference_patient = rp
        self.reference_cents = rv
        self.support, reference_code = np.unique(rv, return_inverse=True)
        self.support = _readonly(self.support)
        self._reference_code = _readonly(reference_code)
        self._point_starts = np.r_[0, np.flatnonzero(np.diff(self.point_patient)) + 1]

        def code(values, name):
            sorted_values = values[order]
            idx = np.searchsorted(self.support, sorted_values)
            in_range = idx < len(self.support)
            if not np.all(in_range) or not np.array_equal(self.support[idx], sorted_values):
                raise ValueError(f"{name} includes a value outside the frozen reference pool")
            return _readonly(idx)

        self._current_code = code(cv, "current_cents")
        self._min7_code = code(m7, "min7_cents")
        self._min48_code = None if m48 is None else code(m48, "min48_cents")
        self._original_labels = self._labels(self.support)
        self._fixed_round_map = _readonly(np.maximum(20, ((self.support + 10) // 20) * 20))
        self._fixed_round_labels = self._labels(self._fixed_round_map)

    def _weights(self, weights=None) -> np.ndarray:
        if weights is None:
            return np.ones(self.n_patients, dtype=np.int64)
        w = _integer_vector(weights, "patient multiplicities")
        if len(w) != self.n_patients:
            raise ValueError("patient multiplicity length differs from frozen cohort")
        total = sum(map(int, w))
        if total < 1 or total > MAX_EXACT_COUNT:
            raise ValueError("patient multiplicities must have a positive bounded total")
        return w

    def reference_frequency(self, weights=None) -> np.ndarray:
        w = self._weights(weights)
        # This bound also makes np.int64 products/sums below exact and safe.
        max_observations = int(w.max()) * len(self.reference_patient)
        if max_observations > MAX_EXACT_COUNT:
            raise ValueError("weighted reference count exceeds exact integer bincount range")
        counts = np.bincount(self._reference_code,
                             weights=w[self.reference_patient],
                             minlength=len(self.support)).astype(np.int64)
        if counts.sum() <= 0:
            raise ValueError("empty weighted reference")
        return counts

    def fit_mappings(self, weights=None) -> MappingFit:
        counts = self.reference_frequency(weights)
        cumulative = np.cumsum(counts, dtype=np.int64)
        total = int(cumulative[-1])
        # ceil(p*N), for p=1/40,39/40 and k/20. No floating probabilities.
        lo_rank, hi_rank = (total + 39) // 40, (39 * total + 39) // 40
        q025 = int(self.support[np.searchsorted(cumulative, lo_rank, side="left")])
        q975 = int(self.support[np.searchsorted(cumulative, hi_rank, side="left")])
        ranks = np.asarray([(k * total + 19) // 20 for k in range(1, 20)], dtype=np.int64)
        nodes = self.support[np.searchsorted(cumulative, ranks, side="left")]
        less = cumulative - counts
        twice_midrank_numerator = 2 * less + counts
        # Nearest k to 20*midrank = 10*(2*less+equal)/N; exact halfway
        # chooses the lower k. Subtracting 1 handles integer exact ties.
        nearest_k = np.clip((20 * twice_midrank_numerator + total - 1) // (2 * total), 1, 19)
        bins = nodes[nearest_k - 1].copy()
        bins[self.support <= q025] = q025
        bins[self.support >= q975] = q975
        tail = np.clip(self.support, q025, q975)
        mappings = {"identity": self.support,
                    "tail_only": _readonly(tail),
                    "tail_plus_bins": _readonly(bins),
                    "lower_bounded_round_0p2": self._fixed_round_map}
        if any(np.any(m <= 0) or np.any(np.diff(m) < 0) for m in mappings.values()):
            raise AssertionError("numeric transform violated positivity/monotonicity")
        return MappingFit(_readonly(counts), q025, q975, _readonly(nodes), mappings)

    def _labels(self, mapping: np.ndarray) -> np.ndarray:
        current = mapping[self._current_code]
        minimum7 = mapping[self._min7_code]
        positive = 2 * current >= 3 * minimum7
        if self.mode == "composite":
            positive |= current - mapping[self._min48_code] >= 30
        return _readonly(np.logical_or.reduceat(positive, self._point_starts))

    def evaluate(self, weights=None) -> Evaluation:
        w = self._weights(weights)
        fit = self.fit_mappings(w)
        labels = {"identity": self._original_labels,
                  "tail_only": self._labels(fit.mappings["tail_only"]),
                  "tail_plus_bins": self._labels(fit.mappings["tail_plus_bins"]),
                  "lower_bounded_round_0p2": self._fixed_round_labels}
        counts = np.empty((len(VERSIONS), 4), dtype=np.int64)
        for k, name in enumerate(VERSIONS):
            pair_code = self._original_labels.astype(np.int8) * 2 + labels[name]
            # Integer np.sum avoids introducing floating weighted cell counts.
            counts[k] = [w[pair_code == cell].sum(dtype=np.int64) for cell in range(4)]
        return Evaluation(_readonly(counts), labels, fit)


def _type1_central_95(values: np.ndarray) -> tuple[float, float]:
    ordered = np.sort(values)
    b = len(ordered)
    return float(ordered[(b + 39) // 40 - 1]), float(ordered[(39*b + 39) // 40 - 1])


def run_bootstrap(design: PairedDesign, *, rng: np.random.Generator,
                  n_resamples: int = 5000) -> BootstrapResult:
    """Resample complete patients, refit maps, and retain all paired count cells.

    Exceptions abort the run; no failed replicate is discarded or replaced.
    Caller selects study_rngs()[center][analysis]; actual generator
    states are returned for auditable reproducibility without claiming a seed
    that may not have been used. Repeated calls advance the supplied generator.
    """
    if not isinstance(n_resamples, (int, np.integer)) or isinstance(n_resamples, (bool, np.bool_)) or n_resamples < 1:
        raise ValueError("n_resamples must be a positive integer")
    if not isinstance(rng, np.random.Generator):
        raise ValueError("an explicit NumPy Generator is required")
    if type(rng.bit_generator) is not np.random.PCG64:
        raise ValueError("the frozen study requires the PCG64 bit generator")
    initial_state = deepcopy(rng.bit_generator.state)
    point = design.evaluate()
    result = np.empty((int(n_resamples), len(VERSIONS), 4), dtype=np.int64)
    distinct_mapping_hashes = set()
    distinct_quantile_nodes = set()
    probabilities = np.full(design.n_patients, 1.0 / design.n_patients)
    for b in range(int(n_resamples)):
        w = rng.multinomial(design.n_patients, probabilities)
        evaluated = design.evaluate(w)
        result[b] = evaluated.counts
        fit = evaluated.mapping_fit
        packed = np.concatenate((fit.mappings["tail_only"], fit.mappings["tail_plus_bins"]))
        distinct_mapping_hashes.add(hashlib.sha256(packed.tobytes()).digest())
        distinct_quantile_nodes.add((fit.q025, *map(int, fit.node_values), fit.q975))
    if not np.all(result.sum(axis=-1) == design.n_patients):
        raise AssertionError("paired denominators changed in bootstrap")
    intervals, per_version = {}, {}
    for k, name in enumerate(VERSIONS):
        rates = _metric_rates(result[:, k, :])
        intervals[name] = {metric: _type1_central_95(np.asarray(rates[metric])) for metric in RATE_NAMES}
        discordance = np.asarray(rates["discordance"])
        all_zero = bool(np.all(discordance == 0))
        point_zero = bool(point.counts[k, 1] + point.counts[k, 2] == 0)
        per_version[name] = {
            "point_discordance_zero": point_zero,
            "all_resampled_discordance_zero": all_zero,
            "zero_degenerate_warning": point_zero and all_zero,
            "discordance_interval_degenerate": intervals[name]["discordance"][0] == intervals[name]["discordance"][1],
            "distinct_discordance_rates": int(len(np.unique(discordance))),
        }
    diagnostics = {
        "resampling_unit": "patient; complete reference and point cluster multiplicities",
        "n_patients": design.n_patients, "n_resamples": int(n_resamples), "mode": design.mode,
        "bit_generator": type(rng.bit_generator).__name__, "failed_replicates": 0,
        "interval_method": "type-1 empirical 2.5/97.5 percentiles; stability interval, coverage not established",
        "distinct_tail_and_bin_mapping_combinations": len(distinct_mapping_hashes),
        "distinct_quantile_node_combinations": len(distinct_quantile_nodes),
        "versions": per_version,
    }
    return BootstrapResult(point, _readonly(result), intervals, diagnostics,
                           initial_state, deepcopy(rng.bit_generator.state))
