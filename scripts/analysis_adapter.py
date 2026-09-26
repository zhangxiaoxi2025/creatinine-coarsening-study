"""Frozen DataFrames to a paired analysis, with no file access or resampling.

The caller must verify current access, source hashes, and the immutable design
manifest BEFORE loading participant DataFrames or calling this module. Building
PairedDesign computes its original and fixed-grid patient labels internally;
therefore this API belongs to the authorized formal-analysis phase.

Put the immutable design-freeze scripts directory on sys.path before importing
this module. Importing it neither reads patient data nor executes an analysis.
"""
from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from typing import Any

import numpy as np
import pandas as pd

from aki_core import decimal, prepare_patient
from paired_resampling_exact import ANALYSIS_ORDER, PairedDesign

UNIT_SCALE = 10**15


COHORT_FLAGS = (
    "documented_renal_exclusion", "primary", "secondary7",
    "procedural_only_primary", "conflict_in_window", "extreme_in_window",
    "duration_over24h", "primary_after5min", "secondary7_after5min",
)
OBSERVATION_FLAGS = ("is_pre7", "is_post", "E48", "E7")
COHORT_COLUMNS = ("case_id", "patient_id", "start", "end", "discharge", *COHORT_FLAGS)
MEASUREMENT_COLUMNS = (
    "case_id", "patient_id", "start", "end", "discharge", "time", "value_str",
    *OBSERVATION_FLAGS,
)


@dataclass(frozen=True)
class AdaptedAnalysis:
    design: PairedDesign
    patient_keys: tuple[str, ...]
    diagnostics: dict[str, Any]


def exact_units(value: str | Decimal) -> int:
    """Parse original decimal text exactly; never round to the fast-path unit."""
    if not isinstance(value, (str, Decimal)):
        raise ValueError("value_str must preserve original decimal text")
    parsed = decimal(value)
    if not Decimal(0) < parsed < Decimal(100):
        raise ValueError("Original SCr must satisfy finite 0 < value < 100")
    # Fraction avoids Decimal-context rounding even for long source literals.
    scaled = Fraction(parsed) * UNIT_SCALE
    if scaled.denominator != 1:
        raise ValueError("Original SCr is not exactly representable in units; STOP, do not round")
    return scaled.numerator


def _require_schema(frame: pd.DataFrame, columns, flags, label: str) -> None:
    absent = set(columns) - set(frame.columns)
    if absent:
        raise ValueError(f"{label} is missing columns: {sorted(absent)}")
    for name in ("case_id", "patient_id"):
        if not frame[name].map(lambda x: isinstance(x, str) and bool(x)).all():
            raise ValueError(f"{label}.{name} must contain original nonempty strings")
    for name in flags:
        if not frame[name].map(lambda x: isinstance(x, (bool, np.bool_))).all():
            raise ValueError(f"{label}.{name} must contain only actual booleans")


def _integer_array(values) -> np.ndarray:
    return np.asarray(values, dtype=np.int64)


def build_analysis(cohort: pd.DataFrame, measurements: pd.DataFrame, *,
                   analysis: str) -> AdaptedAnalysis:
    """Apply one of the eight frozen modes, validating observation contracts.

    Full procedural-cohort DataFrames are accepted. Frozen cohort/timing flags
    first select the mode-specific population; its patients are reconstructed
    with prepare_patient. Discrepancies in selected flags, raw precision, linkage
    or timing abort. The input frames are never mutated. E-outside and early postoperative
    records remain in the reference pool and can contribute to later minima.
    """
    if analysis not in ANALYSIS_ORDER:
        raise ValueError(f"Unknown analysis {analysis!r}; frozen modes are {ANALYSIS_ORDER}")
    _require_schema(cohort, COHORT_COLUMNS, COHORT_FLAGS, "cohort")
    _require_schema(measurements, MEASUREMENT_COLUMNS, OBSERVATION_FLAGS, "measurements")
    if cohort.case_id.duplicated().any() or cohort.patient_id.duplicated().any():
        raise ValueError("Cohort must contain one unique case per unique patient")
    if not set(measurements.case_id) <= set(cohort.case_id):
        raise ValueError("Measurement case is absent from the supplied procedural cohort")
    base_flag = {"parallel7": "secondary7", "S1_no_diagnosis": "procedural_only_primary"}.get(analysis, "primary")
    base = cohort.loc[cohort[base_flag]].set_index("case_id", drop=False)
    base_keys = sorted(base.index)
    grouped = {str(key): frame for key, frame in measurements.loc[measurements.case_id.isin(base_keys)].groupby("case_id", sort=False)}
    references_patient, references_units = [], []
    point_patient, current_units, min7_units, min48_units = [], [], [], []
    keys = []
    count_postop = count_e48 = count_e7 = count_unobservable = 0
    patients_with_unobservable = extra_floor_records = extra_floor_patients = 0
    floor_output_records = target_first5min = first5min_reference = 0

    for key in base_keys:
        row = base.loc[key]
        # Use the hash-verified frozen selection first. Raw precision in a person
        # excluded from this mode does not constrain this mode's numeric engine.
        if ((analysis == "S2_target_after5min" and not row.primary_after5min)
                or (analysis == "S3_exclude_extremes" and row.extreme_in_window)
                or (analysis == "S4_exclude_over24h" and row.duration_over24h)
                or (analysis == "S5_exclude_conflicts" and row.conflict_in_window)):
            continue
        start, end, discharge = (decimal(row[name]) for name in ("start", "end", "discharge"))
        if not start < end <= discharge:
            raise ValueError("Invalid frozen anesthesia/discharge order")
        raw = grouped.get(key)
        if raw is None or raw.empty:
            raise ValueError("Frozen eligible patient has no reference observations")
        if analysis == "S6_all_postop_E48" and (raw.is_post & ~raw.E48).any():
            continue
        if not raw.patient_id.eq(row.patient_id).all():
            raise ValueError("Measurement and cohort patient identity disagree")
        for name, expected in (("start", start), ("end", end), ("discharge", discharge)):
            if any(decimal(x) != expected for x in raw[name]):
                raise ValueError(f"Measurement and cohort {name} disagree")
        lower, upper = Fraction(start) - 604800, min(Fraction(end) + 604800, Fraction(discharge))
        raw_times = [decimal(t) for t in raw.time]
        if any(not lower <= Fraction(t) <= upper for t in raw_times):
            raise ValueError("Frozen reference observation lies outside the global window")
        source_units = [exact_units(v) for v in raw.value_str]
        plan = prepare_patient(zip(raw_times, raw.value_str), aneend=end, discharge=discharge)
        if plan.cleaned.exact_duplicates_removed or plan.cleaned.conflicting_times:
            raise ValueError("Frozen input is not deduplicated and conflict-free")
        clean = plan.cleaned.measurements
        units = [exact_units(m.value) for m in clean]
        by_time = {t: (c, f) for t, c, f in zip(raw_times, source_units, raw[list(OBSERVATION_FLAGS)].itertuples(index=False, name=None))}
        p48 = {p.index: p for p in plan.primary_points}
        p7 = {p.index: p for p in plan.secondary_points}
        postop = set(plan.postoperative_indices)
        pre7 = {i for i, m in enumerate(clean) if lower <= Fraction(m.time) < Fraction(start)}
        for i, m in enumerate(clean):
            expected = (i in pre7, i in postop, i in p48, i in p7)
            if by_time[m.time] != (units[i], expected):
                raise ValueError("Frozen observation flags disagree with exact timestamp reconstruction")
        primary = bool(pre7 and p48 and not row.documented_renal_exclusion)
        secondary = bool(pre7 and p7 and not row.documented_renal_exclusion)
        late48 = [p for p in plan.primary_points if Fraction(clean[p.index].time) - Fraction(end) > 300]
        late7 = [p for p in plan.secondary_points if Fraction(clean[p.index].time) - Fraction(end) > 300]
        extreme = any(c < UNIT_SCALE // 10 or c > 30 * UNIT_SCALE for c in units)
        expected_flags = {
            "primary": primary, "secondary7": secondary,
            "procedural_only_primary": bool(pre7 and p48),
            "primary_after5min": bool(primary and late48),
            "secondary7_after5min": bool(secondary and late7),
            "extreme_in_window": extreme,
            "duration_over24h": Fraction(end) - Fraction(start) > 86400,
        }
        if any(bool(row[name]) != expected for name, expected in expected_flags.items()):
            raise ValueError("Frozen cohort flags disagree with exact observation reconstruction")
        if not pre7:
            raise ValueError("Selected frozen patient lacks preoperative 7-day coverage")
        include = {
            "main": primary,
            "parallel7": secondary,
            "S1_no_diagnosis": bool(pre7 and p48),
            "S2_target_after5min": bool(primary and late48),
            "S3_exclude_extremes": primary and not extreme,
            "S4_exclude_over24h": primary and not row.duration_over24h,
            "S5_exclude_conflicts": primary and not row.conflict_in_window,
            "S6_all_postop_E48": primary and postop == set(p48),
        }[analysis]
        if not include:
            continue
        points = (plan.secondary_points if analysis == "parallel7" else
                  late48 if analysis == "S2_target_after5min" else plan.primary_points)
        if not points:
            raise ValueError("Selected patient has no target points")
        index = len(keys)
        keys.append(key)
        references_patient.extend([index] * len(clean))
        references_units.extend(units)
        for point in points:
            point_patient.append(index)
            current_units.append(units[point.index])
            min7_units.append(min(units[i] for i in point.prior7))
            if analysis != "parallel7":
                min48_units.append(min(units[i] for i in point.prior48))
        count_postop += len(postop)
        count_e48 += len(p48)
        count_e7 += len(p7)
        unobservable = len(postop - set(p48))
        count_unobservable += unobservable
        patients_with_unobservable += bool(unobservable)
        extra_floor_records += sum(c < UNIT_SCALE // 10 for c in units)
        extra_floor_patients += any(c < UNIT_SCALE // 10 for c in units)
        floor_output_records += sum(c < 3 * UNIT_SCALE // 10 for c in units)
        first5min_reference += sum(0 < Fraction(m.time) - Fraction(end) <= 300 for m in clean)
        target_first5min += sum(0 < Fraction(clean[p.index].time) - Fraction(end) <= 300 for p in points)

    if not keys:
        raise ValueError("No eligible patients remain; STOP")
    design = PairedDesign(
        n_patients=len(keys), reference_patient=_integer_array(references_patient),
        reference_units=_integer_array(references_units), point_patient=_integer_array(point_patient),
        current_units=_integer_array(current_units), min7_units=_integer_array(min7_units),
        min48_units=None if analysis == "parallel7" else _integer_array(min48_units),
        mode="ratio" if analysis == "parallel7" else "composite", unit_scale=UNIT_SCALE,
    )
    diagnostics = {
        "analysis": analysis, "base_flag": base_flag, "base_patients": len(base_keys),
        "n_patients": len(keys), "mode_excluded_patients": len(base_keys) - len(keys),
        "reference_measurements": len(references_units), "target_points": len(point_patient),
        "all_observed_postop_points": count_postop, "all_E48_points": count_e48,
        "all_E7_points": count_e7, "postop_points_outside_E48": count_unobservable,
        "patients_with_postop_points_outside_E48": patients_with_unobservable,
        "first5min_reference_measurements_retained": first5min_reference,
        "first5min_target_points": target_first5min,
        "grid_extra_positive_floor_reference_measurements": extra_floor_records,
        "grid_extra_positive_floor_patients": extra_floor_patients,
        "grid_output_equals_0p2_reference_measurements": floor_output_records,
        "patient_order": "original case_id strings, lexicographic",
        "numeric_unit": "exact 10^-15 mg/dL; no source rounding", "unit_scale": UNIT_SCALE,
        "reference_pool": "all legal window observations in mode-specific patient cohort",
    }
    return AdaptedAnalysis(design, tuple(keys), diagnostics)
