"""Exact, array-only kernels for the locked post-result exploratory supplement.

This module does not read participant files, choose cohorts, reconstruct times,
fit models, or resample. All input arrays must come from the frozen design.
Fraction arithmetic determines every weighted quantile and midrank decision;
floating-point numbers appear only in descriptive output proportions.
"""
from __future__ import annotations

from bisect import bisect_left
from fractions import Fraction

import numpy as np


MAX_SAFE_UNITS = np.iinfo(np.int64).max // 4


def _positive_integer(value, name):
    if (not isinstance(value, (int, np.integer))
            or isinstance(value, (bool, np.bool_)) or value < 1):
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _integer_vector(value, name, *, minimum=0, maximum=None):
    raw = np.asarray(value)
    if raw.ndim != 1 or raw.dtype.kind not in "iu":
        raise ValueError(f"{name} must be a one-dimensional integer array")
    upper = np.iinfo(np.int64).max if maximum is None else maximum
    if raw.size and (np.any(raw < minimum) or np.any(raw > upper)):
        raise ValueError(f"{name} contains an out-of-range integer")
    result = raw.astype(np.int64, copy=True)
    result.setflags(write=False)
    return result


def _units(value, name):
    return _integer_vector(value, name, minimum=1, maximum=MAX_SAFE_UNITS)


def _readonly(value):
    value.setflags(write=False)
    return value


def patient_equal_mapping(reference_patient, reference_units, n_patients):
    """Fit B once with exact weight 1/m_i for each retained reference record.

    ``support_weights`` are unnormalized masses (total N);
    ``normalized_weights`` are probability masses (total 1). Both are tuples
    of Fraction in sorted support order. ``record_weights_by_patient`` gives
    1/m_i in dense patient order. The fit is made directly on original units.
    """
    n = _positive_integer(n_patients, "n_patients")
    patient = _integer_vector(reference_patient, "reference_patient")
    values = _units(reference_units, "reference_units")
    if not len(patient) or len(patient) != len(values):
        raise ValueError("nonempty reference arrays must have matching lengths")
    if np.any(patient >= n):
        raise ValueError("reference patient index is outside the frozen cohort")
    counts = np.bincount(patient, minlength=n).astype(np.int64)
    if np.any(counts == 0):
        raise ValueError("every frozen patient must have reference records")
    record_weights = tuple(Fraction(1, int(m)) for m in counts)
    support, codes = np.unique(values, return_inverse=True)
    masses = [Fraction(0) for _ in support]
    for patient_id, code in zip(patient, codes):
        masses[int(code)] += record_weights[int(patient_id)]
    cumulative = []
    total = Fraction(0)
    for weight in masses:
        total += weight
        cumulative.append(total)
    if total != n or any(weight <= 0 for weight in masses):
        raise AssertionError("patient-equal reference mass must equal N")

    def quantile(probability):
        return int(support[bisect_left(cumulative, probability * n)])

    q025, q975 = quantile(Fraction(1, 40)), quantile(Fraction(39, 40))
    nodes = [quantile(Fraction(k, 20)) for k in range(1, 20)]
    mapped = np.empty(len(support), dtype=np.int64)
    midranks = []
    node_indices = []
    for j, (value, weight, upper) in enumerate(zip(support, masses, cumulative)):
        midrank = (upper - weight / 2) / n
        midranks.append(midrank)
        # The nearest integer to 20*M is ceil(20*M - 1/2), with an
        # exact halfway point assigned to the lower integer, then clipped.
        shifted = 20 * midrank - Fraction(1, 2)
        k = max(1, min(19, -(-shifted.numerator // shifted.denominator)))
        node_indices.append(k)
        result = nodes[k - 1]
        if value <= q025:
            result = q025
        if value >= q975:
            result = q975
        mapped[j] = result
    if np.any(mapped <= 0) or np.any(mapped[1:] < mapped[:-1]):
        raise AssertionError("patient-equal transform is not positive/monotone")
    return {
        "support": _readonly(support),
        "mapped_units": _readonly(mapped),
        "support_weights": tuple(masses),
        "normalized_weights": tuple(weight / n for weight in masses),
        "record_weights_by_patient": record_weights,
        "reference_counts_per_patient": _readonly(counts),
        "weighted_midranks": tuple(midranks),
        "nearest_node_indices": node_indices,
        "total_weight": total,
        "q025": q025,
        "q975": q975,
        "node_values": nodes,
    }


def _point_units(design, public_name, private_code_name):
    """Accept explicit fixed vectors or the validated original coded design."""
    if hasattr(design, public_name):
        values = getattr(design, public_name)
        if values is None:
            raise ValueError(f"fixed composite design requires {public_name}")
        return _units(values, public_name)
    if not hasattr(design, "support") or not hasattr(design, private_code_name):
        raise ValueError(f"fixed composite design lacks {public_name}")
    codes = getattr(design, private_code_name)
    if codes is None:
        raise ValueError(f"fixed composite design requires {public_name}")
    codes = _integer_vector(codes, private_code_name)
    source_support = _units(design.support, "design.support")
    if not len(source_support) or np.any(codes >= len(source_support)):
        raise ValueError(f"{public_name} code outside the original support")
    return source_support[codes]


def criterion_states(design, mapping_support, mapped_units):
    """Return int8 patient state 2*A+R; each criterion is an E48-time OR.

    A and R can occur at different eligible times for the same patient.
    Nondecreasing mappings permit the original frozen prior minima to be
    transformed directly. The function never recomputes eligibility or minima.
    """
    n = _positive_integer(design.n_patients, "n_patients")
    scale = _positive_integer(design.unit_scale, "unit_scale")
    if scale % 10 or scale > MAX_SAFE_UNITS:
        raise ValueError("unit_scale must be divisible by 10 and int64-safe")
    if getattr(design, "mode", "composite") != "composite":
        raise ValueError("criterion states require the fixed E48 composite design")
    support = _units(mapping_support, "mapping_support")
    mapped = _units(mapped_units, "mapped_units")
    if not len(support) or len(support) != len(mapped):
        raise ValueError("nonempty mapping support and output must have equal lengths")
    if np.any(support[1:] <= support[:-1]):
        raise ValueError("mapping support must be sorted and unique")
    if np.any(mapped[1:] < mapped[:-1]):
        raise ValueError("mapping must be nondecreasing")
    if hasattr(design, "support"):
        if not np.array_equal(support, _units(design.support, "design.support")):
            raise ValueError("mapping support differs from the fixed design support")
    elif hasattr(design, "reference_units"):
        if not np.array_equal(support, np.unique(_units(design.reference_units, "reference_units"))):
            raise ValueError("mapping support differs from the fixed reference support")

    patient = _integer_vector(design.point_patient, "point_patient")
    current = _point_units(design, "current_units", "_current_code")
    minimum48 = _point_units(design, "min48_units", "_min48_code")
    minimum7 = _point_units(design, "min7_units", "_min7_code")
    if not len(patient) or any(len(x) != len(patient) for x in (current, minimum48, minimum7)):
        raise ValueError("nonempty fixed point vectors must have matching lengths")
    if np.any(patient >= n) or np.any(np.bincount(patient, minlength=n) == 0):
        raise ValueError("every frozen patient must have fixed observable points")
    if np.any(minimum7 > minimum48):
        raise ValueError("nested 7-day minimum exceeds the 48-hour minimum")

    def transform(values):
        code = np.searchsorted(support, values)
        if np.any(code >= len(support)) or not np.array_equal(support[code], values):
            raise ValueError("fixed point value is absent from mapping support")
        return mapped[code]

    current, minimum48, minimum7 = map(transform, (current, minimum48, minimum7))
    absolute = current - minimum48 >= 3 * scale // 10
    ratio = 2 * current >= 3 * minimum7
    patient_absolute = np.zeros(n, dtype=bool)
    patient_ratio = np.zeros(n, dtype=bool)
    np.logical_or.at(patient_absolute, patient, absolute)
    np.logical_or.at(patient_ratio, patient, ratio)
    return _readonly((2 * patient_absolute.astype(np.int8) + patient_ratio).astype(np.int8))


def _state_vector(value, name):
    raw = np.asarray(value)
    if raw.dtype.kind == "b" and raw.ndim == 1:
        raw = raw.astype(np.int8)
    return _integer_vector(raw, name, maximum=3)


def four_by_four(original_states, transformed_states):
    """Rows are original 0/1/2/3 and columns are transformed 0/1/2/3."""
    original = _state_vector(original_states, "original_states")
    transformed = _state_vector(transformed_states, "transformed_states")
    if len(original) != len(transformed):
        raise ValueError("paired state vectors must have equal lengths")
    return _readonly(np.bincount(4 * original + transformed, minlength=16).reshape(4, 4))


def paired_table(original_state, transformed_state, mask=None):
    """Fold all states 1/2/3 as positive and return exact paired count cells.

    Rates are descriptive floats; empty subsets have None rates. Mask input
    must be a Boolean vector so integer indices cannot silently change N.
    """
    original = _state_vector(original_state, "original_state")
    transformed = _state_vector(transformed_state, "transformed_state")
    if len(original) != len(transformed):
        raise ValueError("paired state vectors must have equal lengths")
    if mask is not None:
        chosen = np.asarray(mask)
        if chosen.ndim != 1 or chosen.dtype.kind != "b" or len(chosen) != len(original):
            raise ValueError("mask must be a Boolean vector with one value per patient")
        original, transformed = original[chosen], transformed[chosen]
    code = 2 * (original > 0).astype(np.int8) + (transformed > 0)
    n00, n01, n10, n11 = map(int, np.bincount(code, minlength=4))
    n = len(original)
    result = dict(n=n, n00=n00, n01=n01, n10=n10, n11=n11,
                  discordance_count=n01 + n10, upward_count=n01,
                  downward_count=n10, net_change_count=n01 - n10,
                  original_positive_count=n10 + n11,
                  transformed_positive_count=n01 + n11)
    for key, numerator in (("discordance", n01 + n10), ("upward", n01),
                           ("downward", n10), ("net_change", n01 - n10),
                           ("original_positive_rate", n10 + n11),
                           ("transformed_positive_rate", n01 + n11)):
        result[key] = numerator / n if n else None
    return result
