"""
Target-decoy FDR estimation for MALDI-MSI peptide mass fingerprinting.

Two approaches for estimating FDR when identifying peptides from MS1 MALDI
imaging features:

**Option 4 (database decoys):**
  In silico digest of a protein FASTA, generating both target and decoy
  (reversed) peptides. Match MALDI features against both. Decoy matches
  estimate the false match rate. This is the standard proteomics approach
  adapted for MS1 mass matching.

**Option 1 (mass-shifted decoy features):**
  For each real MALDI feature, create synthetic "decoy features" by shifting
  the observed m/z. Search these decoy features against the target database.
  Any match is false by definition. This approach works without needing to
  modify the reference database and is useful when the reference comes from
  LC-MS/MS rather than in silico digest.
"""

import logging
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# In silico digest
# ---------------------------------------------------------------------------

def _reverse_protein_keep_cleavage(seq: str) -> str:
    """
    Reverse a protein sequence while keeping K and R at their original positions.

    This ensures the reversed protein produces the same number of tryptic peptides
    with the same length distribution as the target, so decoy proteins have
    realistic peptide counts for protein-level consistency scoring.
    """
    # Find positions of K and R
    kr_positions = {i: seq[i] for i in range(len(seq)) if seq[i] in "KR"}

    # Extract non-K/R residues, reverse them
    non_kr = [seq[i] for i in range(len(seq)) if seq[i] not in "KR"]
    non_kr.reverse()

    # Reconstruct: place K/R back at original positions, fill rest with reversed
    result = list(seq)
    j = 0
    for i in range(len(result)):
        if i not in kr_positions:
            result[i] = non_kr[j]
            j += 1

    return "".join(result)


def digest_fasta(
    fasta_path: str,
    enzyme: str = "trypsin",
    missed_cleavages: int = 1,
    min_length: int = 7,
    max_length: int = 25,
    generate_decoys: bool = True,
) -> pd.DataFrame:
    """
    In silico digest a FASTA file to generate target and decoy peptide masses.

    Decoy generation: reverses each protein sequence while keeping K/R at their
    original positions, then digests. This produces decoy proteins with the same
    number of tryptic peptides as the target, giving realistic protein_n_features
    values for target-decoy FDR estimation.

    Parameters
    ----------
    fasta_path
        Path to protein FASTA file.
    enzyme
        Enzyme for digestion. Currently only 'trypsin' (cleave after K/R, not before P).
    missed_cleavages
        Maximum number of missed cleavages.
    min_length, max_length
        Peptide length bounds.
    generate_decoys
        If True, reverse each protein (keeping K/R) and digest as decoy.

    Returns
    -------
    DataFrame with columns: peptide, protein, mass, mh_mz, is_decoy
    """
    from pyteomics import fasta, mass, parser

    def _digest_and_collect(seq, protein_id, is_decoy):
        results = []
        cleaved = parser.cleave(seq, parser.expasy_rules.get(enzyme, enzyme),
                                missed_cleavages=missed_cleavages)
        for pep in cleaved:
            if min_length <= len(pep) <= max_length:
                try:
                    comp = mass.Composition(sequence=pep)
                    pep_mass = mass.calculate_mass(composition=comp)
                    mh = pep_mass + 1.00727646677  # [M+H]+
                    results.append({
                        "peptide": pep,
                        "protein": protein_id,
                        "mass": pep_mass,
                        "mh_mz": mh,
                        "is_decoy": is_decoy,
                        "n_C": comp.get("C", 0),
                        "n_H": comp.get("H", 0),
                        "n_N": comp.get("N", 0),
                        "n_O": comp.get("O", 0),
                        "n_S": comp.get("S", 0),
                    })
                except Exception:
                    continue
        return results

    peptides = []

    for desc, seq in fasta.read(fasta_path):
        protein_id = desc.split("|")[1] if "|" in desc else desc.split()[0]

        # Target digest
        peptides.extend(_digest_and_collect(seq, protein_id, is_decoy=False))

        # Decoy: reverse protein keeping K/R in place, then digest
        if generate_decoys:
            decoy_seq = _reverse_protein_keep_cleavage(seq)
            peptides.extend(_digest_and_collect(
                decoy_seq, f"DECOY_{protein_id}", is_decoy=True
            ))

    df = pd.DataFrame(peptides).drop_duplicates(subset=["peptide", "is_decoy"])
    logger.info(
        f"Digested {fasta_path}: {(~df['is_decoy']).sum()} target, "
        f"{df['is_decoy'].sum()} decoy peptides"
    )
    return df


# ---------------------------------------------------------------------------
# Spatial features from imzML
# ---------------------------------------------------------------------------

def extract_spatial_features(
    imzml_path: str,
    feature_mzs: np.ndarray,
    ppm_tolerance: float = 10.0,
) -> pd.DataFrame:
    """
    Extract spatial features for each m/z feature from an imzML dataset.

    For each feature m/z, extracts the intensity image across all pixels
    and computes spatial quality metrics.

    Parameters
    ----------
    imzml_path
        Path to .imzML file (with .ibd in same directory).
    feature_mzs
        Array of m/z values to extract images for.
    ppm_tolerance
        Mass tolerance for peak matching in ppm.

    Returns
    -------
    DataFrame with spatial features per m/z:
    - n_pixels_detected: number of pixels where this m/z was observed
    - fraction_detected: fraction of total pixels with signal
    - mean_intensity: mean intensity across detected pixels
    - spatial_autocorrelation: Moran's I (higher = more spatially structured)
    - intensity_cv: coefficient of variation of intensities
    """
    from pyimzml.ImzMLParser import ImzMLParser

    logger.info(f"Extracting spatial features from {imzml_path}...")
    p = ImzMLParser(str(imzml_path))

    n_pixels = len(p.coordinates)
    coords = np.array([(c[0], c[1]) for c in p.coordinates])

    # Pre-read all spectra for speed
    logger.info(f"Reading {n_pixels} spectra...")
    all_mzs = []
    all_ints = []
    for i in range(n_pixels):
        mzs, ints = p.getspectrum(i)
        all_mzs.append(mzs)
        all_ints.append(ints)

    results = []
    for feat_idx, feat_mz in enumerate(feature_mzs):
        if feat_idx % 100 == 0:
            logger.debug(f"  Processing feature {feat_idx}/{len(feature_mzs)}")

        tol = feat_mz * ppm_tolerance / 1e6
        intensities = np.zeros(n_pixels)

        for i in range(n_pixels):
            mask = (all_mzs[i] >= feat_mz - tol) & (all_mzs[i] <= feat_mz + tol)
            if mask.any():
                intensities[i] = all_ints[i][mask].max()

        detected = intensities > 0
        n_detected = detected.sum()
        fraction = n_detected / n_pixels

        # Spatial autocorrelation (simplified Moran's I)
        spatial_autocorr = 0.0
        if n_detected >= 10:
            spatial_autocorr = _morans_i_fast(intensities, coords)

        # Intensity statistics
        if n_detected > 0:
            det_ints = intensities[detected]
            mean_int = det_ints.mean()
            cv = det_ints.std() / mean_int if mean_int > 0 else 0
        else:
            mean_int = 0
            cv = 0

        results.append({
            "feature_mz": feat_mz,
            "n_pixels_detected": n_detected,
            "fraction_detected": fraction,
            "mean_intensity": mean_int,
            "spatial_autocorrelation": spatial_autocorr,
            "intensity_cv": cv,
        })

    df = pd.DataFrame(results)
    logger.info(
        f"Extracted spatial features for {len(df)} m/z values "
        f"(median {df['n_pixels_detected'].median():.0f} pixels detected)"
    )
    return df


def _morans_i_fast(values: np.ndarray, coords: np.ndarray, k: int = 8) -> float:
    """
    Compute a fast approximation of Moran's I spatial autocorrelation.

    Uses k-nearest neighbors instead of full distance matrix for speed.
    Positive values indicate spatial clustering, negative = dispersal.
    """
    from scipy.spatial import cKDTree

    n = len(values)
    if n < k + 1:
        return 0.0

    mean_val = values.mean()
    deviations = values - mean_val
    var = np.sum(deviations**2)

    if var == 0:
        return 0.0

    tree = cKDTree(coords)
    _, indices = tree.query(coords, k=k + 1)

    # Vectorized cross-product computation
    neighbor_devs = deviations[indices[:, 1:]]  # exclude self
    cross_products = deviations[:, np.newaxis] * neighbor_devs
    cross_sum = cross_products.sum()
    w_sum = n * k

    morans_i = (n / w_sum) * (cross_sum / var)
    return float(morans_i)


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Isotope envelope features
# ---------------------------------------------------------------------------

def theoretical_isotope_distribution(n_C: int, n_H: int, n_N: int,
                                      n_O: int, n_S: int,
                                      n_peaks: int = 4) -> np.ndarray:
    """
    Compute theoretical isotope distribution using the Poisson approximation.

    Parameters
    ----------
    n_C, n_H, n_N, n_O, n_S
        Elemental counts from the peptide composition.
    n_peaks
        Number of isotope peaks to compute (M0, M1, ...).

    Returns
    -------
    np.ndarray
        Normalized isotope distribution [M0, M1, M2, ...].
    """
    from math import factorial, exp

    # Natural heavy-isotope probabilities
    lam = (n_C * 0.01109 + n_H * 0.000115 + n_N * 0.00364
           + n_O * 0.00205 + n_S * 0.04493)

    dist = np.array([exp(-lam) * lam**k / factorial(k) for k in range(n_peaks)])
    total = dist.sum()
    return dist / total if total > 0 else dist


# HIT-MAP isotope mass offsets:
# Natural heavy isotopes are at +delta from the monoisotopic.
# Decoy "reversed" isotopes are at -delta (fictional lighter isotopes).
NEUTRON = 1.003355  # 13C - 12C mass difference

# Per-element isotope mass offsets (heavy - mono) and abundances
_NATURAL_ISOTOPES = {
    # (mass_offset_from_mono, abundance_of_heavy)
    "C": [(NEUTRON, 0.01109)],
    "H": [(1.00628, 0.000115)],
    "N": [(0.99703, 0.00364)],
    "O": [(1.00422, 0.00038), (2.00425, 0.00205)],
    "S": [(0.99939, 0.0075), (1.99584, 0.0425)],
}


def theoretical_isotope_peaks(
    mono_mz: float,
    n_C: int, n_H: int, n_N: int, n_O: int, n_S: int,
    n_peaks: int = 4,
    reversed_isotopes: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute theoretical isotope peak positions and intensities.

    Parameters
    ----------
    mono_mz
        Monoisotopic [M+H]+ m/z.
    n_C, n_H, n_N, n_O, n_S
        Elemental counts.
    n_peaks
        Number of peaks (including M0).
    reversed_isotopes
        If True, place isotope peaks at M-1, M-2, ... instead of M+1, M+2, ...
        (HIT-MAP "isotope" decoy mode).

    Returns
    -------
    mz_positions : np.ndarray
        m/z values for M0, M+1 (or M-1), M+2 (or M-2), ...
    intensities : np.ndarray
        Normalized intensities.
    """
    intensities = theoretical_isotope_distribution(n_C, n_H, n_N, n_O, n_S, n_peaks)
    direction = -1.0 if reversed_isotopes else 1.0
    mz_positions = np.array([mono_mz + i * NEUTRON * direction for i in range(n_peaks)])
    return mz_positions, intensities


def generate_element_decoy_formula(
    n_C: int, n_H: int, n_N: int, n_O: int, n_S: int,
    target_mass: float,
    ppm_tolerance: float = 10.0,
    max_attempts: int = 50,
    rng: np.random.RandomState | None = None,
) -> tuple[int, int, int, int, int]:
    """
    Generate a decoy elemental composition with equivalent mass but different
    atom counts (HIT-MAP "element" decoy mode).

    Randomly perturbs the elemental composition while keeping mass within
    ppm tolerance. The different composition gives a different isotope pattern.

    Parameters
    ----------
    n_C, n_H, n_N, n_O, n_S
        Target elemental counts.
    target_mass
        Target monoisotopic mass.
    ppm_tolerance
        Mass window for decoy.
    max_attempts
        Max random attempts.
    rng
        Random state.

    Returns
    -------
    Tuple of (n_C, n_H, n_N, n_O, n_S) for the decoy.
    """
    from pyteomics import mass as pyteomics_mass

    if rng is None:
        rng = np.random.RandomState()

    tol = target_mass * ppm_tolerance / 1e6

    # Atom masses
    atom_masses = {"C": 12.0, "H": 1.00783, "N": 14.00307, "O": 15.99491, "S": 31.97207}
    original = {"C": n_C, "H": n_H, "N": n_N, "O": n_O, "S": n_S}

    for _ in range(max_attempts):
        decoy = dict(original)
        # Randomly swap atoms: e.g., replace some C+H with N, or O with C+H4, etc.
        # Strategy: perturb 2-3 elements by +/-1-3 and adjust H to compensate mass
        for elem in rng.choice(["C", "N", "O", "S"], size=rng.randint(1, 3), replace=False):
            delta = rng.randint(-3, 4)
            decoy[elem] = max(0, decoy[elem] + delta)

        # Adjust H to match target mass
        mass_no_h = sum(decoy[e] * atom_masses[e] for e in ["C", "N", "O", "S"])
        needed_h = round((target_mass - mass_no_h) / atom_masses["H"])
        if needed_h < 0:
            continue
        decoy["H"] = needed_h

        # Check mass
        decoy_mass = sum(decoy[e] * atom_masses[e] for e in decoy)
        ppm_err = abs(decoy_mass - target_mass) / target_mass * 1e6
        if ppm_err <= ppm_tolerance and decoy != original:
            return decoy["C"], decoy["H"], decoy["N"], decoy["O"], decoy["S"]

    # Fallback: return original (will have identical isotope pattern)
    return n_C, n_H, n_N, n_O, n_S


def score_sqrtp(
    observed_mz: np.ndarray,
    observed_int: np.ndarray,
    theoretical_mz: np.ndarray,
    theoretical_int: np.ndarray,
    ppm_tolerance: float = 10.0,
) -> float:
    """
    HIT-MAP SQRTP scoring function.

    Score = -log(sqrt(sum((u-v)^2) / (sum(u^2)*sum(v^2)))) + log(N_matched/N_theo)
            - |Phi(mean_ppm/ppm_tol) - 0.5|

    Where u = normalized theoretical, v = normalized observed intensities
    for matched isotope peaks.

    Parameters
    ----------
    observed_mz, observed_int
        Observed spectrum (from MALDI).
    theoretical_mz, theoretical_int
        Theoretical isotope pattern (target or decoy).
    ppm_tolerance
        Mass tolerance in ppm.

    Returns
    -------
    float
        SQRTP score. Higher = better match.
    """
    from scipy.stats import norm

    n_theo = len(theoretical_mz)
    if n_theo == 0:
        return -10.0

    # Match theoretical peaks to observed
    matched_theo = []
    matched_obs = []
    ppm_errors = []

    for t_mz, t_int in zip(theoretical_mz, theoretical_int):
        tol = t_mz * ppm_tolerance / 1e6
        mask = (observed_mz >= t_mz - tol) & (observed_mz <= t_mz + tol)
        if mask.any():
            best_idx = np.argmin(np.abs(observed_mz[mask] - t_mz))
            obs_idx = np.where(mask)[0][best_idx]
            matched_theo.append(t_int)
            matched_obs.append(observed_int[obs_idx])
            ppm_errors.append(abs(observed_mz[obs_idx] - t_mz) / t_mz * 1e6)

    n_matched = len(matched_theo)
    if n_matched == 0:
        return -10.0

    # Normalize
    u = np.array(matched_theo, dtype=np.float64)
    v = np.array(matched_obs, dtype=np.float64)
    u_sum = np.sqrt(np.sum(u**2))
    v_sum = np.sqrt(np.sum(v**2))
    if u_sum > 0:
        u = u / u_sum
    if v_sum > 0:
        v = v / v_sum

    # Similarity: -log(sqrt(sum((u-v)^2) / (sum(u^2)*sum(v^2))))
    u2 = np.sum(u**2)
    v2 = np.sum(v**2)
    diff2 = np.sum((u - v)**2)
    denom = u2 * v2
    if denom <= 0:
        similarity = 0.0
    else:
        ratio = diff2 / denom
        if ratio <= 0:
            similarity = 10.0  # perfect match
        else:
            similarity = -np.log(np.sqrt(ratio))

    # Match fraction bonus
    match_bonus = np.log(n_matched / n_theo) if n_theo > 0 else 0.0

    # ppm penalty
    mean_ppm = np.mean(ppm_errors)
    ppm_penalty = abs(norm.cdf(mean_ppm / ppm_tolerance) - 0.5)

    score = similarity + match_bonus - ppm_penalty
    return float(score)


def extract_observed_envelopes(
    imzml_path: str,
    feature_mzs: np.ndarray,
    ppm_tolerance: float = 30.0,
    n_isotope_peaks: int = 3,
) -> pd.DataFrame:
    """
    Extract observed isotope envelopes from imzML for each MALDI feature.

    For each feature m/z, looks for the monoisotopic peak and n-1 isotope
    peaks (M+1.003, M+2.006, ...) across all pixels, averages intensities.

    Parameters
    ----------
    imzml_path
        Path to .imzML file.
    feature_mzs
        Array of monoisotopic m/z values.
    ppm_tolerance
        Tolerance for matching peaks.
    n_isotope_peaks
        Number of isotope peaks to extract (including M0).

    Returns
    -------
    DataFrame with columns: feature_mz, obs_M0, obs_M1, obs_M2, ...
    """
    from pyimzml.ImzMLParser import ImzMLParser

    logger.info(f"Extracting isotope envelopes from {imzml_path}...")
    p = ImzMLParser(str(imzml_path))
    n_pixels = len(p.coordinates)

    # Pre-read all spectra
    logger.info(f"Reading {n_pixels} spectra...")
    all_mzs = []
    all_ints = []
    for i in range(n_pixels):
        mzs, ints = p.getspectrum(i)
        all_mzs.append(mzs)
        all_ints.append(ints)

    NEUTRON = 1.003355  # C13 - C12 mass difference

    results = []
    for feat_idx, feat_mz in enumerate(feature_mzs):
        if feat_idx % 200 == 0:
            logger.info(f"  Envelope {feat_idx}/{len(feature_mzs)}")

        # For each isotope peak (M0, M1, M2, ...)
        envelope_sums = np.zeros(n_isotope_peaks)
        envelope_counts = np.zeros(n_isotope_peaks)

        for i in range(n_pixels):
            for iso_idx in range(n_isotope_peaks):
                target_mz = feat_mz + iso_idx * NEUTRON
                tol = target_mz * ppm_tolerance / 1e6
                mask = (all_mzs[i] >= target_mz - tol) & (all_mzs[i] <= target_mz + tol)
                if mask.any():
                    envelope_sums[iso_idx] += all_ints[i][mask].max()
                    envelope_counts[iso_idx] += 1

        row = {"feature_mz": feat_mz}
        for iso_idx in range(n_isotope_peaks):
            row[f"obs_M{iso_idx}"] = (
                envelope_sums[iso_idx] / envelope_counts[iso_idx]
                if envelope_counts[iso_idx] > 0 else 0.0
            )
            row[f"obs_M{iso_idx}_n_pixels"] = int(envelope_counts[iso_idx])
        results.append(row)

    df = pd.DataFrame(results)
    logger.info(f"Extracted envelopes for {len(df)} features")
    return df


def score_candidates_sqrtp(
    matches: pd.DataFrame,
    observed_envelopes: pd.DataFrame,
    ppm_tolerance: float = 10.0,
    decoy_mode: str = "isotope",
    n_peaks: int = 4,
) -> pd.DataFrame:
    """
    Score each candidate using HIT-MAP's SQRTP method with isotope decoys.

    For each MALDI feature × candidate pair:
    - Compute theoretical isotope pattern from candidate's elemental composition
    - For decoys, use reversed isotope pattern (M-1, M-2 instead of M+1, M+2)
      or element-shuffled composition
    - Score observed spectrum against theoretical using SQRTP

    Parameters
    ----------
    matches
        Match table with n_C, n_H, n_N, n_O, n_S, feature_mz, is_decoy.
    observed_envelopes
        From extract_observed_envelopes() — observed M0, M1, M2 per feature.
    ppm_tolerance
        Mass tolerance for peak matching.
    decoy_mode
        "isotope" (reversed isotope direction) or "element" (shuffled composition).
    n_peaks
        Number of isotope peaks to use.

    Returns
    -------
    matches with added 'sqrtp_score' column.
    """
    df = matches.copy()
    rng = np.random.RandomState(42)

    # Build observed spectra lookup: feature_mz -> (mz_array, int_array)
    obs_lookup = {}
    for _, row in observed_envelopes.iterrows():
        fmz = row["feature_mz"]
        mzs = []
        ints = []
        for i in range(n_peaks):
            if f"obs_M{i}" in row and row[f"obs_M{i}"] > 0:
                mzs.append(fmz + i * NEUTRON)
                ints.append(row[f"obs_M{i}"])
        if mzs:
            obs_lookup[fmz] = (np.array(mzs), np.array(ints))

    scores = []
    for idx, row in df.iterrows():
        fmz = row["feature_mz"]

        if fmz not in obs_lookup:
            scores.append(-10.0)
            continue

        obs_mz, obs_int = obs_lookup[fmz]

        nc = int(row.get("n_C", 0))
        nh = int(row.get("n_H", 0))
        nn = int(row.get("n_N", 0))
        no = int(row.get("n_O", 0))
        ns = int(row.get("n_S", 0))

        if row["is_decoy"]:
            if decoy_mode == "isotope":
                # Reversed isotope pattern
                theo_mz, theo_int = theoretical_isotope_peaks(
                    fmz, nc, nh, nn, no, ns, n_peaks, reversed_isotopes=True
                )
            elif decoy_mode == "element":
                # Shuffled elemental composition
                dc, dh, dn, do, ds = generate_element_decoy_formula(
                    nc, nh, nn, no, ns,
                    target_mass=row.get("mass", fmz - 1.00728),
                    ppm_tolerance=ppm_tolerance, rng=rng,
                )
                theo_mz, theo_int = theoretical_isotope_peaks(
                    fmz, dc, dh, dn, do, ds, n_peaks, reversed_isotopes=False
                )
            else:
                theo_mz, theo_int = theoretical_isotope_peaks(
                    fmz, nc, nh, nn, no, ns, n_peaks
                )
        else:
            # Target: normal isotope pattern
            theo_mz, theo_int = theoretical_isotope_peaks(
                fmz, nc, nh, nn, no, ns, n_peaks
            )

        score = score_sqrtp(obs_mz, obs_int, theo_mz, theo_int, ppm_tolerance)
        scores.append(score)

    df["sqrtp_score"] = scores

    # Stats
    t = df[~df["is_decoy"]]["sqrtp_score"]
    d = df[df["is_decoy"]]["sqrtp_score"]
    logger.info(
        f"SQRTP scores ({decoy_mode} decoys): "
        f"target mean={t.mean():.3f}, decoy mean={d.mean():.3f}, "
        f"separation={t.mean()-d.mean():.3f}"
    )
    return df


def compute_envelope_quality(
    observed_envelopes: pd.DataFrame,
    n_peaks: int = 3,
) -> pd.DataFrame:
    """
    Compute a feature-level envelope quality score.

    Checks whether the observed isotope pattern looks like a real peptide:
    - M0 must be present
    - M1 must be present with M1/M0 in a reasonable range for peptides
    - M2 should be present for larger peptides

    This is NOT candidate-specific — it tells us which MALDI features are
    likely real peptide signals vs noise/matrix.

    Returns
    -------
    DataFrame with added columns: envelope_quality (0-1), M1_M0_ratio.
    """
    df = observed_envelopes.copy()

    df["M1_M0_ratio"] = np.where(
        df["obs_M0"] > 0, df["obs_M1"] / df["obs_M0"], 0.0
    )
    df["M2_M0_ratio"] = np.where(
        df["obs_M0"] > 0, df["obs_M2"] / df["obs_M0"], 0.0
    )

    # Quality score: M0 present, M1/M0 in peptide range [0.2, 1.5]
    has_m0 = df["obs_M0"] > 0
    has_m1 = df["obs_M1"] > 0
    m1_ok = (df["M1_M0_ratio"] > 0.2) & (df["M1_M0_ratio"] < 1.5)

    # Continuous quality: how peptide-like is the M1/M0 ratio?
    # Ideal range for tryptic peptides: M1/M0 ~ 0.35 (800 Da) to 1.0 (2000 Da)
    # Penalize values far outside this range
    ideal_center = df["feature_mz"] * 0.0005  # rough scaling
    ideal_center = ideal_center.clip(lower=0.3, upper=1.1)
    ratio_deviation = np.abs(df["M1_M0_ratio"] - ideal_center)
    df["envelope_quality"] = np.where(
        has_m0 & has_m1,
        np.exp(-ratio_deviation ** 2 / 0.1),  # Gaussian penalty
        0.0,
    )

    return df


def compute_isotope_scores(
    matches: pd.DataFrame,
    observed_envelopes: pd.DataFrame,
    n_peaks: int = 3,
) -> pd.DataFrame:
    """
    Compute isotope envelope match score for each candidate.

    For each candidate peptide, computes the theoretical isotope distribution
    from its elemental composition and compares it to the observed envelope
    using cosine similarity.

    Parameters
    ----------
    matches
        Match table with columns: feature_mz, n_C, n_H, n_N, n_O, n_S.
    observed_envelopes
        DataFrame from extract_observed_envelopes().
    n_peaks
        Number of isotope peaks to use for scoring.

    Returns
    -------
    matches DataFrame with added 'isotope_score' column.
    """
    df = matches.copy()

    # Merge observed envelopes
    env_map = observed_envelopes.set_index("feature_mz")
    obs_cols = [f"obs_M{i}" for i in range(n_peaks)]

    for col in obs_cols:
        df[col] = df["feature_mz"].map(env_map[col]).fillna(0)

    # Compute theoretical distribution and cosine similarity per candidate
    scores = []
    for _, row in df.iterrows():
        # Theoretical envelope from elemental composition
        if all(col in row.index for col in ["n_C", "n_H", "n_N", "n_O", "n_S"]):
            theo = theoretical_isotope_distribution(
                int(row.get("n_C", 0)), int(row.get("n_H", 0)),
                int(row.get("n_N", 0)), int(row.get("n_O", 0)),
                int(row.get("n_S", 0)), n_peaks=n_peaks,
            )
        else:
            # Fallback: averagine approximation from mass
            n_c_approx = int(row["mass"] * 0.0444)
            theo = theoretical_isotope_distribution(
                n_c_approx, n_c_approx * 2, int(n_c_approx * 0.3),
                int(n_c_approx * 0.4), 0, n_peaks=n_peaks,
            )

        # Observed envelope
        obs = np.array([row[f"obs_M{i}"] for i in range(n_peaks)])

        # Cosine similarity
        if np.linalg.norm(obs) > 0 and np.linalg.norm(theo) > 0:
            cosine = np.dot(obs, theo) / (np.linalg.norm(obs) * np.linalg.norm(theo))
            scores.append(float(cosine))
        else:
            scores.append(0.0)

    df["isotope_score"] = scores

    n_scored = sum(1 for s in scores if s > 0)
    logger.info(
        f"Computed isotope scores for {n_scored}/{len(df)} candidates "
        f"(mean={np.mean([s for s in scores if s > 0]):.3f})"
    )
    return df


# ---------------------------------------------------------------------------
# Option 4: Match MALDI features to target+decoy database
# ---------------------------------------------------------------------------

def match_features_to_database(
    features_mz: np.ndarray,
    database: pd.DataFrame,
    ppm_tolerance: float = 10.0,
    features_intensity: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """
    Match MALDI m/z features against a target+decoy peptide database.

    Parameters
    ----------
    features_mz
        Array of MALDI feature m/z values.
    database
        Peptide database from digest_fasta() with 'mh_mz' and 'is_decoy' columns.
    ppm_tolerance
        Mass tolerance in ppm.
    features_intensity
        Optional array of MALDI feature intensities.

    Returns
    -------
    DataFrame with all matches, including target/decoy labels.
    """
    matches = []
    for i, mz in enumerate(features_mz):
        tol = mz * ppm_tolerance / 1e6
        candidates = database[
            (database["mh_mz"] >= mz - tol) &
            (database["mh_mz"] <= mz + tol)
        ].copy()

        if len(candidates) == 0:
            continue

        candidates["feature_mz"] = mz
        candidates["feature_idx"] = i
        candidates["ppm_error"] = (mz - candidates["mh_mz"]) / candidates["mh_mz"] * 1e6
        candidates["ppm_error_abs"] = candidates["ppm_error"].abs()
        if features_intensity is not None:
            candidates["feature_intensity"] = features_intensity[i]

        matches.append(candidates)

    if not matches:
        return pd.DataFrame()

    result = pd.concat(matches, ignore_index=True)

    # Add protein-level consistency (computed separately for targets and decoys
    # so that decoy proteins get realistic protein_n_features values)
    protein_feature_count = result.groupby("protein")["feature_mz"].nunique()
    result["protein_n_features"] = result["protein"].map(protein_feature_count).fillna(0)

    # Count candidates per feature
    result["n_candidates"] = result.groupby("feature_mz")["feature_mz"].transform("count")

    logger.info(
        f"Matched {result['feature_mz'].nunique()}/{len(features_mz)} features to "
        f"{(~result['is_decoy']).sum()} target + {result['is_decoy'].sum()} decoy candidates"
    )
    return result


# ---------------------------------------------------------------------------
# Option 1: Mass-shifted decoy features
# ---------------------------------------------------------------------------

def create_decoy_features(
    features_mz: np.ndarray,
    shift_da: float = 10.5,
    n_shifts: int = 1,
) -> np.ndarray:
    """
    Create decoy MALDI features by shifting observed m/z values.

    Parameters
    ----------
    features_mz
        Array of real MALDI feature m/z values.
    shift_da
        Fixed mass shift in Da. 10.5 Da avoids overlap with real peptides.
    n_shifts
        Number of shifted versions to create.

    Returns
    -------
    Array of decoy m/z values.
    """
    decoy_mz = []
    for i in range(1, n_shifts + 1):
        decoy_mz.append(features_mz + i * shift_da)
    return np.concatenate(decoy_mz)


def match_with_decoy_features(
    features_mz: np.ndarray,
    reference: pd.DataFrame,
    ppm_tolerance: float = 10.0,
    shift_da: float = 10.5,
    features_intensity: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """
    Match real and decoy MALDI features against a reference peptide list.

    The reference can be LC-MS/MS identifications (no need for decoy peptides
    in the database — the decoy features provide the null distribution).

    Parameters
    ----------
    features_mz
        Real MALDI feature m/z values.
    reference
        LC-MS/MS reference with 'mh_mz' and 'peptide' columns.
    ppm_tolerance
        Mass tolerance in ppm.
    shift_da
        Mass shift for decoy features in Da.
    features_intensity
        Optional intensities for real features.

    Returns
    -------
    DataFrame with matches, labeled as target or decoy feature.
    """
    decoy_mz = create_decoy_features(features_mz, shift_da=shift_da)

    # Ensure reference has is_decoy=False (these are real peptides)
    ref = reference.copy()
    if "is_decoy" not in ref.columns:
        ref["is_decoy"] = False

    # Match real features
    target_matches = match_features_to_database(
        features_mz, ref, ppm_tolerance, features_intensity
    )
    if len(target_matches) > 0:
        target_matches["is_decoy_feature"] = False

    # Match decoy features (any match is false by definition)
    decoy_matches = match_features_to_database(
        decoy_mz, ref, ppm_tolerance
    )
    if len(decoy_matches) > 0:
        decoy_matches["is_decoy_feature"] = True

    result = pd.concat([target_matches, decoy_matches], ignore_index=True)

    n_target = target_matches["feature_mz"].nunique() if len(target_matches) > 0 else 0
    n_decoy = decoy_matches["feature_mz"].nunique() if len(decoy_matches) > 0 else 0
    logger.info(
        f"Matched {n_target} target features, {n_decoy} decoy features "
        f"(shift={shift_da} Da)"
    )
    return result


# ---------------------------------------------------------------------------
# Scoring and FDR calculation
# ---------------------------------------------------------------------------

def compute_match_score(
    matches: pd.DataFrame,
    weights: Optional[dict] = None,
    spatial_features: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """
    Compute a composite score for each candidate match.

    Score components (all normalized to [0, 1]):
    - Mass accuracy: 1 - ppm_error_abs / max_ppm
    - Protein consistency: log(1 + protein_n_features) / max_log
    - Candidate uniqueness: 1 / n_candidates
    - Spatial quality: fraction_detected * spatial_autocorrelation (when available)

    Parameters
    ----------
    matches
        Match table from match_features_to_database or match_with_decoy_features.
    weights
        Optional weights for score components. Default: equal weights.
    spatial_features
        Optional DataFrame from extract_spatial_features(), keyed by feature_mz.

    Returns
    -------
    matches with added 'score' column.
    """
    df = matches.copy()

    if weights is None:
        weights = {
            "mass_accuracy": 1.0,
            "protein_consistency": 1.0,
            "uniqueness": 0.5,
            "spatial": 1.0,
        }

    # Mass accuracy score (higher = better)
    max_ppm = df["ppm_error_abs"].max()
    if max_ppm > 0:
        df["score_mass"] = 1.0 - df["ppm_error_abs"] / max_ppm
    else:
        df["score_mass"] = 1.0

    # Protein consistency score
    log_prot = np.log1p(df["protein_n_features"])
    max_log = log_prot.max()
    df["score_protein"] = log_prot / max_log if max_log > 0 else 0.0

    # Uniqueness (fewer candidates = more confident)
    df["score_uniqueness"] = 1.0 / df["n_candidates"]

    # Spatial quality score
    if spatial_features is not None and "spatial" in weights and weights["spatial"] > 0:
        # Merge spatial features by feature_mz
        spatial_map = spatial_features.set_index("feature_mz")
        df["fraction_detected"] = df["feature_mz"].map(
            spatial_map["fraction_detected"]
        ).fillna(0)
        df["spatial_autocorrelation"] = df["feature_mz"].map(
            spatial_map["spatial_autocorrelation"]
        ).fillna(0)

        # Spatial score: features that are detected in many pixels AND show
        # spatial structure (high Moran's I) are more likely real peptides.
        # Clip Moran's I to [0, 1] range
        moran_clipped = df["spatial_autocorrelation"].clip(lower=0)
        max_moran = moran_clipped.max()
        if max_moran > 0:
            df["score_spatial"] = (
                df["fraction_detected"] * 0.3
                + (moran_clipped / max_moran) * 0.7
            )
        else:
            df["score_spatial"] = df["fraction_detected"]
    else:
        df["score_spatial"] = 0.0
        if "spatial" in weights:
            weights = {k: v for k, v in weights.items() if k != "spatial"}

    # Isotope envelope score (candidate-specific!)
    if "isotope_score" in df.columns and "isotope" in weights and weights["isotope"] > 0:
        df["score_isotope"] = df["isotope_score"].clip(lower=0)
    else:
        df["score_isotope"] = 0.0
        if "isotope" in weights:
            weights = {k: v for k, v in weights.items() if k != "isotope"}

    # Composite score
    all_components = ["mass", "protein", "uniqueness", "spatial", "isotope"]
    total_weight = sum(weights.get(k, 0) for k in all_components)
    if total_weight == 0:
        total_weight = 1.0
    df["score"] = sum(
        weights.get(k, 0) * df[f"score_{k}"]
        for k in all_components
        if k in weights
    ) / total_weight

    return df


def compute_fdr(
    matches: pd.DataFrame,
    score_col: str = "score",
    decoy_col: str = "is_decoy",
) -> pd.DataFrame:
    """
    Compute q-values using target-decoy competition.

    For each score threshold, FDR = n_decoy_above / n_target_above.

    Parameters
    ----------
    matches
        Match table with score and decoy label columns.
        For Option 4: decoy_col = 'is_decoy' (decoy peptides in database).
        For Option 1: decoy_col = 'is_decoy_feature' (decoy m/z features).
    score_col
        Column name for the match score (higher = better).
    decoy_col
        Column name for the decoy label.

    Returns
    -------
    matches with added 'qvalue' column.
    """
    df = matches.copy()

    # Take best candidate per feature (highest score)
    best = df.loc[df.groupby("feature_mz")[score_col].idxmax()].copy()
    best = best.sort_values(score_col, ascending=False).reset_index(drop=True)

    # Cumulative target and decoy counts
    is_decoy = best[decoy_col].values.astype(bool)
    cum_decoy = np.cumsum(is_decoy)
    cum_target = np.cumsum(~is_decoy)

    # FDR at each threshold
    fdr = np.where(cum_target > 0, cum_decoy / cum_target, 0.0)

    # Convert to q-values (monotonize: each q-value >= the minimum FDR at lower ranks)
    qvalues = np.minimum.accumulate(fdr[::-1])[::-1]
    best["qvalue"] = qvalues

    n_1pct = ((best["qvalue"] <= 0.01) & (~best[decoy_col])).sum()
    n_5pct = ((best["qvalue"] <= 0.05) & (~best[decoy_col])).sum()
    logger.info(
        f"FDR results: {n_1pct} target hits at 1% FDR, "
        f"{n_5pct} at 5% FDR (of {(~best[decoy_col]).sum()} total target)"
    )

    return best


# ---------------------------------------------------------------------------
# Approach A: Generate decoy peptides from LC-MS/MS reference
# ---------------------------------------------------------------------------

def generate_decoy_peptides_from_reference(
    reference_df: pd.DataFrame,
    peptide_col: str = "Sequence",
    protein_col: str = "Accession",
    mass_col: str = "Mass",
) -> pd.DataFrame:
    """
    Generate decoy peptides from an LC-MS/MS reference by reversing at the
    protein level (K/R-preserving).

    For each protein, all its peptide sequences are reversed while keeping
    K and R at their original positions. This produces decoy proteins with
    the same number of peptides and similar length distributions.

    Parameters
    ----------
    reference_df
        LC-MS/MS peptide identifications with sequence, protein, and mass.
    peptide_col, protein_col, mass_col
        Column names.

    Returns
    -------
    Combined target+decoy DataFrame with 'is_decoy' flag and 'mh_mz' column.
    """
    from pyteomics import mass as pyteomics_mass

    # Build target entries
    targets = reference_df[[peptide_col, protein_col, mass_col]].copy()
    targets = targets.rename(columns={peptide_col: "peptide", protein_col: "protein", mass_col: "mass"})
    targets["mh_mz"] = targets["mass"] + 1.00728
    targets["is_decoy"] = False

    # Generate decoys: reverse each peptide (K/R preserving)
    # Reversed peptides have the SAME amino acid composition → SAME mass.
    # We copy the target mass directly (not recomputing from sequence, which
    # would give the theoretical monoisotopic mass instead of the observed mass).
    decoys = []
    for _, row in targets.drop_duplicates(subset=["peptide", "protein"]).iterrows():
        rev_seq = _reverse_peptide_keep_terminal(row["peptide"])
        if rev_seq == row["peptide"]:
            continue  # skip palindromes
        decoys.append({
            "peptide": rev_seq,
            "protein": f"DECOY_{row['protein']}",
            "mass": row["mass"],       # same mass as target
            "mh_mz": row["mh_mz"],     # same [M+H]+ as target
            "is_decoy": True,
        })

    decoy_df = pd.DataFrame(decoys)
    combined = pd.concat([targets, decoy_df], ignore_index=True)
    combined = combined.drop_duplicates(subset=["peptide", "is_decoy"])

    logger.info(
        f"Generated reference database: {(~combined['is_decoy']).sum()} target, "
        f"{combined['is_decoy'].sum()} decoy peptides"
    )
    return combined


def _reverse_peptide_keep_terminal(seq: str) -> str:
    """Reverse a peptide sequence, keeping the C-terminal residue (K/R) in place."""
    if len(seq) <= 1:
        return seq
    if seq[-1] in "KR":
        return seq[-2::-1] + seq[-1]
    return seq[::-1]


# ---------------------------------------------------------------------------
# Approach B: Load unfiltered ProteomeDiscoverer results
# ---------------------------------------------------------------------------

def load_pd_unfiltered(msf_path: str) -> pd.DataFrame:
    """
    Load ALL target and decoy PSMs from a ProteomeDiscoverer .msf file.

    Returns a unified DataFrame with both targets and decoys, unfiltered,
    including all search engine scores. Each row is a peptide identification
    that can be matched to MALDI features by m/z.

    Parameters
    ----------
    msf_path
        Path to .msf file (SQLite database).

    Returns
    -------
    DataFrame with columns: peptide, protein, mass, mh_mz, is_decoy,
    xcorr, percolator_score, charge, and any available scores.
    """
    import sqlite3

    conn = sqlite3.connect(msf_path)

    # Target PSMs with protein mapping
    targets = pd.read_sql("""
        SELECT tp.Sequence as peptide, tp.Mass as mass, tp.Charge as charge,
               tp.XCorr as xcorr, tp.PercolatorqValue as percolator_qvalue,
               tp.PercolatorSVMScore as percolator_score,
               tp.DeltaMassInPPM as delta_ppm,
               tp.MissedCleavages as missed_cleavages,
               tp.Modifications as modifications,
               tprot.Accession as protein, tprot.Description as protein_desc
        FROM TargetPsms tp
        LEFT JOIN TargetProteinsTargetPsms tptp
            ON tp.PeptideID = tptp.TargetPsmsPeptideID
        LEFT JOIN TargetProteins tprot
            ON tptp.TargetProteinsUniqueSequenceID = tprot.UniqueSequenceID
    """, conn)
    targets["is_decoy"] = False

    # Decoy PSMs with protein mapping
    decoys = pd.read_sql("""
        SELECT dp.Sequence as peptide, dp.Mass as mass, dp.Charge as charge,
               dp.XCorr as xcorr, dp.PercolatorqValue as percolator_qvalue,
               dp.PercolatorSVMScore as percolator_score,
               dp.DeltaMassInPPM as delta_ppm,
               dp.MissedCleavages as missed_cleavages,
               dp.Modifications as modifications,
               dprot.Accession as protein, dprot.Description as protein_desc
        FROM DecoyPsms dp
        LEFT JOIN DecoyProteinsDecoyPsms dptp
            ON dp.PeptideID = dptp.DecoyPsmsPeptideID
        LEFT JOIN DecoyProteins dprot
            ON dptp.DecoyProteinsUniqueSequenceID = dprot.UniqueSequenceID
    """, conn)
    decoys["is_decoy"] = True

    conn.close()

    combined = pd.concat([targets, decoys], ignore_index=True)
    combined["mh_mz"] = combined["mass"] + 1.00728

    # Compute spectral count per peptide before deduplication
    spec_counts = combined.groupby(["peptide", "is_decoy"]).size().reset_index(name="spectral_count")
    # Aggregate: max xcorr, min PEP, mean percolator_score
    agg = combined.groupby(["peptide", "is_decoy"]).agg(
        max_xcorr=("xcorr", "max"),
        mean_xcorr=("xcorr", "mean"),
        min_pep=("percolator_qvalue", "min"),
    ).reset_index()

    # Deduplicate: keep best-scoring PSM per peptide
    combined = combined.sort_values("xcorr", ascending=False).drop_duplicates(
        subset=["peptide", "is_decoy"], keep="first"
    )

    # Merge spectral count and aggregated scores
    combined = combined.merge(spec_counts, on=["peptide", "is_decoy"], how="left")
    combined = combined.merge(agg, on=["peptide", "is_decoy"], how="left")

    logger.info(
        f"Loaded PD results: {(~combined['is_decoy']).sum()} target, "
        f"{combined['is_decoy'].sum()} decoy peptides "
        f"(XCorr range: {combined['xcorr'].min():.2f} - {combined['xcorr'].max():.2f})"
    )
    return combined


# ---------------------------------------------------------------------------
# Isotope envelope extraction for MALDI vs LC-MS/MS comparison
# ---------------------------------------------------------------------------

def extract_lcms_envelopes(
    msf_path: str,
    n_peaks: int = 3,
) -> dict[str, np.ndarray]:
    """
    Extract per-peptide isotope envelopes from a ProteomeDiscoverer .msf file.

    For each peptide with an LC-MS/MS feature, extracts the isotope peak
    heights (M0, M1, M2, ...) and normalizes them to sum to 1.

    Parameters
    ----------
    msf_path
        Path to ProteomeDiscoverer .msf SQLite file.
    n_peaks
        Number of isotope peaks to extract.

    Returns
    -------
    dict mapping peptide sequence -> normalized envelope array [M0, M1, M2, ...].
    Only peptides with LC-MS/MS features are included.
    """
    import sqlite3

    conn = sqlite3.connect(msf_path)

    envelopes_raw = pd.read_sql("""
        SELECT tp.Sequence as peptide,
               lf.MonoisotopicMassOverCharge as mono_mz,
               lf.ChargeState as charge,
               lp.MassOverCharge as peak_mz,
               lp.PeakHeight as height
        FROM TargetPsms tp
        JOIN LcmsFeaturesTargetPsms lft ON tp.PeptideID = lft.TargetPsmsPeptideID
        JOIN LcmsFeatures lf ON lft.LcmsFeaturesId = lf.Id
        JOIN LcmsFeaturesLcmsPeaks lflp ON lf.Id = lflp.LcmsFeaturesId
        JOIN LcmsPeaks lp ON lflp.LcmsPeaksId = lp.Id
        ORDER BY tp.Sequence, lp.MassOverCharge
    """, conn)
    conn.close()

    NEUTRON = 1.003355
    envelopes = {}

    for peptide, group in envelopes_raw.groupby("peptide"):
        mono_mz = group.iloc[0]["mono_mz"]
        charge = group.iloc[0]["charge"]

        peaks = np.zeros(n_peaks)
        for _, row in group.iterrows():
            delta_mass = (row["peak_mz"] - mono_mz) * charge
            iso_idx = round(delta_mass / NEUTRON)
            if 0 <= iso_idx < n_peaks:
                peaks[iso_idx] = max(peaks[iso_idx], row["height"])

        total = peaks.sum()
        if total > 0:
            envelopes[peptide] = peaks / total

    logger.info(f"Extracted LC-MS/MS isotope envelopes for {len(envelopes)} peptides")
    return envelopes


def extract_maldi_envelopes_from_tsf(
    tsf_path: str,
    feature_mzs: np.ndarray,
    n_peaks: int = 3,
    ppm_tolerance: float = 20.0,
    n_sample_pixels: int = 100,
    ion_images: Optional[np.ndarray] = None,
    ion_image_mzs: Optional[np.ndarray] = None,
) -> dict[float, np.ndarray]:
    """
    Extract isotope envelopes from MALDI TSF data for each feature m/z.

    For each feature, selects pixels WHERE the feature is detected (from
    ion images) and averages isotope peak intensities across those pixels.
    This avoids diluting signal with noise from pixels where the peptide
    is absent.

    Parameters
    ----------
    tsf_path
        Path to Bruker .d directory containing analysis.tsf.
    feature_mzs
        Array of monoisotopic m/z values.
    n_peaks
        Number of isotope peaks to extract.
    ppm_tolerance
        Tolerance for peak matching.
    n_sample_pixels
        Max pixels to sample per feature (from detected pixels).
    ion_images
        3D array (n_features, height, width) of ion images. If provided,
        only pixels where the feature is detected are sampled.
    ion_image_mzs
        Array of m/z values corresponding to ion_images.

    Returns
    -------
    dict mapping feature_mz -> normalized envelope array.
    """
    import imzy

    reader = imzy.get_reader(tsf_path)
    n_pixels = reader.n_pixels
    NEUTRON = 1.003355
    rng = np.random.RandomState(42)

    # Build m/z → image index lookup for smart pixel selection
    mz_to_img_idx = {}
    if ion_images is not None and ion_image_mzs is not None:
        mz_to_img_idx = {mz: i for i, mz in enumerate(ion_image_mzs)}

    # Pre-determine which pixels to read for each feature
    feature_pixel_map = {}
    all_pixels_needed = set()

    for feat_mz in feature_mzs:
        if feat_mz in mz_to_img_idx:
            # Select pixels where this feature is detected
            img = ion_images[mz_to_img_idx[feat_mz]]
            detected_flat = np.where(img.flatten() > 0)[0]
            if len(detected_flat) > n_sample_pixels:
                # Sample brightest pixels preferentially
                flat_ints = img.flatten()[detected_flat]
                top_idx = np.argsort(flat_ints)[-n_sample_pixels:]
                pixel_indices = detected_flat[top_idx]
            elif len(detected_flat) > 0:
                pixel_indices = detected_flat
            else:
                pixel_indices = rng.choice(n_pixels, size=min(n_sample_pixels, n_pixels), replace=False)
        else:
            # No ion image — fall back to random sampling
            pixel_indices = rng.choice(n_pixels, size=min(n_sample_pixels, n_pixels), replace=False)

        feature_pixel_map[feat_mz] = pixel_indices
        all_pixels_needed.update(pixel_indices)

    # Read spectra for all needed pixels (cache to avoid re-reading)
    logger.info(f"Reading {len(all_pixels_needed)} unique pixels for envelope extraction...")
    pixel_spectra = {}
    for px in all_pixels_needed:
        mzs, ints = reader.get_spectrum(px)
        pixel_spectra[px] = (mzs, ints)
    reader.close()

    # Extract envelopes per feature from its selected pixels
    envelopes = {}
    for feat_mz in feature_mzs:
        pixel_indices = feature_pixel_map[feat_mz]
        envelope_sums = np.zeros(n_peaks)
        envelope_counts = np.zeros(n_peaks)

        for px in pixel_indices:
            if px not in pixel_spectra:
                continue
            mzs, ints = pixel_spectra[px]
            for iso_idx in range(n_peaks):
                target_mz = feat_mz + iso_idx * NEUTRON
                tol = target_mz * ppm_tolerance / 1e6
                mask = (mzs >= target_mz - tol) & (mzs <= target_mz + tol)
                if mask.any():
                    envelope_sums[iso_idx] += ints[mask].max()
                    envelope_counts[iso_idx] += 1

        avg = np.where(envelope_counts > 0, envelope_sums / envelope_counts, 0.0)
        total = avg.sum()
        if total > 0:
            envelopes[feat_mz] = avg / total

    n_valid = sum(1 for v in envelopes.values() if v is not None)
    n_smart = sum(1 for mz in feature_mzs if mz in mz_to_img_idx)
    logger.info(
        f"Extracted MALDI isotope envelopes for {n_valid}/{len(feature_mzs)} features "
        f"({n_smart} with smart pixel selection, {len(feature_mzs)-n_smart} random)"
    )
    return envelopes
