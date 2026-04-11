"""
MALDI-MSI rescoring pipeline using MS2Rescore/mokapot.

Treats each MALDI m/z feature as a "spectrum" and each candidate peptide
(target or decoy) as a "PSM". Computes per-candidate features (mass accuracy,
spatial structure, protein consistency) and uses mokapot for semi-supervised
target-decoy FDR estimation.

Two decoy strategies are supported:
- Approach A: reverse LC-MS/MS reference peptides (K/R-preserving)
- Approach B: pass through unfiltered search engine target+decoy PSMs
"""

import logging
from typing import Optional

import numpy as np
import pandas as pd
from psm_utils import PSM, PSMList

from ms2rescore.maldi_fdr import (
    generate_decoy_peptides_from_reference,
    load_pd_unfiltered,
    match_features_to_database,
)

logger = logging.getLogger(__name__)


def candidates_to_psm_list(
    matches: pd.DataFrame,
    charge: int = 1,
) -> PSMList:
    """
    Convert a MALDI feature × candidate match table to a psm_utils.PSMList.

    Each row becomes a PSM where:
    - spectrum_id groups candidates by MALDI feature
    - peptidoform is the candidate peptide in ProForma notation
    - precursor_mz is the MALDI observed m/z
    - score is initialized from ppm error (for rank assignment)
    - is_decoy comes from the database label

    Parameters
    ----------
    matches
        Match table from match_features_to_database(). Must have columns:
        feature_mz, feature_idx, peptide, protein, is_decoy, ppm_error_abs,
        and optionally: xcorr, spatial_autocorrelation, fraction_detected, etc.
    charge
        Default charge state for MALDI [M+H]+ ions.

    Returns
    -------
    PSMList ready for feature generation and mokapot rescoring.
    """
    psms = []
    for i, (_, row) in enumerate(matches.iterrows()):
        metadata = {}
        # Store all match data in metadata for feature generators to use
        for col in matches.columns:
            if col not in ("peptide", "protein", "is_decoy"):
                val = row[col]
                if pd.notna(val):
                    metadata[col] = val

        psms.append(PSM(
            peptidoform=f"{row['peptide']}/{charge}",
            spectrum_id=f"maldi_{int(row.get('feature_idx', i))}",
            run="maldi",
            precursor_mz=float(row["feature_mz"]),
            retention_time=0.0,
            score=float(-row.get("ppm_error_abs", 0)),  # lower ppm = higher score
            is_decoy=bool(row["is_decoy"]),
            protein_list=[str(row.get("protein", ""))],
            metadata=metadata,
        ))

    psm_list = PSMList(psm_list=psms)
    logger.info(
        f"Built PSMList: {len(psm_list)} candidates "
        f"({sum(not p.is_decoy for p in psm_list)} target, "
        f"{sum(p.is_decoy for p in psm_list)} decoy) "
        f"from {matches['feature_mz'].nunique()} MALDI features"
    )
    return psm_list


def maldi_rescore(
    maldi_features_mz: np.ndarray,
    reference_peptides: Optional[pd.DataFrame] = None,
    msf_path: Optional[str] = None,
    spatial_features: Optional[pd.DataFrame] = None,
    ppm_tolerance: float = 20.0,
    output_path: str = "maldi_rescoring",
    decoy_approach: str = "A",
    train_fdr: float = 0.05,
) -> PSMList:
    """
    Run the full MALDI-MSI rescoring pipeline.

    Parameters
    ----------
    maldi_features_mz
        Array of MALDI feature m/z values.
    reference_peptides
        LC-MS/MS peptide identifications (Approach A). Must have
        Sequence, Accession, Mass columns.
    msf_path
        Path to ProteomeDiscoverer .msf file (Approach B).
    spatial_features
        Pre-computed spatial features DataFrame (optional).
    ppm_tolerance
        Mass tolerance for matching MALDI to database (ppm).
    output_path
        Prefix for output files.
    decoy_approach
        "A" (reverse targets) or "B" (pass-through PD decoys).
    train_fdr
        FDR threshold for mokapot training.

    Returns
    -------
    Rescored PSMList with q-values.
    """
    from ms2rescore.maldi_rescoring import MALDISpatialFeatureGenerator
    from ms2rescore.rescoring_engines import mokapot

    # Step 1: Build candidate database
    if decoy_approach == "A":
        if reference_peptides is None:
            raise ValueError("reference_peptides required for Approach A")
        database = generate_decoy_peptides_from_reference(reference_peptides)
    elif decoy_approach == "B":
        if msf_path is None:
            raise ValueError("msf_path required for Approach B")
        database = load_pd_unfiltered(msf_path)
    else:
        raise ValueError(f"Unknown decoy_approach: {decoy_approach}")

    # Step 2: Match MALDI features to database
    matches = match_features_to_database(
        maldi_features_mz, database, ppm_tolerance=ppm_tolerance
    )
    if len(matches) == 0:
        raise ValueError("No matches found. Check ppm_tolerance or m/z ranges.")

    # Merge spatial features if available
    if spatial_features is not None:
        matches = matches.merge(
            spatial_features, on="feature_mz", how="left"
        ).fillna(0)

    # Step 3: Convert to PSMList
    psm_list = candidates_to_psm_list(matches)

    # Step 4: Add rescoring features
    fgen = MALDISpatialFeatureGenerator(
        spatial_features=spatial_features,
        use_im=False,
    )
    fgen.add_features(psm_list)

    # Step 5: Rescore with mokapot
    all_feature_names = fgen.feature_names
    logger.info(f"Rescoring with {len(all_feature_names)} features using mokapot")

    mokapot.rescore(
        psm_list,
        output_file_root=output_path,
        train_fdr=train_fdr,
    )

    # Log results
    n_1pct = sum(1 for p in psm_list if p.qvalue is not None
                 and p.qvalue <= 0.01 and not p.is_decoy)
    n_5pct = sum(1 for p in psm_list if p.qvalue is not None
                 and p.qvalue <= 0.05 and not p.is_decoy)
    logger.info(f"Results: {n_1pct} at 1% FDR, {n_5pct} at 5% FDR")

    return psm_list
