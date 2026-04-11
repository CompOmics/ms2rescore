"""
MALDI-MSI rescoring module.

Provides feature generation and scoring for MALDI-MSI peptide identification
from MS1 feature matching. Works with or without ion mobility data.

Features used for rescoring:
- Mass accuracy (ppm error between MALDI m/z and theoretical peptide mass)
- CCS agreement (predicted vs observed, when ion mobility is available)
- LC-MS/MS search quality (probability, expectation, spectral count)
- Peptide properties (length, number of candidates)
- MALDI signal quality (intensity)

The framework is designed to work within the ms2rescore philosophy:
each candidate gets a feature vector, and a semi-supervised classifier
(e.g., mokapot/Percolator) separates true from false matches.
"""

import logging
from typing import Optional

import numpy as np
import pandas as pd
from psm_utils import PSMList

from ms2rescore.feature_generators.base import FeatureGeneratorBase

logger = logging.getLogger(__name__)


class MALDIFeatureGenerator:
    """
    Generate rescoring features for MALDI-MSI peptide candidates.

    Each MALDI feature (m/z) may match multiple LC-MS/MS peptide candidates.
    This generator computes per-candidate features that help distinguish
    the correct assignment from incorrect ones.

    Parameters
    ----------
    use_im : bool
        Whether to compute ion mobility-based features (requires CCS data).
    """

    def __init__(self, use_im: bool = True):
        self.use_im = use_im

    @property
    def feature_names(self) -> list[str]:
        names = [
            # Mass accuracy features
            "ppm_error",
            "ppm_error_abs",
            "ppm_rank",               # Rank of this candidate by ppm within group
            "ppm_best_ratio",         # This candidate's ppm / best ppm in group
            # Candidate ambiguity
            "n_candidates",
            "log_n_candidates",
            # LC-MS/MS search quality
            "probability",
            "log_expectation",
            "spectral_count",
            "log_spectral_count",
            # Peptide properties
            "peptide_length",
            "n_missed_cleavages",
            "has_modifications",
            # MALDI signal
            "log_maldi_intensity",
            # Protein-level consistency
            "protein_n_features",     # How many MALDI features match peptides from same protein
            "log_protein_n_features",
        ]
        if self.use_im:
            names.extend([
                "ccs_error",              # predicted CCS - MALDI observed CCS
                "ccs_error_abs",
                "ccs_error_pct",
                "ccs_rank",              # Rank by CCS agreement within group
                "inv_k0_error_abs",       # |MALDI 1/K0 - LC-MS/MS 1/K0|
            ])
        return names

    def compute_features(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Compute rescoring features for a MALDI-MSI matched candidate table.

        Parameters
        ----------
        df : pd.DataFrame
            Matched candidates with columns:
            - maldi_mz: MALDI observed m/z
            - maldi_ccs: MALDI observed CCS (optional if use_im=False)
            - maldi_inv_k0: MALDI observed 1/K0 (optional)
            - maldi_intensity: MALDI feature intensity
            - peptide: peptide sequence
            - mz_ppm_error: ppm error from matching
            - probability: LC-MS/MS search probability
            - expectation: LC-MS/MS search expectation value
            - spectral_count (or Spectral Count): LC-MS/MS spectral count
            - peptide_length (or Peptide Length): peptide length
            - n_candidates: number of candidates for this MALDI feature
            - predicted_ccs: IM2Deep predicted CCS (optional)
            - lcms_inv_k0: LC-MS/MS observed 1/K0 (optional)

        Returns
        -------
        pd.DataFrame
            Input dataframe with added feature columns.
        """
        df = df.copy()

        # --- Mass accuracy features ---
        ppm_col = _find_col(df, ["mz_ppm_error", "ppm"])
        df["ppm_error"] = df[ppm_col]
        df["ppm_error_abs"] = df[ppm_col].abs()

        # Rank within group (1 = best mass match)
        df["ppm_rank"] = df.groupby("maldi_mz")["ppm_error_abs"].rank(method="min")

        # Ratio to best ppm in group
        best_ppm = df.groupby("maldi_mz")["ppm_error_abs"].transform("min")
        df["ppm_best_ratio"] = df["ppm_error_abs"] / best_ppm.clip(lower=0.1)

        # --- Candidate ambiguity ---
        n_cand_col = _find_col(df, ["n_candidates"])
        if n_cand_col is None:
            df["n_candidates"] = df.groupby("maldi_mz")["maldi_mz"].transform("count")
        df["log_n_candidates"] = np.log1p(df["n_candidates"])

        # --- LC-MS/MS search quality ---
        prob_col = _find_col(df, ["probability", "Probability"])
        df["probability"] = df[prob_col] if prob_col else 0.0

        exp_col = _find_col(df, ["expectation", "Expectation"])
        df["log_expectation"] = np.log10(df[exp_col].clip(lower=1e-15)) if exp_col else 0.0

        sc_col = _find_col(df, ["spectral_count", "Spectral Count"])
        df["spectral_count"] = df[sc_col] if sc_col else 0
        df["log_spectral_count"] = np.log1p(df["spectral_count"])

        # --- Peptide properties ---
        pl_col = _find_col(df, ["peptide_length", "Peptide Length"])
        df["peptide_length"] = df[pl_col] if pl_col else df["peptide"].str.len()

        mod_col = _find_col(df, ["assigned_modifications", "Assigned Modifications"])
        if mod_col:
            df["has_modifications"] = df[mod_col].notna().astype(int)
        else:
            df["has_modifications"] = df["peptidoform"].str.contains(r"\[", regex=True).astype(int)

        # Missed cleavages: count internal K/R not followed by P
        df["n_missed_cleavages"] = df["peptide"].apply(
            lambda s: sum(1 for i, c in enumerate(s[:-1]) if c in "KR" and s[i+1] != "P")
        )

        # --- MALDI signal ---
        int_col = _find_col(df, ["maldi_intensity", "Intensity"])
        df["log_maldi_intensity"] = np.log1p(df[int_col]) if int_col else 0.0

        # --- Protein-level consistency ---
        prot_col = _find_col(df, ["protein", "Protein"])
        if prot_col:
            # Count distinct MALDI features that match peptides from the same protein
            protein_feature_count = df.groupby(prot_col)["maldi_mz"].nunique()
            df["protein_n_features"] = df[prot_col].map(protein_feature_count).fillna(0)
        else:
            df["protein_n_features"] = 0
        df["log_protein_n_features"] = np.log1p(df["protein_n_features"])

        # --- Ion mobility features ---
        if self.use_im:
            ccs_obs_col = _find_col(df, ["maldi_ccs", "CCS "])
            ccs_pred_col = _find_col(df, ["predicted_ccs", "calibrated_ccs_spline"])

            if ccs_obs_col and ccs_pred_col:
                obs = df[ccs_obs_col].values.astype(float)
                pred = df[ccs_pred_col].values.astype(float)
                df["ccs_error"] = pred - obs
                df["ccs_error_abs"] = np.abs(df["ccs_error"])
                df["ccs_error_pct"] = df["ccs_error_abs"] / np.clip(obs, 1, None) * 100
                df["ccs_rank"] = df.groupby("maldi_mz")["ccs_error_abs"].rank(method="min")
            else:
                logger.warning("CCS columns not found, setting CCS features to 0")
                for c in ["ccs_error", "ccs_error_abs", "ccs_error_pct", "ccs_rank"]:
                    df[c] = 0.0

            im_maldi_col = _find_col(df, ["maldi_inv_k0", "1/K0"])
            im_lcms_col = _find_col(df, ["lcms_inv_k0", "Apex Ion Mobility"])
            if im_maldi_col and im_lcms_col:
                df["inv_k0_error_abs"] = np.abs(
                    df[im_maldi_col].astype(float) - df[im_lcms_col].astype(float)
                )
            else:
                df["inv_k0_error_abs"] = 0.0

        logger.info(f"Computed {len(self.feature_names)} features for {len(df)} candidates")
        return df

    def build_training_data(
        self,
        df: pd.DataFrame,
        label_by: str = "ppm",
    ) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
        """
        Build labeled training data for rescoring from ambiguous matches.

        Uses only ambiguous groups (2+ candidates per MALDI feature). Within
        each group, the candidate with the best value for `label_by` is labeled
        positive (1), the rest negative (0).

        Features derived from the label criterion are excluded from the
        returned feature set to avoid circularity.

        Parameters
        ----------
        df : pd.DataFrame
            Feature-annotated candidates (output of compute_features).
        label_by : str
            Which criterion defines "true" for labeling:
            - "ppm": best mass accuracy is true (tests if CCS/protein/quality
              features can recover mass accuracy ranking)
            - "ccs": best CCS match is true (tests if mass/protein/quality
              features can recover CCS ranking)

        Returns
        -------
        features : pd.DataFrame
            Feature matrix (only ambiguous candidates).
        labels : np.ndarray
            Binary labels.
        feature_names_used : list[str]
            Feature names actually used (after excluding label-derived features).
        """
        amb = df[df["n_candidates"] > 1].copy()

        if len(amb) == 0:
            raise ValueError("No ambiguous candidates to build training data from.")

        # Assign labels based on rank within each group
        if label_by == "ppm":
            amb["label"] = (amb["ppm_rank"] == 1).astype(int)
            exclude = {"ppm_rank", "ppm_best_ratio", "ppm_error", "ppm_error_abs"}
        elif label_by == "ccs":
            if "ccs_rank" not in amb.columns:
                raise ValueError("CCS features not available for 'ccs' labeling.")
            amb["label"] = (amb["ccs_rank"] == 1).astype(int)
            exclude = {"ccs_rank", "ccs_error", "ccs_error_abs", "ccs_error_pct"}
        else:
            raise ValueError(f"Unknown label_by: {label_by}")

        use_features = [f for f in self.feature_names if f not in exclude]
        features = amb[use_features].copy()
        labels = amb["label"].values

        n_pos = labels.sum()
        n_neg = len(labels) - n_pos
        logger.info(
            f"Training data from {amb['maldi_mz'].nunique()} ambiguous features: "
            f"{n_pos} positive, {n_neg} negative. "
            f"Labels by: {label_by}. "
            f"Features: {len(use_features)} (excluded {exclude})."
        )

        return features, labels, use_features


def _find_col(df: pd.DataFrame, candidates: list[str]) -> Optional[str]:
    """Find the first column name that exists in the DataFrame."""
    for c in candidates:
        if c in df.columns:
            return c
    return None


# ---------------------------------------------------------------------------
# FeatureGeneratorBase implementation for MS2Rescore integration
# ---------------------------------------------------------------------------


class MALDISpatialFeatureGenerator(FeatureGeneratorBase):
    """
    MALDI-MSI feature generator for MS2Rescore.

    Computes per-candidate features from PSMList metadata, including mass
    accuracy, protein consistency, spatial quality, and peptide properties.
    Integrates into the MS2Rescore pipeline via the FeatureGeneratorBase
    interface.

    Parameters
    ----------
    spatial_features : pd.DataFrame, optional
        Pre-computed spatial features (feature_mz, spatial_autocorrelation, etc.).
    use_im : bool
        Whether to include ion mobility features.
    """

    required_ms_data = set()  # No spectrum files needed

    def __init__(
        self,
        *args,
        spatial_features: Optional[pd.DataFrame] = None,
        ion_images: Optional[np.ndarray] = None,
        ion_image_mzs: Optional[np.ndarray] = None,
        use_im: bool = False,
        **kwargs,
    ):
        """
        Parameters
        ----------
        spatial_features
            Pre-computed spatial features DataFrame.
        ion_images
            3D array (n_features, height, width) of ion images. Required for
            co-localization feature.
        ion_image_mzs
            Array of m/z values corresponding to ion_images. Required for
            co-localization feature.
        use_im
            Whether to include ion mobility features.
        """
        super().__init__(*args, **kwargs)
        self.spatial_features = spatial_features
        self.ion_images = ion_images
        self.ion_image_mzs = ion_image_mzs
        self.use_im = use_im
        self.lcms_envelopes = kwargs.pop("lcms_envelopes", None)
        self.maldi_envelopes = kwargs.pop("maldi_envelopes", None)

        # Build m/z → image index lookup
        self._mz_to_img_idx = {}
        if ion_images is not None and ion_image_mzs is not None:
            self._mz_to_img_idx = {mz: i for i, mz in enumerate(ion_image_mzs)}

    @property
    def feature_names(self) -> list[str]:
        names = [
            # Mass accuracy (multiple views)
            "ppm_error_abs",
            "ppm_rank",
            "ppm_best_ratio",
            # Candidate ambiguity
            "n_candidates",
            "log_n_candidates",
            # Protein consistency (multiple views)
            "protein_n_features",
            "log_protein_n_features",
            "protein_coverage",
            "protein_rank",
            "protein_best_ratio",
            # Peptide properties
            "peptide_length",
            "n_missed_cleavages",
            "has_modifications",
            # MALDI signal
            "log_maldi_intensity",
            # LC-MS/MS confidence + quantification (no Percolator scores — they leak T/D info)
            "lcms_xcorr",
            "lcms_spectral_count",
            "lcms_log_spectral_count",
            # Sequence-specific theoretical isotope comparison (MALDI obs vs theoretical)
            "theo_isotope_cosine",
            "theo_isotope_chi2",
            "theo_isotope_kl",
            "theo_has_sulfur",
            # Averagine deviation
            "averagine_deviation",
            "averagine_deviation_sulfur",
            # Individual isotope ratio features vs theoretical
            "theo_m1_ratio_diff",
            "theo_m2_ratio_diff",
            "theo_m1_ratio_diff_lcms",
            "theo_m2_ratio_diff_lcms",
            # MALDI-specific ionization features
            "n_arginine",
            "n_basic_residues",
            "n_phenylalanine",
            "n_aromatic",
            "gravy_score",
            "charge_proxy",
        ]
        if self.spatial_features is not None:
            names.extend([
                "spatial_autocorrelation",
                "fraction_detected",
                "intensity_cv",
                "log_mean_intensity",
                "spatial_entropy",
            ])
        if self.ion_images is not None:
            names.extend([
                "protein_colocalization",
                "protein_colocalization_max",
                "protein_colocalization_median",
                "protein_colocalization_n_partners",
            ])
        if self.lcms_envelopes is not None and self.maldi_envelopes is not None:
            names.extend([
                "isotope_envelope_cosine",
                "isotope_envelope_pearson",
                "isotope_envelope_mse",
                "isotope_m1_ratio_diff",
                "isotope_m2_ratio_diff",
                "isotope_n_matched",
            ])
        if self.use_im:
            names.extend([
                "ccs_error_abs",
                "ccs_error_pct",
                "ccs_rank",
                "inv_k0_error_abs",
            ])
        return names

    def add_features(self, psm_list: PSMList) -> None:
        """
        Add MALDI rescoring features to each PSM in the list.

        Extracts match metadata from PSM.metadata, computes group-level
        features (ranks, protein consistency), and writes features into
        PSM.rescoring_features.
        """
        logger.info("Adding MALDI spatial rescoring features")

        # Extract metadata into a working DataFrame
        records = []
        for i, psm in enumerate(psm_list):
            meta = psm.metadata or {}
            records.append({
                "idx": i,
                "feature_mz": meta.get("feature_mz", psm.precursor_mz),
                "ppm_error_abs": meta.get("ppm_error_abs", 0.0),
                "peptide": psm.peptidoform.sequence,
                "protein": psm.protein_list[0] if psm.protein_list else "",
                "is_decoy": psm.is_decoy,
                "xcorr": meta.get("xcorr", 0.0),
                "maldi_intensity": meta.get("feature_intensity", meta.get("mean_intensity", 0.0)),
                "spectral_count": meta.get("spectral_count", meta.get("Spectral Count", 0)),
                "n_pixels_detected": meta.get("n_pixels_detected", 0),
                "fraction_detected": meta.get("fraction_detected", 0.0),
                "spatial_autocorrelation": meta.get("spatial_autocorrelation", 0.0),
                "intensity_cv": meta.get("intensity_cv", 0.0),
            })

        df = pd.DataFrame(records)

        # Ensure numeric types
        for col in ["feature_mz", "ppm_error_abs", "xcorr", "maldi_intensity",
                     "spectral_count", "n_pixels_detected", "fraction_detected",
                     "spatial_autocorrelation", "intensity_cv"]:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

        # --- Mass accuracy features ---
        df["ppm_rank"] = df.groupby("feature_mz")["ppm_error_abs"].rank(method="min")
        best_ppm = df.groupby("feature_mz")["ppm_error_abs"].transform("min")
        df["ppm_best_ratio"] = df["ppm_error_abs"] / best_ppm.clip(lower=0.01)

        # --- Candidate ambiguity ---
        df["n_candidates"] = df.groupby("feature_mz")["feature_mz"].transform("count")
        df["log_n_candidates"] = np.log1p(df["n_candidates"])

        # --- Protein consistency (multiple views) ---
        protein_feature_count = df.groupby("protein")["feature_mz"].nunique()
        protein_total_peptides = df.groupby("protein")["peptide"].nunique()
        df["protein_n_features"] = df["protein"].map(protein_feature_count).fillna(0)
        df["log_protein_n_features"] = np.log1p(df["protein_n_features"])
        # Coverage: fraction of protein's peptides that match MALDI features
        total_peps = df["protein"].map(protein_total_peptides).fillna(1).clip(lower=1)
        df["protein_coverage"] = df["protein_n_features"] / total_peps
        # Rank: within each feature group, rank proteins by n_features (best = 1)
        df["protein_rank"] = df.groupby("feature_mz")["protein_n_features"].rank(
            method="min", ascending=False
        )
        # Ratio: this protein's count / best protein's count in the group
        best_prot = df.groupby("feature_mz")["protein_n_features"].transform("max")
        df["protein_best_ratio"] = df["protein_n_features"] / best_prot.clip(lower=1)

        # --- Peptide properties ---
        df["peptide_length"] = df["peptide"].str.len()
        df["n_missed_cleavages"] = df["peptide"].apply(
            lambda s: sum(1 for i, c in enumerate(s[:-1]) if c in "KR" and s[i+1] != "P")
            if len(s) > 1 else 0
        )
        df["has_modifications"] = 0  # Can be updated if modification info is in metadata

        # --- MALDI signal ---
        df["log_maldi_intensity"] = np.log1p(df["maldi_intensity"])

        # --- LC-MS/MS confidence + quantification ---
        # Note: Percolator score and PEP are excluded — they leak target-decoy info
        df["lcms_xcorr"] = pd.to_numeric(df["xcorr"], errors="coerce").fillna(0)
        sc_vals = [float(psm_list[int(r["idx"])].metadata.get("spectral_count", 0) or 0)
                   for _, r in df.iterrows()]
        df["lcms_spectral_count"] = sc_vals
        df["lcms_log_spectral_count"] = np.log1p(sc_vals)

        # --- Sequence-specific theoretical isotope features ---
        df = self._compute_theoretical_isotope_features(df)

        # --- MALDI-specific ionization features ---
        df = self._compute_maldi_ionization_features(df)

        # --- Spatial features (from pre-computed data or metadata) ---
        if self.spatial_features is not None:
            spatial_map = self.spatial_features.set_index("feature_mz")
            for col in ["spatial_autocorrelation", "fraction_detected", "intensity_cv"]:
                if col in spatial_map.columns:
                    mapped = df["feature_mz"].map(spatial_map[col])
                    df[col] = mapped.fillna(df[col])
            # Additional spatial views
            if "mean_intensity" in spatial_map.columns:
                df["log_mean_intensity"] = np.log1p(
                    df["feature_mz"].map(spatial_map["mean_intensity"]).fillna(0)
                )
            else:
                df["log_mean_intensity"] = df["log_maldi_intensity"]
            # Spatial entropy: -sum(p * log(p)) of pixel intensities
            # Approximated from fraction_detected and CV
            frac = df["fraction_detected"].clip(lower=0.001, upper=0.999)
            df["spatial_entropy"] = -(frac * np.log(frac) + (1 - frac) * np.log(1 - frac))

        # --- Protein co-localization (candidate-specific, multiple views) ---
        if self.ion_images is not None:
            coloc = self._compute_colocalization_detailed(df)
            df["protein_colocalization"] = coloc["mean"]
            df["protein_colocalization_max"] = coloc["max"]
            df["protein_colocalization_median"] = coloc["median"]
            df["protein_colocalization_n_partners"] = coloc["n_partners"]

        # --- Isotope envelope similarity (candidate-specific, multiple views) ---
        if self.lcms_envelopes is not None and self.maldi_envelopes is not None:
            iso = self._compute_envelope_similarity_detailed(df)
            df["isotope_envelope_cosine"] = iso["cosine"]
            df["isotope_envelope_pearson"] = iso["pearson"]
            df["isotope_envelope_mse"] = iso["mse"]
            df["isotope_m1_ratio_diff"] = iso["m1_ratio_diff"]
            df["isotope_m2_ratio_diff"] = iso["m2_ratio_diff"]
            df["isotope_n_matched"] = iso["n_matched"]

        # --- Write features back to PSMs ---
        for _, row in df.iterrows():
            psm = psm_list[int(row["idx"])]
            if psm.rescoring_features is None:
                psm.rescoring_features = {}
            for feat in self.feature_names:
                psm.rescoring_features[feat] = float(row.get(feat, 0.0))

        logger.info(f"Added {len(self.feature_names)} MALDI features to {len(psm_list)} PSMs")

    def _compute_colocalization_detailed(self, df: pd.DataFrame) -> dict[str, pd.Series]:
        """
        Compute multiple views of protein co-localization per candidate.

        Returns dict with keys: mean, max, median, n_partners.
        """
        from scipy.stats import pearsonr

        protein_to_mzs = df.groupby("protein")["feature_mz"].apply(
            lambda x: list(x.unique())
        ).to_dict()

        mean_scores = np.zeros(len(df))
        max_scores = np.zeros(len(df))
        median_scores = np.zeros(len(df))
        n_partners = np.zeros(len(df))

        for i, (_, row) in enumerate(df.iterrows()):
            protein = row["protein"]
            this_mz = row["feature_mz"]

            other_mzs = [
                mz for mz in protein_to_mzs.get(protein, [])
                if mz != this_mz and mz in self._mz_to_img_idx
            ]

            if not other_mzs or this_mz not in self._mz_to_img_idx:
                continue

            this_img = self.ion_images[self._mz_to_img_idx[this_mz]].flatten().astype(float)
            if this_img.std() == 0:
                continue

            correlations = []
            for other_mz in other_mzs:
                other_img = self.ion_images[self._mz_to_img_idx[other_mz]].flatten().astype(float)
                if other_img.std() > 0:
                    r, _ = pearsonr(this_img, other_img)
                    correlations.append(r)

            if correlations:
                mean_scores[i] = float(np.mean(correlations))
                max_scores[i] = float(np.max(correlations))
                median_scores[i] = float(np.median(correlations))
                n_partners[i] = len(correlations)

        n_scored = np.count_nonzero(mean_scores)
        logger.info(
            f"Protein co-localization: {n_scored}/{len(df)} candidates scored "
            f"(mean={mean_scores[mean_scores != 0].mean():.3f})"
            if n_scored > 0 else
            "Protein co-localization: no candidates had same-protein features"
        )

        return {
            "mean": pd.Series(mean_scores, index=df.index),
            "max": pd.Series(max_scores, index=df.index),
            "median": pd.Series(median_scores, index=df.index),
            "n_partners": pd.Series(n_partners, index=df.index),
        }

    def _compute_envelope_similarity_detailed(self, df: pd.DataFrame) -> dict[str, pd.Series]:
        """
        Compute multiple views of MALDI vs LC-MS/MS isotope envelope similarity.

        Returns dict with keys: cosine, pearson, mse, m1_ratio_diff, m2_ratio_diff, n_matched.
        Candidates without envelopes get median fill values.
        """
        from scipy.stats import pearsonr

        n = len(df)
        cosine = np.full(n, np.nan)
        pearson_r = np.full(n, np.nan)
        mse = np.full(n, np.nan)
        m1_diff = np.full(n, np.nan)
        m2_diff = np.full(n, np.nan)
        n_matched = np.zeros(n)

        for i, (_, row) in enumerate(df.iterrows()):
            lcms_env = self.lcms_envelopes.get(row["peptide"])
            maldi_env = self.maldi_envelopes.get(row["feature_mz"])

            if lcms_env is None or maldi_env is None:
                continue

            k = min(len(lcms_env), len(maldi_env))
            a = np.array(lcms_env[:k], dtype=np.float64)
            b = np.array(maldi_env[:k], dtype=np.float64)

            # Count matched (both non-zero)
            matched = np.sum((a > 0) & (b > 0))
            n_matched[i] = matched

            na, nb = np.linalg.norm(a), np.linalg.norm(b)

            # Cosine similarity
            if na > 0 and nb > 0:
                cosine[i] = float(np.dot(a, b) / (na * nb))

            # Pearson correlation
            if a.std() > 0 and b.std() > 0 and k >= 2:
                r, _ = pearsonr(a, b)
                pearson_r[i] = float(r)

            # Mean squared error (lower = better)
            mse[i] = float(np.mean((a - b) ** 2))

            # M1/M0 ratio difference
            if a[0] > 0 and b[0] > 0 and k >= 2:
                m1_diff[i] = abs(a[1] / a[0] - b[1] / b[0])
            if a[0] > 0 and b[0] > 0 and k >= 3:
                m2_diff[i] = abs(a[2] / a[0] - b[2] / b[0])

        # Fill NaN with median (fair to candidates without envelopes)
        result = {}
        for name, arr in [("cosine", cosine), ("pearson", pearson_r), ("mse", mse),
                          ("m1_ratio_diff", m1_diff), ("m2_ratio_diff", m2_diff)]:
            valid = arr[~np.isnan(arr)]
            fill = float(np.median(valid)) if len(valid) > 0 else 0.0
            result[name] = pd.Series(np.where(np.isnan(arr), fill, arr), index=df.index)

        result["n_matched"] = pd.Series(n_matched, index=df.index)

        n_scored = np.sum(~np.isnan(cosine))
        logger.info(f"Isotope envelope: {n_scored}/{n} candidates scored")

        return result

    def _compute_theoretical_isotope_features(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Compute sequence-specific theoretical isotope features.

        For each candidate peptide:
        - Compute exact elemental composition → theoretical isotope distribution
        - Compare MALDI observed envelope against theoretical (cosine, chi2, KL)
        - Compare LC-MS/MS observed against theoretical (ratio diffs)
        - Compute averagine deviation
        """
        from pyteomics import mass as pyteomics_mass
        from math import factorial, exp

        n = len(df)
        theo_cosine = np.zeros(n)
        theo_chi2 = np.zeros(n)
        theo_kl = np.zeros(n)
        theo_has_sulfur = np.zeros(n)
        avg_dev = np.zeros(n)
        avg_dev_sulfur = np.zeros(n)
        theo_m1_diff = np.zeros(n)
        theo_m2_diff = np.zeros(n)
        theo_m1_diff_lcms = np.zeros(n)
        theo_m2_diff_lcms = np.zeros(n)

        for i, (_, row) in enumerate(df.iterrows()):
            peptide = row["peptide"]

            # Get elemental composition
            try:
                comp = pyteomics_mass.Composition(sequence=peptide)
                nc = comp.get("C", 0)
                nh = comp.get("H", 0)
                nn = comp.get("N", 0)
                no = comp.get("O", 0)
                ns = comp.get("S", 0)
            except Exception:
                continue

            theo_has_sulfur[i] = 1.0 if ns > 0 else 0.0

            # Sequence-specific theoretical distribution (Poisson)
            lam = nc * 0.01109 + nh * 0.000115 + nn * 0.00364 + no * 0.00205 + ns * 0.04493
            theo = np.array([exp(-lam) * lam**k / factorial(k) for k in range(3)])
            theo_total = theo.sum()
            if theo_total > 0:
                theo = theo / theo_total

            # Averagine theoretical for the same mass
            pep_mass = pyteomics_mass.calculate_mass(composition=comp)
            nc_avg = int(pep_mass * 0.0444)
            nh_avg = int(pep_mass * 0.0698)
            nn_avg = int(pep_mass * 0.0123)
            no_avg = int(pep_mass * 0.0133)
            # averagine has negligible S
            lam_avg = nc_avg * 0.01109 + nh_avg * 0.000115 + nn_avg * 0.00364 + no_avg * 0.00205
            avg_theo = np.array([exp(-lam_avg) * lam_avg**k / factorial(k) for k in range(3)])
            avg_total = avg_theo.sum()
            if avg_total > 0:
                avg_theo = avg_theo / avg_total

            # Averagine deviation
            na_t, na_a = np.linalg.norm(theo), np.linalg.norm(avg_theo)
            if na_t > 0 and na_a > 0:
                avg_dev[i] = 1.0 - float(np.dot(theo, avg_theo) / (na_t * na_a))
            # Sulfur-specific: difference in M+2/M+0
            if theo[0] > 0 and avg_theo[0] > 0:
                avg_dev_sulfur[i] = abs(theo[2] / theo[0] - avg_theo[2] / avg_theo[0])

            # Compare MALDI observed vs sequence-specific theoretical
            maldi_env = self.maldi_envelopes.get(row["feature_mz"]) if self.maldi_envelopes else None
            if maldi_env is not None and len(maldi_env) >= 3:
                obs = np.array(maldi_env[:3], dtype=np.float64)
                na_o = np.linalg.norm(obs)
                if na_o > 0 and na_t > 0:
                    theo_cosine[i] = float(np.dot(obs, theo) / (na_o * na_t))
                # Chi-squared
                expected = theo * obs.sum()
                mask = expected > 0
                if mask.any():
                    theo_chi2[i] = float(np.sum((obs[mask] - expected[mask]) ** 2 / expected[mask]))
                # KL divergence (obs || theo)
                obs_norm = obs / obs.sum() if obs.sum() > 0 else obs
                theo_safe = np.clip(theo, 1e-10, None)
                obs_safe = np.clip(obs_norm, 1e-10, None)
                theo_kl[i] = float(np.sum(obs_safe * np.log(obs_safe / theo_safe)))
                # Ratio diffs
                if obs[0] > 0 and theo[0] > 0:
                    theo_m1_diff[i] = abs(obs[1] / obs[0] - theo[1] / theo[0])
                    if len(obs) >= 3 and len(theo) >= 3:
                        theo_m2_diff[i] = abs(obs[2] / obs[0] - theo[2] / theo[0])

            # Compare LC-MS/MS observed vs theoretical
            lcms_env = self.lcms_envelopes.get(peptide) if self.lcms_envelopes else None
            if lcms_env is not None and len(lcms_env) >= 3 and theo[0] > 0:
                lobs = np.array(lcms_env[:3], dtype=np.float64)
                if lobs[0] > 0:
                    theo_m1_diff_lcms[i] = abs(lobs[1] / lobs[0] - theo[1] / theo[0])
                    theo_m2_diff_lcms[i] = abs(lobs[2] / lobs[0] - theo[2] / theo[0])

        df["theo_isotope_cosine"] = theo_cosine
        df["theo_isotope_chi2"] = theo_chi2
        df["theo_isotope_kl"] = theo_kl
        df["theo_has_sulfur"] = theo_has_sulfur
        df["averagine_deviation"] = avg_dev
        df["averagine_deviation_sulfur"] = avg_dev_sulfur
        df["theo_m1_ratio_diff"] = theo_m1_diff
        df["theo_m2_ratio_diff"] = theo_m2_diff
        df["theo_m1_ratio_diff_lcms"] = theo_m1_diff_lcms
        df["theo_m2_ratio_diff_lcms"] = theo_m2_diff_lcms

        logger.info(f"Theoretical isotope features: {(theo_cosine > 0).sum()}/{n} candidates scored")
        return df

    def _compute_maldi_ionization_features(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Compute MALDI-specific ionization sequence features.

        Based on Giese et al. (BMC Bioinformatics 2008) and known MALDI
        ionization properties.
        """
        # Kyte-Doolittle hydropathy scale
        kd = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
              "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
              "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
              "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}

        n_arg = []
        n_basic = []
        n_phe = []
        n_aromatic = []
        gravy = []
        charge_proxy = []

        for _, row in df.iterrows():
            seq = row["peptide"]
            n_arg.append(seq.count("R"))
            n_basic.append(seq.count("R") + seq.count("K") + seq.count("H"))
            n_phe.append(seq.count("F"))
            n_aromatic.append(seq.count("F") + seq.count("W") + seq.count("Y"))
            # GRAVY: average hydropathy
            if len(seq) > 0:
                gravy.append(sum(kd.get(aa, 0) for aa in seq) / len(seq))
            else:
                gravy.append(0)
            # Charge proxy at ~pH 2 (MALDI conditions):
            # R, K, H all protonated (+1 each), D, E protonated (neutral)
            # N-term +1
            charge_proxy.append(seq.count("R") + seq.count("K") + seq.count("H") + 1)

        df["n_arginine"] = n_arg
        df["n_basic_residues"] = n_basic
        df["n_phenylalanine"] = n_phe
        df["n_aromatic"] = n_aromatic
        df["gravy_score"] = gravy
        df["charge_proxy"] = charge_proxy

        logger.info(f"MALDI ionization features computed for {len(df)} candidates")
        return df
