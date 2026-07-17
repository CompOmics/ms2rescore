import json
import logging
from multiprocessing import cpu_count
from typing import Dict, Optional

import psm_utils.io
from psm_utils import PSMList

from ms2rescore import _ristretto_utils, _utils, exceptions, rescoring
from ms2rescore.feature_generators import FEATURE_GENERATORS
from ms2rescore.parse_psms import parse_psms
from ms2rescore.parse_spectra import MSDataType, add_precursor_values, annotate_spectra
from ms2rescore.report import generate
from ms2rescore.report.data import ReportData

logger = logging.getLogger(__name__)


def rescore(configuration: Dict, psm_list: Optional[PSMList] = None) -> None:
    """
    Run full MS²Rescore workflow with passed configuration.

    Parameters
    ----------
    configuration
        Dictionary containing ms2rescore configuration.
    psm_list
        PSMList object containing PSMs. If None, PSMs will be read from configuration ``psm_file``.

    """
    logger.debug(
        f"Running MS²Rescore with following configuration: {json.dumps(configuration, indent=4)}"
    )
    config = configuration["ms2rescore"]
    output_file_root = config["output_path"].split(".intermediate.")[
        0
    ]  # if no intermediate, takes full name

    # Write full configuration including defaults to file
    with open(output_file_root + ".full-config.json", "w") as f:
        json.dump(configuration, f, indent=4)

    logger.debug("Using %i of %i available CPUs.", int(config["processes"]), int(cpu_count()))

    # Parse PSMs
    psm_list = parse_psms(config, psm_list)

    # Evaluate PSMs' pre-rescoring score with ristretto, for report baselines
    before_result = _ristretto_utils.evaluate_before(psm_list, config)
    usi_by_native_id = None
    if config["rename_to_usi"]:
        # Build the (run, native spectrum_id) -> USI lookup now, while psm_list still has its
        # native spectrum_id (the actual rename, further below, must stay after spectrum-file
        # matching, which needs that native ID). Reused below for the actual rename, and here
        # to remap before_result's spectrum_id into the same representation psm_list will end
        # up with, so the report can match "before" and "after" rows by spectrum. The USI
        # omits the peptidoform, so every candidate PSM of a given spectrum maps to the same
        # USI, making (run, spectrum_id) a safe, collision-free lookup key.
        usi_by_native_id = {
            (psm.run, psm.spectrum_id): f"mzspec:{psm.collection}:{psm.run}:scan:{psm.spectrum_id}"
            for psm in psm_list
        }
        before_result.psms["spectrum_id"] = [
            usi_by_native_id[(run, spectrum_id)]
            for run, spectrum_id in zip(before_result.psms["run"], before_result.psms["spectrum_id"])
        ]
    n_id_before = (
        (before_result.psms["qvalue"] <= config["report_fdr"]) & ~before_result.psms["is_decoy"]
    ).sum()
    logger.info(
        f"Found {n_id_before} identified PSMs at {config['report_fdr']:.2%} FDR before rescoring."
    )

    # Define feature names; get existing feature names from PSM file
    feature_names = dict()
    psm_list_feature_names = {
        feature_name
        for psm_list_features in psm_list["rescoring_features"]
        for feature_name in psm_list_features.keys()
    }
    feature_names["psm_file"] = psm_list_feature_names
    logger.debug(
        f"PSMs already contain the following rescoring features: {psm_list_feature_names}"
    )
    # Check if all features are already present; collect generators to skip
    skip_fgens = set()
    for fgen_name, fgen_config in config["feature_generators"].items():
        fgen_features = FEATURE_GENERATORS[fgen_name]().feature_names
        if set(fgen_features).issubset(psm_list_feature_names):
            logger.debug(
                f"Skipping feature generator {fgen_name} because all features are already "
                "present in the PSM file."
            )
            feature_names[fgen_name] = set(fgen_features)
            feature_names["psm_file"] = psm_list_feature_names - set(fgen_features)
            skip_fgens.add(fgen_name)

    # Add missing precursor info from spectrum file if needed
    required_ms_data = {
        ms_data
        for fgen_name in config["feature_generators"].keys()
        if fgen_name not in skip_fgens
        for ms_data in FEATURE_GENERATORS[fgen_name].required_ms_data
    }
    available_ms_data = add_precursor_values(
        psm_list,
        required_ms_data,
        spectrum_path=config["spectrum_path"],
        spectrum_id_pattern=config["spectrum_id_pattern"],
    )

    # Annotate spectra once so MS2/MS2PIP can both consume preloaded fragment annotations.
    if MSDataType.ms2_spectra in available_ms_data:
        annotate_spectra(
            psm_list,
            fragmentation_model=config.get("fragmentation_model", "cidhcd"),
            ms2_tolerance=config.get("tolerance_value", 0.02),
            ms2_tolerance_mode=config.get("tolerance_mode", "Da"),
        )

    # Add rescoring features
    for fgen_name, fgen_config in config["feature_generators"].items():
        if fgen_name in skip_fgens:
            continue
        # Compile configuration
        conf = config.copy()
        conf.update(fgen_config)
        fgen = FEATURE_GENERATORS[fgen_name](**conf)

        # Check if required MS data is available
        missing_ms_data = fgen.required_ms_data - available_ms_data
        if missing_ms_data:
            logger.warning(
                f"Skipping feature generator {fgen_name} because required MS data is missing: "
                f"{missing_ms_data}. Ensure that the required MS data is present in the input "
                "files or disable the feature generator."
            )
            continue
        try:
            fgen.add_features(psm_list)
        except (
            Exception,
            KeyboardInterrupt,
        ) as e:  # Intentionally broad to save intermediate output before re-raising
            logger.error(
                f"Error while adding features from {fgen_name}: {e}, writing intermediary output..."
            )
            # Write intermediate TSV
            psm_utils.io.write_file(
                psm_list, output_file_root + ".intermediate.tsv", filetype="tsv"
            )
            raise
        logger.debug(f"Adding features from {fgen_name}: {set(fgen.feature_names)}")
        feature_names[fgen_name] = set(fgen.feature_names)

        # Remove overlapping features from psm_file to avoid duplicates
        # (e.g., hyperscore can be in both psm_file and ms2pip)
        overlap = feature_names.get("psm_file", set()) & feature_names[fgen_name]
        if overlap:
            feature_names["psm_file"] = feature_names["psm_file"] - overlap

    # Release the annotated MS2 spectra now that all feature generators have consumed them.
    psm_list["spectrum"] = [None] * len(psm_list)

    # Filter out psms that do not have all added features
    all_feature_names = {f for fgen in feature_names.values() for f in fgen}
    psms_with_features = [
        (set(psm.rescoring_features.keys()) == all_feature_names) for psm in psm_list
    ]

    if psms_with_features.count(False) > 0:
        removed_psms = psm_list[[not psm for psm in psms_with_features]]
        missing_features = {
            feature_name
            for psm in removed_psms
            for feature_name in all_feature_names - set(psm.rescoring_features.keys())
        }
        logger.warning(
            f"Removed {psms_with_features.count(False)} PSMs that were missing one or more "
            f"rescoring feature(s), {missing_features}."
        )
    psm_list = psm_list[psms_with_features]

    if "mumble" in config["psm_generator"]:
        # Remove PSMs where matched_ions_pct drops 25% below the original hit
        psm_list = _utils.filter_mumble_psms(psm_list, threshold=0.50)

        if config["max_psm_rank_output"] == 1:
            logger.warning(
                "Mumble adds multiple candidate PSMs per spectrum, some of which can end up "
                "with an identical rescoring score. With `max_psm_rank_output` set to 1, only "
                "one PSM per spectrum is kept, and which candidate 'wins' a tie is not "
                "deterministic. See the Mumble user guide for details."
            )

    # Write feature names to file
    _write_feature_names(feature_names, output_file_root)

    # Rename PSMs to USIs if requested, reusing the lookup built above
    if config["rename_to_usi"]:
        logging.debug(f"Creating USIs for {len(psm_list)} PSMs")
        psm_list["spectrum_id"] = [usi_by_native_id[(psm.run, psm.spectrum_id)] for psm in psm_list]

    # Rescore PSMs
    logger.info(f"Rescoring {len(psm_list)} PSMs with ristretto...")
    try:
        psm_list, after_result = rescoring.rescore(psm_list, config, output_file_root)
    except (
        Exception,
        KeyboardInterrupt,
    ):  # Intentionally broad to save intermediate output before re-raising
        # Write output
        logger.info(f"Writing intermediary output to {output_file_root}.intermediate.tsv...")
        psm_utils.io.write_file(psm_list, output_file_root + ".intermediate.tsv", filetype="tsv")

        # Reraise exception
        raise

    # Post-rescoring processing. before_result was already trimmed to max_psm_rank_output,
    # same as after_result, so this comparison stays fair regardless of its value.
    n_after = (
        (after_result.psms["qvalue"] <= config["report_fdr"]) & ~after_result.psms["is_decoy"]
    ).sum()
    n_before = n_id_before
    diff = n_after - n_before
    diff_perc = f" ({diff / n_before:.2%})" if n_before > 0 else ""
    logger.info(
        f"Identified {diff:+d}{diff_perc} PSMs at {config['report_fdr']:.2%} FDR after "
        "rescoring, compared to before."
    )

    _ristretto_utils.write_rescoring_tables(after_result, output_file_root)

    # Write output
    logger.info(f"Writing output to {output_file_root}.tsv...")
    psm_utils.io.write_file(psm_list, output_file_root + ".tsv", filetype="tsv")

    if config["write_flashlfq"]:
        logger.info(f"Writing output to {output_file_root}.flashlfq.tsv...")
        psm_utils.io.write_file(
            psm_list,
            output_file_root + ".flashlfq.tsv",
            filetype="flashlfq",
            fdr_threshold=config["report_fdr"],
            only_target=True,
        )

    # Write report
    if config["write_report"]:
        try:
            report_data = ReportData.from_run(
                psm_list,
                feature_names,
                configuration,
                before_result,
                after_result,
                fdr_threshold=config["report_fdr"],
            )
            generate.generate_report(output_file_root, report_data)
        except exceptions.ReportGenerationError as e:
            logger.exception(e)


def _write_feature_names(feature_names, output_file_root):
    """Write feature names to file."""
    with open(output_file_root + ".feature_names.tsv", "w") as f:
        f.write("feature_generator\tfeature_name\n")
        for fgen, fgen_features in feature_names.items():
            for feature in fgen_features:
                f.write(f"{fgen}\t{feature}\n")
