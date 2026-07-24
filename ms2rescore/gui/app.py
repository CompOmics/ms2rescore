"""Graphical user interface for MS²Rescore using CustomTkinter."""

import importlib.resources
import logging
import multiprocessing
import os
import platform
import sys
import webbrowser
from pathlib import Path
from typing import Any

import customtkinter as ctk
from joblib import parallel_backend
from ms2pip.constants import MODELS as ms2pip_models
from PIL import Image
from psm_utils.io import FILETYPES
from rich.console import Console
from rich.logging import RichHandler

import ms2rescore.package_data.img as pkg_data_img
from ms2rescore import __version__ as ms2rescore_version
from ms2rescore._version import check_for_update
from ms2rescore.config_parser import parse_configurations
from ms2rescore.core import rescore
from ms2rescore.exceptions import MS2RescoreConfigurationError
from ms2rescore.gui import widgets
from ms2rescore.gui.function2ctk import Function2CTk, PopupWindow

_IMG_DIR = Path(str(importlib.resources.files(pkg_data_img)))

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

UPDATE_CHECK_TIMEOUT = 0.8  # seconds; keep small to avoid noticeable blocking
UPDATE_CHECK_ENABLED = os.environ.get("MS2RESCORE_SKIP_UPDATE_CHECK", "0") not in (
    "1",
    "true",
    "True",
)

try:
    import matplotlib.pyplot as plt

    plt.set_loglevel("warning")
except ImportError:
    pass

# TODO Does this disable multiprocessing everywhere?
parallel_backend("threading")

CONFIG_WIDTH = 600
CITATIONS = [
    (
        "MS²Rescore: Declercq et al. JPR (2024)",
        "https://doi.org/10.1021/acs.jproteome.3c00785",
    ),
    (
        "MS²PIP: Declercq et al. NAR (2023)",
        "https://doi.org/10.1093/nar/gkad335",
    ),
    (
        "IM2Deep: Declercq, Devreese et al. JPR (2025)",
        "https://doi.org/10.1021/acs.jproteome.4c00609",
    ),
    (
        "DeepLC: Bouwmeester et al. Nat Methods (2021)",
        "https://doi.org/10.1038/s41592-021-01301-5",
    ),
    (
        "Semi-supervised learning: Käll et al. Nat Methods (2007)",
        "https://doi.org/10.1038/nmeth1113",
    ),
]
LINKS = [
    (
        "User guide",
        "https://ms2rescore.readthedocs.io/en/stable/userguide/configuration/",
        "docs",
    ),
    (
        "Discussion forum",
        "https://github.com/compomics/ms2rescore/discussions/categories/q-a",
        "comments",
    ),
    (
        "CompOmics/ms2rescore",
        "https://github.com/compomics/ms2rescore",
        "github",
    ),
]


class SideBar(ctk.CTkFrame):
    def __init__(self, *args, **kwargs):
        """Create the UI sidebar"""
        super().__init__(*args, **kwargs)

        # Configure layout (three rows, one column)
        self.grid_rowconfigure(0, weight=1)
        row_count = 0

        # Top row: logo
        self.logo = ctk.CTkImage(
            light_image=Image.open(os.path.join(str(_IMG_DIR), "ms2rescore_logo.png")),
            size=(130, 130),
        )
        self.logo_label = ctk.CTkLabel(self, text="", image=self.logo)
        self.logo_label.grid(row=row_count, column=0, padx=0, pady=(20, 50), sticky="n")
        row_count += 1

        # Links
        self.links = LinkFrame(self, LINKS)
        self.links.configure(fg_color="transparent")
        self.links.grid(row=row_count, column=0, padx=20, pady=(0, 10), sticky="nsew")
        row_count += 1

        # Citations
        self.citations = CitationFrame(self, CITATIONS)
        self.citations.configure(fg_color="transparent")
        self.citations.grid(row=row_count, column=0, padx=20, pady=(0, 10), sticky="nsew")
        row_count += 1

        # Bottom row: Appearance and UI scaling
        self.ui_control = widgets.UIControl(self)
        self.ui_control.configure(fg_color="transparent")
        self.ui_control.grid(row=row_count, column=0, padx=20, pady=(0, 10), sticky="nsew")
        row_count += 1

        # Bottom row: version
        self.version_label = ctk.CTkLabel(self, text=f"v{ms2rescore_version}")
        self.version_label.grid(row=row_count, column=0, padx=20, pady=(10, 10))
        row_count += 1


class LinkFrame(ctk.CTkFrame):
    def __init__(self, master, links: list[tuple[str, str, str]], *args, **kwargs):
        super().__init__(master, *args, **kwargs)
        self.heading = ctk.CTkLabel(
            self, text="Useful links", font=ctk.CTkFont(weight="bold"), anchor="w"
        )
        self.heading.grid(row=0, column=0, padx=0, pady=0, sticky="ew")

        for i, (ref, url, icon) in enumerate(links):
            button = ctk.CTkButton(
                self,
                text=ref,
                text_color=("#000000", "#fefdff"),
                hover_color=("#3a7ebf", "#1f538d"),
                fg_color="transparent",
                anchor="w",
                image=ctk.CTkImage(
                    dark_image=Image.open(os.path.join(str(_IMG_DIR), f"{icon}_icon_white.png")),
                    light_image=Image.open(os.path.join(str(_IMG_DIR), f"{icon}_icon_black.png")),
                    size=(20, 20),
                ),
                command=lambda x=url: webbrowser.open_new(x),
            )
            button.grid(row=i + 1, column=0, padx=0, pady=(0, 5), sticky="ew")


class CitationFrame(ctk.CTkFrame):
    def __init__(self, master, citations: list[tuple[str]], *args, **kwargs):
        super().__init__(master, *args, **kwargs)
        self.heading = ctk.CTkLabel(
            self, text="Please cite", font=ctk.CTkFont(weight="bold"), anchor="w"
        )
        self.heading.grid(row=0, column=0, padx=0, pady=0, sticky="ew")

        self.buttons = []
        for i, (ref, url) in enumerate(citations):
            button = ctk.CTkButton(
                self,
                text=ref,
                text_color=("#000000", "#fefdff"),
                hover_color=("#3a7ebf", "#1f538d"),
                fg_color="transparent",
                anchor="w",
                height=8,
                command=lambda x=url: webbrowser.open_new(x),
            )
            button.grid(row=i + 1, column=0, padx=0, pady=0, sticky="ew")
            self.buttons.append(button)


class ConfigFrame(ctk.CTkTabview):
    def __init__(self, *args, **kwargs):
        """MS²Rescore configuration frame."""
        super().__init__(*args, **kwargs)

        self.configure(width=CONFIG_WIDTH)

        for tab in ["Main", "Advanced", "Feature generators", "Rescoring"]:
            self.add(tab)
            self.tab(tab).grid_columnconfigure(0, weight=1)
            self.tab(tab).grid_rowconfigure(0, weight=1)
        self.set("Main")

        self.main_config = MainConfiguration(self.tab("Main"))
        self.main_config.grid(row=0, column=0, padx=5, sticky="nsew")

        self.advanced_config = AdvancedConfiguration(self.tab("Advanced"))
        self.advanced_config.grid(row=0, column=0, padx=5, sticky="nsew")

        self.fgen_config = FeatureGeneratorConfig(self.tab("Feature generators"))
        self.fgen_config.grid(row=0, column=0, padx=5, sticky="nsew")

        self.rescoring_config = RescoringConfiguration(self.tab("Rescoring"))
        self.rescoring_config.grid(row=0, column=0, padx=5, sticky="nsew")

    def get(self):
        """Create MS²Rescore config file"""
        main_config = self.main_config.get()
        advanced_config = self.advanced_config.get()

        config = {"ms2rescore": main_config}
        config["ms2rescore"].update(advanced_config)
        feature_generators, annotation_config = self.fgen_config.get()
        config["ms2rescore"]["feature_generators"] = feature_generators
        config["ms2rescore"].update(annotation_config)
        config["ms2rescore"].update(self.rescoring_config.get())

        args = (config,)  # Comma required to wrap in tuple
        kwargs = {}

        return args, kwargs


class MainConfiguration(ctk.CTkFrame):
    def __init__(self, *args, **kwargs):
        """Main MS²Rescore configuration frame."""
        super().__init__(*args, **kwargs)

        self.configure(fg_color="transparent")
        self.grid_columnconfigure(0, weight=1)
        row_n = 0

        self.psm_file = widgets.LabeledFileSelect(
            self,
            label="Identification file",
            description=(
                "Select the PSM file generated by the search engine. This file should contain "
                "all unfiltered target and decoy identifications. Multiple files can be selected."
            ),
            wraplength=CONFIG_WIDTH - 20,
            file_option="openfiles",
        )
        self.psm_file.grid(row=row_n, column=0, pady=(0, 10), sticky="nsew")
        row_n += 1

        self.psm_file_type = widgets.LabeledOptionMenu(
            self,
            vertical=False,
            label="Identification file type",
            description=(
                "Select the file type of the PSM file. The 'infer' option will attempt to infer "
                "the file type from the file extension. If your file type is not listed, do not "
                "hesitate to open a feature request on the discussion forum."
            ),
            wraplength=CONFIG_WIDTH - 150,
            values=["infer"] + list(FILETYPES.keys()),
        )
        self.psm_file_type.grid(row=row_n, column=0, pady=(0, 10), sticky="nsew")
        row_n += 1

        self.spectrum_path = widgets.LabeledFileSelect(
            self,
            label="Spectrum file or directory",
            description=(
                "Select the MGF, mzML, or Bruker raw file(s) containing the spectra. If the "
                "search engine wrote the file names to the PSM file, select the directory "
                "containing the files. The file path should not contain spaces. "
            ),
            wraplength=CONFIG_WIDTH - 20,
            file_option="file/dir",
        )
        self.spectrum_path.grid(row=row_n, column=0, pady=(0, 10), sticky="nsew")
        row_n += 1

        self.modification_mapping = widgets.TableInput(
            self,
            label="Modification mapping",
            description=(
                "Map search engine modification labels to ProForma labels (PSI-MOD, UniMod, "
                "formula, or mass shift). This is required for correct modification parsing. "
                "If this field is left empty, the search engine labels will be used as is, "
                "which may lead to incorrect feature generation for modified peptides. "
                "Check out the user guide for more information."
            ),
            columns=2,
            header_labels=["Search engine label", "ProForma label"],
            wraplength=CONFIG_WIDTH - 20,  # width of the frame minus padding; hardcoded for now
        )
        self.modification_mapping.grid(row=row_n, column=0, pady=(0, 10), sticky="new")
        row_n += 1

        self.fixed_modifications = widgets.TableInput(
            self,
            label="Fixed modifications (MaxQuant only)",
            description=(
                "Add fixed modifications that are not included in the PSM file by the search "
                "engine. If the search engine writes fixed modifications to the PSM file (as most "
                "do), leave this field empty. However, if you are using MaxQuant, which does not "
                "write fixed modifications to the PSM file, you should add them here."
            ),
            columns=2,
            header_labels=["ProForma label", "Amino acids (comma-separated)"],
            wraplength=CONFIG_WIDTH - 20,  # width of the frame minus padding; hardcoded for now
        )
        self.fixed_modifications.grid(row=row_n, column=0, pady=(0, 10), sticky="nsew")
        row_n += 1

    def get(self) -> dict:
        """Get the configured values as a dictionary."""
        try:
            # there cannot be spaces in the file path
            # TODO: Fix this in widgets.LabeledFileSelect
            psm_files = self.psm_file.get().split(" ")
        except AttributeError:
            raise MS2RescoreConfigurationError("No PSM file provided. Please select a file.")
        return {
            "psm_file": psm_files,
            "psm_file_type": self.psm_file_type.get(),
            "spectrum_path": self.spectrum_path.get(),
            "modification_mapping": self._parse_modification_mapping(
                self.modification_mapping.get()
            ),
            "fixed_modifications": self._parse_fixed_modifications(self.fixed_modifications.get()),
        }

    @staticmethod
    def _parse_modification_mapping(table_output):
        """Parse text input modifications mapping"""
        modification_map = {}
        for mod in table_output:
            if mod[0] and mod[1]:
                modification_map[mod[0].strip()] = mod[1].strip()
        return modification_map or None

    @staticmethod
    def _parse_fixed_modifications(table_output):
        """Parse text input fixed modifications"""
        fixed_modifications = {}
        for mod in table_output:
            if mod[0] and mod[1]:
                amino_acids = [aa.upper() for aa in mod[1].strip().split(",")]
                fixed_modifications[mod[0]] = amino_acids
        return fixed_modifications or None


class AdvancedConfiguration(ctk.CTkFrame):
    def __init__(self, *args, **kwargs):
        """Advanced MS²Rescore configuration frame."""
        super().__init__(*args, **kwargs)

        self.configure(fg_color="transparent")
        self.grid_columnconfigure(0, weight=1)

        self.usi = widgets.LabeledSwitch(
            self,
            label="Rename spectrum IDs to USIs",
            description=(
                "Rename the spectrum identifiers to Universal Spectrum Identifiers "
                "(USIs) for full provenance tracking."
            ),
            wraplength=CONFIG_WIDTH - 180,
        )
        self.usi.grid(row=0, column=0, pady=(0, 10), sticky="nsew")

        self.write_flashlfq = widgets.LabeledSwitch(
            self,
            label="Write FlashLFQ input file",
            description=(
                "Write a file that can be used as input for FlashLFQ. This file only contains "
                "target PSMs that pass the FDR threshold."
            ),
            wraplength=CONFIG_WIDTH - 180,
        )
        self.write_flashlfq.grid(row=1, column=0, pady=(0, 10), sticky="nsew")

        self.generate_report = widgets.LabeledSwitch(
            self,
            label="Generate interactive report",
            description=(
                "Generate an interactive report with quality control charts about the MS²Rescore "
                "run. This report can be viewed in a web browser."
            ),
            wraplength=CONFIG_WIDTH - 180,
            default=True,
        )
        self.generate_report.grid(row=2, column=0, pady=(0, 10), sticky="nsew")

        self.id_decoy_pattern = widgets.LabeledEntry(
            self,
            label="Decoy protein regex pattern",
            description=(
                "A regular expression pattern to identify decoy PSMs by the associated protein "
                "names. Most PSM file types contain a dedicated field indicating decoy PSMs, in "
                "which case this field can be left empty."
            ),
            wraplength=CONFIG_WIDTH - 180,
        )
        self.id_decoy_pattern.grid(row=3, column=0, pady=(0, 10), sticky="nsew")

        self.psm_id_pattern = widgets.LabeledEntry(
            self,
            label="PSM ID regex pattern",
            description=(
                "A regular expression pattern to extract the spectrum ID from the IDs in the PSM "
                "file. In most cases, this field can be left empty. Check the user guide for more "
                "information."
            ),
            wraplength=CONFIG_WIDTH - 180,
        )
        self.psm_id_pattern.grid(row=4, column=0, pady=(0, 10), sticky="nsew")

        self.spectrum_id_pattern = widgets.LabeledEntry(
            self,
            label="Spectrum ID regex pattern",
            description=(
                "Similar to the PSM ID regex pattern, but for the IDs in the spectrum file."
            ),
            wraplength=CONFIG_WIDTH - 180,
        )
        self.spectrum_id_pattern.grid(row=5, column=0, pady=(0, 10), sticky="nsew")

        self.processes = widgets.LabeledOptionMenu(
            self,
            label="Number of parallel processes",
            description=(
                "Choose higher values for faster processing, and lower values for less memory "
                "usage."
            ),
            wraplength=CONFIG_WIDTH - 180,
            # Limit to 16 processes to avoid memory overhead
            values=[str(x) for x in list(range(1, min(16, multiprocessing.cpu_count()) + 1))],
            default_value=str(min(16, multiprocessing.cpu_count())),
        )
        self.processes.grid(row=6, column=0, pady=(0, 10), sticky="nsew")

        self.file_prefix = widgets.LabeledFileSelect(
            self,
            label="Filename for output files",
            file_option="savefile",
            description=(
                "Select the output file prefix. The output files will be saved in the selected "
                "directory with the selected filename plus a suffix denoting the file type. If "
                "left empty, this will be based on the PSM file."
            ),
            wraplength=CONFIG_WIDTH - 20,
        )
        self.file_prefix.grid(row=7, column=0, columnspan=2, sticky="nsew")

        self.config_file = widgets.LabeledFileSelect(
            self,
            label="Configuration file",
            file_option="openfile",
            description=(
                "Select a configuration file. Any options that are left empty in the application "
                "will be filled in with values from this file."
            ),
            wraplength=CONFIG_WIDTH - 20,
        )
        self.config_file.grid(row=8, column=0, columnspan=2, sticky="nsew")

    def get(self) -> dict:
        """Get the configured values as a dictionary."""
        return {
            "rename_to_usi": self.usi.get(),
            "write_flashlfq": self.write_flashlfq.get(),
            "write_report": self.generate_report.get(),
            "id_decoy_pattern": self.id_decoy_pattern.get(),
            "psm_id_pattern": self.psm_id_pattern.get(),
            "spectrum_id_pattern": self.spectrum_id_pattern.get(),
            "processes": int(self.processes.get()),
            "output_path": self.file_prefix.get(),
            "config_file": self.config_file.get(),
        }


class FeatureGeneratorConfig(ctk.CTkFrame):
    def __init__(self, *args, **kwargs):
        """Feature generator configuration frame."""
        super().__init__(*args, **kwargs)

        self.configure(fg_color="transparent")
        self.grid_columnconfigure(0, weight=1)

        self.basic_config = BasicFeatureConfiguration(self)
        self.basic_config.grid(row=0, column=0, pady=(0, 20), sticky="nsew")

        self.ms2pip_config = MS2PIPConfiguration(self)
        self.ms2pip_config.grid(row=1, column=0, pady=(0, 20), sticky="nsew")

        self.deeplc_config = DeepLCConfiguration(self)
        self.deeplc_config.grid(row=2, column=0, pady=(0, 20), sticky="nsew")

        self.im2deep_config = Im2DeepConfiguration(self)
        self.im2deep_config.grid(row=3, column=0, pady=(0, 20), sticky="nsew")

    def get(self) -> dict:
        """Return the configuration as a dictionary."""
        basic_enabled, basic_config = self.basic_config.get()
        ms2pip_enabled, ms2pip_config, annotation_config = self.ms2pip_config.get()
        deeplc_enabled, deeplc_config = self.deeplc_config.get()
        im2deep_enabled, im2deep_config = self.im2deep_config.get()

        config = {}
        if basic_enabled:
            config["basic"] = basic_config
        if ms2pip_enabled:
            config["ms2pip"] = ms2pip_config
        if deeplc_enabled:
            config["deeplc"] = deeplc_config
        if im2deep_enabled:
            config["im2deep"] = im2deep_config

        return config, annotation_config


class BasicFeatureConfiguration(ctk.CTkFrame):
    def __init__(self, *args, **kwargs):
        """Basic configuration frame."""
        super().__init__(*args, **kwargs)

        self.configure(fg_color="transparent")
        self.grid_columnconfigure(0, weight=1)

        self.title = widgets._Heading(self, text="Basic features")
        self.title.grid(row=0, column=0, columnspan=2, pady=(0, 5), sticky="ew")

        self.enabled = widgets.LabeledSwitch(self, label="Enable Basic features", default=True)
        self.enabled.grid(row=1, column=0, pady=(0, 10), sticky="nsew")

    def get(self) -> dict:
        """Return the configuration as a dictionary."""
        enabled = self.enabled.get()
        config = {}
        return enabled, config


class MS2PIPConfiguration(ctk.CTkFrame):
    def __init__(self, *args, **kwargs):
        """MS²PIP configuration frame."""
        super().__init__(*args, **kwargs)

        self.configure(fg_color="transparent")
        self.grid_columnconfigure(0, weight=1)

        self.title = widgets._Heading(self, text="MS²PIP (spectrum intensity prediction)")
        self.title.grid(row=0, column=0, columnspan=2, pady=(0, 5), sticky="ew")

        self.enabled = widgets.LabeledSwitch(self, label="Enable MS²PIP", default=True)
        self.enabled.grid(row=1, column=0, pady=(0, 10), sticky="nsew")

        self.model = widgets.LabeledOptionMenu(
            self, label="MS²PIP model", values=list(ms2pip_models.keys()), default_value="HCD2021"
        )
        self.model.grid(row=2, column=0, pady=(0, 10), sticky="nsew")

        self.fragmentation_model = widgets.LabeledOptionMenu(
            self,
            label="Fragmentation model",
            values=["cidhcd", "etd", "ethcd", "all"],
            default_value="cidhcd",
        )
        self.fragmentation_model.grid(row=3, column=0, pady=(0, 10), sticky="nsew")

        self.ms2_tolerance = widgets.LabeledFloatSpinbox(
            self,
            label="MS² error tolerance in Da",
            step_size=0.01,
            initial_value=0.02,
        )
        self.ms2_tolerance.grid(row=4, column=0, pady=(0, 10), sticky="nsew")

        self.tolerance_mode = widgets.LabeledOptionMenu(
            self,
            label="MS² error tolerance mode",
            values=["Da", "ppm"],
            default_value="Da",
        )
        self.tolerance_mode.grid(row=5, column=0, pady=(0, 10), sticky="nsew")

    def get(self) -> dict:
        """Return the configuration as a dictionary."""
        enabled = self.enabled.get()
        feature_config = {"model": self.model.get()}
        annotation_config = {
            "fragmentation_model": self.fragmentation_model.get(),
            "tolerance_value": self.ms2_tolerance.get(),
            "tolerance_mode": self.tolerance_mode.get(),
        }
        return enabled, feature_config, annotation_config


class DeepLCConfiguration(ctk.CTkFrame):
    def __init__(self, *args, **kwargs):
        """DeepLC configuration frame."""
        super().__init__(*args, **kwargs)

        self.configure(fg_color="transparent")
        self.grid_columnconfigure(0, weight=1)

        self.title = widgets._Heading(self, text="DeepLC (retention time prediction)")
        self.title.grid(row=0, column=0, columnspan=2, pady=(0, 5), sticky="ew")

        self.enabled = widgets.LabeledSwitch(self, label="Enable DeepLC", default=True)
        self.enabled.grid(row=1, column=0, pady=(0, 10), sticky="nsew")

        self.transfer_learning = widgets.LabeledSwitch(self, label="Use transfer learning")
        self.transfer_learning.grid(row=2, column=0, pady=(0, 10), sticky="nsew")

        self.num_epochs = widgets.LabeledFloatSpinbox(
            self,
            label="Number of transfer learning epochs",
            step_size=5,
            initial_value=20,
        )  # way to remove float in spinbox label?
        self.num_epochs.grid(row=3, column=0, pady=(0, 10), sticky="nsew")

        self.calibration_set_size = widgets.LabeledEntry(
            self,
            label="Set calibration set size (fraction or number of PSMs)",
            placeholder_text="0.15",
        )
        self.calibration_set_size.grid(row=4, column=0, pady=(0, 10), sticky="nsew")

    def get(self) -> dict:
        """Return the configuration as a dictionary."""
        if self.calibration_set_size.get() == "":
            calibration_set_size = 0.15
        elif not self.calibration_set_size.get().replace(".", "", 1).isdigit():
            raise MS2RescoreConfigurationError(
                f"Error parsing {self.calibration_set_size.get()}\nMake sure calibration set size "
                f"is a number or percentage"
            )
        elif "." in self.calibration_set_size.get():
            calibration_set_size = float(self.calibration_set_size.get())
        else:
            calibration_set_size = int(self.calibration_set_size.get())

        enabled = self.enabled.get()
        config = {
            "finetune": self.transfer_learning.get(),
            "epochs": int(self.num_epochs.get()),
            "calibration_set_size": calibration_set_size,
        }
        return enabled, config


class Im2DeepConfiguration(ctk.CTkFrame):
    def __init__(self, *args, **kwargs):
        """IM2Deep configuration frame."""
        super().__init__(*args, **kwargs)

        self.configure(fg_color="transparent")
        self.grid_columnconfigure(0, weight=1)

        self.title = widgets._Heading(self, text="IM2Deep (ion mobility prediction)")
        self.title.grid(row=0, column=0, columnspan=2, pady=(0, 5), sticky="ew")

        self.enabled = widgets.LabeledSwitch(self, label="Enable im2deep", default=False)
        self.enabled.grid(row=1, column=0, pady=(0, 10), sticky="nsew")

    def get(self) -> tuple[bool, dict[str, Any]]:
        """Return the configuration as a dictionary."""
        enabled = self.enabled.get()
        config: dict[str, Any] = {}
        return enabled, config


class RescoringConfiguration(ctk.CTkFrame):
    def __init__(self, *args, **kwargs):
        """Ristretto rescoring configuration frame."""
        super().__init__(*args, **kwargs)

        self.configure(fg_color="transparent")
        self.grid_columnconfigure(0, weight=1)
        row_n = 0

        self.title = widgets._Heading(self, text="Ristretto configuration")
        self.title.grid(row=row_n, column=0, columnspan=2, pady=(0, 5), sticky="ew")
        row_n += 1

        self.train_fdr = widgets.LabeledEntry(
            self,
            label="Training FDR threshold",
            description="FDR threshold used for positive selection during ristretto's training.",
            placeholder_text="0.01",
        )
        self.train_fdr.grid(row=row_n, column=0, pady=(0, 10), sticky="nsew")
        row_n += 1

        self.model = widgets.LabeledOptionMenu(
            self,
            label="Rescoring model",
            description=(
                "'svm': iterative Percolator-style SVM (default). 'lda': single-pass Fisher "
                "LDA, faster but less powerful on hard-to-separate datasets."
            ),
            wraplength=CONFIG_WIDTH - 180,
            values=["svm", "lda"],
            default_value="svm",
        )
        self.model.grid(row=row_n, column=0, pady=(0, 10), sticky="nsew")
        row_n += 1

    def get(self) -> dict:
        """Return the configuration as a dictionary."""
        train_fdr_str = self.train_fdr.get()
        if train_fdr_str == "":
            train_fdr = 0.01
        elif not train_fdr_str.replace(".", "", 1).isdigit():
            raise MS2RescoreConfigurationError(
                f"Error parsing {train_fdr_str}\nMake sure the training FDR is a number."
            )
        else:
            train_fdr = float(train_fdr_str)

        return {"rescoring": {"train_fdr": train_fdr, "model": self.model.get()}}


class UpdateDialog(PopupWindow):
    def __init__(self, master, current: str, latest: str, url: str):
        msg = (
            f"A new version of MS²Rescore is available.\n\n"
            f"Current version: {current}\n"
            f"Latest version:  {latest}"
        )
        super().__init__(master, "Update available", msg, action_button=True)

        # Add an "Open release page" button alongside the existing close button.
        def open_release():
            if url:
                webbrowser.open_new_tab(url)
            self.destroy()

        self.action_button.configure(text="Open download page", command=open_release)
        self.close_button.configure(text="Ignore")


def _check_updates_sync(root):
    """Run the update check synchronously with a very short timeout."""
    if not UPDATE_CHECK_ENABLED:
        return
    try:
        update_info = check_for_update(
            timeout_seconds=UPDATE_CHECK_TIMEOUT,
        )
        if update_info.get("update_available", False):
            latest = update_info.get("latest_version") or "unknown"
            url = (
                update_info.get("html_url")
                or "https://github.com/CompOmics/ms2rescore/releases/latest"
            )
            UpdateDialog(root, ms2rescore_version, latest, url)
        # If not ok / offline / rate-limited, we do nothing (no errors to user)
    except Exception:
        # Fully silent to the user; still log for debugging
        logger.exception("Update check failed")


def _setup_logging(log_level: str, log_file: str) -> Console:
    """Set up file logging plus a recording Rich console for the HTML log export."""
    log_level_map = {
        "critical": logging.CRITICAL,
        "error": logging.ERROR,
        "warning": logging.WARNING,
        "info": logging.INFO,
        "debug": logging.DEBUG,
    }
    level = log_level_map.get(log_level, logging.INFO)

    file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    )
    file_handler.setLevel(level)

    console = Console(record=True)
    rich_handler = RichHandler(console=console, rich_tracebacks=True, show_path=False)
    rich_handler.setLevel(level)

    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    root_logger.addHandler(file_handler)
    root_logger.addHandler(rich_handler)

    return console


def function(config):
    """Function to be executed in a separate process."""
    config = config.copy()
    if config["ms2rescore"]["config_file"]:
        config_list = [[config["ms2rescore"]["config_file"], config]]
    else:
        config_list = [config]
    config = parse_configurations(config_list)

    # Set up file logging for GUI
    console = _setup_logging(
        config["ms2rescore"]["log_level"], config["ms2rescore"]["output_path"] + ".log.txt"
    )

    try:
        rescore(configuration=config)
    finally:
        console.save_html(config["ms2rescore"]["output_path"] + ".log.html")

    if config["ms2rescore"]["write_report"]:
        webbrowser.open_new_tab(config["ms2rescore"]["output_path"] + ".report.html")


def app():
    """Start the application."""
    root = Function2CTk(
        sidebar_frame=SideBar,
        config_frame=ConfigFrame,
        function=function,
    )
    root.protocol("WM_DELETE_WINDOW", sys.exit)
    dpi = root.winfo_fpixels("1i")
    root.geometry(f"{int(15 * dpi)}x{int(10 * dpi)}")
    root.minsize(int(13 * dpi), int(9 * dpi))
    root.title("MS²Rescore")
    if platform.system() != "Linux":
        root.wm_iconbitmap(os.path.join(str(_IMG_DIR), "program_icon.ico"))

    # Schedule a synchronous, fast update check shortly after the window paints.
    root.after(150, lambda: _check_updates_sync(root))

    root.mainloop()
