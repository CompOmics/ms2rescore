"""Graphical user interface for MS²Rescore using CustomTkinter."""

import importlib.resources
import logging
import multiprocessing
import os
import platform
import sys
import webbrowser
from pathlib import Path
from typing import Dict, List, Tuple

import customtkinter as ctk
from joblib import parallel_backend
from ms2pip.constants import MODELS as ms2pip_models
from PIL import Image
from psm_utils.io import FILETYPES

import ms2rescore.gui.widgets as widgets
import ms2rescore.package_data.img as pkg_data_img
from ms2rescore import __version__ as ms2rescore_version
from ms2rescore.config_parser import parse_configurations
from ms2rescore.core import rescore
from ms2rescore.exceptions import MS2RescoreConfigurationError
from ms2rescore.gui.function2ctk import Function2CTk

with importlib.resources.path(pkg_data_img, "config_icon.png") as resource:
    _IMG_DIR = Path(resource).parent

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

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
        "ionmob: Teschner et al. Bioinformatics (2023)",
        "https://doi.org/10.1093/bioinformatics/btad486",
    ),
    (
        "Mokapot: Fondrie et al. JPR (2021)",
        "https://doi.org/10.1021/acs.jproteome.0c01010",
    ),
    (
        "Percolator: Käll et al. Nat Methods (2007)",
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
    def __init__(self, master, links: List[Tuple[str, str, str]], *args, **kwargs):
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
    def __init__(self, master, citations: List[Tuple[str]], *args, **kwargs):
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

        for tab in ["Main", "Advanced", "Feature generators", "Rescoring engine"]:
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

        self.rescoring_engine_config = RescoringEngineConfig(self.tab("Rescoring engine"))
        self.rescoring_engine_config.grid(row=0, column=0, padx=5, sticky="nsew")

    def get(self):
        """Create MS²Rescore config file"""
        main_config = self.main_config.get()
        advanced_config = self.advanced_config.get()

        config = {"ms2rescore": main_config}
        config["ms2rescore"].update(advanced_config)
        config["ms2rescore"]["feature_generators"] = self.fgen_config.get()
        config["ms2rescore"]["rescoring_engine"] = self.rescoring_engine_config.get()

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

    def get(self) -> Dict:
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

        self.lower_score = widgets.LabeledSwitch(
            self,
            label="Lower score is better",
            description=(
                "When enabled, a lower search engine score is considered to denote a better PSM."
            ),
            wraplength=CONFIG_WIDTH - 180,
        )
        self.lower_score.grid(row=0, column=0, pady=(0, 10), sticky="nsew")

        self.usi = widgets.LabeledSwitch(
            self,
            label="Rename spectrum IDs to USIs",
            description=(
                "Rename the spectrum identifiers to Universal Spectrum Identifiers "
                "(USIs) for full provenance tracking."
            ),
            wraplength=CONFIG_WIDTH - 180,
        )
        self.usi.grid(row=1, column=0, pady=(0, 10), sticky="nsew")

        self.write_flashlfq = widgets.LabeledSwitch(
            self,
            label="Write FlashLFQ input file",
            description=(
                "Write a file that can be used as input for FlashLFQ. This file only contains "
                "target PSMs that pass the FDR threshold."
            ),
            wraplength=CONFIG_WIDTH - 180,
        )
        self.write_flashlfq.grid(row=2, column=0, pady=(0, 10), sticky="nsew")

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
        self.generate_report.grid(row=3, column=0, pady=(0, 10), sticky="nsew")

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
        self.id_decoy_pattern.grid(row=4, column=0, pady=(0, 10), sticky="nsew")

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
        self.psm_id_pattern.grid(row=5, column=0, pady=(0, 10), sticky="nsew")

        self.spectrum_id_pattern = widgets.LabeledEntry(
            self,
            label="Spectrum ID regex pattern",
            description=(
                "Similar to the PSM ID regex pattern, but for the IDs in the spectrum file."
            ),
            wraplength=CONFIG_WIDTH - 180,
        )
        self.spectrum_id_pattern.grid(row=6, column=0, pady=(0, 10), sticky="nsew")

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
        self.processes.grid(row=7, column=0, pady=(0, 10), sticky="nsew")

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
        self.file_prefix.grid(row=8, column=0, columnspan=2, sticky="nsew")

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
        self.config_file.grid(row=9, column=0, columnspan=2, sticky="nsew")

    def get(self) -> Dict:
        """Get the configured values as a dictionary."""
        return {
            "lower_score_is_better": bool(int(self.lower_score.get())),  # str repr of 0 or 1
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

        self.ionmob_config = IonmobConfiguration(self)
        self.ionmob_config.grid(row=4, column=0, pady=(0, 20), sticky="nsew")

    def get(self) -> Dict:
        """Return the configuration as a dictionary."""
        basic_enabled, basic_config = self.basic_config.get()
        ms2pip_enabled, ms2pip_config = self.ms2pip_config.get()
        deeplc_enabled, deeplc_config = self.deeplc_config.get()
        im2deep_enabled, im2deep_config = self.im2deep_config.get()
        ionmob_enabled, ionmob_config = self.ionmob_config.get()

        config = {}
        if basic_enabled:
            config["basic"] = basic_config
        if ms2pip_enabled:
            config["ms2pip"] = ms2pip_config
        if deeplc_enabled:
            config["deeplc"] = deeplc_config
        if ionmob_enabled:
            config["ionmob"] = ionmob_config
        if im2deep_enabled:
            config["im2deep"] = im2deep_config

        return config


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

    def get(self) -> Dict:
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

        self.ms2_tolerance = widgets.LabeledFloatSpinbox(
            self,
            label="MS² error tolerance in Da",
            step_size=0.01,
            initial_value=0.02,
        )
        self.ms2_tolerance.grid(row=3, column=0, pady=(0, 10), sticky="nsew")

    def get(self) -> Dict:
        """Return the configuration as a dictionary."""
        enabled = self.enabled.get()
        config = {
            "model": self.model.get(),
            "ms2_tolerance": self.ms2_tolerance.get(),
        }
        return enabled, config


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

    def get(self) -> Dict:
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
            "deeplc_retrain": self.transfer_learning.get(),
            "n_epochs": int(self.num_epochs.get()),
            "calibration_set_size": calibration_set_size,
        }
        return enabled, config


class IonmobConfiguration(ctk.CTkFrame):
    def __init__(self, *args, **kwargs):
        """IonMob configuration frame."""
        super().__init__(*args, **kwargs)

        self.configure(fg_color="transparent")
        self.grid_columnconfigure(0, weight=1)

        self.title = widgets._Heading(self, text="Ionmob (ion mobility prediction)")
        self.title.grid(row=0, column=0, columnspan=2, pady=(0, 5), sticky="ew")

        self.enabled = widgets.LabeledSwitch(self, label="Enable Ionmob", default=False)
        self.enabled.grid(row=1, column=0, pady=(0, 10), sticky="nsew")

        self.model = widgets.LabeledEntry(
            self,
            label="Name of built-in model or path to custom model",
            placeholder_text="GRUPredictor",
            default_value="GRUPredictor",
        )
        self.model.grid(row=3, column=0, pady=(0, 10), sticky="nsew")

    def get(self) -> Dict:
        """Return the configuration as a dictionary."""
        enabled = self.enabled.get()
        config = {"ionmob_model": self.model.get()}
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

    def get(self) -> Dict:
        """Return the configuration as a dictionary."""
        enabled = self.enabled.get()
        config = {}
        return enabled, config


class RescoringEngineConfig(ctk.CTkFrame):
    def __init__(self, *args, **kwargs):
        """Rescoring engine configuration frame."""
        super().__init__(*args, **kwargs)

        self.configure(fg_color="transparent")
        self.grid_columnconfigure(0, weight=1)

        self.radio_button = widgets.LabeledRadioButtons(
            self,
            label="Rescoring engine",
            options=["Mokapot", "Percolator"],
            default_value="Mokapot",
        )
        self.radio_button.grid(row=0, column=0, pady=(0, 10), sticky="nsew")

        self.mokapot_config = MokapotRescoringConfiguration(self)
        self.mokapot_config.grid(row=1, column=0, pady=(0, 10), sticky="nsew")

        self.percolator_config = PercolatorRescoringConfiguration(self)
        self.percolator_config.grid(row=2, column=0, pady=(0, 10), sticky="nsew")

    def get(self) -> Dict:
        """Return the configuration as a dictionary."""
        if self.radio_button.get().lower() == "mokapot":
            return {self.radio_button.get().lower(): self.mokapot_config.get()}
        elif self.radio_button.get().lower() == "percolator":
            return {self.radio_button.get().lower(): self.percolator_config.get()}


class MokapotRescoringConfiguration(ctk.CTkFrame):
    def __init__(self, *args, **kwargs):
        """Rescoring engine configuration frame."""
        super().__init__(*args, **kwargs)

        self.configure(fg_color="transparent")
        self.grid_columnconfigure(0, weight=1)
        row_n = 0

        self.title = widgets._Heading(self, text="Mokapot coffeeguration")
        self.title.grid(row=row_n, column=0, columnspan=2, pady=(0, 5), sticky="ew")
        row_n += 1

        self.write_weights = widgets.LabeledSwitch(
            self, label="Write model weights to file", default=True
        )
        self.write_weights.grid(row=row_n, column=0, pady=(0, 10), sticky="nsew")
        row_n += 1

        self.write_txt = widgets.LabeledSwitch(self, label="Write TXT output files", default=True)
        self.write_txt.grid(row=row_n, column=0, pady=(0, 10), sticky="nsew")
        row_n += 1

        self.fasta_file = widgets.LabeledFileSelect(
            self,
            label="Select FASTA file (optional, required for protein inference)",
            file_option="openfile",
        )
        self.fasta_file.grid(row=row_n, column=0, pady=(0, 10), sticky="nsew")
        row_n += 1

        self.protein_kwargs = widgets.TableInput(
            self,
            label="`mokapot.read_fasta` options (see Mokapot documentation)",
            columns=2,
            header_labels=["Parameter", "Value"],
        )
        self.protein_kwargs.grid(row=row_n, column=0, sticky="nsew")
        row_n += 1

    def get(self) -> Dict:
        """Return the configuration as a dictionary."""
        config = {
            "write_weights": self.write_weights.get(),
            "write_txt": self.write_txt.get(),
            "fasta_file": self.fasta_file.get(),
            "protein_kwargs": self._parse_protein_kwargs(self.protein_kwargs.get()),
        }
        return config

    @staticmethod
    def _parse_protein_kwargs(table_output):
        """Parse text input modifications mapping"""
        protein_kwargs = {}
        for mod in table_output:
            if mod[0] and mod[1]:
                protein_kwargs[mod[0].strip()] = mod[1].strip()
        return protein_kwargs


class PercolatorRescoringConfiguration(ctk.CTkFrame):
    def __init__(self, *args, **kwargs):
        """Rescoring engine configuration frame."""
        super().__init__(*args, **kwargs)

        self.configure(fg_color="transparent")
        self.grid_columnconfigure(0, weight=1)

        self.title = widgets._Heading(self, text="Percolator coffeeguration")
        self.title.grid(row=0, column=0, columnspan=2, pady=(0, 5), sticky="ew")

        self.weights_file = widgets.LabeledFileSelect(
            self, label="Pretrained Percolator model weights", file_option="openfile"
        )
        self.weights_file.grid(row=1, column=0, columnspan=2, sticky="nsew")

    def get(self) -> Dict:
        """Return the configuration as a dictionary."""
        config = {"init-weights": self.weights_file.get()}
        return config


def function(config):
    """Function to be executed in a separate process."""
    config = config.copy()
    if config["ms2rescore"]["config_file"]:
        config_list = [[config["ms2rescore"]["config_file"], config]]
    else:
        config_list = [config]
    config = parse_configurations(config_list)
    rescore(configuration=config)
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
    root.geometry(f"{int(15*dpi)}x{int(10*dpi)}")
    root.minsize(int(13 * dpi), int(9 * dpi))
    root.title("MS²Rescore")
    if platform.system() != "Linux":
        root.wm_iconbitmap(os.path.join(str(_IMG_DIR), "program_icon.ico"))

    root.mainloop()
