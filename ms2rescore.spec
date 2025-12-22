import importlib.metadata
import os
import re

from PyInstaller.building.build_main import COLLECT, EXE, PYZ, Analysis
from PyInstaller.utils.hooks import collect_all

from ms2rescore import __version__

# Package info
exe_name = "ms2rescore"
script_name = "ms2rescore/gui/__main__.py"
icon = "./img/ms2rescore.ico"
location = os.getcwd()
project = "ms2rescore"
bundle_name = "ms2rescore"
bundle_identifier = f"{bundle_name}.{__version__}"

extra_requirements = {"ionmob"}

# Requirements config
skip_requirements_regex = r"^(?:.*\..*)"


# Collect hidden imports and data for all requirements
requirements = importlib.metadata.requires(project)
requirements = {
    re.match(r"^[\w\-]+", req)[0]  # Remove version specifiers
    for req in requirements
    if "; extra ==" not in req  # Exclude optional dependencies
}
requirements.update([project, "xgboost"])
requirements.update(extra_requirements)

hidden_imports = set()
datas = []
binaries = []
checked = set()
while requirements:
    requirement = requirements.pop()
    if re.match(skip_requirements_regex, requirement):
        continue
    if requirement in ["tomli"]:
        continue
    checked.add(requirement)
    module_version = importlib.metadata.version(re.match(r"^[\w\-]+", requirement)[0])
    try:
        # Use filter to exclude problematic xgboost.testing module
        filter_func = lambda name: not name.startswith("xgboost.testing") if requirement == "xgboost" else True
        datas_, binaries_, hidden_imports_ = collect_all(
            requirement,
            include_py_files=True,
            filter_submodules=filter_func
        )
    except (ImportError, RuntimeError) as e:
        # Skip packages that fail to collect (e.g., xgboost.testing requires hypothesis)
        print(f"Warning: Failed to collect {requirement}: {e}")
        continue
    datas += datas_
    hidden_imports_ = set(hidden_imports_)
    if "" in hidden_imports_:
        hidden_imports_.remove("")
    if None in hidden_imports_:
        hidden_imports_.remove(None)
    requirements |= hidden_imports_ - checked
    hidden_imports |= hidden_imports_

hidden_imports = sorted([h for h in hidden_imports if "tests" not in h.split(".")])
hidden_imports = [h for h in hidden_imports if "__pycache__" not in h]

# Add hdf5plugin imports to fix runtime import issues
hidden_imports.extend([
    "hdf5plugin.plugins.bshuf",
    "hdf5plugin.plugins.blosc",
    "hdf5plugin.plugins.blosc2",
    "hdf5plugin.plugins.lz4",
    "hdf5plugin.plugins.fcidecomp",
    "hdf5plugin.plugins.zfp",
    "hdf5plugin.plugins.zstd",
])

datas = [
    d
    for d in datas
    if ("__pycache__" not in d[0]) and (d[1] not in [".", "build", "dist", "Output"])
]
datas += [("ms2rescore\package_data", "package_data")]

block_cipher = None
# Build package
a = Analysis(
    [script_name],
    pathex=[location],
    binaries=binaries,
    datas=datas,
    hiddenimports=hidden_imports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name=exe_name,
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    windowed=True,
    disable_windowed_traceback=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon="./img/ms2rescore.ico",
)

coll = COLLECT(
    exe, a.binaries, a.zipfiles, a.datas, strip=False, upx=True, upx_exclude=[], name=exe_name
)
