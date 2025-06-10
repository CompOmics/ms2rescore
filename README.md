<img src="https://github.com/compomics/ms2rescore/raw/main/img/ms2rescore_logo.png" width="150" height="150" alt="MS²Rescore"/>
<br/><br/>

[![GitHub release](https://img.shields.io/github/release-pre/compomics/ms2rescore.svg?style=flat-square)](https://github.com/compomics/ms2rescore/releases)
[![PyPI](https://flat.badgen.net/pypi/v/ms2rescore)](https://pypi.org/project/ms2rescore/)
[![GitHub Workflow Status](https://flat.badgen.net/github/checks/compomics/ms2rescore/main)](https://github.com/compomics/ms2rescore/actions/)
[![GitHub issues](https://img.shields.io/github/issues/compomics/ms2rescore?style=flat-square)](https://github.com/compomics/ms2rescore/issues)
[![GitHub](https://img.shields.io/github/license/compomics/ms2rescore.svg?style=flat-square)](https://www.apache.org/licenses/LICENSE-2.0)
[![Last commit](https://flat.badgen.net/github/last-commit/compomics/ms2rescore)](https://github.com/compomics/ms2rescore/commits/)

Modular and user-friendly platform for AI-assisted rescoring of peptide identifications

## About MS²Rescore

MS²Rescore performs ultra-sensitive peptide identification rescoring with LC-MS predictors such as
[MS²PIP][ms2pip] and [DeepLC][deeplc], and with ML-driven rescoring engines
[Percolator][percolator] or [Mokapot][mokapot]. This results in more confident peptide
identifications, which allows you to get **more peptide IDs** at the same false discovery rate
(FDR) threshold, or to set a **more stringent FDR threshold** while still retaining a similar
number of peptide IDs. MS²Rescore is **ideal for challenging proteomics identification workflows**,
such as proteogenomics, metaproteomics, or immunopeptidomics.

![MS²Rescore overview](https://raw.githubusercontent.com/compomics/ms2rescore/main/docs/source/_static/img/ms2rescore-overview.png)

MS²Rescore can read peptide identifications in any format supported by [psm_utils][psm_utils]
(see [Supported file formats][file-formats]) and has been tested with various search engines output
files:

- [MS Amanda](http://ms.imp.ac.at/?goto=msamanda) `.csv`
- [Sage](https://github.com/lazear/sage) `.sage.tsv`
- [PeptideShaker](https://compomics.github.io/projects/peptide-shaker.html) `.mzid`
- [ProteomeDiscoverer](#)`.msf`
- [MSGFPlus](https://omics.pnl.gov/software/ms-gf) `.mzid`
- [Mascot](https://www.matrixscience.com/) `.mzid`
- [MaxQuant](https://www.maxquant.org/) `msms.txt`
- [X!Tandem](https://www.thegpm.org/tandem/) `.xml`
- [PEAKS](https://www.bioinfor.com/peaksdb/) `.mzid`

MS²Rescore is available as a [desktop application][desktop], a [command line tool][cli], and a
[modular Python API][python-package].

## TIMS²Rescore: Direct support for DDA-PASEF data

MS²Rescore v3.1+ includes TIMS²Rescore, a usage mode with specialized default configurations for
DDA-PASEF data from timsTOF instruments. TIMS²Rescore makes use of new MS²PIP prediction models for
timsTOF fragmentation and IM2Deep for ion mobility separation. Bruker .d and miniTDF spectrum
files are directly supported through the [timsrust](https://github.com/MannLabs/timsrust) library.

Checkout our [paper](https://doi.org/10.1021/acs.jproteome.4c00609) for more information and the
[TIMS²Rescore documentation][tims2rescore] to get started.

## Citing

**Latest MS²Rescore publication:**

> **MS²Rescore 3.0 is a modular, flexible, and user-friendly platform to boost peptide identifications, as showcased with MS Amanda 3.0.**
> Louise Marie Buur*, Arthur Declercq*, Marina Strobl, Robbin Bouwmeester, Sven Degroeve, Lennart Martens, Viktoria Dorfer*, and Ralf Gabriels*.
> _Journal of Proteome Research_ (2024) [doi:10.1021/acs.jproteome.3c00785](https://doi.org/10.1021/acs.jproteome.3c00785) <br/> \*contributed equally <span class="__dimensions_badge_embed__" data-doi="10.1021/acs.jproteome.3c00785" data-hide-zero-citations="true" data-style="small_rectangle"></span>

**MS²Rescore for immunopeptidomics:**

> **MS²Rescore: Data-driven rescoring dramatically boosts immunopeptide identification rates.**
> Arthur Declercq, Robbin Bouwmeester, Aurélie Hirschler, Christine Carapito, Sven Degroeve, Lennart Martens, and Ralf Gabriels.
> _Molecular & Cellular Proteomics_ (2021) [doi:10.1016/j.mcpro.2022.100266](https://doi.org/10.1016/j.mcpro.2022.100266) <span class="__dimensions_badge_embed__" data-doi="10.1016/j.mcpro.2022.100266" data-hide-zero-citations="true" data-style="small_rectangle"></span>

**MS²Rescore for timsTOF DDA-PASEF data:**

> **TIMS²Rescore: A DDA-PASEF optimized data-driven rescoring pipeline based on MS²Rescore.**
> Arthur Declercq*, Robbe Devreese*, Jonas Scheid, Caroline Jachmann, Tim Van Den Bossche, Annica Preikschat, David Gomez-Zepeda, Jeewan Babu Rijal, Aurélie Hirschler, Jonathan R Krieger, Tharan Srikumar, George Rosenberger, Dennis Trede, Christine Carapito, Stefan Tenzer, Juliane S Walz, Sven Degroeve, Robbin Bouwmeester, Lennart Martens, and Ralf Gabriels.
> _Journal of Proteome Research_ (2025) [doi:10.1021/acs.jproteome.4c00609](https://doi.org/10.1021/acs.jproteome.4c00609) <span class="__dimensions_badge_embed__" data-doi="10.1021/acs.jproteome.4c00609" data-hide-zero-citations="true" data-style="small_rectangle"></span>

**Original publication describing the concept of rescoring with predicted spectra:**

> **Accurate peptide fragmentation predictions allow data driven approaches to replace and improve upon proteomics search engine scoring functions.**
> Ana S C Silva, Robbin Bouwmeester, Lennart Martens, and Sven Degroeve.
> _Bioinformatics_ (2019) [doi:10.1093/bioinformatics/btz383](https://doi.org/10.1093/bioinformatics/btz383) <span class="__dimensions_badge_embed__" data-doi="10.1093/bioinformatics/btz383" data-hide-zero-citations="true" data-style="small_rectangle"></span>

To replicate the experiments described in this article, check out the
[publication branch][publication-branch] of the repository.

## Getting started

The desktop application can be installed on Windows with a [one-click installer][desktop-installer].
The Python package and command line interface can be installed with `pip`, `conda`, or `docker`.
Check out the [full documentation][docs] to get started.

## Questions or issues?

Have questions on how to apply MS²Rescore on your data? Or ran into issues while using MS²Rescore?
Post your questions on the [GitHub Discussions][discussions] forum and we are happy to help!

## How to contribute

Bugs, questions or suggestions? Feel free to post an issue in the [issue tracker][issues] or to
make a [pull request][pr]!

[docs]: https://ms2rescore.readthedocs.io/
[issues]: https://github.com/compomics/ms2rescore/issues/
[discussions]: https://github.com/compomics/ms2rescore/discussions/
[pr]: https://github.com/compomics/ms2rescore/pulls/
[desktop]: https://ms2rescore.readthedocs.io/en/stable/gui/
[desktop-installer]: https://github.com/compomics/ms2rescore/releases/latest
[cli]: https://ms2rescore.readthedocs.io/en/stable/cli/
[python-package]: https://ms2rescore.readthedocs.io/en/stable/api/ms2rescore/
[docker]: https://ms2rescore.readthedocs.io/en/stable/installation#docker-container
[publication-branch]: https://github.com/compomics/ms2rescore/tree/pub
[ms2pip]: https://github.com/compomics/ms2pip
[deeplc]: https://github.com/compomics/deeplc
[percolator]: https://github.com/percolator/percolator/
[mokapot]: https://mokapot.readthedocs.io/
[psm_utils]: https://github.com/compomics/psm_utils
[file-formats]: https://psm-utils.readthedocs.io/en/stable/#supported-file-formats
[tims2rescore]: https://ms2rescore.readthedocs.io/en/stable/userguide/tims2Rescore
