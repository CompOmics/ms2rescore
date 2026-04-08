import logging
from pathlib import Path

import click
import psm_utils.io
from rich.logging import RichHandler

from ms2rescore.report.generate import generate_report

logger = logging.getLogger(__name__)


@click.command()
@click.argument("psm_file", type=click.Path(exists=True))
@click.option(
    "--output",
    "-o",
    type=click.Path(),
    default=None,
    help="Output path for the report HTML file. If not provided, will be based on PSM file name.",
)
def main(psm_file, output):
    """Generate MS²Rescore report from a PSM TSV file.

    PSM_FILE: Path to the PSM TSV file (e.g., output.psms.tsv)
    """
    logging.getLogger("mokapot").setLevel(logging.WARNING)
    logging.basicConfig(
        level=logging.INFO,
        handlers=[RichHandler(rich_tracebacks=True)],
        format="%(message)s",
    )

    try:
        psm_file_path = Path(psm_file)

        # Infer output prefix from PSM file name
        if ".ms2rescore.psms.tsv" in psm_file_path.name:
            output_prefix = str(psm_file_path).replace(".psms.tsv", "")
        else:
            output_prefix = str(psm_file_path.with_suffix(""))

        # Determine output path
        if output:
            output_path = Path(output)
        else:
            output_path = Path(output_prefix + ".report.html")

        logger.info(f"Reading PSMs from {psm_file_path}...")
        psm_list = psm_utils.io.read_file(psm_file_path, filetype="tsv", show_progressbar=True)

        logger.info("Generating report...")
        generate_report(
            output_path_prefix=output_prefix,
            psm_list=psm_list,
            output_file=output_path,
        )

        logger.info(f"✓ Report generated: {output_path}")

    except Exception as e:
        logger.exception(e)
        exit(1)


if __name__ == "__main__":
    main()
