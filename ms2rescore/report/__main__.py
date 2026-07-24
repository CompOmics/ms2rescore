import logging
import sys
from pathlib import Path

import click
from rich.logging import RichHandler

from ms2rescore.report.data import ReportData
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
@click.option(
    "--fdr",
    type=click.FloatRange(0, 1),
    default=None,
    help=(
        "FDR threshold for the report's identification stats/charts. Defaults to the "
        "report_fdr used in the original run (from its saved full-config.json)."
    ),
)
def main(psm_file, output, fdr):
    """Generate MS²Rescore report from a PSM TSV file.

    PSM_FILE: Path to the main PSM TSV file (e.g., output.ms2rescore.tsv)
    """
    logging.basicConfig(
        level=logging.INFO,
        handlers=[RichHandler(rich_tracebacks=True)],
        format="%(message)s",
    )

    try:
        # The main PSM list is always written to "<output_prefix>.tsv"
        output_prefix = str(Path(psm_file).with_suffix(""))

        # Determine output path
        if output:
            output_path = Path(output)
        else:
            output_path = Path(output_prefix + ".report.html")

        logger.info("Generating report...")
        report_data = ReportData.from_files(output_prefix, fdr_threshold=fdr)
        generate_report(output_prefix, report_data, output_file=output_path)

        logger.info(f"✓ Report generated: {output_path}")

    except Exception as e:
        logger.exception(e)
        sys.exit(1)


if __name__ == "__main__":
    main()
