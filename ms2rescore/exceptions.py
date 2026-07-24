"""MS²Rescore exceptions."""


class MS2RescoreError(Exception):
    """Generic MS2Rescore error."""



class MS2RescoreConfigurationError(MS2RescoreError):
    """Invalid MS2Rescore configuration."""



class IDFileParsingError(MS2RescoreError):
    """Identification file parsing error."""



class ModificationParsingError(IDFileParsingError):
    """Identification file parsing error."""



class MissingValuesError(MS2RescoreError):
    """Missing values in PSMs and/or spectra."""



class ReportGenerationError(MS2RescoreError):
    """Error while generating report."""



class RescoringError(MS2RescoreError):
    """Error while rescoring PSMs."""



class ParseSpectrumError(MS2RescoreError):
    """Error while parsing spectrum files."""

