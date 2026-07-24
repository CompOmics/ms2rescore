from abc import ABC, abstractmethod

from psm_utils import PSMList

from ms2rescore.parse_spectra import MSDataType


class FeatureGeneratorBase(ABC):
    """Base class from which all feature generators must inherit."""

    # List of required MS data types for feature generation
    required_ms_data: set[MSDataType] = set()

    def __init__(self, *args, **kwargs) -> None:
        super().__init__()

    @property
    @abstractmethod
    def feature_names(self) -> list[str]:
        pass

    @abstractmethod
    def add_features(self, psm_list: PSMList) -> None:
        pass


class FeatureGeneratorException(Exception):
    """Base class for exceptions raised by feature generators."""

