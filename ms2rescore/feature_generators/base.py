from abc import ABC, abstractmethod
from typing import Set

from psm_utils import PSMList

from ms2rescore.parse_spectra import MSDataType


class FeatureGeneratorBase(ABC):
    """Base class from which all feature generators must inherit."""

    # List of required MS data types for feature generation
    required_ms_data: Set[MSDataType] = set()

    def __init__(self, *args, **kwargs) -> None:
        super().__init__()

    @property
    @abstractmethod
    def feature_names(self):
        pass

    @abstractmethod
    def add_features(self, psm_list: PSMList):
        pass


class FeatureGeneratorException(Exception):
    """Base class for exceptions raised by feature generators."""

    pass
