"""
Feature generators to add rescoring features to PSMs from various (re)sources and prediction tools.
"""

from ms2rescore.feature_generators.base import FeatureGeneratorBase
from ms2rescore.feature_generators.basic import BasicFeatureGenerator
from ms2rescore.feature_generators.deeplc import DeepLCFeatureGenerator
from ms2rescore.feature_generators.im2deep import IM2DeepFeatureGenerator
from ms2rescore.feature_generators.ms2 import MS2FeatureGenerator
from ms2rescore.feature_generators.ms2pip import MS2PIPFeatureGenerator

FEATURE_GENERATORS: dict[str, type[FeatureGeneratorBase]] = {
    "basic": BasicFeatureGenerator,
    "ms2pip": MS2PIPFeatureGenerator,
    "deeplc": DeepLCFeatureGenerator,
    "im2deep": IM2DeepFeatureGenerator,
    "ms2": MS2FeatureGenerator,
}
