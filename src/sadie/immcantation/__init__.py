"""
SADIE Immcantation Pipeline Module

This module provides a comprehensive wrapper for the Immcantation framework
(pRESTO, Change-O, Dowser) for VDJ antibody analysis from Illumina reads.

Pipeline workflow:
1. pRESTO: Quality control, UMI processing, primer identification, consensus building
2. Change-O: VDJ alignment with IgBLAST, clonal clustering
3. Lineage Analysis: Germline reconstruction, phylogenetic tree building

Author: SADIE Team
License: MIT
"""

from sadie.immcantation.pipeline import ImmcantationPipeline
from sadie.immcantation.presto_wrapper import PrestoWrapper
from sadie.immcantation.changeo_wrapper import ChangeoWrapper
from sadie.immcantation.lineage import LineageAnalyzer

__all__ = [
    "ImmcantationPipeline",
    "PrestoWrapper",
    "ChangeoWrapper",
    "LineageAnalyzer",
]

__version__ = "1.0.0"
