"""
Builders - BLAST Database, Auxiliary File, and Gapping
=======================================================

Utilities for building IgBLAST-compatible databases and auxiliary files.

- BlastDBBuilder: Creates BLAST databases from FASTA files
- AuxFileBuilder: Generates IgBLAST auxiliary files from gapped sequences
- GapperService: Gaps ungapped sequences to IMGT numbering using BioPython
"""

from .blast import BlastDBBuilder
from .aux import AuxFileBuilder
from .gapper import GapperService, gap_sequences_batch

__all__ = ["BlastDBBuilder", "AuxFileBuilder", "GapperService", "gap_sequences_batch"]
