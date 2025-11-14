"""
Configuration module for Immcantation pipeline

Provides default parameters and configuration options for the pipeline.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Dict, List


@dataclass
class PrestoConfig:
    """Configuration for pRESTO processing"""

    # Quality filtering
    min_quality: int = 20
    min_length: int = 200
    max_error: float = 0.1
    max_missing: int = 10

    # UMI parameters
    umi_length: int = 12
    umi_start: int = 0

    # Primer parameters
    max_primer_error: float = 0.2
    primer_start: int = 0
    primer_gap: int = 10

    # Consensus building
    min_consensus_freq: float = 0.6
    min_consensus_count: int = 2
    max_consensus_diversity: float = 0.1

    # Assembly (for paired-end)
    min_overlap: int = 10
    max_assembly_error: float = 0.3
    alpha: float = 1e-5

    # Deduplication
    max_duplicate_distance: int = 0

    # Performance
    nproc: int = 4

    def to_dict(self) -> Dict:
        """Convert to dictionary"""
        return {
            "min_quality": self.min_quality,
            "min_length": self.min_length,
            "max_error": self.max_error,
            "max_missing": self.max_missing,
            "umi_length": self.umi_length,
            "umi_start": self.umi_start,
            "max_primer_error": self.max_primer_error,
            "primer_start": self.primer_start,
            "primer_gap": self.primer_gap,
            "min_consensus_freq": self.min_consensus_freq,
            "min_consensus_count": self.min_consensus_count,
            "max_consensus_diversity": self.max_consensus_diversity,
            "min_overlap": self.min_overlap,
            "max_assembly_error": self.max_assembly_error,
            "alpha": self.alpha,
            "max_duplicate_distance": self.max_duplicate_distance,
            "nproc": self.nproc,
        }


@dataclass
class ChangeoConfig:
    """Configuration for Change-O processing"""

    # IgBLAST parameters
    organism: str = "human"
    loci: str = "ig"  # ig or tr
    igdata: Optional[Path] = None

    # Germline database paths (will be auto-detected if None)
    germline_dir: Optional[Path] = None
    v_field: str = "v_call"
    d_field: str = "d_call"
    j_field: str = "j_call"

    # Clonal clustering
    distance_threshold: float = 0.15
    distance_metric: str = "ham"  # ham, aa, hh_s1f, hh_s5f, hh_s7f
    linkage_method: str = "single"  # single, average, complete

    # Germline reconstruction
    germline_min_call: int = 3
    cloned: bool = False

    # Quality filters
    min_v_call_score: float = 200.0
    min_sequence_length: int = 150

    # Performance
    nproc: int = 4

    def to_dict(self) -> Dict:
        """Convert to dictionary"""
        return {
            "organism": self.organism,
            "loci": self.loci,
            "igdata": str(self.igdata) if self.igdata else None,
            "germline_dir": str(self.germline_dir) if self.germline_dir else None,
            "v_field": self.v_field,
            "d_field": self.d_field,
            "j_field": self.j_field,
            "distance_threshold": self.distance_threshold,
            "distance_metric": self.distance_metric,
            "linkage_method": self.linkage_method,
            "germline_min_call": self.germline_min_call,
            "cloned": self.cloned,
            "min_v_call_score": self.min_v_call_score,
            "min_sequence_length": self.min_sequence_length,
            "nproc": self.nproc,
        }


@dataclass
class LineageConfig:
    """Configuration for lineage tree analysis"""

    # Tree building method
    method: str = "igphyml"  # igphyml, dnapars, raxml

    # IgPhyML parameters
    igphyml_mode: str = "ml"  # ml, gtr
    optimize: str = "lr"  # l=length, r=rate, t=tree

    # Filtering
    min_clone_size: int = 5
    max_clone_size: int = 1000

    # Visualization
    plot_trees: bool = True
    tree_format: str = "pdf"

    # Performance
    nproc: int = 2

    def to_dict(self) -> Dict:
        """Convert to dictionary"""
        return {
            "method": self.method,
            "igphyml_mode": self.igphyml_mode,
            "optimize": self.optimize,
            "min_clone_size": self.min_clone_size,
            "max_clone_size": self.max_clone_size,
            "plot_trees": self.plot_trees,
            "tree_format": self.tree_format,
            "nproc": self.nproc,
        }


@dataclass
class PipelineConfig:
    """Main pipeline configuration"""

    # Sub-configurations
    presto: PrestoConfig = field(default_factory=PrestoConfig)
    changeo: ChangeoConfig = field(default_factory=ChangeoConfig)
    lineage: LineageConfig = field(default_factory=LineageConfig)

    # Input/Output
    input_format: str = "fasta"  # fasta, fastq, fastq-pe
    output_dir: Path = field(default_factory=lambda: Path("./immcantation_output"))

    # Pipeline steps to run
    run_presto: bool = False  # Set to False for pre-processed data
    run_changeo: bool = True
    run_lineage: bool = True

    # Logging
    verbose: bool = True
    log_file: Optional[Path] = None

    # Cleanup
    keep_intermediate: bool = True

    def to_dict(self) -> Dict:
        """Convert to dictionary"""
        return {
            "presto": self.presto.to_dict(),
            "changeo": self.changeo.to_dict(),
            "lineage": self.lineage.to_dict(),
            "input_format": self.input_format,
            "output_dir": str(self.output_dir),
            "run_presto": self.run_presto,
            "run_changeo": self.run_changeo,
            "run_lineage": self.run_lineage,
            "verbose": self.verbose,
            "log_file": str(self.log_file) if self.log_file else None,
            "keep_intermediate": self.keep_intermediate,
        }


# Default primer sets for common protocols
PRIMER_SETS = {
    "human_heavy_5race": {
        "name": "Human IgH 5'RACE",
        "primers": [
            {"name": "IGHV1", "seq": "ACAGGTGCCCACTCCCAGGTGCAG"},
            {"name": "IGHV3", "seq": "AGGTGCAGCTGGTGGAGTCTGG"},
            {"name": "IGHV4", "seq": "AGGTGCAGCTACAGTCAGTGG"},
        ],
    },
    "human_kappa_5race": {
        "name": "Human IgK 5'RACE",
        "primers": [
            {"name": "IGKV", "seq": "ATGAGGSTCCCYGCTCAGCTGCTGG"},
        ],
    },
    "human_lambda_5race": {
        "name": "Human IgL 5'RACE",
        "primers": [
            {"name": "IGLV", "seq": "ATGAGCTACKCCTGGGCTCCTGG"},
        ],
    },
}
