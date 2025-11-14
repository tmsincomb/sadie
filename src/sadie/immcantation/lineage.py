"""
Lineage Analysis Module

Provides tools for B-cell lineage tree reconstruction and analysis.
Integrates with Dowser and other phylogenetic tools.
"""

import logging
import subprocess
from pathlib import Path
from typing import Optional, List, Dict, Tuple
import pandas as pd
import json

from sadie.immcantation.config import LineageConfig

logger = logging.getLogger(__name__)


class LineageAnalyzer:
    """
    Analyzer for B-cell lineage trees.

    Provides methods for filtering clones, building trees, and analyzing
    phylogenetic relationships within clonal families.
    """

    def __init__(self, config: Optional[LineageConfig] = None):
        """
        Initialize lineage analyzer.

        Parameters
        ----------
        config : LineageConfig, optional
            Configuration object for lineage analysis
        """
        self.config = config or LineageConfig()
        self.logger = logging.getLogger(self.__class__.__name__)

    def _run_command(
        self,
        cmd: List[str],
        capture_output: bool = True,
        check: bool = True,
    ) -> Tuple[int, str, str]:
        """
        Run a command-line tool.

        Parameters
        ----------
        cmd : List[str]
            Command and arguments
        capture_output : bool
            Whether to capture stdout/stderr
        check : bool
            Whether to raise exception on error

        Returns
        -------
        Tuple[int, str, str]
            Return code, stdout, stderr
        """
        self.logger.debug(f"Running command: {' '.join(cmd)}")

        result = subprocess.run(
            cmd,
            capture_output=capture_output,
            text=True,
            check=False,
        )

        if check and result.returncode != 0:
            raise RuntimeError(
                f"Command failed with return code {result.returncode}:\n"
                f"Command: {' '.join(cmd)}\n"
                f"Stderr: {result.stderr}"
            )

        return result.returncode, result.stdout, result.stderr

    def filter_clones_by_size(
        self,
        df: pd.DataFrame,
        clone_field: str = "clone_id",
        min_size: Optional[int] = None,
        max_size: Optional[int] = None,
    ) -> pd.DataFrame:
        """
        Filter database to keep only clones within size range.

        Parameters
        ----------
        df : pd.DataFrame
            Input database
        clone_field : str
            Column name for clone ID
        min_size : int, optional
            Minimum clone size
        max_size : int, optional
            Maximum clone size

        Returns
        -------
        pd.DataFrame
            Filtered database
        """
        min_size = min_size or self.config.min_clone_size
        max_size = max_size or self.config.max_clone_size

        if clone_field not in df.columns:
            self.logger.warning(f"Clone field '{clone_field}' not found")
            return df

        # Calculate clone sizes
        clone_sizes = df[clone_field].value_counts()

        # Filter by size
        valid_clones = clone_sizes[
            (clone_sizes >= min_size) & (clone_sizes <= max_size)
        ].index

        filtered_df = df[df[clone_field].isin(valid_clones)].copy()

        self.logger.info(
            f"Filtered {len(df)} -> {len(filtered_df)} sequences "
            f"({len(valid_clones)} clones with size {min_size}-{max_size})"
        )

        return filtered_df

    def split_clones(
        self,
        df: pd.DataFrame,
        output_dir: Path,
        clone_field: str = "clone_id",
    ) -> Dict[str, Path]:
        """
        Split database into separate files per clone.

        Parameters
        ----------
        df : pd.DataFrame
            Input database
        output_dir : Path
            Output directory for clone files
        clone_field : str
            Column name for clone ID

        Returns
        -------
        Dict[str, Path]
            Mapping of clone ID to file path
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        if clone_field not in df.columns:
            self.logger.error(f"Clone field '{clone_field}' not found")
            return {}

        clone_files = {}

        for clone_id, clone_df in df.groupby(clone_field):
            if pd.isna(clone_id):
                continue

            output_file = output_dir / f"clone_{clone_id}.tsv"
            clone_df.to_csv(output_file, sep="\t", index=False)
            clone_files[str(clone_id)] = output_file

        self.logger.info(f"Split {len(df)} sequences into {len(clone_files)} clone files")

        return clone_files

    def build_trees_simple(
        self,
        clone_files: Dict[str, Path],
        output_dir: Path,
        method: str = "nj",
    ) -> Dict[str, Path]:
        """
        Build simple phylogenetic trees using neighbor-joining.

        This is a fallback method that doesn't require R/Dowser.

        Parameters
        ----------
        clone_files : Dict[str, Path]
            Mapping of clone ID to file path
        output_dir : Path
            Output directory for trees
        method : str
            Tree building method (nj, upgma)

        Returns
        -------
        Dict[str, Path]
            Mapping of clone ID to tree file
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        tree_files = {}

        for clone_id, clone_file in clone_files.items():
            try:
                # Read clone sequences
                clone_df = pd.read_csv(clone_file, sep="\t")

                if len(clone_df) < 2:
                    self.logger.debug(f"Clone {clone_id} has <2 sequences, skipping")
                    continue

                # Build tree using Python (Bio.Phylo)
                tree_file = output_dir / f"clone_{clone_id}_tree.nwk"
                self._build_tree_biopython(clone_df, tree_file, method)

                tree_files[clone_id] = tree_file

            except Exception as e:
                self.logger.warning(f"Failed to build tree for clone {clone_id}: {e}")

        self.logger.info(f"Built {len(tree_files)} trees")

        return tree_files

    def _build_tree_biopython(
        self,
        clone_df: pd.DataFrame,
        output_file: Path,
        method: str = "nj",
    ) -> None:
        """
        Build tree using Biopython.

        Parameters
        ----------
        clone_df : pd.DataFrame
            Clone sequences
        output_file : Path
            Output tree file
        method : str
            Tree building method
        """
        try:
            from Bio import Phylo
            from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
            from Bio.Align import MultipleSeqAlignment
            from Bio.SeqRecord import SeqRecord
            from Bio.Seq import Seq

            # Create alignment from sequences
            sequences = []
            for idx, row in clone_df.iterrows():
                if "sequence_alignment" in row:
                    seq = row["sequence_alignment"]
                elif "sequence" in row:
                    seq = row["sequence"]
                else:
                    continue

                seq_id = row.get("sequence_id", f"seq_{idx}")
                sequences.append(SeqRecord(Seq(seq), id=str(seq_id)))

            if len(sequences) < 2:
                self.logger.debug("Not enough sequences for tree building")
                return

            # Pad sequences to same length
            max_len = max(len(s.seq) for s in sequences)
            for seq_rec in sequences:
                if len(seq_rec.seq) < max_len:
                    seq_rec.seq = seq_rec.seq + Seq("-" * (max_len - len(seq_rec.seq)))

            alignment = MultipleSeqAlignment(sequences)

            # Calculate distance matrix
            calculator = DistanceCalculator("identity")
            dm = calculator.get_distance(alignment)

            # Build tree
            constructor = DistanceTreeConstructor()
            if method == "nj":
                tree = constructor.nj(dm)
            else:  # upgma
                tree = constructor.upgma(dm)

            # Write tree
            Phylo.write(tree, output_file, "newick")

        except ImportError:
            self.logger.warning("Biopython Phylo not available for tree building")
        except Exception as e:
            self.logger.warning(f"Failed to build tree with Biopython: {e}")

    def calculate_tree_stats(
        self,
        tree_files: Dict[str, Path],
    ) -> pd.DataFrame:
        """
        Calculate statistics for phylogenetic trees.

        Parameters
        ----------
        tree_files : Dict[str, Path]
            Mapping of clone ID to tree file

        Returns
        -------
        pd.DataFrame
            Tree statistics
        """
        stats = []

        for clone_id, tree_file in tree_files.items():
            try:
                from Bio import Phylo

                tree = Phylo.read(tree_file, "newick")

                stats.append(
                    {
                        "clone_id": clone_id,
                        "num_terminals": len(tree.get_terminals()),
                        "num_nonterminals": len(tree.get_nonterminals()),
                        "tree_depth": tree.depths().get(tree.root, 0)
                        if tree.root
                        else 0,
                        "total_branch_length": tree.total_branch_length(),
                    }
                )

            except Exception as e:
                self.logger.warning(f"Failed to calculate stats for clone {clone_id}: {e}")

        return pd.DataFrame(stats) if stats else pd.DataFrame()

    def summarize_lineages(
        self,
        clone_df: pd.DataFrame,
        tree_files: Dict[str, Path],
        output_file: Path,
    ) -> pd.DataFrame:
        """
        Create summary of lineage analysis.

        Parameters
        ----------
        clone_df : pd.DataFrame
            Database with clone assignments
        tree_files : Dict[str, Path]
            Tree files
        output_file : Path
            Output summary file

        Returns
        -------
        pd.DataFrame
            Summary statistics
        """
        summary = []

        # Clone size distribution
        clone_sizes = clone_df["clone_id"].value_counts()

        # V gene usage per clone
        v_gene_usage = (
            clone_df.groupby("clone_id")["v_call"]
            .apply(lambda x: x.mode()[0] if len(x.mode()) > 0 else None)
            .to_dict()
        )

        # Junction length per clone
        junction_lengths = clone_df.groupby("clone_id")["junction_aa_length"].mean().to_dict()

        for clone_id, size in clone_sizes.items():
            if pd.isna(clone_id):
                continue

            summary.append(
                {
                    "clone_id": clone_id,
                    "clone_size": size,
                    "v_gene": v_gene_usage.get(clone_id, "NA"),
                    "mean_junction_length": junction_lengths.get(clone_id, 0),
                    "has_tree": str(clone_id) in tree_files,
                }
            )

        summary_df = pd.DataFrame(summary)

        if not summary_df.empty:
            summary_df = summary_df.sort_values("clone_size", ascending=False)
            summary_df.to_csv(output_file, sep="\t", index=False)
            self.logger.info(f"Saved lineage summary: {output_file}")

        return summary_df
