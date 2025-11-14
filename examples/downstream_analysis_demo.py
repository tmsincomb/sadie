#!/usr/bin/env python3
"""
Downstream VDJ Analysis Demo

Demonstrates clonal clustering and lineage tree analysis using
pre-annotated AIRR data. This workflow represents the downstream
analysis portion of the Immcantation pipeline.

Usage:
    python downstream_analysis_demo.py

Author: SADIE Team
Date: 2025-11-14
"""

import sys
from pathlib import Path
import pandas as pd
import logging

# Add sadie to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from sadie.airr.airrtable import AirrTable
from sadie.airr.methods import run_mutational_analysis, run_igl_assignment
from sadie.cluster import Cluster
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def build_clone_tree(clone_df, output_file):
    """Build a neighbor-joining tree for a clonal lineage"""
    try:
        sequences = []

        for idx, row in clone_df.iterrows():
            # Use sequence_alignment if available, otherwise sequence
            seq = row.get("sequence_alignment", row.get("sequence", ""))
            if not seq or pd.isna(seq):
                continue

            seq_id = row.get("sequence_id", f"seq_{idx}")
            sequences.append(SeqRecord(Seq(str(seq)), id=str(seq_id)[:20]))  # Limit ID length

        if len(sequences) < 2:
            logger.debug(f"    Not enough sequences ({len(sequences)}) for tree building")
            return False

        # Pad sequences to same length
        max_len = max(len(s.seq) for s in sequences)
        for seq_rec in sequences:
            if len(seq_rec.seq) < max_len:
                gap_seq = "-" * (max_len - len(seq_rec.seq))
                seq_rec.seq = seq_rec.seq + Seq(gap_seq)

        # Build alignment and tree
        alignment = MultipleSeqAlignment(sequences)
        calculator = DistanceCalculator("identity")
        dm = calculator.get_distance(alignment)

        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)

        # Save tree
        Phylo.write(tree, output_file, "newick")

        return True

    except Exception as e:
        logger.warning(f"    Failed to build tree: {e}")
        return False


def main():
    output_dir = Path("/tmp/downstream_analysis_demo")
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 80)
    logger.info("DOWNSTREAM VDJ ANALYSIS DEMONSTRATION")
    logger.info("=" * 80)
    logger.info(f"Output: {output_dir}")
    logger.info("=" * 80)

    # Load pre-annotated AIRR data
    logger.info("\n[STEP 1] Loading pre-annotated AIRR data...")
    test_data = Path(__file__).parent.parent / "tests/data/fixtures/airr_tables"

    # Try to load heavy chain data
    heavy_file = test_data / "heavy_airrtable_with_mutational.feather"

    if heavy_file.exists():
        airr_table = pd.read_feather(str(heavy_file))
        airr_table = AirrTable(airr_table)
        logger.info(f"  Loaded {len(airr_table)} heavy chain sequences")
    else:
        # Fallback to catnap data
        heavy_file = test_data / "catnap_heavy_airrtable.feather"
        airr_table = pd.read_feather(str(heavy_file))
        airr_table = AirrTable(airr_table)
        logger.info(f"  Loaded {len(airr_table)} sequences from CAT-NAP dataset")

    # Show basic stats
    logger.info(f"  Loci: {airr_table['locus'].value_counts().to_dict()}")
    if "productive" in airr_table.columns:
        n_productive = (airr_table["productive"] == True).sum()
        logger.info(f"  Productive: {n_productive} / {len(airr_table)}")

    # Step 2: Filter to productive sequences
    logger.info("\n[STEP 2] Filtering productive sequences...")
    if "productive" in airr_table.columns:
        productive = airr_table[airr_table["productive"] == True].copy()
    else:
        productive = airr_table.copy()

    logger.info(f"  Retained {len(productive)} sequences")

    # Save filtered
    filtered_file = output_dir / "01_productive.tsv"
    productive.to_csv(filtered_file, sep="\t", index=False)

    # Step 3: Clonal Clustering
    logger.info("\n[STEP 3] Performing clonal clustering...")
    logger.info("  Using CDR1, CDR2, CDR3 amino acid sequences")

    cluster_api = Cluster(
        productive,
        linkage="single",
        lookup=["cdr1_aa", "cdr2_aa", "cdr3_aa"],
    )

    clustered = cluster_api.cluster(0.15)  # distance threshold

    # Save clustered data
    cluster_file = output_dir / "02_clustered.tsv"
    clustered.to_csv(cluster_file, sep="\t", index=False)

    # Analyze clusters
    clone_sizes = clustered["cluster"].value_counts()
    logger.info(f"  Identified {len(clone_sizes)} clonal groups")
    logger.info(f"  Clone size range: {clone_sizes.min()} - {clone_sizes.max()}")
    logger.info(f"  Mean clone size: {clone_sizes.mean():.1f}")

    # Step 4: Clone Statistics
    logger.info("\n[STEP 4] Calculating clone statistics...")

    clone_stats = []
    for cluster_id in clone_sizes.index:
        cluster_df = clustered[clustered["cluster"] == cluster_id]

        # Get modal V and J genes
        v_gene = cluster_df["v_call"].mode()[0] if len(cluster_df) > 0 else "NA"
        j_gene = cluster_df["j_call"].mode()[0] if len(cluster_df) > 0 else "NA"

        # Get mean CDR3 length
        if "cdr3_aa" in cluster_df.columns:
            cdr3_lens = cluster_df["cdr3_aa"].str.len()
            mean_cdr3_len = cdr3_lens.mean()
        else:
            mean_cdr3_len = 0

        clone_stats.append({
            "clone_id": cluster_id,
            "size": len(cluster_df),
            "v_gene": v_gene,
            "j_gene": j_gene,
            "mean_cdr3_length": mean_cdr3_len,
        })

    stats_df = pd.DataFrame(clone_stats).sort_values("size", ascending=False)

    # Save stats
    stats_file = output_dir / "03_clone_statistics.tsv"
    stats_df.to_csv(stats_file, sep="\t", index=False)
    logger.info(f"  Saved clone statistics: {stats_file}")

    # Show top clones
    logger.info("\n  Top 10 clones by size:")
    logger.info(stats_df.head(10).to_string(index=False))

    # Step 5: Build Lineage Trees
    logger.info("\n[STEP 5] Building lineage trees for large clones...")

    min_clone_size = 5
    max_clone_size = 100

    large_clones = stats_df[
        (stats_df["size"] >= min_clone_size) &
        (stats_df["size"] <= max_clone_size)
    ]

    logger.info(f"  Building trees for {len(large_clones)} clones (size {min_clone_size}-{max_clone_size})")

    trees_dir = output_dir / "trees"
    trees_dir.mkdir(exist_ok=True)

    trees_built = 0
    trees_failed = 0

    for idx, clone_info in large_clones.iterrows():
        clone_id = clone_info["clone_id"]
        clone_seqs = clustered[clustered["cluster"] == clone_id]

        tree_file = trees_dir / f"clone_{clone_id}.nwk"

        if build_clone_tree(clone_seqs, tree_file):
            trees_built += 1
        else:
            trees_failed += 1

    logger.info(f"  Successfully built {trees_built} trees")
    if trees_failed > 0:
        logger.info(f"  Failed to build {trees_failed} trees")

    # Step 6: Tree Statistics
    if trees_built > 0:
        logger.info("\n[STEP 6] Calculating tree statistics...")

        tree_stats = []
        for tree_file in trees_dir.glob("*.nwk"):
            try:
                tree = Phylo.read(tree_file, "newick")
                clone_id = tree_file.stem.replace("clone_", "")

                tree_stats.append({
                    "clone_id": clone_id,
                    "terminals": len(tree.get_terminals()),
                    "internal_nodes": len(tree.get_nonterminals()),
                    "total_branch_length": tree.total_branch_length(),
                    "tree_depth": max(tree.depths().values()) if tree.depths() else 0,
                })

            except Exception as e:
                logger.warning(f"  Failed to analyze tree {tree_file}: {e}")

        if tree_stats:
            tree_stats_df = pd.DataFrame(tree_stats)
            tree_stats_file = output_dir / "04_tree_statistics.tsv"
            tree_stats_df.to_csv(tree_stats_file, sep="\t", index=False)
            logger.info(f"  Saved tree statistics: {tree_stats_file}")

            logger.info(f"\n  Mean tree depth: {tree_stats_df['tree_depth'].mean():.3f}")
            logger.info(f"  Mean branch length: {tree_stats_df['total_branch_length'].mean():.3f}")

    # Final Summary
    logger.info("\n" + "=" * 80)
    logger.info("ANALYSIS COMPLETE")
    logger.info("=" * 80)
    logger.info(f"Total sequences:       {len(airr_table)}")
    logger.info(f"Productive sequences:  {len(productive)}")
    logger.info(f"Clonal groups:         {len(clone_sizes)}")
    logger.info(f"Large clones:          {len(large_clones)}")
    logger.info(f"Lineage trees built:   {trees_built}")
    logger.info(f"\nOutput directory:      {output_dir}")
    logger.info("\nOutput files:")
    for f in sorted(output_dir.glob("*.tsv")):
        size_kb = f.stat().st_size / 1024
        logger.info(f"  - {f.name:40s} ({size_kb:,.1f} KB)")
    logger.info(f"  - trees/ ({trees_built} Newick files)")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
