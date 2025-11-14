#!/usr/bin/env python3
"""
Simplified VDJ analysis pipeline using SADIE's existing tools.

This bypasses Change-O and uses SADIE's native IgBLAST integration,
then adds clonal clustering and basic lineage analysis.

Usage:
    python run_simple_vdj_pipeline.py INPUT.fasta OUTPUT_DIR

Author: SADIE Team
Date: 2025-11-14
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import logging

# Add sadie to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from sadie.airr import Airr
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


def main():
    parser = argparse.ArgumentParser(description="Simple VDJ analysis pipeline")
    parser.add_argument("input", type=Path, help="Input FASTA file")
    parser.add_argument("output", type=Path, help="Output directory")
    parser.add_argument("--organism", default="human", help="Organism (default: human)")
    parser.add_argument("--distance", type=float, default=0.15, help="Clone distance threshold")
    parser.add_argument("--min-clone-size", type=int, default=3, help="Minimum clone size")

    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 80)
    logger.info("SIMPLIFIED VDJ ANALYSIS PIPELINE")
    logger.info("=" * 80)
    logger.info(f"Input:    {args.input}")
    logger.info(f"Output:   {args.output}")
    logger.info(f"Organism: {args.organism}")
    logger.info("=" * 80)

    try:
        # Step 1: VDJ Annotation with SADIE
        logger.info("\n[STEP 1] Running VDJ annotation with SADIE IgBLAST...")
        airr_api = Airr(args.organism, num_cpus=2)
        airr_table = airr_api.run_fasta(str(args.input))

        logger.info(f"  Annotated {len(airr_table)} sequences")

        # Save annotation
        annotation_file = args.output / "01_annotation.tsv"
        airr_table.to_airr(str(annotation_file))
        logger.info(f"  Saved: {annotation_file}")

        # Step 2: Filter productive
        logger.info("\n[STEP 2] Filtering productive sequences...")
        productive = airr_table[airr_table["productive"] == True]
        logger.info(f"  {len(productive)} / {len(airr_table)} productive sequences")

        productive_file = args.output / "02_productive.tsv"
        productive.to_airr(str(productive_file))

        # Step 3: Mutational Analysis
        logger.info("\n[STEP 3] Running mutational analysis...")
        productive = run_mutational_analysis(productive, scheme="imgt")
        productive = run_igl_assignment(productive)

        analysis_file = args.output / "03_analyzed.tsv"
        productive.to_airr(str(analysis_file))
        logger.info(f"  Saved: {analysis_file}")

        # Step 4: Clonal Clustering
        logger.info("\n[STEP 4] Clonal clustering...")
        cluster_api = Cluster(
            productive,
            linkage="single",
            lookup=["cdr1_aa", "cdr2_aa", "cdr3_aa"],
        )

        clustered = cluster_api.cluster(distance=args.distance)

        cluster_file = args.output / "04_clustered.tsv"
        clustered.to_airr(str(cluster_file))

        # Count clones
        clone_sizes = clustered["cluster"].value_counts()
        logger.info(f"  Found {len(clone_sizes)} clonal groups")
        logger.info(f"  Size range: {clone_sizes.min()} - {clone_sizes.max()}")

        # Step 5: Clone Statistics
        logger.info("\n[STEP 5] Calculating clone statistics...")
        clone_stats = []

        for cluster_id, cluster_df in clustered.groupby("cluster"):
            clone_stats.append({
                "clone_id": cluster_id,
                "size": len(cluster_df),
                "v_gene": cluster_df["v_call"].mode()[0] if len(cluster_df) > 0 else "NA",
                "j_gene": cluster_df["j_call"].mode()[0] if len(cluster_df) > 0 else "NA",
                "mean_mutations": cluster_df["mutations"].mean() if "mutations" in cluster_df.columns else 0,
            })

        stats_df = pd.DataFrame(clone_stats).sort_values("size", ascending=False)
        stats_file = args.output / "05_clone_stats.tsv"
        stats_df.to_csv(stats_file, sep="\t", index=False)
        logger.info(f"  Saved: {stats_file}")

        # Step 6: Build trees for large clones
        logger.info("\n[STEP 6] Building lineage trees...")
        trees_dir = args.output / "trees"
        trees_dir.mkdir(exist_ok=True)

        large_clones = stats_df[stats_df["size"] >= args.min_clone_size]
        logger.info(f"  Building trees for {len(large_clones)} clones")

        trees_built = 0
        for _, clone_info in large_clones.iterrows():
            clone_id = clone_info["clone_id"]
            clone_seqs = clustered[clustered["cluster"] == clone_id]

            if len(clone_seqs) >= 2:
                try:
                    tree_file = trees_dir / f"clone_{clone_id}.nwk"
                    build_simple_tree(clone_seqs, tree_file)
                    trees_built += 1
                except Exception as e:
                    logger.warning(f"    Failed to build tree for clone {clone_id}: {e}")

        logger.info(f"  Built {trees_built} trees")

        # Final Summary
        logger.info("\n" + "=" * 80)
        logger.info("PIPELINE COMPLETE")
        logger.info("=" * 80)
        logger.info(f"Total sequences:      {len(airr_table)}")
        logger.info(f"Productive:           {len(productive)}")
        logger.info(f"Clonal groups:        {len(clone_sizes)}")
        logger.info(f"Trees built:          {trees_built}")
        logger.info(f"\nResults in: {args.output}")
        logger.info("=" * 80)

    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        sys.exit(1)


def build_simple_tree(clone_df, output_file):
    """Build a simple NJ tree from clone sequences"""
    sequences = []

    for idx, row in clone_df.iterrows():
        seq = row.get("sequence_alignment", row.get("sequence", ""))
        seq_id = row.get("sequence_id", f"seq_{idx}")
        sequences.append(SeqRecord(Seq(seq), id=str(seq_id)))

    if len(sequences) < 2:
        return

    # Pad to same length
    max_len = max(len(s.seq) for s in sequences)
    for seq_rec in sequences:
        if len(seq_rec.seq) < max_len:
            seq_rec.seq = seq_rec.seq + Seq("-" * (max_len - len(seq_rec.seq)))

    # Build tree
    alignment = MultipleSeqAlignment(sequences)
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # Save
    Phylo.write(tree, output_file, "newick")


if __name__ == "__main__":
    main()
