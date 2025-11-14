#!/usr/bin/env python3
"""
Example script for running the Immcantation pipeline on VDJ antibody data.

This script demonstrates how to use the SADIE Immcantation pipeline for
comprehensive analysis of human VDJ antibody repertoires from Illumina reads.

Usage:
    python run_immcantation_pipeline.py --input INPUT.fasta --output OUTPUT_DIR

Author: SADIE Team
Date: 2025-11-14
"""

import argparse
import logging
from pathlib import Path
import sys

# Add sadie to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from sadie.immcantation.pipeline import ImmcantationPipeline
from sadie.immcantation.config import PipelineConfig, ChangeoConfig, LineageConfig


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Run Immcantation VDJ antibody analysis pipeline"
    )

    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        required=True,
        help="Input FASTA file with VDJ sequences",
    )

    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("./immcantation_results"),
        help="Output directory (default: ./immcantation_results)",
    )

    parser.add_argument(
        "--sample-name",
        "-n",
        type=str,
        help="Sample name (default: input filename)",
    )

    parser.add_argument(
        "--organism",
        type=str,
        default="human",
        choices=["human", "mouse", "rat", "rabbit"],
        help="Organism (default: human)",
    )

    parser.add_argument(
        "--distance-threshold",
        type=float,
        default=0.15,
        help="Distance threshold for clonal clustering (default: 0.15)",
    )

    parser.add_argument(
        "--min-clone-size",
        type=int,
        default=5,
        help="Minimum clone size for lineage analysis (default: 5)",
    )

    parser.add_argument(
        "--max-clone-size",
        type=int,
        default=1000,
        help="Maximum clone size for lineage analysis (default: 1000)",
    )

    parser.add_argument(
        "--skip-presto",
        action="store_true",
        help="Skip pRESTO preprocessing (use for pre-processed data)",
    )

    parser.add_argument(
        "--skip-lineage",
        action="store_true",
        help="Skip lineage tree analysis",
    )

    parser.add_argument(
        "--nproc",
        type=int,
        default=4,
        help="Number of parallel processes (default: 4)",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Verbose logging",
    )

    args = parser.parse_args()

    # Validate input
    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)

    # Configure pipeline
    config = PipelineConfig(
        output_dir=args.output,
        run_presto=not args.skip_presto,
        run_changeo=True,
        run_lineage=not args.skip_lineage,
        verbose=args.verbose,
        log_file=args.output / "pipeline.log",
    )

    # Configure Change-O
    config.changeo = ChangeoConfig(
        organism=args.organism,
        distance_threshold=args.distance_threshold,
        nproc=args.nproc,
    )

    # Configure Lineage
    config.lineage = LineageConfig(
        min_clone_size=args.min_clone_size,
        max_clone_size=args.max_clone_size,
    )

    # Print configuration
    print("=" * 80)
    print("IMMCANTATION PIPELINE")
    print("=" * 80)
    print(f"Input:          {args.input}")
    print(f"Output:         {args.output}")
    print(f"Sample:         {args.sample_name or args.input.stem}")
    print(f"Organism:       {args.organism}")
    print(f"Processes:      {args.nproc}")
    print(f"Distance:       {args.distance_threshold}")
    print(f"Clone size:     {args.min_clone_size}-{args.max_clone_size}")
    print(f"Run pRESTO:     {config.run_presto}")
    print(f"Run lineage:    {config.run_lineage}")
    print("=" * 80)
    print()

    # Create and run pipeline
    try:
        pipeline = ImmcantationPipeline(config)

        outputs = pipeline.run(
            input_file=args.input,
            sample_name=args.sample_name,
        )

        # Print results
        print("\n" + "=" * 80)
        print("PIPELINE COMPLETED SUCCESSFULLY")
        print("=" * 80)
        print("\nOutput files:")
        for key, path in outputs.items():
            if path.exists():
                print(f"  {key:20s} : {path}")

        print(f"\nFull results in: {args.output}")

        if "report" in outputs:
            print(f"\nSee report: {outputs['report']}")

    except Exception as e:
        print(f"\n ERROR: Pipeline failed: {e}")
        if args.verbose:
            import traceback

            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
