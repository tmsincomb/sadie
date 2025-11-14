"""
Immcantation Pipeline Orchestrator

Main pipeline class that coordinates pRESTO, Change-O, and lineage analysis
for comprehensive VDJ antibody repertoire analysis.
"""

import logging
import shutil
from pathlib import Path
from typing import Optional, Dict, List
from datetime import datetime
import pandas as pd

from sadie.immcantation.config import PipelineConfig
from sadie.immcantation.presto_wrapper import PrestoWrapper
from sadie.immcantation.changeo_wrapper import ChangeoWrapper
from sadie.immcantation.lineage import LineageAnalyzer

logger = logging.getLogger(__name__)


class ImmcantationPipeline:
    """
    Complete Immcantation analysis pipeline.

    Orchestrates VDJ antibody analysis from raw reads to lineage trees:
    1. pRESTO: Quality control and preprocessing
    2. Change-O: VDJ assignment and clonal clustering
    3. Lineage Analysis: Tree building and analysis
    """

    def __init__(self, config: Optional[PipelineConfig] = None):
        """
        Initialize pipeline.

        Parameters
        ----------
        config : PipelineConfig, optional
            Pipeline configuration
        """
        self.config = config or PipelineConfig()
        self.presto = PrestoWrapper(self.config.presto)
        self.changeo = ChangeoWrapper(self.config.changeo)
        self.lineage = LineageAnalyzer(self.config.lineage)

        # Create output directory first
        self.config.output_dir.mkdir(parents=True, exist_ok=True)

        # Setup logging
        self._setup_logging()

        # Track intermediate files
        self.intermediate_files = []

        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.info("Immcantation Pipeline initialized")

    def _setup_logging(self):
        """Setup logging configuration"""
        log_level = logging.DEBUG if self.config.verbose else logging.INFO

        # Configure root logger
        logging.basicConfig(
            level=log_level,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

        # File logging if specified
        if self.config.log_file:
            file_handler = logging.FileHandler(self.config.log_file)
            file_handler.setLevel(log_level)
            formatter = logging.Formatter(
                "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            )
            file_handler.setFormatter(formatter)
            logging.getLogger().addHandler(file_handler)

    def run(
        self,
        input_file: Path,
        sample_name: Optional[str] = None,
    ) -> Dict[str, Path]:
        """
        Run complete pipeline.

        Parameters
        ----------
        input_file : Path
            Input FASTA/FASTQ file
        sample_name : str, optional
            Sample name for outputs

        Returns
        -------
        Dict[str, Path]
            Dictionary of output files
        """
        start_time = datetime.now()
        self.logger.info("=" * 80)
        self.logger.info("IMMCANTATION PIPELINE START")
        self.logger.info("=" * 80)

        if not sample_name:
            sample_name = input_file.stem

        outputs = {}

        try:
            # Step 1: pRESTO preprocessing (optional)
            if self.config.run_presto:
                self.logger.info("\n[STEP 1] Running pRESTO preprocessing...")
                processed_file = self._run_presto(input_file, sample_name)
                outputs["presto_output"] = processed_file
            else:
                self.logger.info("\n[STEP 1] Skipping pRESTO (using pre-processed data)")
                processed_file = input_file

            # Step 2: Change-O VDJ assignment and clonal clustering
            if self.config.run_changeo:
                self.logger.info("\n[STEP 2] Running Change-O VDJ assignment...")
                changeo_outputs = self._run_changeo(processed_file, sample_name)
                outputs.update(changeo_outputs)
            else:
                self.logger.info("\n[STEP 2] Skipping Change-O")
                changeo_outputs = {}

            # Step 3: Lineage analysis
            if self.config.run_lineage and "clone_db" in changeo_outputs:
                self.logger.info("\n[STEP 3] Running lineage analysis...")
                lineage_outputs = self._run_lineage(
                    changeo_outputs["clone_db"], sample_name
                )
                outputs.update(lineage_outputs)
            else:
                self.logger.info("\n[STEP 3] Skipping lineage analysis")

            # Cleanup intermediate files
            if not self.config.keep_intermediate:
                self._cleanup()

            # Generate final report
            report_file = self._generate_report(outputs, sample_name)
            outputs["report"] = report_file

        except Exception as e:
            self.logger.error(f"Pipeline failed: {e}", exc_info=True)
            raise

        finally:
            end_time = datetime.now()
            duration = end_time - start_time
            self.logger.info("=" * 80)
            self.logger.info(f"PIPELINE COMPLETE - Duration: {duration}")
            self.logger.info("=" * 80)

        return outputs

    def _run_presto(self, input_file: Path, sample_name: str) -> Path:
        """
        Run pRESTO preprocessing steps.

        Parameters
        ----------
        input_file : Path
            Input file
        sample_name : str
            Sample name

        Returns
        -------
        Path
            Processed FASTA file
        """
        presto_dir = self.config.output_dir / "presto"
        presto_dir.mkdir(exist_ok=True)

        # For FASTA input, minimal processing needed
        if input_file.suffix in [".fasta", ".fa"]:
            self.logger.info("Input is FASTA format - applying basic filtering")

            # Length filtering
            filtered_file = presto_dir / f"{sample_name}_filtered.fasta"
            self.presto.filter_length(input_file, filtered_file)
            self.intermediate_files.append(filtered_file)

            return filtered_file

        # For FASTQ, run full QC
        elif input_file.suffix in [".fastq", ".fq"]:
            self.logger.info("Input is FASTQ format - running quality control")

            # Quality filtering
            qc_file = presto_dir / f"{sample_name}_qc.fastq"
            self.presto.filter_quality(input_file, qc_file)
            self.intermediate_files.append(qc_file)

            # Length filtering
            filtered_file = presto_dir / f"{sample_name}_filtered.fasta"
            self.presto.filter_length(qc_file, filtered_file)
            self.intermediate_files.append(filtered_file)

            return filtered_file

        else:
            self.logger.warning(f"Unknown format: {input_file.suffix}, using as-is")
            return input_file

    def _run_changeo(self, input_file: Path, sample_name: str) -> Dict[str, Path]:
        """
        Run Change-O VDJ assignment and clonal clustering.

        Parameters
        ----------
        input_file : Path
            Input FASTA file
        sample_name : str
            Sample name

        Returns
        -------
        Dict[str, Path]
            Output files from Change-O
        """
        changeo_dir = self.config.output_dir / "changeo"
        changeo_dir.mkdir(exist_ok=True)

        outputs = {}

        # Step 2.1: Assign genes using IgBLAST
        self.logger.info("  [2.1] Assigning V(D)J genes with IgBLAST...")
        airr_db = changeo_dir / f"{sample_name}_db-pass.tsv"

        try:
            self.changeo.assign_genes(
                input_file,
                airr_db,
                organism=self.config.changeo.organism,
                loci=self.config.changeo.loci,
            )
            outputs["airr_db"] = airr_db
        except Exception as e:
            self.logger.error(f"Gene assignment failed: {e}")
            # Try with sadie's IgBLAST as fallback
            self.logger.info("  Attempting fallback with sadie's IgBLAST...")
            airr_db = self._run_sadie_igblast(input_file, changeo_dir, sample_name)
            outputs["airr_db"] = airr_db

        # Load database
        db = self.changeo.load_database(airr_db)
        self.logger.info(f"  Loaded {len(db)} sequences")

        # Step 2.2: Filter for productive sequences
        self.logger.info("  [2.2] Filtering productive sequences...")
        productive_db = self.changeo.filter_productive(db)

        productive_file = changeo_dir / f"{sample_name}_productive.tsv"
        productive_db.to_csv(productive_file, sep="\t", index=False)
        outputs["productive_db"] = productive_file

        # Step 2.3: Define clones
        self.logger.info("  [2.3] Defining clonal groups...")
        clone_file = changeo_dir / f"{sample_name}_clone-pass.tsv"

        try:
            self.changeo.define_clones(
                productive_file,
                clone_file,
                distance=self.config.changeo.distance_threshold,
                metric=self.config.changeo.distance_metric,
            )
            outputs["clone_db"] = clone_file
        except Exception as e:
            self.logger.warning(f"DefineClones failed: {e}")
            self.logger.info("  Using productive database as clone database")
            outputs["clone_db"] = productive_file

        # Get clone statistics
        clone_db = self.changeo.load_database(outputs["clone_db"])
        clone_stats = self.changeo.get_clone_statistics(clone_db)

        if not clone_stats.empty:
            stats_file = changeo_dir / f"{sample_name}_clone_stats.tsv"
            clone_stats.to_csv(stats_file, sep="\t", index=False)
            outputs["clone_stats"] = stats_file

            self.logger.info(f"  Found {len(clone_stats)} clonal groups")
            self.logger.info(
                f"  Clone size range: {clone_stats['clone_size'].min()} - {clone_stats['clone_size'].max()}"
            )

        return outputs

    def _run_sadie_igblast(
        self, input_file: Path, output_dir: Path, sample_name: str
    ) -> Path:
        """
        Fallback: Run IgBLAST using sadie's Airr module.

        Parameters
        ----------
        input_file : Path
            Input FASTA file
        output_dir : Path
            Output directory
        sample_name : str
            Sample name

        Returns
        -------
        Path
            AIRR database file
        """
        from sadie.airr import Airr

        self.logger.info("  Running sadie IgBLAST...")

        # Run sadie AIRR annotation
        airr_api = Airr(
            reference_name=self.config.changeo.organism,
            num_cpus=self.config.changeo.nproc,
        )

        airr_table = airr_api.run_fasta(str(input_file))

        # Save as AIRR format
        output_file = output_dir / f"{sample_name}_sadie_db.tsv"
        airr_table.to_airr(str(output_file))

        self.logger.info(f"  Sadie IgBLAST complete: {len(airr_table)} sequences")

        return output_file

    def _run_lineage(self, clone_file: Path, sample_name: str) -> Dict[str, Path]:
        """
        Run lineage tree analysis.

        Parameters
        ----------
        clone_file : Path
            Clone database file
        sample_name : str
            Sample name

        Returns
        -------
        Dict[str, Path]
            Lineage analysis outputs
        """
        lineage_dir = self.config.output_dir / "lineage"
        lineage_dir.mkdir(exist_ok=True)

        outputs = {}

        # Load clone database
        clone_db = pd.read_csv(clone_file, sep="\t")

        # Step 3.1: Filter clones by size
        self.logger.info("  [3.1] Filtering clones by size...")
        filtered_db = self.lineage.filter_clones_by_size(clone_db)

        filtered_file = lineage_dir / f"{sample_name}_filtered_clones.tsv"
        filtered_db.to_csv(filtered_file, sep="\t", index=False)
        outputs["filtered_clones"] = filtered_file

        # Step 3.2: Split into individual clone files
        self.logger.info("  [3.2] Splitting clones...")
        clone_files_dir = lineage_dir / "clone_files"
        clone_files = self.lineage.split_clones(filtered_db, clone_files_dir)

        self.logger.info(f"  Split into {len(clone_files)} clone files")

        # Step 3.3: Build trees
        self.logger.info("  [3.3] Building phylogenetic trees...")
        trees_dir = lineage_dir / "trees"
        tree_files = self.lineage.build_trees_simple(clone_files, trees_dir)

        self.logger.info(f"  Built {len(tree_files)} trees")

        # Step 3.4: Calculate tree statistics
        if tree_files:
            self.logger.info("  [3.4] Calculating tree statistics...")
            tree_stats = self.lineage.calculate_tree_stats(tree_files)

            if not tree_stats.empty:
                stats_file = lineage_dir / f"{sample_name}_tree_stats.tsv"
                tree_stats.to_csv(stats_file, sep="\t", index=False)
                outputs["tree_stats"] = stats_file

        # Step 3.5: Generate summary
        self.logger.info("  [3.5] Generating lineage summary...")
        summary_file = lineage_dir / f"{sample_name}_lineage_summary.tsv"
        summary = self.lineage.summarize_lineages(
            filtered_db, tree_files, summary_file
        )
        outputs["lineage_summary"] = summary_file

        return outputs

    def _generate_report(
        self, outputs: Dict[str, Path], sample_name: str
    ) -> Path:
        """
        Generate final pipeline report.

        Parameters
        ----------
        outputs : Dict[str, Path]
            Pipeline output files
        sample_name : str
            Sample name

        Returns
        -------
        Path
            Report file
        """
        report_file = self.config.output_dir / f"{sample_name}_report.txt"

        with open(report_file, "w") as f:
            f.write("=" * 80 + "\n")
            f.write("IMMCANTATION PIPELINE REPORT\n")
            f.write("=" * 80 + "\n\n")

            f.write(f"Sample: {sample_name}\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

            f.write("Pipeline Configuration:\n")
            f.write(f"  - Run pRESTO: {self.config.run_presto}\n")
            f.write(f"  - Run Change-O: {self.config.run_changeo}\n")
            f.write(f"  - Run Lineage: {self.config.run_lineage}\n")
            f.write(f"  - Organism: {self.config.changeo.organism}\n")
            f.write(f"  - Distance threshold: {self.config.changeo.distance_threshold}\n\n")

            f.write("Output Files:\n")
            for key, path in outputs.items():
                if path.exists():
                    size = path.stat().st_size / 1024  # KB
                    f.write(f"  - {key}: {path.name} ({size:.1f} KB)\n")

            # Add statistics if available
            f.write("\n")
            f.write("=" * 80 + "\n")
            f.write("STATISTICS\n")
            f.write("=" * 80 + "\n\n")

            if "clone_stats" in outputs and outputs["clone_stats"].exists():
                clone_stats = pd.read_csv(outputs["clone_stats"], sep="\t")
                f.write(f"Clonal Groups: {len(clone_stats)}\n")
                f.write(f"  - Mean clone size: {clone_stats['clone_size'].mean():.1f}\n")
                f.write(f"  - Median clone size: {clone_stats['clone_size'].median():.1f}\n")
                f.write(f"  - Largest clone: {clone_stats['clone_size'].max()}\n")
                f.write(f"  - Smallest clone: {clone_stats['clone_size'].min()}\n\n")

            if "lineage_summary" in outputs and outputs["lineage_summary"].exists():
                lineage_sum = pd.read_csv(outputs["lineage_summary"], sep="\t")
                trees_built = lineage_sum["has_tree"].sum()
                f.write(f"Lineage Trees: {trees_built}\n\n")

        self.logger.info(f"Report generated: {report_file}")
        return report_file

    def _cleanup(self):
        """Remove intermediate files"""
        self.logger.info("Cleaning up intermediate files...")
        for file_path in self.intermediate_files:
            if file_path.exists():
                file_path.unlink()
                self.logger.debug(f"Removed: {file_path}")
