"""
Change-O Wrapper Module

Provides Python wrappers for Change-O command-line tools for
VDJ alignment, clonal clustering, and germline reconstruction.
"""

import logging
import subprocess
from pathlib import Path
from typing import Optional, List, Dict, Tuple
import pandas as pd

from sadie.immcantation.config import ChangeoConfig

logger = logging.getLogger(__name__)


class ChangeoWrapper:
    """
    Wrapper class for Change-O toolkit commands.

    Change-O is a suite of tools for B-cell and T-cell receptor repertoire analysis,
    including VDJ assignment, clonal clustering, and germline reconstruction.
    """

    def __init__(self, config: Optional[ChangeoConfig] = None):
        """
        Initialize Change-O wrapper.

        Parameters
        ----------
        config : ChangeoConfig, optional
            Configuration object for Change-O parameters
        """
        self.config = config or ChangeoConfig()
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

    def make_db(
        self,
        input_file: Path,
        output_file: Path,
        format: str = "fasta",
    ) -> Path:
        """
        Convert FASTA/FASTQ to Change-O database format.

        Uses MakeDb.py fasta

        Parameters
        ----------
        input_file : Path
            Input FASTA or FASTQ file
        output_file : Path
            Output TSV file
        format : str
            Input format ('fasta' or 'fastq')

        Returns
        -------
        Path
            Output database file
        """
        cmd = [
            "MakeDb.py",
            format,
            "-s",
            str(input_file),
            "-o",
            str(output_file),
            "--outname",
            str(output_file.stem),
        ]

        self._run_command(cmd)
        self.logger.info(f"Database creation complete: {output_file}")
        return output_file

    def assign_genes(
        self,
        input_file: Path,
        output_file: Path,
        blast_output: Optional[Path] = None,
        igdata: Optional[Path] = None,
        organism: Optional[str] = None,
        loci: Optional[str] = None,
    ) -> Path:
        """
        Assign V(D)J genes using IgBLAST output.

        Uses AssignGenes.py igblast

        Parameters
        ----------
        input_file : Path
            Input FASTA file
        output_file : Path
            Output AIRR TSV file
        blast_output : Path, optional
            Pre-computed IgBLAST output (fmt 7)
        igdata : Path, optional
            IgBLAST database directory
        organism : str, optional
            Organism name
        loci : str, optional
            Loci (ig or tr)

        Returns
        -------
        Path
            Output AIRR database file
        """
        organism = organism or self.config.organism
        loci = loci or self.config.loci
        igdata = igdata or self.config.igdata

        if blast_output:
            # Use pre-computed IgBLAST output
            cmd = [
                "AssignGenes.py",
                "igblast",
                "-s",
                str(input_file),
                "-b",
                str(blast_output),
                "--organism",
                organism,
                "--loci",
                loci,
                "--format",
                "blast",
                "--outname",
                str(output_file.stem),
                "--outdir",
                str(output_file.parent),
            ]
        else:
            # Run IgBLAST
            cmd = [
                "AssignGenes.py",
                "igblast",
                "-s",
                str(input_file),
                "--organism",
                organism,
                "--loci",
                loci,
                "--format",
                "airr",
                "--outname",
                str(output_file.stem),
                "--outdir",
                str(output_file.parent),
                "--nproc",
                str(self.config.nproc),
            ]

        if igdata:
            cmd.extend(["--igdata", str(igdata)])

        self._run_command(cmd)
        self.logger.info(f"Gene assignment complete: {output_file}")
        return output_file

    def parse_db(
        self,
        input_file: Path,
        output_file: Path,
        select: Optional[List[str]] = None,
        filter: Optional[str] = None,
    ) -> Path:
        """
        Filter and select records from Change-O database.

        Uses ParseDb.py select/split

        Parameters
        ----------
        input_file : Path
            Input database file
        output_file : Path
            Output database file
        select : List[str], optional
            Fields to select
        filter : str, optional
            Filter expression

        Returns
        -------
        Path
            Output database file
        """
        cmd = [
            "ParseDb.py",
            "select",
            "-d",
            str(input_file),
            "-o",
            str(output_file),
            "--outname",
            str(output_file.stem),
        ]

        if select:
            cmd.extend(["-f"] + select)

        if filter:
            cmd.extend(["--regex", filter])

        self._run_command(cmd)
        self.logger.info(f"Database parsing complete: {output_file}")
        return output_file

    def create_germlines(
        self,
        input_file: Path,
        output_file: Path,
        germline_dir: Optional[Path] = None,
        format: str = "airr",
        v_field: Optional[str] = None,
        d_field: Optional[str] = None,
        j_field: Optional[str] = None,
        cloned: bool = False,
    ) -> Path:
        """
        Create germline sequences for database records.

        Uses CreateGermlines.py

        Parameters
        ----------
        input_file : Path
            Input AIRR database
        output_file : Path
            Output database with germlines
        germline_dir : Path, optional
            Directory containing germline FASTA files
        format : str
            Database format (airr or changeo)
        v_field : str, optional
            V gene field name
        d_field : str, optional
            D gene field name
        j_field : str, optional
            J gene field name
        cloned : bool
            Whether sequences are clonally grouped

        Returns
        -------
        Path
            Output database file
        """
        germline_dir = germline_dir or self.config.germline_dir
        v_field = v_field or self.config.v_field
        d_field = d_field or self.config.d_field
        j_field = j_field or self.config.j_field

        cmd = [
            "CreateGermlines.py",
            "-d",
            str(input_file),
            "-g",
            "dmask",
            "-r",
            str(germline_dir / "imgt_human_IGHV.fasta"),
            str(germline_dir / "imgt_human_IGHD.fasta"),
            str(germline_dir / "imgt_human_IGHJ.fasta"),
            "--format",
            format,
            "--vf",
            v_field,
            "--df",
            d_field,
            "--jf",
            j_field,
            "--outname",
            str(output_file.stem),
            "--outdir",
            str(output_file.parent),
        ]

        if cloned:
            cmd.append("--cloned")

        if germline_dir and germline_dir.exists():
            self._run_command(cmd)
            self.logger.info(f"Germline creation complete: {output_file}")
        else:
            self.logger.warning(f"Germline directory not found: {germline_dir}")
            self.logger.info("Skipping germline creation")

        return output_file

    def define_clones(
        self,
        input_file: Path,
        output_file: Path,
        mode: str = "gene",
        distance: Optional[float] = None,
        metric: Optional[str] = None,
        linkage: Optional[str] = None,
        vf: Optional[str] = None,
        jf: Optional[str] = None,
    ) -> Path:
        """
        Assign sequences to clonal groups.

        Uses DefineClones.py

        Parameters
        ----------
        input_file : Path
            Input AIRR database
        output_file : Path
            Output database with clone assignments
        mode : str
            Cloning mode ('gene', 'allele', 'aa')
        distance : float, optional
            Distance threshold for clustering
        metric : str, optional
            Distance metric (ham, aa, etc.)
        linkage : str, optional
            Linkage method (single, average, complete)
        vf : str, optional
            V gene field
        jf : str, optional
            J gene field

        Returns
        -------
        Path
            Output database file with clone_id
        """
        distance = distance or self.config.distance_threshold
        metric = metric or self.config.distance_metric
        linkage = linkage or self.config.linkage_method
        vf = vf or self.config.v_field
        jf = jf or self.config.j_field

        cmd = [
            "DefineClones.py",
            "-d",
            str(input_file),
            "--mode",
            mode,
            "--act",
            "set",
            "--model",
            metric,
            "--dist",
            str(distance),
            "--link",
            linkage,
            "--vf",
            vf,
            "--jf",
            jf,
            "--outname",
            str(output_file.stem),
            "--outdir",
            str(output_file.parent),
            "--nproc",
            str(self.config.nproc),
        ]

        self._run_command(cmd)
        self.logger.info(f"Clone definition complete: {output_file}")
        return output_file

    def build_trees(
        self,
        input_file: Path,
        output_dir: Path,
        germline_field: str = "germline_alignment",
        clone_field: str = "clone_id",
        nproc: Optional[int] = None,
    ) -> Path:
        """
        Build phylogenetic trees for clonal lineages.

        Uses BuildTrees.py

        Parameters
        ----------
        input_file : Path
            Input AIRR database with clones
        output_dir : Path
            Output directory for trees
        germline_field : str
            Field containing germline sequence
        clone_field : str
            Field containing clone ID
        nproc : int, optional
            Number of processes

        Returns
        -------
        Path
            Output directory
        """
        nproc = nproc or self.config.nproc

        output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            "BuildTrees.py",
            "-d",
            str(input_file),
            "--outdir",
            str(output_dir),
            "--outname",
            "trees",
            "--gf",
            germline_field,
            "--cf",
            clone_field,
            "--nproc",
            str(nproc),
        ]

        try:
            self._run_command(cmd)
            self.logger.info(f"Tree building complete: {output_dir}")
        except RuntimeError as e:
            self.logger.warning(f"BuildTrees.py may not be available: {e}")
            self.logger.info("Continuing without tree building")

        return output_dir

    def load_database(self, file_path: Path) -> pd.DataFrame:
        """
        Load Change-O database into pandas DataFrame.

        Parameters
        ----------
        file_path : Path
            Path to TSV database file

        Returns
        -------
        pd.DataFrame
            Database content
        """
        try:
            df = pd.read_csv(file_path, sep="\t")
            self.logger.info(f"Loaded database with {len(df)} records")
            return df
        except Exception as e:
            self.logger.error(f"Failed to load database: {e}")
            raise

    def filter_productive(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Filter database for productive sequences.

        Parameters
        ----------
        df : pd.DataFrame
            Input database

        Returns
        -------
        pd.DataFrame
            Filtered database
        """
        if "productive" in df.columns:
            productive_df = df[df["productive"] == "T"].copy()
            self.logger.info(
                f"Filtered {len(df)} -> {len(productive_df)} productive sequences"
            )
            return productive_df
        else:
            self.logger.warning("'productive' column not found, returning all records")
            return df

    def get_clone_statistics(self, df: pd.DataFrame, clone_field: str = "clone_id") -> pd.DataFrame:
        """
        Calculate statistics for each clone.

        Parameters
        ----------
        df : pd.DataFrame
            Input database with clone assignments
        clone_field : str
            Column name for clone ID

        Returns
        -------
        pd.DataFrame
            Clone statistics
        """
        if clone_field not in df.columns:
            self.logger.warning(f"Clone field '{clone_field}' not found")
            return pd.DataFrame()

        clone_stats = (
            df.groupby(clone_field)
            .agg(
                {
                    "sequence_id": "count",
                    "v_call": lambda x: x.mode()[0] if len(x.mode()) > 0 else None,
                    "j_call": lambda x: x.mode()[0] if len(x.mode()) > 0 else None,
                    "junction_aa_length": "mean",
                }
            )
            .rename(columns={"sequence_id": "clone_size"})
            .reset_index()
        )

        self.logger.info(f"Calculated statistics for {len(clone_stats)} clones")
        return clone_stats
