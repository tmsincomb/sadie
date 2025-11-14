"""
pRESTO Wrapper Module

Provides Python wrappers for pRESTO command-line tools.
"""

import logging
import subprocess
from pathlib import Path
from typing import Optional, List, Dict, Tuple
import pandas as pd

from sadie.immcantation.config import PrestoConfig

logger = logging.getLogger(__name__)


class PrestoWrapper:
    """
    Wrapper class for pRESTO toolkit commands.

    pRESTO (The REpertoire Sequencing TOolkit) processes high-throughput
    sequencing data for adaptive immune receptor repertoires.
    """

    def __init__(self, config: Optional[PrestoConfig] = None):
        """
        Initialize pRESTO wrapper.

        Parameters
        ----------
        config : PrestoConfig, optional
            Configuration object for pRESTO parameters
        """
        self.config = config or PrestoConfig()
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

    def filter_quality(
        self,
        input_file: Path,
        output_file: Path,
        min_quality: Optional[int] = None,
        min_length: Optional[int] = None,
    ) -> Path:
        """
        Filter sequences by quality score.

        Uses FilterSeq.py quality

        Parameters
        ----------
        input_file : Path
            Input FASTQ file
        output_file : Path
            Output FASTQ file
        min_quality : int, optional
            Minimum quality score
        min_length : int, optional
            Minimum sequence length

        Returns
        -------
        Path
            Output file path
        """
        min_quality = min_quality or self.config.min_quality
        min_length = min_length or self.config.min_length

        cmd = [
            "FilterSeq.py",
            "quality",
            "-s",
            str(input_file),
            "-q",
            str(min_quality),
            "--outname",
            str(output_file.stem),
            "--outdir",
            str(output_file.parent),
            "--nproc",
            str(self.config.nproc),
        ]

        self._run_command(cmd)
        self.logger.info(f"Quality filtering complete: {output_file}")
        return output_file

    def filter_length(
        self,
        input_file: Path,
        output_file: Path,
        min_length: Optional[int] = None,
        max_length: Optional[int] = None,
    ) -> Path:
        """
        Filter sequences by length.

        Uses FilterSeq.py length

        Parameters
        ----------
        input_file : Path
            Input FASTA/FASTQ file
        output_file : Path
            Output file
        min_length : int, optional
            Minimum sequence length
        max_length : int, optional
            Maximum sequence length

        Returns
        -------
        Path
            Output file path
        """
        min_length = min_length or self.config.min_length

        cmd = [
            "FilterSeq.py",
            "length",
            "-s",
            str(input_file),
            "-n",
            str(min_length),
            "--outname",
            str(output_file.stem),
            "--outdir",
            str(output_file.parent),
            "--nproc",
            str(self.config.nproc),
        ]

        if max_length:
            cmd.extend(["-m", str(max_length)])

        self._run_command(cmd)
        self.logger.info(f"Length filtering complete: {output_file}")
        return output_file

    def mask_primers(
        self,
        input_file: Path,
        output_file: Path,
        primer_file: Path,
        mode: str = "mask",
        max_error: Optional[float] = None,
    ) -> Path:
        """
        Identify and mask/cut primers from sequences.

        Uses MaskPrimers.py align

        Parameters
        ----------
        input_file : Path
            Input FASTA file
        output_file : Path
            Output file
        primer_file : Path
            FASTA file containing primer sequences
        mode : str
            Either 'mask' or 'cut'
        max_error : float, optional
            Maximum error rate for primer matching

        Returns
        -------
        Path
            Output file path
        """
        max_error = max_error or self.config.max_primer_error

        cmd = [
            "MaskPrimers.py",
            "align",
            "-s",
            str(input_file),
            "-p",
            str(primer_file),
            "--mode",
            mode,
            "--maxerror",
            str(max_error),
            "--outname",
            str(output_file.stem),
            "--outdir",
            str(output_file.parent),
            "--nproc",
            str(self.config.nproc),
            "--log",
            str(output_file.parent / f"{output_file.stem}.log"),
        ]

        self._run_command(cmd)
        self.logger.info(f"Primer masking complete: {output_file}")
        return output_file

    def pair_sequences(
        self,
        forward_file: Path,
        reverse_file: Path,
        output_file: Path,
        coord_type: str = "illumina",
    ) -> Path:
        """
        Match forward and reverse reads by sequence ID.

        Uses PairSeq.py

        Parameters
        ----------
        forward_file : Path
            Forward read file
        reverse_file : Path
            Reverse read file
        output_file : Path
            Output file base name
        coord_type : str
            Coordinate format (illumina, sra, solexa, 454, presto)

        Returns
        -------
        Path
            Output file path
        """
        cmd = [
            "PairSeq.py",
            "-1",
            str(forward_file),
            "-2",
            str(reverse_file),
            "--coord",
            coord_type,
            "--outname",
            str(output_file.stem),
            "--outdir",
            str(output_file.parent),
        ]

        self._run_command(cmd)
        self.logger.info(f"Sequence pairing complete: {output_file}")
        return output_file

    def build_consensus(
        self,
        input_file: Path,
        output_file: Path,
        barcode_field: str = "BARCODE",
        min_freq: Optional[float] = None,
        min_count: Optional[int] = None,
    ) -> Path:
        """
        Build consensus sequences from UMI groups.

        Uses BuildConsensus.py

        Parameters
        ----------
        input_file : Path
            Input file with UMI annotations
        output_file : Path
            Output file
        barcode_field : str
            Header field containing barcode/UMI
        min_freq : float, optional
            Minimum frequency for consensus
        min_count : int, optional
            Minimum count for consensus

        Returns
        -------
        Path
            Output file path
        """
        min_freq = min_freq or self.config.min_consensus_freq
        min_count = min_count or self.config.min_consensus_count

        cmd = [
            "BuildConsensus.py",
            "-s",
            str(input_file),
            "-n",
            str(min_count),
            "-q",
            str(min_freq),
            "--bf",
            barcode_field,
            "--outname",
            str(output_file.stem),
            "--outdir",
            str(output_file.parent),
            "--nproc",
            str(self.config.nproc),
            "--log",
            str(output_file.parent / f"{output_file.stem}.log"),
        ]

        self._run_command(cmd)
        self.logger.info(f"Consensus building complete: {output_file}")
        return output_file

    def assemble_pairs(
        self,
        input_file: Path,
        output_file: Path,
        min_overlap: Optional[int] = None,
        max_error: Optional[float] = None,
        alpha: Optional[float] = None,
    ) -> Path:
        """
        Assemble paired-end reads.

        Uses AssemblePairs.py align

        Parameters
        ----------
        input_file : Path
            Input file with paired reads
        output_file : Path
            Output file
        min_overlap : int, optional
            Minimum overlap length
        max_error : float, optional
            Maximum error rate
        alpha : float, optional
            Alpha parameter for assembly

        Returns
        -------
        Path
            Output file path
        """
        min_overlap = min_overlap or self.config.min_overlap
        max_error = max_error or self.config.max_assembly_error
        alpha = alpha or self.config.alpha

        cmd = [
            "AssemblePairs.py",
            "align",
            "-1",
            str(input_file),
            "--minlen",
            str(min_overlap),
            "--maxerror",
            str(max_error),
            "--alpha",
            str(alpha),
            "--outname",
            str(output_file.stem),
            "--outdir",
            str(output_file.parent),
            "--nproc",
            str(self.config.nproc),
            "--log",
            str(output_file.parent / f"{output_file.stem}.log"),
        ]

        self._run_command(cmd)
        self.logger.info(f"Pair assembly complete: {output_file}")
        return output_file

    def collapse_duplicates(
        self,
        input_file: Path,
        output_file: Path,
        max_distance: Optional[int] = None,
        inner: bool = False,
    ) -> Path:
        """
        Remove duplicate sequences.

        Uses CollapseSeq.py

        Parameters
        ----------
        input_file : Path
            Input file
        output_file : Path
            Output file
        max_distance : int, optional
            Maximum distance for duplicates
        inner : bool
            Use inner field for grouping

        Returns
        -------
        Path
            Output file path
        """
        max_distance = max_distance or self.config.max_duplicate_distance

        cmd = [
            "CollapseSeq.py",
            "-s",
            str(input_file),
            "-n",
            str(max_distance),
            "--outname",
            str(output_file.stem),
            "--outdir",
            str(output_file.parent),
        ]

        if inner:
            cmd.append("--inner")

        self._run_command(cmd)
        self.logger.info(f"Duplicate collapse complete: {output_file}")
        return output_file

    def parse_log(self, log_file: Path) -> pd.DataFrame:
        """
        Parse pRESTO log file into DataFrame.

        Uses ParseLog.py

        Parameters
        ----------
        log_file : Path
            pRESTO log file

        Returns
        -------
        pd.DataFrame
            Parsed log data
        """
        cmd = [
            "ParseLog.py",
            "-l",
            str(log_file),
            "-f",
            "ID",
            "PRIMER",
            "ERROR",
        ]

        returncode, stdout, stderr = self._run_command(cmd)

        # Read the output tab file
        output_file = log_file.parent / f"{log_file.stem}_table.tab"
        if output_file.exists():
            return pd.read_csv(output_file, sep="\t")
        else:
            self.logger.warning(f"Could not find parsed log: {output_file}")
            return pd.DataFrame()
