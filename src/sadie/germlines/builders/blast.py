"""
BLAST Database Builder
======================

Builds BLAST databases for IgBLAST from germline FASTA files.

IgBLAST requires:
- Separate BLAST databases for V, D, J segments
- Standard BLAST database format (created by makeblastdb)
- Specific naming convention: <species>_<segment>

Design Principles:
- Use subprocess for makeblastdb (explicit over implicit)
- Validate inputs before building
- Provide clear error messages
"""

import logging
import subprocess
from pathlib import Path
from typing import Optional
from Bio import SeqIO

logger = logging.getLogger(__name__)


# Constants
SEGMENTS = ["V", "D", "J"]
CHAINS = ["H", "K", "L"]


class BlastDBBuilder:
    """
    Build BLAST databases for IgBLAST.

    Creates separate BLAST databases for each segment type (V, D, J)
    from normalized ungapped FASTA files.

    Examples
    --------
    >>> builder = BlastDBBuilder()
    >>> builder.build_for_species(
    ...     "human",
    ...     source_dir=Path("normalized/human/ungapped"),
    ...     output_dir=Path("igblast/database/human")
    ... )
    """

    def build_for_species(
        self,
        species: str,
        source_dir: Path,
        output_dir: Path
    ) -> None:
        """
        Build all BLAST databases for a species.

        Parameters
        ----------
        species : str
            Species name
        source_dir : Path
            Directory containing ungapped FASTA files (normalized/)
        output_dir : Path
            Output directory for BLAST databases
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Building BLAST databases for {species}")

        # Combine all chains for each segment
        for segment in SEGMENTS:
            self._build_segment_database(
                species,
                segment,
                source_dir,
                output_dir
            )

    def _build_segment_database(
        self,
        species: str,
        segment: str,
        source_dir: Path,
        output_dir: Path
    ) -> None:
        """
        Build BLAST database for single segment.

        Combines all chains (H, K, L) for this segment.

        Parameters
        ----------
        species : str
            Species name
        segment : str
            Segment type (V, D, or J)
        source_dir : Path
            Source directory with FASTA files
        output_dir : Path
            Output directory
        """
        # Collect sequences from all chains
        combined_sequences = []

        for chain in CHAINS:
            fasta_path = source_dir / f"IG{chain}{segment}.fasta"

            if not fasta_path.exists():
                logger.debug(f"No file: {fasta_path}")
                continue

            try:
                records = list(SeqIO.parse(fasta_path, "fasta"))
                combined_sequences.extend(records)
                logger.info(
                    f"Added {len(records)} sequences from {fasta_path.name}"
                )
            except Exception as e:
                logger.error(f"Failed to read {fasta_path}: {e}")

        # Guard: no sequences found
        if not combined_sequences:
            logger.warning(f"No sequences for {species} {segment}")
            return

        # Write combined FASTA
        combined_fasta = output_dir / f"{species}_{segment}.fasta"
        SeqIO.write(combined_sequences, combined_fasta, "fasta")
        logger.info(
            f"Wrote {len(combined_sequences)} sequences to {combined_fasta}"
        )

        # Build BLAST database
        self._run_makeblastdb(combined_fasta, segment)

    def _run_makeblastdb(self, fasta_path: Path, segment: str) -> None:
        """
        Run makeblastdb to create BLAST database.

        Parameters
        ----------
        fasta_path : Path
            Input FASTA file
        segment : str
            Segment type (for logging)
        """
        db_name = fasta_path.stem

        try:
            cmd = [
                "makeblastdb",
                "-dbtype", "nucl",
                "-in", str(fasta_path),
                "-out", str(fasta_path.with_suffix("")),
                "-parse_seqids",
                "-hash_index",
            ]

            logger.debug(f"Running: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            logger.info(f"Built BLAST database: {db_name}")
            logger.debug(result.stdout)

        except subprocess.CalledProcessError as e:
            logger.error(
                f"makeblastdb failed for {db_name}: {e.stderr}"
            )
            raise

        except FileNotFoundError:
            logger.error(
                "makeblastdb not found. "
                "Ensure BLAST+ is installed and in PATH."
            )
            raise

    def validate_database(self, db_path: Path) -> bool:
        """
        Validate BLAST database was created successfully.

        Parameters
        ----------
        db_path : Path
            Path to BLAST database (without extension)

        Returns
        -------
        bool
            True if valid
        """
        # Check for required BLAST database files
        required_extensions = [".nhr", ".nin", ".nsq"]

        for ext in required_extensions:
            file_path = db_path.with_suffix(ext)
            if not file_path.exists():
                logger.error(f"Missing BLAST file: {file_path}")
                return False

        return True
