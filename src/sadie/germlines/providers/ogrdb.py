"""
OGRDB Provider
==============

Provider for OGRDB (Open Germline Receptor Database).

OGRDB is a community-curated database providing:
- Novel alleles discovered through repertoire sequencing
- Inferred germline sequences with supporting evidence
- AIRR Community governance and review

Data Source: https://ogrdb.airr-community.org/
API Documentation: https://ogrdb.airr-community.org/api/docs

TODO: Implement OGRDB API client and data parsing
Current Status: Stub implementation - reads from pre-downloaded FASTA files
"""

import logging
from pathlib import Path
from typing import List, Optional
from datetime import datetime
from Bio import SeqIO

from .base import GermlineProvider
from ..models import GermlineGene, ProviderMetadata


logger = logging.getLogger(__name__)


class OGRDBProvider(GermlineProvider):
    """
    Provider for OGRDB germline database.

    OGRDB provides novel alleles and inferred germline sequences
    discovered through repertoire sequencing.

    Expected Directory Structure:
        sources/ogrdb/
        ├── human/
        │   ├── IGHV.fasta    # May be ungapped
        │   ├── IGHJ.fasta
        │   ├── IGKV.fasta
        │   ├── IGKJ.fasta
        │   ├── IGLV.fasta
        │   └── IGLJ.fasta
        └── mouse/
            └── ...

    Examples
    --------
    >>> provider = OGRDBProvider()
    >>> genes = provider.fetch_genes("human", "V", "H")
    >>> print(f"Found {len(genes)} OGRDB IGHV genes")
    """

    def fetch_genes(
        self,
        species: str,
        segment: str,
        chain: str
    ) -> List[GermlineGene]:
        """
        Fetch OGRDB genes from pre-downloaded FASTA files.

        TODO: Implement automatic download from OGRDB API if files missing

        Parameters
        ----------
        species : str
            Species name
        segment : str
            Segment type
        chain : str
            Chain type

        Returns
        -------
        List[GermlineGene]
            OGRDB genes
        """
        fasta_path = self.get_fasta_path(species, segment, chain)

        # Guard: file doesn't exist
        if not fasta_path.exists():
            logger.debug(f"No OGRDB file: {fasta_path}")
            logger.info(
                "Run download script or add FASTA manually. "
                "See sources/ogrdb/README.md"
            )
            return []

        logger.info(f"Loading OGRDB FASTA: {fasta_path}")

        genes = self._parse_ogrdb_fasta(fasta_path, species, segment, chain)

        logger.info(f"Loaded {len(genes)} OGRDB genes")

        return genes

    def _parse_ogrdb_fasta(
        self,
        fasta_path: Path,
        species: str,
        segment: str,
        chain: str
    ) -> List[GermlineGene]:
        """
        Parse OGRDB FASTA file.

        TODO: Parse OGRDB-specific metadata

        Parameters
        ----------
        fasta_path : Path
            Path to FASTA file
        species : str
            Species name
        segment : str
            Segment type
        chain : str
            Chain type

        Returns
        -------
        List[GermlineGene]
            Parsed genes
        """
        genes = []

        try:
            records = list(SeqIO.parse(fasta_path, "fasta"))
        except Exception as e:
            logger.error(f"Failed to parse {fasta_path}: {e}")
            return []

        for record in records:
            gene = self._create_ogrdb_gene(record, species, segment, chain)
            if gene:
                genes.append(gene)

        return genes

    def _create_ogrdb_gene(
        self,
        record,
        species: str,
        segment: str,
        chain: str
    ) -> Optional[GermlineGene]:
        """
        Create GermlineGene from OGRDB SeqRecord.

        OGRDB sequences may be ungapped or gapped.
        TODO: Parse OGRDB metadata from headers

        Parameters
        ----------
        record : SeqRecord
            BioPython sequence record
        species : str
            Species name
        segment : str
            Segment type
        chain : str
            Chain type

        Returns
        -------
        GermlineGene or None
            Gene object if successful
        """
        gene_name = record.id.split("|")[0]

        sequence = str(record.seq).upper()
        is_gapped = "." in sequence or "-" in sequence

        if is_gapped:
            sequence_gapped = sequence
            sequence_ungapped = sequence.replace(".", "").replace("-", "")
        else:
            sequence_ungapped = sequence
            sequence_gapped = None  # TODO: Gap using aligner

        try:
            gene = GermlineGene(
                name=gene_name,
                species=species,
                segment=segment,
                chain=chain,
                sequence=sequence_ungapped,
                sequence_gapped=sequence_gapped,
                is_functional=True,  # OGRDB sequences are typically functional
                functionality="F",
                source="ogrdb",
            )
            return gene

        except Exception as e:
            logger.error(f"Failed to create OGRDB gene {gene_name}: {e}")
            return None

    def fetch_gene_by_name(
        self,
        name: str,
        species: str
    ) -> Optional[GermlineGene]:
        """
        Fetch specific OGRDB gene by name.

        Parameters
        ----------
        name : str
            Gene name
        species : str
            Species name

        Returns
        -------
        GermlineGene or None
            Gene if found
        """
        for segment in ["V", "D", "J"]:
            for chain in ["H", "K", "L"]:
                genes = self.fetch_genes(species, segment, chain)
                for gene in genes:
                    if gene.name == name:
                        return gene

        return None

    def is_available(self, species: str) -> bool:
        """
        Check if OGRDB data is available for species.

        Parameters
        ----------
        species : str
            Species name

        Returns
        -------
        bool
            True if data available
        """
        species_dir = self.data_dir / species

        if not species_dir.exists():
            return False

        fasta_files = list(species_dir.glob("IG*.fasta"))
        return len(fasta_files) > 0

    def get_metadata(self) -> ProviderMetadata:
        """
        Get OGRDB provider metadata.

        Returns
        -------
        ProviderMetadata
            Provider metadata
        """
        species = []
        if self.data_dir.exists():
            for species_dir in self.data_dir.iterdir():
                if species_dir.is_dir():
                    if list(species_dir.glob("IG*.fasta")):
                        species.append(species_dir.name)

        return ProviderMetadata(
            name="ogrdb",
            version="latest",  # TODO: Extract from OGRDB API
            last_updated=datetime.now(),
            species_available=species,
            url="https://ogrdb.airr-community.org/",
        )

    def download(self, species: List[str]) -> None:
        """
        Download OGRDB data for species.

        TODO: Implement OGRDB API download automation
        See scripts/download_ogrdb.py for manual download instructions

        Parameters
        ----------
        species : List[str]
            Species to download

        Raises
        ------
        NotImplementedError
            Automatic download not yet implemented
        """
        raise NotImplementedError(
            "OGRDB automatic download not yet implemented. "
            "See scripts/download_ogrdb.py for manual download instructions."
        )
