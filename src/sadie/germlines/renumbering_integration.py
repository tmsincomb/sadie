"""
Renumbering Integration for Germlines Module
=============================================

Provides HMM generation from local germlines database for use with
sadie.renumbering.aligners.hmmer module.

This replaces G3 API calls with local database queries.
"""

import logging
from pathlib import Path
from typing import Optional, List, Tuple
from functools import lru_cache

import pyhmmer
from sadie.typing import Chain, Source, Species

logger = logging.getLogger(__name__)


class LocalHMMBuilder:
    """
    Build HMM models from local germlines database.

    This class replaces G3 client for HMM generation in renumbering workflows.
    It queries the germlines manager and builds HMM models using pyhmmer.

    Examples
    --------
    >>> builder = LocalHMMBuilder()
    >>> hmm = builder.get_hmm(species="human", chain="H")
    >>> # Use in renumbering
    >>> from sadie.renumbering.aligners.hmmer import HMMER
    >>> aligner = HMMER(species="human", chains="H")
    """

    def __init__(self):
        """Initialize HMM builder with pyhmmer components."""
        self.alphabet = pyhmmer.easel.Alphabet.amino()
        self.builder = pyhmmer.plan7.Builder(self.alphabet)
        self.background = pyhmmer.plan7.Background(self.alphabet)

        # Cache directory for HMM files
        from sadie.germlines import get_germlines_base_dir
        self.hmm_dir = get_germlines_base_dir() / "hmms"
        self.hmm_dir.mkdir(parents=True, exist_ok=True)

    @lru_cache(maxsize=None)
    def get_hmm(
        self,
        species: Species,
        chain: Chain,
        source: Source = "imgt"
    ) -> pyhmmer.plan7.HMM:
        """
        Get or build HMM model for species/chain combination.

        Parameters
        ----------
        species : Species
            Species name (e.g., "human", "mouse")
        chain : Chain
            Chain type ("H", "K", or "L")
        source : Source
            Data source (default: "imgt")

        Returns
        -------
        pyhmmer.plan7.HMM
            Compiled HMM model
        """
        hmm_path = self.hmm_dir / f"{species}_{chain}.hmm"

        # Return cached HMM if exists
        if hmm_path.exists():
            with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
                return next(hmm_file)

        # Build new HMM
        logger.info(f"Building HMM for {species} {chain}")
        return self._build_hmm(species, chain, source)

    def _build_hmm(
        self,
        species: str,
        chain: str,
        source: str
    ) -> pyhmmer.plan7.HMM:
        """
        Build HMM from germlines database.

        Parameters
        ----------
        species : str
            Species name
        chain : str
            Chain type
        source : str
            Data source

        Returns
        -------
        pyhmmer.plan7.HMM
            Compiled HMM model
        """
        # Get V and J sequences from germlines
        vj_pairs = self._get_vj_alignment_pairs(species, chain, source)

        if not vj_pairs:
            raise ValueError(f"No sequences found for {species} {chain}")

        # Create Stockholm alignment
        sto_path = self._write_stockholm(vj_pairs, species, chain)

        # Build HMM using pyhmmer
        hmm_path = self.hmm_dir / f"{species}_{chain}.hmm"

        with pyhmmer.easel.MSAFile(
            sto_path,
            digital=True,
            alphabet=self.alphabet,
            format="stockholm"
        ) as msa_file:
            msa = next(msa_file)
            hmm, _, _ = self.builder.build_msa(msa, self.background)

        # Save HMM to disk
        with open(hmm_path, "wb") as f:
            hmm.write(f)

        logger.info(f"Built HMM: {hmm_path}")
        return hmm

    def _get_vj_alignment_pairs(
        self,
        species: str,
        chain: str,
        source: str
    ) -> List[Tuple[str, str]]:
        """
        Get V-J alignment pairs from germlines database.

        Parameters
        ----------
        species : str
            Species name
        chain : str
            Chain type
        source : str
            Data source

        Returns
        -------
        List[Tuple[str, str]]
            List of (gene_name, gapped_aa_sequence) tuples
        """
        from sadie.germlines import get_manager

        manager = get_manager()
        pairs = []

        # Get V and J genes
        for segment in ["V", "J"]:
            genes = manager.get_genes(species, segment, chain)

            for gene in genes:
                if gene.sequence_aa_gapped:
                    # Use gapped amino acid sequence for alignment
                    pairs.append((gene.name, gene.sequence_aa_gapped))

        return pairs

    def _write_stockholm(
        self,
        pairs: List[Tuple[str, str]],
        species: str,
        chain: str
    ) -> Path:
        """
        Write Stockholm alignment file.

        Parameters
        ----------
        pairs : List[Tuple[str, str]]
            List of (name, gapped_sequence) tuples
        species : str
            Species name
        chain : str
            Chain type

        Returns
        -------
        Path
            Path to Stockholm file
        """
        sto_dir = self.hmm_dir.parent / "stockholms"
        sto_dir.mkdir(parents=True, exist_ok=True)

        sto_path = sto_dir / f"{species}_{chain}.sto"

        lines = [
            "# STOCKHOLM 1.0",
            f"#=GF ID {species}_{chain}",
            ""
        ]

        # Find max name length for alignment
        max_len = max(len(name) for name, _ in pairs)

        for name, seq in pairs:
            lines.append(f"{name.ljust(max_len)}  {seq}")

        # Add reference line and terminator
        lines.append("")
        lines.append(f"#=GC RF{''.ljust(max_len)}  {'x' * 128}")
        lines.append("//")

        sto_path.write_text("\n".join(lines))

        return sto_path


def use_local_hmm_builder() -> bool:
    """
    Check if renumbering should use local HMM builder.

    Returns True if germlines module should be used (default),
    False if G3 API should be used (legacy).

    Returns
    -------
    bool
        True to use local HMM builder, False to use G3
    """
    from sadie.germlines.utils import use_germlines_module
    return use_germlines_module()
