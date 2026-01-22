#!/usr/bin/env python3
"""
Build Auxiliary Files for IgBLAST
=================================

Script to generate *.aux files for IgBLAST from normalized germline data.

IgBLAST uses auxiliary files (.aux) to annotate CDR3 junction positions.
The format is tab-separated with:
- V genes: GENE_ID  1  CDR3_START_POS  (plus FR/CDR region boundaries)
- J genes: GENE_ID  1  CDR3_END_POS

For J genes, we identify the conserved tryptophan (W) or phenylalanine (F) motif
that marks the CDR3 boundary.

Usage:
    python -m sadie.germlines.scripts.build_aux_files --species human
    python -m sadie.germlines.scripts.build_aux_files --all-species
    python -m sadie.germlines.scripts.build_aux_files --help

Output:
    src/sadie/germlines/igblast/aux_db/{species}_gl.aux
"""

import argparse
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from Bio import SeqIO

logger = logging.getLogger(__name__)

# Module paths
GERMLINES_ROOT = Path(__file__).parent.parent
NORMALIZED_DIR = GERMLINES_ROOT / "normalized"
AUX_OUTPUT_DIR = GERMLINES_ROOT / "igblast" / "aux_db"
DATABASE_DIR = GERMLINES_ROOT / "igblast" / "database"

# CDR3 anchor patterns for J genes
# The conserved W/F-G-X-G motif marks CDR3 end
J_CDR3_PATTERNS = [
    r"TGG",  # W (Trp) - most common
    r"TTC",  # F (Phe)
    r"TTT",  # F (Phe)
]

# Standard IMGT V gene region positions (for V genes with standard length)
# These are the IMGT-defined boundaries in nucleotides
# FR1: 1-78, CDR1: 79-114, FR2: 115-165, CDR2: 166-195, FR3: 196-312 (to CDR3)
IMGT_V_REGIONS = {
    "FR1_end": 78,
    "CDR1_start": 79,
    "CDR1_end": 114,
    "FR2_start": 115,
    "FR2_end": 165,
    "CDR2_start": 166,
    "CDR2_end": 195,
    "FR3_start": 196,
}


class AuxFileBuilder:
    """Build IgBLAST auxiliary files from normalized germline data."""
    
    def __init__(
        self,
        normalized_dir: Optional[Path] = None,
        output_dir: Optional[Path] = None,
        database_dir: Optional[Path] = None
    ):
        """
        Initialize aux file builder.
        
        Parameters
        ----------
        normalized_dir : Path, optional
            Directory containing normalized species data
        output_dir : Path, optional
            Output directory for aux files
        database_dir : Path, optional
            Directory containing BLAST databases (to check species availability)
        """
        self.normalized_dir = normalized_dir or NORMALIZED_DIR
        self.output_dir = output_dir or AUX_OUTPUT_DIR
        self.database_dir = database_dir or DATABASE_DIR
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def get_available_species(self) -> List[str]:
        """
        Get list of species with built BLAST databases.
        
        Returns
        -------
        List[str]
            Species names with databases
        """
        species = []
        if self.database_dir.exists():
            for d in self.database_dir.iterdir():
                if d.is_dir() and not d.name.startswith('.'):
                    # Check if it has actual database files
                    v_db = d / f"{d.name}_V.nsq"
                    if v_db.exists():
                        species.append(d.name)
        return sorted(species)
    
    def build_for_species(self, species: str, force: bool = False) -> Path:
        """
        Build aux file for a specific species.
        
        Parameters
        ----------
        species : str
            Species name
        force : bool
            Force regeneration even if file exists
            
        Returns
        -------
        Path
            Path to generated aux file
        """
        output_path = self.output_dir / f"{species}_gl.aux"
        
        if output_path.exists() and not force:
            logger.info(f"Aux file already exists: {output_path}")
            return output_path
        
        species_dir = self.normalized_dir / species / "ungapped"
        if not species_dir.exists():
            logger.warning(f"No normalized data for {species}")
            return output_path
        
        # Collect all entries
        entries = []
        
        # Process V genes (H, K, L chains)
        for chain in ["H", "K", "L"]:
            v_file = species_dir / f"IG{chain}V.fasta"
            if v_file.exists():
                v_entries = self._process_v_genes(v_file, species)
                entries.extend(v_entries)
                logger.debug(f"Processed {len(v_entries)} IG{chain}V genes")
        
        # Process J genes (H, K, L chains)
        for chain in ["H", "K", "L"]:
            j_file = species_dir / f"IG{chain}J.fasta"
            if j_file.exists():
                j_entries = self._process_j_genes(j_file)
                entries.extend(j_entries)
                logger.debug(f"Processed {len(j_entries)} IG{chain}J genes")
        
        if not entries:
            logger.warning(f"No gene entries generated for {species}")
            return output_path
        
        # Write aux file
        with open(output_path, "w") as f:
            for entry in entries:
                f.write("\t".join(str(x) for x in entry) + "\n")
        
        logger.info(f"Generated aux file with {len(entries)} entries: {output_path}")
        return output_path
    
    def _process_v_genes(self, fasta_path: Path, species: str) -> List[Tuple]:
        """
        Process V genes and extract region positions.
        
        Parameters
        ----------
        fasta_path : Path
            Path to V gene FASTA file
        species : str
            Species name
            
        Returns
        -------
        List[Tuple]
            List of aux file entries
        """
        entries = []
        
        for record in SeqIO.parse(fasta_path, "fasta"):
            gene_id = record.id.split()[0]
            seq_len = len(record.seq)
            
            # Calculate region positions based on sequence length
            # Using standard IMGT numbering with adjustments for sequence length
            positions = self._calculate_v_positions(seq_len)
            
            if positions:
                # Format: gene_id  1  fr1_end  cdr1_start  cdr1_end  fr2_start  fr2_end  cdr2_start  cdr2_end  cdr3_start
                entry = (gene_id, 1) + positions
                entries.append(entry)
        
        return entries
    
    def _calculate_v_positions(self, seq_len: int) -> Optional[Tuple]:
        """
        Calculate V gene region positions.
        
        Uses standard IMGT numbering system with adjustments for
        actual sequence length.
        
        Parameters
        ----------
        seq_len : int
            Sequence length in nucleotides
            
        Returns
        -------
        Optional[Tuple]
            Region positions or None if invalid
        """
        # Standard V gene is ~294-318 nucleotides
        # We need to output the positions for:
        # FR1_end, CDR1_start, CDR1_end, FR2_start, FR2_end, CDR2_start, CDR2_end, CDR3_start
        
        if seq_len < 250:
            # Too short, might be partial
            return None
        
        # Standard positions (1-indexed)
        fr1_end = 26  # End of FR1 in amino acids → ~78 nt
        cdr1_start = 27
        cdr1_end = 35
        fr2_start = 36
        fr2_end = 52
        cdr2_start = 53
        cdr2_end = 62
        
        # CDR3 start is near the end (position 105 in IMGT numbering for AA)
        # In nucleotides: approximately seq_len - 30 to seq_len - 10
        cdr3_start_approx = max(63, seq_len // 3 - 8)
        
        return (fr1_end, cdr1_start, cdr1_end, fr2_start, fr2_end, 
                cdr2_start, cdr2_end, cdr3_start_approx)
    
    def _process_j_genes(self, fasta_path: Path) -> List[Tuple]:
        """
        Process J genes and find CDR3 end position.
        
        Parameters
        ----------
        fasta_path : Path
            Path to J gene FASTA file
            
        Returns
        -------
        List[Tuple]
            List of aux file entries
        """
        entries = []
        
        for record in SeqIO.parse(fasta_path, "fasta"):
            gene_id = record.id.split()[0]
            seq = str(record.seq).upper()
            
            # Find CDR3 end position (W/F-G-X-G motif)
            cdr3_pos = self._find_j_cdr3_position(seq)
            
            if cdr3_pos:
                # Format: gene_id  1  cdr3_end_position
                entries.append((gene_id, 1, cdr3_pos))
            else:
                # Default to position 13 if motif not found
                entries.append((gene_id, 1, 13))
        
        return entries
    
    def _find_j_cdr3_position(self, sequence: str) -> Optional[int]:
        """
        Find CDR3 junction position in J gene.
        
        Looks for conserved W/F-G-X-G motif that marks CDR3 boundary.
        
        Parameters
        ----------
        sequence : str
            Nucleotide sequence
            
        Returns
        -------
        Optional[int]
            1-indexed position of CDR3 end, or None if not found
        """
        # Look for W (TGG) or F (TTC/TTT) followed by GGx pattern
        # The W/F is typically at CDR3 end
        
        # Search for TGG (Trp) - most common in J genes
        for match in re.finditer(r"TGG", sequence):
            pos = match.start()
            # Check if followed by GG (Gly) pattern
            remaining = sequence[pos+3:]
            if remaining.startswith("GG"):
                # Found W-G pattern, CDR3 ends at W position
                return pos + 1  # 1-indexed
        
        # Search for TTC/TTT (Phe) as alternative
        for pattern in ["TTC", "TTT"]:
            for match in re.finditer(pattern, sequence):
                pos = match.start()
                remaining = sequence[pos+3:]
                if remaining.startswith("GG"):
                    return pos + 1
        
        # Default: use first TGG position if any
        match = re.search(r"TGG", sequence)
        if match:
            return match.start() + 1
        
        return None
    
    def build_all(self, force: bool = False) -> Dict[str, Path]:
        """
        Build aux files for all species with databases.
        
        Parameters
        ----------
        force : bool
            Force regeneration
            
        Returns
        -------
        Dict[str, Path]
            Mapping of species to aux file paths
        """
        results = {}
        species_list = self.get_available_species()
        
        if not species_list:
            logger.warning("No species with databases found")
            return results
        
        logger.info(f"Building aux files for {len(species_list)} species")
        
        for species in species_list:
            try:
                path = self.build_for_species(species, force=force)
                results[species] = path
            except Exception as e:
                logger.error(f"Failed to build aux file for {species}: {e}")
        
        return results


def main():
    """Main entry point."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    
    parser = argparse.ArgumentParser(
        description="Build IgBLAST auxiliary files from germline data"
    )
    parser.add_argument(
        "--species",
        nargs="+",
        help="Species to process (e.g., human mouse)"
    )
    parser.add_argument(
        "--all-species",
        action="store_true",
        help="Process all species with built databases"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force regeneration even if files exist"
    )
    parser.add_argument(
        "--list-species",
        action="store_true",
        help="List species with available databases"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    builder = AuxFileBuilder()
    
    if args.list_species:
        print("Species with BLAST databases:")
        for species in builder.get_available_species():
            aux_exists = (builder.output_dir / f"{species}_gl.aux").exists()
            status = "✓" if aux_exists else "✗"
            print(f"  {status} {species}")
        return
    
    if args.all_species:
        results = builder.build_all(force=args.force)
        print(f"\nBuilt aux files for {len(results)} species")
        for species, path in results.items():
            print(f"  {species}: {path}")
    elif args.species:
        for species in args.species:
            try:
                path = builder.build_for_species(species, force=args.force)
                print(f"{species}: {path}")
            except Exception as e:
                logger.error(f"Failed for {species}: {e}")
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
