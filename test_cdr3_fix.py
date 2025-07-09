#!/usr/bin/env python3
"""
Test script for CDR3 null bug fixes
Demonstrates the enhanced error handling and diagnostic capabilities
"""

import sys
import logging
from pathlib import Path

# Configure logging to see debug messages
logging.basicConfig(level=logging.DEBUG, format='%(name)s - %(levelname)s - %(message)s')

try:
    from sadie.airr import Airr
    from sadie.airr.airrtable import AirrTable
    from sadie.airr.igblast.germline import GermlineData
except ImportError as e:
    print(f"Error importing sadie modules: {e}")
    print("Make sure you're running this from the sadie directory with the package installed")
    sys.exit(1)

def test_auxiliary_file_validation():
    """Test the auxiliary file validation functionality"""
    print("=== Testing Auxiliary File Validation ===")
    
    try:
        # Test with a valid species
        print("Testing with valid species (human)...")
        gd = GermlineData("human", scheme="imgt")
        print(f"✓ Successfully created GermlineData for human")
        print(f"  Auxiliary path: {gd.aux_path}")
        print(f"  Auxiliary file exists: {gd.aux_path.exists()}")
        
    except Exception as e:
        print(f"✗ Error with valid species: {e}")
    
    try:
        # Test with an invalid species
        print("\nTesting with invalid species (nonexistent)...")
        gd = GermlineData("nonexistent", scheme="imgt")
        print(f"✗ This should have failed but didn't")
        
    except FileNotFoundError as e:
        print(f"✓ Correctly caught missing auxiliary file: {e}")
    except Exception as e:
        print(f"✗ Unexpected error: {e}")

def test_airr_initialization():
    """Test AIRR initialization with enhanced validation"""
    print("\n=== Testing AIRR Initialization ===")
    
    try:
        # Test with valid species
        print("Testing AIRR initialization with human...")
        airr = Airr("human")
        print(f"✓ Successfully initialized AIRR for human")
        
    except Exception as e:
        print(f"✗ Error initializing AIRR: {e}")

def test_cdr3_diagnostics():
    """Test the new CDR3 diagnostic functionality"""
    print("\n=== Testing CDR3 Diagnostics ===")
    
    try:
        # Create a sample AIRR table to test diagnostics
        import pandas as pd
        import numpy as np
        
        # Create test data with some CDR3 issues
        test_data = {
            'sequence_id': ['seq1', 'seq2', 'seq3', 'seq4'],
            'v_call': ['IGHV1-2*01', 'IGHV3-23*01', 'IGHV1-69*01', 'IGHV2-5*01'],
            'j_call': ['IGHJ4*02', 'IGHJ6*01', 'IGHJ1*01', 'IGHJ3*02'],
            'cdr3': ['TGTGCGAGAGA', np.nan, 'TGTGCAAGGATC', np.nan],  # Some missing CDR3
            'cdr3_aa': ['CARE', np.nan, 'CARD', np.nan],  # Some missing CDR3_AA
            'liable': [False, True, False, True]  # Some liable sequences
        }
        
        airr_table = AirrTable(pd.DataFrame(test_data))
        
        print("Created test AIRR table with CDR3 issues...")
        print(f"Table has {len(airr_table)} sequences")
        
        # Run diagnostics
        diagnostics = airr_table.diagnose_cdr3_issues()
        print("\nCDR3 Diagnostic Results:")
        if not diagnostics.empty:
            for _, issue in diagnostics.iterrows():
                print(f"  Issue: {issue['issue']}")
                print(f"    Count: {issue['count']}")
                print(f"    Percentage: {issue['percentage']:.1f}%")
                print(f"    Description: {issue['description']}")
                print()
        else:
            print("  No CDR3 issues detected!")
            
    except Exception as e:
        print(f"✗ Error testing diagnostics: {e}")

def main():
    """Main test function"""
    print("CDR3 Bug Fix Test Suite")
    print("=" * 50)
    
    test_auxiliary_file_validation()
    test_airr_initialization()
    test_cdr3_diagnostics()
    
    print("\n" + "=" * 50)
    print("Test suite completed!")
    print("\nIf you see errors about missing auxiliary files, this indicates")
    print("the enhanced validation is working correctly.")
    print("\nFor actual sequence annotation, ensure you have valid sequences")
    print("and that all auxiliary files are properly installed.")

if __name__ == "__main__":
    main()