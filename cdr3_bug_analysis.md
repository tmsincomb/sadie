# CDR3 Null Bug Analysis - Issue #241

## Issue Summary
CDR3 fields (both `cdr3` and `cdr3_aa`) are returning null/NaN values even when VDJ regions are properly populated in the AIRR annotation results.

## Root Cause Analysis

After investigating the codebase, I've identified the primary cause of this bug:

### 1. **Auxiliary Data Path Construction Issue**

The CDR3 detection in IgBLAST depends on auxiliary data files (`.aux` files) that contain information about CDR3 boundaries. The path to these files is constructed in the `GermlineData` class:

**File**: `src/sadie/airr/igblast/germline.py:53`
```python
self.aux_path = self.base_dir / f"aux_db/{scheme}/{name}_gl.aux"
```

**Problem**: The path construction uses a `scheme` parameter that may not match the available auxiliary data directory structure.

### 2. **IgBLAST Auxiliary Data Dependency**

**File**: `src/sadie/airr/igblast/igblast.py:616-621`
```python
@property
def aux_path(self) -> ArgumentType:
    """Auxilary data path. This is needed to lookup the J genes and tell them when the CDR3 stops."""

@aux_path.setter
def aux_path(self, aux_path: Path | str) -> None:
    if not aux_path.exists():
        raise BadIgBLASTArgument(aux_path, "valid path to Auxilary database")
```

The auxiliary path is critical for CDR3 detection. If this path is incorrect or the file doesn't exist, IgBLAST will not be able to determine CDR3 boundaries, resulting in null values.

### 3. **Scheme Parameter Mismatch**

**File**: `src/sadie/airr/igblast/germline.py:24-30`
```python
def __init__(
    self,
    name: str,
    receptor: str = "Ig",
    database_dir: Optional[str | Path] = None,
    scheme: str = "imgt",  # Default scheme
):
```

The default scheme is "imgt", but there may be a mismatch between the scheme being used and the actual directory structure of auxiliary files.

## Potential Causes

1. **Missing auxiliary files**: The required `.aux` files for the specified scheme and species combination don't exist
2. **Incorrect path construction**: The scheme parameter doesn't match the actual directory structure
3. **File permissions**: The auxiliary files exist but are not readable
4. **Case sensitivity**: File name case mismatch (e.g., `human_gl.aux` vs `Human_gl.aux`)

## Suggested Fixes

### 1. **Add Auxiliary File Validation**

Add validation in the `GermlineData.__init__` method to check if the auxiliary file exists:

```python
def __init__(self, name: str, receptor: str = "Ig", database_dir: Optional[str | Path] = None, scheme: str = "imgt"):
    # ... existing code ...
    self.aux_path = self.base_dir / f"aux_db/{scheme}/{name}_gl.aux"
    
    # Add validation
    if not self.aux_path.exists():
        # Try alternative schemes or naming conventions
        alternative_schemes = ["imgt", "kabat"]
        found = False
        for alt_scheme in alternative_schemes:
            alt_path = self.base_dir / f"aux_db/{alt_scheme}/{name}_gl.aux"
            if alt_path.exists():
                self.aux_path = alt_path
                found = True
                warnings.warn(f"Using auxiliary file from {alt_scheme} scheme instead of {scheme}")
                break
        
        if not found:
            raise FileNotFoundError(f"Auxiliary file not found: {self.aux_path}")
```

### 2. **Improve Error Handling in AIRR Class**

Add better error handling and debugging information in the `Airr` class:

```python
def __init__(self, reference_name: str, ...):
    # ... existing code ...
    
    # Validate auxiliary path explicitly
    if not self.germline_data.aux_path.exists():
        raise FileNotFoundError(
            f"Critical auxiliary file missing: {self.germline_data.aux_path}. "
            f"This file is required for CDR3 detection. Please check the installation."
        )
```

### 3. **Add Debugging Output**

Add logging to help diagnose path issues:

```python
def __init__(self, reference_name: str, ...):
    # ... existing code ...
    
    logger.debug(f"Using auxiliary file: {self.germline_data.aux_path}")
    logger.debug(f"Auxiliary file exists: {self.germline_data.aux_path.exists()}")
```

### 4. **Directory Structure Validation**

Add a method to validate the expected directory structure:

```python
def validate_data_structure(self) -> bool:
    """Validate that all required data files exist"""
    required_files = [
        self.aux_path,
        # Add other critical files
    ]
    
    missing_files = [f for f in required_files if not f.exists()]
    if missing_files:
        raise FileNotFoundError(f"Missing required files: {missing_files}")
    return True
```

## Testing Steps

1. **Check auxiliary file existence**: Verify that the auxiliary files exist for the species and scheme being used
2. **Test with different schemes**: Try different schemes (imgt, kabat) to see if auxiliary files exist
3. **Add debug logging**: Enable debug logging to see what paths are being constructed
4. **Test with known working sequences**: Use sequences that are known to work to isolate the issue

## Implementation Priority

1. **High Priority**: Add auxiliary file validation and better error messages
2. **Medium Priority**: Add fallback scheme logic
3. **Low Priority**: Add comprehensive debugging and logging

This fix should resolve the CDR3 null bug by ensuring that the auxiliary data files required for CDR3 detection are properly located and accessible to IgBLAST.