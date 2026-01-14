"""
Germline Providers - Data Source Abstraction
=============================================

Providers abstract different germline database sources:
- IMGT: International ImMunoGeneTics information system
- OGRDB: Open Germline Receptor Database
- Custom: User-supplied sequences

All providers implement a common interface defined in base.py
"""

from .base import GermlineProvider

__all__ = ["GermlineProvider"]
