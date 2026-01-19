"""
Feature flag utilities for germlines module migration.

This module provides feature flag controls for the gradual migration from G3 API
to the local germlines module, following Constitution Principle V (Integration Compatibility).
"""

import os
import logging
from typing import Literal

logger = logging.getLogger(__name__)


FeatureFlagValue = Literal["true", "false"]


def use_germlines_module() -> bool:
    """
    Determine whether to use the germlines module or fall back to G3 API.

    This feature flag controls the migration strategy from G3 API to the local
    germlines module. During the validation period, users can toggle between
    implementations for testing and rollback capability.

    Environment Variable:
        SADIE_USE_GERMLINES_MODULE: Set to "true" (default) or "false"

    Returns:
        bool: True to use germlines module, False to use G3 API

    Behavior:
        - Default: "true" (use germlines module)
        - "true" -> Use germlines module exclusively (no G3 API calls)
        - "false" -> Use G3 API exclusively (for rollback during validation period only)

    Examples:
        >>> # Default behavior (germlines module)
        >>> use_germlines_module()
        True

        >>> # Explicit germlines module
        >>> os.environ["SADIE_USE_GERMLINES_MODULE"] = "true"
        >>> use_germlines_module()
        True

        >>> # Rollback to G3 (validation period only)
        >>> os.environ["SADIE_USE_GERMLINES_MODULE"] = "false"
        >>> use_germlines_module()
        False

    Notes:
        - No automatic fallback implemented
        - G3 fallback is only for validation period
        - After validation period, G3 dependencies will be removed
    """
    flag_value = os.getenv("SADIE_USE_GERMLINES_MODULE", "true").lower()

    # Log deprecation warning when using G3 fallback
    if flag_value == "false":
        logger.warning(
            "G3 API mode is active (SADIE_USE_GERMLINES_MODULE=false). "
            "This is only supported during the validation period. "
            "G3 API will be deprecated and removed in a future release. "
            "Please report any issues with germlines module and prepare to migrate."
        )

    # Parse flag value (case-insensitive)
    use_germlines = flag_value in ("true", "1", "yes", "on")

    logger.debug(
        f"Feature flag SADIE_USE_GERMLINES_MODULE={flag_value} -> "
        f"use_germlines={use_germlines}"
    )

    return use_germlines
