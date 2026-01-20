import os
import logging
from functools import lru_cache

logger = logging.getLogger(__name__)

SADIE_USE_GERMLINES_MODULE = "SADIE_USE_GERMLINES_MODULE"


@lru_cache(maxsize=1)
def use_germlines_module() -> bool:
    env_value = os.environ.get(SADIE_USE_GERMLINES_MODULE, "true").lower()
    use_germlines = env_value in ("true", "1", "yes")
    if not use_germlines:
        logger.warning(
            "G3 API is deprecated. Set SADIE_USE_GERMLINES_MODULE=true. "
            "G3 will be removed after 2026-06-01."
        )
    return use_germlines


def clear_feature_flag_cache():
    use_germlines_module.cache_clear()
