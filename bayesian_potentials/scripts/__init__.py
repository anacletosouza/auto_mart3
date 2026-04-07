"""Scripts for mapping and topology generation."""

from .map_aa_to_cg import main as map_main
from .generate_cg_top import main as top_main

__all__ = ["map_main", "top_main"]
