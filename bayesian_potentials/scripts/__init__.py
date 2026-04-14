"""Scripts for mapping and topology generation."""

from .map_aa_to_cg import main as map_main
from .generate_cg_top import main as top_main
from .generate_bonds_angles_dihedrals import main as params_main
from .bp_distributions import main as bp_main

__all__ = ["map_main", "top_main", "params_main", "bp_main"]
