"""Auto_Mart3 - Tools for coarse-grained molecular dynamics."""

__version__ = "1.0.0"
__author__ = "Anacleto Silva de Souza"

from . import scripts
from .cli import main

__all__ = ["main", "scripts"]
