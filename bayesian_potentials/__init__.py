"""Bayesian Potentials - Tools for coarse-grained molecular dynamics."""

__version__ = "0.1.0"
__author__ = "Anacleto Silva de Souza"

from . import scripts
from .cli import main

__all__ = ["main", "scripts"]
