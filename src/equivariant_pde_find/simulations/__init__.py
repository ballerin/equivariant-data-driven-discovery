"""Simulation utilities for generating ground truth for experiments."""

from ._imex import _imex_step, _run
from ._etdrk4 import _etdrk4_coeffs, _etdrk4_run
from .allen_cahn import allen_cahn, allen_cahn_nd
from .swift_hohenberg import swift_hohenberg, swift_hohenberg_nd
from .kuramoto_sivashinsky import kuramoto_sivashinsky
from .cahn_hilliard import cahn_hilliard, cahn_hilliard_nd
from .porous_medium import porous_medium, porous_medium_nd
from .regression_data import make_regression_data
from .weak_data import (
    make_weak_regression_data,
    weak_form_from_snapshots,
    rescale_coefficients,
)
