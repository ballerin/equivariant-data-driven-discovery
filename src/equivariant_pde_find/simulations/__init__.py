"""Simulation utilities for generating ground truth for experiments."""

from ._imex import _imex_step, _run
from ._etdrk4 import _etdrk4_coeffs, _etdrk4_run
from .allen_cahn import allen_cahn
from .swift_hohenberg import swift_hohenberg
from .kuramoto_sivashinsky import kuramoto_sivashinsky
from .cahn_hilliard import cahn_hilliard
from .porous_medium import porous_medium
from .regression_data import make_regression_data
from .weak_data import (
    make_weak_regression_data,
    weak_form_from_snapshots,
    rescale_coefficients,
)
