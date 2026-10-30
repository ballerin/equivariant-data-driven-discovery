"""
Library construction for PDE discovery — Public APIs facade

Builds single-snapshot libraries in the pointwise or space-weak formulation,
selected by the ``formulation`` parameter (see ``library.spacetime`` for weak
forms in time).  Actual implementations are in the following modules:

=====================  ===========================================
Module                  Description
=====================  ===========================================
``multigraph.py``       O(2)-equivariant multigraph operators
``mod_multigraph.py``   Modified multigraph (IBP-friendly) basis
``pdefind.py``          Coordinate-basis (PDE-FIND) library
``wsindy.py``           WSINDy library
``polyjet.py``          Polynomial-jet monomial library
``weak_form.py``        Weak-form infrastructure
``conversions.py``      Basis change-of-basis & library→jet mappings
=====================  ===========================================
"""


# -----------------------------------------------------------------------
# Formulation names
#
# Every builder here works on a single snapshot, so its weak form can only
# be weak in *space*: ``formulation="space_weak"``, with ``space_ibp`` to move
# spatial derivatives onto the test function.  Weak forms in time are
# assembled on top of these builders by ``library.spacetime``.  The names
# ``"weak"`` and ``use_ibp`` predate the time-weak extension and are accepted
# as aliases.
# -----------------------------------------------------------------------

FORMULATION_ALIASES = {"weak": "space_weak"}
WEAK_KWARG_ALIASES = {"use_ibp": "space_ibp"}


def _resolve_formulation(formulation, weak_kwargs):
    """Canonical formulation name and weak kwargs (aliases translated)."""
    formulation = FORMULATION_ALIASES.get(formulation, formulation)
    for old, new in WEAK_KWARG_ALIASES.items():
        if old in weak_kwargs:
            if new in weak_kwargs:
                raise TypeError(f"Pass either {new!r} or its alias {old!r}, "
                                f"not both.")
            weak_kwargs[new] = weak_kwargs.pop(old)
    if formulation not in ("pointwise", "space_weak"):
        raise ValueError(
            f"Unknown formulation {formulation!r}. "
            f"Expected 'pointwise' or 'space_weak'."
        )
    return formulation, weak_kwargs


# -----------------------------------------------------------------------
# Unified library constructors
#
# Extra kwargs (**weak_kwargs) are passed to the weak-form library
# constructors
# -----------------------------------------------------------------------


def build_multigraph_library(f, Lx, Ly, max_p=2, max_n=3,
                             formulation="pointwise", Lz=None, **weak_kwargs):
    """Build the equivariant (multigraph) library matrix.

    Supports 2-D and 3-D fields.

    Parameters
    ----------
    f     : ndarray (Nx, Ny) or (Nx, Ny, Nz)
    Lx, Ly, Lz : float
    max_p : int
    max_n : int
    formulation : {"pointwise", "space_weak"} ("weak" is an alias)
    """
    formulation, weak_kwargs = _resolve_formulation(formulation, weak_kwargs)
    if formulation == "pointwise":
        from .multigraph import _build_multigraph_pointwise
        return _build_multigraph_pointwise(f, Lx, Ly, max_p, max_n, Lz=Lz)
    elif formulation == "space_weak":
        from .multigraph import _build_multigraph_space_weak
        return _build_multigraph_space_weak(f, Lx, Ly, max_p, max_n,
                                            Lz=Lz, **weak_kwargs)


def build_mod_multigraph_library(f, Lx, Ly, max_p=2, max_n=3,
                                 formulation="pointwise", Lz=None, **weak_kwargs):
    """Build the modified multigraph (IBP-friendly) library matrix.

    Supports 2-D and 3-D fields.
    """
    formulation, weak_kwargs = _resolve_formulation(formulation, weak_kwargs)
    if formulation == "pointwise":
        from .mod_multigraph import _build_mod_multigraph_pointwise
        return _build_mod_multigraph_pointwise(f, Lx, Ly, max_p, max_n, Lz=Lz)
    elif formulation == "space_weak":
        from .mod_multigraph import _build_mod_multigraph_space_weak
        return _build_mod_multigraph_space_weak(f, Lx, Ly, max_p, max_n,
                                                Lz=Lz, **weak_kwargs)


def build_pdefind_library(f, Lx, Ly, max_order=4, max_poly=3,
                          formulation="pointwise", Lz=None, **weak_kwargs):
    """Build a coordinate-based (PDE-FIND) library.

    Supports 2-D and 3-D fields.
    For max_order=4, max_poly=3: 60 columns (2-D) or 140 columns (3-D).
    """
    formulation, weak_kwargs = _resolve_formulation(formulation, weak_kwargs)
    if formulation == "pointwise":
        from .pdefind import _build_pdefind_pointwise
        return _build_pdefind_pointwise(f, Lx, Ly, max_order, max_poly, Lz=Lz)
    elif formulation == "space_weak":
        from .pdefind import _build_pdefind_space_weak
        return _build_pdefind_space_weak(f, Lx, Ly, max_order, max_poly,
                                         Lz=Lz, **weak_kwargs)


def build_wsindy_library(f, Lx, Ly, max_order=4, max_poly=3,
                         formulation="pointwise", Lz=None, **weak_kwargs):
    """Build the WSINDy library (∂^α(f^q)).

    Supports any ndim (isotropic domain beyond ndim=3; see
    derivatives/spectral.py's module docstring).

    For max_order=4 and max_poly=3 this gives 46 columns (2-D).

    Parameters
    ----------
    f         : ndarray, shape (N_1, ..., N_ndim) - data
    Lx, Ly, Lz : float - domain lengths (isotropic beyond ndim=3)
    max_order : int - highest partial derivative order (default 4)
    max_poly  : int - highest monomial power q in ∂^α(f^q) (default 3)
    formulation : {"pointwise", "space_weak"} ("weak" is an alias)

    Returns
    -------
    Theta        : ndarray, shape (f.size, n_cols)  [pointwise]  or
                   ndarray, shape (N_valid, n_cols) [weak]
    descriptions : list of n_cols strings
    """
    formulation, weak_kwargs = _resolve_formulation(formulation, weak_kwargs)
    if formulation == "pointwise":
        from .wsindy import _build_wsindy_pointwise
        return _build_wsindy_pointwise(f, Lx, Ly, max_order, max_poly, Lz=Lz)
    elif formulation == "space_weak":
        from .wsindy import _build_wsindy_space_weak
        return _build_wsindy_space_weak(f, Lx, Ly, max_order, max_poly,
                                        Lz=Lz, **weak_kwargs)


def build_polyjet_library(f, Lx, Ly, max_order=4, max_degree=3,
                          formulation="pointwise", Lz=None, **weak_kwargs):
    """Build the polynomial jet library: all monomials in jet variables.

    Supports any ndim (isotropic domain beyond ndim=3; see
    derivatives/spectral.py's module docstring).  The *jet variables* are
    all partial derivatives ∂^α f for 0 <= |α| <= max_order (including f
    itself).  A library column is a product (monomial) of these jet
    variables.

    For max_order=4 and max_degree=3 this gives 105 columns (2-D).

    .. note::
       This enumerates *all* monomials satisfying the degree/order caps.
       For cross-library evaluation where only the monomials needed to
       represent specific operators are required, use
       :func:`build_polyjet_from_descriptions` together with
       ``collect_needed_jet_monomials`` from ``library.conversions``.

    Parameters
    ----------
    f          : ndarray, shape (N_1, ..., N_ndim) - data
    Lx, Ly, Lz : float - domain lengths (isotropic beyond ndim=3)
    max_order  : int - maximum total differential order (default 4)
    max_degree : int - maximum total polynomial degree (default 3)
    formulation : {"pointwise", "space_weak"} ("weak" is an alias)

    Returns
    -------
    Theta        : ndarray, shape (f.size, n_cols)  [pointwise]  or
                   ndarray, shape (N_valid, n_cols) [weak]
    descriptions : list of n_cols strings
    """
    formulation, weak_kwargs = _resolve_formulation(formulation, weak_kwargs)
    if formulation == "pointwise":
        from .polyjet import _build_polyjet_pointwise
        return _build_polyjet_pointwise(f, Lx, Ly, max_order, max_degree, Lz=Lz)
    elif formulation == "space_weak":
        from .polyjet import _build_polyjet_space_weak
        return _build_polyjet_space_weak(f, Lx, Ly, max_order,
                                         max_degree, Lz=Lz, **weak_kwargs)


def build_polyjet_from_descriptions(f, Lx, Ly, descriptions,
                                    formulation="pointwise", Lz=None,
                                    **weak_kwargs):
    """Build a polynomial-jet library from an explicit list of monomial names.

    Unlike ``build_polyjet_library``, which enumerates all monomials under
    ``max_order``/``max_degree`` caps, this builds exactly the columns named
    in *descriptions*.  The set of monomials is typically produced by
    ``collect_needed_jet_monomials`` from ``library.conversions``.

    This is the builder used by ``evaluation.build_eval_theta_jet`` when
    ``jet_monomials`` is provided — the recommended path for cross-library
    evaluation.

    Parameters
    ----------
    f : ndarray, shape (N_1, ..., N_ndim)
    Lx, Ly, Lz : float — domain lengths
    descriptions : list of str
        Monomial description strings (canonical format from
        ``_describe_exponents`` / ``describe_monomial``).
    formulation : {"pointwise", "space_weak"} ("weak" is an alias)
        Only ``"pointwise"`` is currently supported.

    Returns
    -------
    Theta : ndarray, shape (f.size, len(descriptions))  [pointwise]
    descriptions : list of str (echoed back for API consistency)
    """
    formulation, weak_kwargs = _resolve_formulation(formulation, weak_kwargs)
    if formulation == "pointwise":
        from .polyjet import _build_polyjet_from_descriptions
        return _build_polyjet_from_descriptions(f, Lx, Ly, descriptions, Lz=Lz)
    elif formulation == "space_weak":
        raise NotImplementedError(
            "Weak-form dynamic polyjet is not yet implemented."
        )


def build_pdefind_library_1d(u, L, max_deriv=3, max_poly=3):
    """Build a coordinate-based library for a 1D periodic scalar field.

    Parameters
    ----------
    u         : ndarray, shape (N,) — field values on a uniform periodic grid
    L         : float — domain length
    max_deriv : int — highest derivative order (default 3)
    max_poly  : int — highest polynomial degree factor (default 3)

    Returns
    -------
    Theta        : ndarray, shape (N, n_cols)
    descriptions : list of n_cols strings
    """
    from .pdefind import _build_pdefind_pointwise_1d
    return _build_pdefind_pointwise_1d(u, L, max_deriv, max_poly)
