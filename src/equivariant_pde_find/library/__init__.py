"""
Library construction for PDE discovery — Public APIs facade

Builds libraries in pointwise/weak-form formulations based on the
``formulation`` parameter.  Actual implementations are in the following
modules:

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
# Unified library constructors
#
# Extra kwargs (**weak_kwargs) are passed to the weak-form library
# constructors
# -----------------------------------------------------------------------


def build_multigraph_library(f, Lx, Ly, max_p=2, max_n=3,
                             formulation="pointwise", **weak_kwargs):
    """Build the equivariant (multigraph) library matrix for a 2-D scalar field.

    Parameters
    ----------
    f     : ndarray (Nx, Ny) - data
    Lx    : float - domain length in x
    Ly    : float - domain length in y
    max_p : int - maximum multigraph edge count p (differential order = 2*max_p)
    max_n : int - maximum polynomial degree
    formulation : {"pointwise", "weak"}

    Returns
    -------
    Theta        : ndarray, shape (Nx*Ny, K)  [pointwise]  or
                   ndarray, shape (N_valid, K) [weak]
    descriptions : list of K strings
    """
    if formulation == "pointwise":
        from .multigraph import _build_multigraph_pointwise
        return _build_multigraph_pointwise(f, Lx, Ly, max_p, max_n)
    elif formulation == "weak":
        from .multigraph import _build_multigraph_weak
        return _build_multigraph_weak(f, Lx, Ly, max_p, max_n,
                                      **weak_kwargs)
    else:
        raise ValueError(
            f"Unknown formulation {formulation!r}. "
            f"Expected 'pointwise' or 'weak'."
        )


def build_mod_multigraph_library(f, Lx, Ly, max_p=2, max_n=3,
                                 formulation="pointwise", **weak_kwargs):
    """Build the modified multigraph (IBP-friendly) library matrix.

    Parameters
    ----------
    f     : ndarray, shape (Nx, Ny)
    Lx, Ly : float
    max_p : int - maximum edge count (0, 1, or 2)
    max_n : int - maximum polynomial degree (0..3)
    formulation : {"pointwise", "weak"}

    Returns
    -------
    Theta        : ndarray, shape (Nx*Ny, K)  [pointwise]  or
                   ndarray, shape (N_valid, K) [weak]
    descriptions : list of K strings
    """
    if formulation == "pointwise":
        from .mod_multigraph import _build_mod_multigraph_pointwise
        return _build_mod_multigraph_pointwise(f, Lx, Ly, max_p, max_n)
    elif formulation == "weak":
        from .mod_multigraph import _build_mod_multigraph_weak
        return _build_mod_multigraph_weak(f, Lx, Ly, max_p, max_n,
                                          **weak_kwargs)
    else:
        raise ValueError(
            f"Unknown formulation {formulation!r}. "
            f"Expected 'pointwise' or 'weak'."
        )


def build_pdefind_library(f, Lx, Ly, max_order=4, max_poly=3,
                          formulation="pointwise", **weak_kwargs):
    """Build a coordinate-based (PDE-FIND) library for a 2-D scalar field.

    For max_order=4 and max_poly=3 this gives 60 columns.

    Parameters
    ----------
    f         : ndarray, shape (Nx, Ny) - data
    Lx, Ly    : float - domain lengths
    max_order : int - highest partial derivative order (default 4)
    max_poly  : int - highest polynomial degree factor (default 3)
    formulation : {"pointwise", "weak"}

    Returns
    -------
    Theta        : ndarray, shape (Nx*Ny, n_cols)  [pointwise]  or
                   ndarray, shape (N_valid, n_cols) [weak]
    descriptions : list of strings
    """
    if formulation == "pointwise":
        from .pdefind import _build_pdefind_pointwise
        return _build_pdefind_pointwise(f, Lx, Ly, max_order, max_poly)
    elif formulation == "weak":
        from .pdefind import _build_pdefind_weak
        return _build_pdefind_weak(f, Lx, Ly, max_order, max_poly,
                                   **weak_kwargs)
    else:
        raise ValueError(
            f"Unknown formulation {formulation!r}. "
            f"Expected 'pointwise' or 'weak'."
        )


def build_wsindy_library(f, Lx, Ly, max_order=4, max_poly=3,
                         formulation="pointwise", **weak_kwargs):
    """Build the WSINDy library (∂^α(f^q)) for a 2-D scalar field.

    For max_order=4 and max_poly=3 this gives 46 columns.

    Parameters
    ----------
    f         : ndarray shape (Nx, Ny) - data
    Lx, Ly    : float - domain lengths
    max_order : int - highest partial derivative order (default 4)
    max_poly  : int - highest monomial power q in ∂^α(f^q) (default 3)
    formulation : {"pointwise", "weak"}

    Returns
    -------
    Theta        : ndarray, shape (Nx*Ny, n_cols)  [pointwise]  or
                   ndarray, shape (N_valid, n_cols) [weak]
    descriptions : list of n_cols strings
    """
    if formulation == "pointwise":
        from .wsindy import _build_wsindy_pointwise
        return _build_wsindy_pointwise(f, Lx, Ly, max_order, max_poly)
    elif formulation == "weak":
        from .wsindy import _build_wsindy_weak
        return _build_wsindy_weak(f, Lx, Ly, max_order, max_poly,
                                  **weak_kwargs)
    else:
        raise ValueError(
            f"Unknown formulation {formulation!r}. "
            f"Expected 'pointwise' or 'weak'."
        )


def build_polyjet_library(f, Lx, Ly, max_order=4, max_degree=3,
                          formulation="pointwise", **weak_kwargs):
    """Build the polynomial jet library: all monomials in jet variables.

    The *jet variables* are all partial derivatives ∂_x^a ∂_y^b f for
    0 ≤ a+b ≤ max_order (including f itself).  A library column is a
    product (monomial) of these jet variables.

    For max_order=4 and max_degree=3 this gives 105 columns.

    Parameters
    ----------
    f          : ndarray, shape (Nx, Ny) - data
    Lx, Ly     : float - domain lengths
    max_order  : int - maximum total differential order (default 4)
    max_degree : int - maximum total polynomial degree (default 3)
    formulation : {"pointwise", "weak"}

    Returns
    -------
    Theta        : ndarray, shape (Nx*Ny, n_cols)  [pointwise]  or
                   ndarray, shape (N_valid, n_cols) [weak]
    descriptions : list of n_cols strings
    """
    if formulation == "pointwise":
        from .polyjet import _build_polyjet_pointwise
        return _build_polyjet_pointwise(f, Lx, Ly, max_order, max_degree)
    elif formulation == "weak":
        from .polyjet import _build_polyjet_weak
        return _build_polyjet_weak(f, Lx, Ly, max_order,
                                   max_degree, **weak_kwargs)
    else:
        raise ValueError(
            f"Unknown formulation {formulation!r}. "
            f"Expected 'pointwise' or 'weak'."
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
