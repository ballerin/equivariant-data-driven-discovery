"""Shared helpers for the derivative backends.

Kept dependency-free (numpy only, no intra-package imports) so both
``spectral.py`` and ``polynomial.py`` can import it without creating a cycle
with ``derivatives/__init__.py``, which imports them.
"""


def axis_lengths(f, Lx, Ly, Lz=None):
    """Per-axis domain lengths for ``f.ndim`` axes.

    ``ndim <= 3``: ``(Lx,)``, ``(Lx, Ly)`` or ``(Lx, Ly, Lz)`` — independent
    per-axis lengths, so rectangular/anisotropic domains keep working exactly
    as they did before the N-D generalization.

    ``ndim > 3``: **isotropic** — ``Lx`` is reused for every axis.  Every
    simulator in this codebase that goes beyond 3-D (``allen_cahn_nd``,
    ``cahn_hilliard_nd``, ...) produces isotropic hypercube domains, so there
    is no anisotropic N-D caller to support; handling one would require a
    signature change (a full list of N lengths).

    This is the single source of truth for the convention — ``spectral.py``,
    ``polynomial.py`` and ``derivatives/__init__.py`` all delegate here, so
    the three cannot drift apart.

    Raises
    ------
    ValueError
        If ``f`` is 3-D and ``Lz`` was not supplied.  (Without this check the
        missing length surfaces much later as a confusing ``TypeError:
        unsupported operand type(s) for /: 'float' and 'NoneType'`` inside
        the FFT wavenumber construction.)
    """
    ndim = f.ndim
    if ndim == 1:
        return (Lx,)
    if ndim == 2:
        return (Lx, Ly)
    if ndim == 3:
        if Lz is None:
            raise ValueError(
                "Lz is required for 3-D fields (got Lz=None). Pass the third "
                "domain length explicitly; it is only optional for ndim<=2 "
                "and for ndim>3, where the domain is isotropic and Lx is "
                "reused for every axis."
            )
        return (Lx, Ly, Lz)
    return (Lx,) * ndim
