`circlify.py` is vendored from circlify 0.15.1 (PyPI, MIT), unmodified, so
`circle_pack.py` needs no new package in the shared environment. It is pure
stdlib apart from an optional matplotlib import used only by `bubbles()`,
which we do not call.

Note: circlify packs with the A1.0 hole-degree algorithm, not d3-pack, so
coordinates will not match a ggraph `circlepack` layout. That is fine here —
the Python path lays out a different clustering (diffusion), so there is no
cached R layout to reproduce.
