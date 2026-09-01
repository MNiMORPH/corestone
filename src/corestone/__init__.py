"""
corestone: ONE LINE saying what this evolves, and of what.

Say what the name means, if it means something. Then, in a sentence or two,
the physics: what is conserved, what drives the change, and what sets the rate.

The pieces:

- :class:`corestone.Model` -- the state and its time evolution.

Units are SI throughout unless a docstring says otherwise; see ``CLAUDE.md``
for the sign and axis conventions this model assumes.
"""

from ._version import __version__
from .corestone import Model

__all__ = ["__version__", "Model"]
