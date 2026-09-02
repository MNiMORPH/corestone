"""
corestone: fracture-controlled chemical weathering of granite.

Named for what the model leaves behind. A corestone is not tougher rock -- it
is the same granite, the same minerals, the same temperature. It survives
because the water never reached it, or reached it already saturated. The rock
around it, which the water did reach, falls apart into grus.

The physics, in one line: dissolution runs at an Arrhenius rate constant
multiplied by how far the pore water is from equilibrium,

    R = k(T) * A * (1 - C / C_eq)

so water that has equilibrated stops weathering rock however soluble the rock
and however warm the water. Fresh water descends the joints; the joints
therefore decide where weathering happens.

The pieces:

- :class:`corestone.FractureNetwork` -- seed the joints that carry the water,
  or supply your own with :meth:`~corestone.FractureNetwork.from_masks`.
- :class:`corestone.JointSet` -- one family of subparallel joints.
- :func:`corestone.conjugate_sets` -- a throughgoing set and one that
  abuts it, 90 degrees apart by default.
- :func:`corestone.orthogonal_grid` -- the perfectly regular case: exactly
  vertical and horizontal joints, evenly spaced, identical square blocks.
- :class:`corestone.Weathering` -- dissolve the rock along those joints.

Units are SI throughout unless a docstring says otherwise; see ``CLAUDE.md``
for the sign and axis conventions this model assumes -- in particular that the
vertical coordinate is depth, positive downward.
"""

from ._version import __version__
from .fractures import (FractureNetwork, JointSet, conjugate_sets,
                        orthogonal_grid, uniform_grid_shape,
                        periodic_grid_shape,
                        GRANITE_SETS)
from .weathering import Weathering, YEAR

__all__ = ["__version__", "FractureNetwork", "JointSet",
           "conjugate_sets", "orthogonal_grid", "uniform_grid_shape",
           "periodic_grid_shape",
           "GRANITE_SETS",
           "Weathering", "YEAR"]
