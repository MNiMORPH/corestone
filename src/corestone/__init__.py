"""
corestone: fracture-controlled chemical weathering of granite.

Named for what the model leaves behind. A corestone is not tougher rock -- it
is the same granite, the same minerals, the same temperature. It survives
because the water never reached it, or reached it spent. The rock around it,
which the water did reach, falls apart into grus.

The physics, in one line, and the model carries two versions of it:

    R = k_ox * A * C                    oxidising, the default: biotite Fe(II)
                                        by dissolved O2, driven by how much
                                        oxygen is there
    R = k(T) * A * (1 - C / C_eq)       dissolving: plagioclase into water
                                        approaching quartz saturation, driven
                                        by how far it is from the ceiling

A REACTANT runs out; a PRODUCT accumulates until it stops the reaction. Either
way the water can stop doing work, and rock the water never reached -- or
reached spent -- survives. Fresh water descends the joints, so the joints
decide where weathering happens; ``Weathering.driver`` decides which reaction
it is. See ``weathering.py`` for why the default changed.

The pieces:

- :class:`corestone.FractureNetwork` -- seed the joints that carry the water,
  or supply your own with :meth:`~corestone.FractureNetwork.from_masks`.
- :class:`corestone.JointSet` -- one family of subparallel joints.
- :func:`corestone.conjugate_sets` -- a throughgoing set and one that
  abuts it, 90 degrees apart by default.
- :func:`corestone.orthogonal_grid` -- the perfectly regular case: exactly
  vertical and horizontal joints, evenly spaced, identical square blocks.
- :class:`corestone.Weathering` -- weather the rock along those joints,
  by oxidation or by dissolution.

Units are SI throughout unless a docstring says otherwise; see ``CLAUDE.md``
for the sign and axis conventions this model assumes -- in particular that the
vertical coordinate is depth, positive downward.
"""

from ._version import __version__
from .fractures import (FractureNetwork, JointSet, conjugate_sets,
                        orthogonal_grid, uniform_grid_shape,
                        periodic_grid_shape, rotated_grid_shape,
                        tiling_angles, tiling_spacing, tiling_spacings,
                        GRANITE_SETS)
from .weathering import Weathering, YEAR

__all__ = ["__version__", "FractureNetwork", "JointSet",
           "conjugate_sets", "orthogonal_grid", "uniform_grid_shape",
           "periodic_grid_shape", "rotated_grid_shape",
           "tiling_angles", "tiling_spacing", "tiling_spacings",
           "GRANITE_SETS",
           "Weathering", "YEAR"]
