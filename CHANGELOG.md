# Changelog

All notable changes to this project are documented here, following
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Nothing released yet. The model runs, but the browser front end does not exist
and every physical parameter is a placeholder.

### Added

- `corestone.FractureNetwork` and `corestone.JointSet`: seed a joint network
  into a 2D vertical cross section. Joints are placed by *set*, following the
  outcrop hierarchy -- a primary set of throughgoing joints clipped from
  infinite lines so both tips land on the domain boundary, and a secondary set
  cut back to run from one primary joint to the next and terminate there.
- `corestone.conjugate_sets()`: a throughgoing set and one that abuts it, 90
  degrees apart by default (vertical joints cut by horizontal ones).
  `conjugate_sets(45, -45)` gives the symmetric shear pair instead.
- `FractureNetwork.distance_to_fracture()`: distance from every cell to the
  nearest joint, the quantity that decides which rock can weather at all.
- `corestone.Weathering`: steady gravity-driven flow routed down the joints in
  one sweep, and dissolution at an Arrhenius rate constant modulated by the
  chemical affinity of the pore water. Two solid phases, one soluble and one
  inert, so the rock disaggregates into grus rather than dissolving to a
  cavity.
- `corestone.orthogonal_grid()`: a perfectly regular vertical/horizontal joint
  network, with no orientation scatter and exact spacing.
- `examples/seed_a_joint_network.py` and `examples/figure_three_panel.py`.
- Design documents recording each decision, the probe that settled it, and the
  parameters it introduced, in `design/`.

### Notes

- Default resolution is 5 cm over a 20 x 15 m section, 120,000 cells. This is
  an accuracy choice: at 40 cm a joint was smeared across a whole cell and grus
  came out 4 percentage points high. The result is converged by 10 cm.
- The joint network is validated against
  [fractopo](https://github.com/nialov/fractopo) in
  `prototypes/probe_c_topology.py`. fractopo is not a dependency: it is a check
  run by hand.
