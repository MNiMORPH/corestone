"""
Probe C: is the generated network geologically realistic, by an outside measure?

Checked with `fractopo` (Ovaskainen 2023, JOSS), which is an established tool
for characterising real outcrop trace maps. The measure is the Sanderson &
Nixon node classification:

    X  crossing            two joints cut through each other
    Y  abutting (T)        a joint terminates against another -- CONNECTED
    I  isolated tip        a joint ends in intact rock -- NOT connected
    E  edge                a tip at the domain boundary, uninformative

A well-connected outcrop network is Y-dominated with few I nodes. A scatter of
free-floating segments is I-dominated. That is the difference the throughgoing
plus abutting construction is supposed to make, and this probe measures it
rather than asserting it.

Needs fractopo, which is NOT a dependency of corestone -- it is a check run by
hand, not part of the model. Run it in an environment that has it:

    <venv>/bin/python prototypes/probe_c_topology.py
"""
import numpy as np
import geopandas as gpd
from shapely.geometry import LineString, Polygon
from fractopo import Network

from corestone import FractureNetwork, JointSet, conjugate_sets

NZ, NX, DX = 38, 50, 0.40
LX, LZ = NX * DX, NZ * DX


def topology(fn, name):
    traces = [LineString([tuple(p0), tuple(p1)]) for p0, p1 in fn.segments
              if np.hypot(*(p1 - p0)) > 1e-9]
    gdf = gpd.GeoDataFrame(geometry=traces)
    area = gpd.GeoDataFrame(
        geometry=[Polygon([(0, 0), (LX, 0), (LX, LZ), (0, LZ)])])
    # Abutting traces terminate exactly on their host, so many endpoints are
    # coincident. fractopo's snapping needs a threshold that is neither so
    # tight that it loops nor so loose that it merges distinct nodes.
    d = fn.distance_to_fracture()
    try:
        net = Network(gdf, area, name=name, determine_branches_nodes=True,
                      snap_threshold=0.05 * DX, truncate_traces=True)
        n = net.node_counts
    except RecursionError:
        # fractopo's snapping does not converge when very many traces
        # terminate on exactly the same host points. Report it; do not hide it.
        print(f"{name:26s} traces {len(traces):3d}  P21 {fn.p21:5.2f}   "
              f"fractopo: snapping did not converge   "
              f"median d {np.median(d):.2f} m  p90 {np.percentile(d, 90):.2f} m")
        return
    print(f"{name:26s} traces {len(traces):3d}  P21 {fn.p21:5.2f}   "
          f"X {n['X']:3d}  Y {n['Y']:3d}  I {n['I']:3d}   "
          f"Y/(Y+I) {n['Y'] / max(n['Y'] + n['I'], 1):5.2f}   "
          f"CpB {net.parameters['Connections per Branch']:4.2f}   "
          f"median d {np.median(d):.2f} m  p90 {np.percentile(d, 90):.2f} m")


rng = lambda: np.random.default_rng(12345)

# Three sets with NO abutting: every joint spans the domain and simply crosses
# the others. An X-dominated network, and far denser than the intensity asked
# for -- the contrast that shows what the abutting rule buys.
allcross = [JointSet("F1", dip_deg=75., spacing=1.5),
            JointSet("F2", dip_deg=-75., spacing=1.5),
            JointSet("F3", dip_deg=5., spacing=1.5)]
print("Sanderson & Nixon node topology, measured by fractopo\n")
topology(FractureNetwork(NZ, NX, DX).seed(sets=allcross, rng=rng()),
         "3 sets, all throughgoing")
topology(FractureNetwork(NZ, NX, DX).seed(sets=conjugate_sets(90., 0.), rng=rng()),
         "conjugate 90/0, abutting")
topology(FractureNetwork(NZ, NX, DX).seed(sets=conjugate_sets(45., -45.), rng=rng()),
         "conjugate +/-45, abutting")
topology(FractureNetwork(NZ, NX, DX).seed(
             sets=conjugate_sets(90., 0., density=0.6), rng=rng()),
         "conjugate 90/0, density 0.6")
topology(FractureNetwork(NZ, NX, DX).seed(
             sets=conjugate_sets(90., 0., density=0.25), rng=rng()),
         "conjugate 90/0, density 0.25")
