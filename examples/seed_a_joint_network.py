#! /usr/bin/python3
"""
Seed a conjugate joint network and look at how far the rock is from water.

Distance to the nearest joint is the quantity that decides corestones: rock
further from a joint than the equilibration length never sees undersaturated
water, whatever its mineralogy. This example seeds a network, measures that
distance, and plots both.

    python3 examples/seed_a_joint_network.py
"""
import numpy as np
from matplotlib import pyplot as plt

from corestone import FractureNetwork, conjugate_sets

NZ, NX, DX = 75, 100, 0.20                     # a 20 x 15 m section
L_EQ = 0.50                                    # equilibration length [m]

net = FractureNetwork(NZ, NX, DX).seed(
    sets=conjugate_sets(dip_primary=90.0, dip_secondary=0.0, spacing=1.5),
    rng=np.random.default_rng(12345))

d = net.distance_to_fracture()

print("%d traces, P21 = %.2f m/m2" % (len(net.segments), net.p21))
print("distance to the nearest joint: median %.2f m, p90 %.2f m, max %.2f m"
      % (np.median(d), np.percentile(d, 90), d.max()))
print("%.1f %% of the rock lies within L_eq = %.2f m of a joint"
      % ((d <= L_EQ).mean() * 100, L_EQ))

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(13, 4.2), sharey=True,
                               constrained_layout=True)
extent = [0.0, net.lx, net.lz, 0.0]            # depth increases downward

for (p0, p1), name in zip(net.segments, net.segment_set):
    through = name == "J1"
    ax0.plot([p0[0], p1[0]], [p0[1], p1[1]], color="#1a1a1a",
             lw=1.0 if through else 0.7,
             ls="-" if through else (0, (3.5, 2.2)))
ax0.set_xlim(0, net.lx)
ax0.set_ylim(net.lz, 0)
ax0.set_aspect("equal")
ax0.set_title("the joints", loc="left")
ax0.set_ylabel("Depth [m]")
ax0.set_xlabel("Distance [m]")

im = ax1.imshow(d, extent=extent, origin="upper", cmap="Purples")
ax1.contour(np.linspace(DX / 2, net.lx - DX / 2, NX),
            np.linspace(DX / 2, net.lz - DX / 2, NZ),
            d, levels=[L_EQ], colors="#4a1d6a", linewidths=1.4)
ax1.set_aspect("equal")
ax1.set_title("distance to the nearest joint; the line is $L_{eq}$", loc="left")
ax1.set_xlabel("Distance [m]")
fig.colorbar(im, ax=ax1, fraction=0.046, label="distance [m]")

plt.show()
