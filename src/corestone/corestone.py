#! /usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt


class Model(object):
    """
    ONE LINE: what this evolves.
    """

    def __init__(self):
        # State
        self.z = None           # the evolving field [m]
        self.x = None           # node positions [m]
        self.dx = None          # node spacing [m]

        # Parameters
        self.k = None           # rate coefficient [m2/s]

        # Time
        self.t = 0.             # model time [s]
        self.dt = None          # time step [s]

        # Output
        self.t_out = []         # times at which state was recorded [s]
        self.z_out = []         # recorded state [m]

    # ---- Parameter setters (one per parameter; keep the units in the comment)

    def set_k(self, value):
        """Rate coefficient [m2/s]."""
        self.k = value

    def set_dt(self, value):
        """Time step [s]."""
        self.dt = value

    # ---- Lifecycle

    def initialize(self, x, z):
        """
        Set the grid and the initial state. Run once, before update().
        """
        self.x = np.asarray(x, dtype=float)
        self.z = np.asarray(z, dtype=float).copy()
        self.dx = self.x[1] - self.x[0]
        self.record()

    def update(self, dt=None):
        """
        Advance the state by one time step.

        REPLACE the body: this is an explicit-Euler linear diffusion step, here
        only so that the skeleton runs and the test has something to check.
        """
        if dt is None:
            dt = self.dt
        d2z_dx2 = np.zeros_like(self.z)
        d2z_dx2[1:-1] = (self.z[2:] - 2. * self.z[1:-1] + self.z[:-2]) \
                        / self.dx**2
        self.z = self.z + self.k * d2z_dx2 * dt
        self.t += dt
        self.record()

    def run(self, nt):
        """
        Advance nt time steps.
        """
        for i in range(nt):
            self.update()

    def finalize(self):
        """
        Convert recorded output to arrays. Run once, after the last update().
        """
        self.t_out = np.array(self.t_out)
        self.z_out = np.array(self.z_out)

    # ---- Output

    def record(self):
        """Append the current state to the output lists."""
        self.t_out.append(self.t)
        self.z_out.append(self.z.copy())

    def plot(self, show=True):
        """Plot the initial and final state."""
        plt.figure()
        plt.plot(self.x, self.z_out[0], '--', label='initial')
        plt.plot(self.x, self.z_out[-1], '-', label='final')
        plt.xlabel('Distance [m]')
        plt.ylabel('z [m]')
        plt.legend()
        if show:
            plt.show()
