import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from Grid import Grid
from Momentum import Momentum


class Simulation:

    def __init__(self):
        self.grid = Grid()

        self.wind1 = Momentum(flux_choice="center", TI ="AB")
        self.wind2 = Momentum(flux_choice="JS", TI="AB")
        self.wind3 = Momentum(flux_choice="center", TI="RK3")
        self.wind4 = Momentum(flux_choice="JS", TI="RK3")

    def initialize(self):
        self.wind1.initialize(self.grid)
        self.wind2.initialize(self.grid)
        self.wind3.initialize(self.grid)
        self.wind4.initialize(self.grid)



    def run(self, ntimestep, plot=False):
        self.wind1.initial_stepAB()

        self.wind2.initial_stepAB()


        for _ in range(2):
            self.wind3.update()
            self.wind4.update()

        for _ in range(ntimestep):

            self.wind1.update()

            self.wind2.update()

            self.wind3.update()

            self.wind4.update()

        if plot == True:
            fig = plt.figure(3)
            plt.plot(self.grid.x, self.wind1.u)
            plt.show()

    def plot_energy(self):
        Efig = plt.figure(2)
        Eax = plt.axes()

        plt.xlabel(r"$t$", size=19)
        plt.ylabel(r"$Energy\quad(\;\frac{1}{2}\;u^2\;)$", size=19)

        Eax.grid()

        line1, = plt.plot(self.wind1.timesteps, self.wind1.energy, "b", label="AB center")
        line2, = plt.plot(self.wind2.timesteps, self.wind2.energy, "g", label="AB WENO")
        line3, = plt.plot(self.wind4.timesteps, self.wind4.energy, "y", label="RK WENO")
        line4, = plt.plot(self.wind3.timesteps, self.wind3.energy, "r", label="RK center")

        first_legend = plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                   ncol=2, mode="expand", borderaxespad=0.)

        plt.show()

    def plot_movie(self):
        fig = plt.figure(1)
        ax  = plt.axes(xlim=(0, (self.grid.nx + 5)), ylim=((np.amin(self.wind1.u)-.5), (np.amax(self.wind1.u)+.5)))
        ax.grid()

        line, = ax.plot([], [])

        time_text   = ax.text(8, (np.amax(self.wind1.u) + .1), '', size=15)
        energy_text = ax.text(8, (np.amax(self.wind1.u)), '', size=15)

        # called between frames to clear the figure
        def init():
            line.set_data([], [])
            time_text.set_text('')
            energy_text.set_text('')
            return line, time_text, energy_text

        self.wind1.initial_stepAB()

        # animation function.  This is called sequentially
        def animate(i):
            self.wind1.update()
            line.set_data(self.grid.x, self.wind1.u)
            time_text.set_text('time = %.1f' % self.wind1.time_elapsed)
            energy_text.set_text('energy = %.3f J' % self.wind1.Energy())
            return line, time_text, energy_text

        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=10, interval=7, blit=True)

        plt.show()
