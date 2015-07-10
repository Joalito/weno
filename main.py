from Simulation import Simulation


def main():

    sim = Simulation()

    sim.initialize()
    #sim.run(5000, False)
    #sim.plot_energy()
    sim.plot_movie()


if __name__ == "__main__":
    main()
