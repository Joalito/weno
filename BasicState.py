import numpy as np

class BasicState:

    def __init__(self, grid):
        #self.u_momentum = np.sin(2*np.pi*grid.x/grid.lx)**3
        #self.u_momentum = np.sin(2*np.pi*grid.x/grid.lx)**3 + 3
        self.u_momentum = np.sin(2*np.pi*grid.x/grid.lx)**2
        #self.u_momentum = np.asarray(delta)

nx = 150
delta = []
for _ in range(nx/3):
    delta.append(0)
for _ in range(nx/3):
    delta.append(1)
for _ in range(nx/3):
    delta.append(0)
