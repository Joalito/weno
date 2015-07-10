import numpy as np
from Values import Namelist

class Grid:
    def __init__(self):
        self.nml = Namelist()
        self.nx = self.nml.nx
        self.dx = self.nml.dx
        self.x = np.asarray([i*self.dx for i in range(self.nx)]) # uniform grid
        self.lx = self.nx * self.dx # physical length

    def __str__(self):
        return "{0}".format(self.x)
