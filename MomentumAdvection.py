import numpy as np
from Values import Namelist
from wenoInterpolation import JSWENO_interpolation, WENOM_interpolation, WENOZ_interpolation, WENOSG_interpolation

class MomentumAdvection:

    def __init__(self, wind):
        self.nml = Namelist()
        self.um2 = np.roll(wind, 2)
        self.um1 = np.roll(wind, 1)
        self.u0  = np.roll(wind, 0)
        self.up1 = np.roll(wind, -1)
        self.up2 = np.roll(wind, -2)
        self.up3 = np.roll(wind, -3)
        self.tendency = None
        self.dx = self.nml.dx
        self.dy = self.nml.dy

    def __mul__(self, other):
        return self.tendency * other

    __rmul__ = __mul__


    def fill(self,MA):

            if MA == 'WENO-JS':
                    up = JSWENO_interpolation(self.um2,self.um1,self.u0,self.up1,self.up2)
                    um = JSWENO_interpolation(self.up3,self.up2,self.up1,self.u0,self.um1)
            elif MA == 'WENO-M':
                    up = WENOM_interpolation(self.um2,self.um1,self.u0,self.up1,self.up2)
                    um = WENOM_interpolation(self.up3,self.up2,self.up1,self.u0,self.um1)
            elif MA == 'WENO-Z':
                    up = WENOZ_interpolation(self.um2,self.um1,self.u0,self.up1,self.up2)
                    um = WENOZ_interpolation(self.up3,self.up2,self.up1,self.u0,self.um1)       
            elif MA == 'WENO-SG':
                    up = WENOSG_interpolation(self.um2,self.um1,self.u0,self.up1,self.up2)
                    um = WENOSG_interpolation(self.up3,self.up2,self.up1,self.u0,self.um1) 

            advecting_velocity = .5 * (up + um)
            flux = .5*(advecting_velocity + abs(advecting_velocity))*up + .5*(advecting_velocity - abs(advecting_velocity))*um

            self.tendency = - (flux - np.roll(flux, 1)) / self.dx

