import numpy as np
import matplotlib.pyplot as plt
from Grid import Grid
from Values import Namelist
from MomentumAdvection import MomentumAdvection
from BasicState import BasicState

class Momentum:

    def __init__(self, flux_choice="JS", TI="AB"):
        self.nml = Namelist()
        self.u = None
        self.flux_choice = flux_choice
        self.time_elapsed = 0
        self.timeintegration = TI
        self.timesteps = []
        self.energy = []

    def initialize(self, grid):
            initial = BasicState(grid)
            self.u = initial.u_momentum
            self.dt = self.nml.cfl * grid.dx / np.amax(abs(self.u))

    def initial_stepAB(self):
        global dudt
        global dudtm1
        global JSdudt
        global JSdudtm1
        global Mdudt
        global Mdudtm1
        global Zdudt
        global Zdudtm1

        if self.flux_choice == "center":
            rflux = .25*(self.u + np.roll(self.u, -1))**2
            dudt = -(rflux - np.roll(rflux, 1)) / 1.

            
            dudtm1 = dudt

            self.u += self.dt*dudt
            self.time_elapsed += self.dt
            self.timesteps.append(self.time_elapsed)
            self.energy.append(self.Energy())

            rflux = .25*(self.u + np.roll(self.u, -1))**2
            dudt = -(rflux - np.roll(rflux, 1)) / 1.
            
            self.u = self.u + (1.5*dudt - .5*dudtm1) * self.dt
            self.time_elapsed += self.dt
            self.timesteps.append(self.time_elapsed)
            self.energy.append(self.Energy())

        elif self.flux_choice == "JS":
            JSdudt = MomentumAdvection(self.u)
            JSdudt.fill("WENO-JS")
            
            JSdudtm1 = JSdudt

            self.u += self.dt*JSdudt
            self.time_elapsed += self.dt
            self.timesteps.append(self.time_elapsed)
            self.energy.append(self.Energy())

            JSdudt.fill("WENO-JS")
            
            self.u = self.u + (1.5*JSdudt - .5*JSdudtm1) * self.dt
            self.time_elapsed += self.dt
            self.timesteps.append(self.time_elapsed)
            self.energy.append(self.Energy())

        elif self.flux_choice == "M":
            Mdudt = MomentumAdvection(self.u)
            Mdudt.fill("WENO-M")
            
            Mdudtm1 = Mdudt

            self.u += self.dt*Mdudt
            self.time_elapsed += self.dt
            self.timesteps.append(self.time_elapsed)
            self.energy.append(self.Energy())

            Mdudt.fill("WENO-M")
            
            self.u = self.u + (1.5*Mdudt - .5*Mdudtm1) * self.dt
            self.time_elapsed += self.dt
            self.timesteps.append(self.time_elapsed)
            self.energy.append(self.Energy())

        elif self.flux_choice == "Z":
            Zdudt = MomentumAdvection(self.u)
            Zdudt.fill("WENO-Z")
            
            Zdudtm1 = Zdudt

            self.u += self.dt*Zdudt
            self.time_elapsed += self.dt
            self.timesteps.append(self.time_elapsed)
            self.energy.append(self.Energy())

            Zdudt.fill("WENO-Z")
            
            self.u = self.u + (1.5*Zdudt - .5*Zdudtm1) * self.dt
            self.time_elapsed += self.dt
            self.timesteps.append(self.time_elapsed)
            self.energy.append(self.Energy())

    def update(self):

        global dudt
        global dudtm1
        global dudtm2
        global JSdudt
        global JSdudtm1
        global JSdudtm2
        global Mdudt
        global Mdudtm1
        global Mdudtm2
        global Zdudt
        global Zdudtm1
        global Zdudtm2

        if self.timeintegration == "AB":

            if self.flux_choice == "center":
                dudtm2 = dudtm1
                dudtm1 = dudt

                rflux = .25*(self.u + np.roll(self.u, -1))**2
                dudt = -(rflux - np.roll(rflux, 1)) / 1.

                self.u = self.u + ((23./12.)*dudt - (4./3.)*dudtm1 + (5./12.)*dudtm2) * self.dt
                self.time_elapsed += self.dt
                self.timesteps.append(self.time_elapsed)
                self.energy.append(self.Energy())

            elif self.flux_choice == "JS":
                JSdudtm2 = JSdudtm1
                JSdudtm1 = JSdudt

                JSdudt = MomentumAdvection(self.u)
                JSdudt.fill("WENO-JS")

                self.u = self.u + ((23./12.)*JSdudt - (4./3.)*JSdudtm1 + (5./12.)*JSdudtm2) * self.dt
                self.time_elapsed += self.dt
                self.timesteps.append(self.time_elapsed)
                self.energy.append(self.Energy())

            elif self.flux_choice == "M":
                Mdudtm2 = Mdudtm1
                Mdudtm1 = Mdudt

                Mdudt = MomentumAdvection(self.u)
                Mdudt.fill("WENO-M")

                self.u = self.u + ((23./12.)*Mdudt - (4./3.)*Mdudtm1 + (5./12.)*Mdudtm2) * self.dt
                self.time_elapsed += self.dt
                self.timesteps.append(self.time_elapsed)
                self.energy.append(self.Energy())

            elif self.flux_choice == "Z":
                Zdudtm2 = Zdudtm1
                Zdudtm1 = Zdudt

                Zdudt = MomentumAdvection(self.u)
                Zdudt.fill("WENO-Z")

                self.u = self.u + ((23./12.)*Zdudt - (4./3.)*Zdudtm1 + (5./12.)*Zdudtm2) * self.dt
                self.time_elapsed += self.dt
                self.timesteps.append(self.time_elapsed)
                self.energy.append(self.Energy())

        elif self.timeintegration == "RK3":

            if self.flux_choice == "center":
                rflux = .25*(self.u + np.roll(self.u, -1))**2
                dudtR = -(rflux - np.roll(rflux, 1)) / 1.

                temp1 = self.u + self.dt*dudtR
                
                rflux = .25*(temp1 + np.roll(temp1, -1))**2
                dudtR = -(rflux - np.roll(rflux, 1)) / 1.

                temp2 = .75*self.u + .25*temp1 +.25*self.dt*dudtR
     
                rflux = .25*(temp2 + np.roll(temp2, -1))**2
                dudtR = -(rflux - np.roll(rflux, 1)) / 1.

                self.u = (1./3.)*self.u + (2./3.)*temp2 + (2./3.)*self.dt*dudtR
                self.time_elapsed += self.dt
                self.timesteps.append(self.time_elapsed)
                self.energy.append(self.Energy())

            elif self.flux_choice == "JS":
                dudtR = MomentumAdvection(self.u)
                dudtR.fill("WENO-JS")

                temp1 = self.u + self.dt*dudtR.tendency
                
                dudtR = MomentumAdvection(temp1)
                dudtR.fill("WENO-JS")

                temp2 = .75*self.u + .25*temp1 +.25*self.dt*dudtR.tendency
            
                dudtR = MomentumAdvection(temp2)
                dudtR.fill("WENO-JS")

                self.u = (1./3.)*self.u + (2./3.)*temp2 + (2./3.)*self.dt*dudtR.tendency
                self.time_elapsed += self.dt
                self.timesteps.append(self.time_elapsed)
                self.energy.append(self.Energy())

            elif self.flux_choice == "M":
                dudtR = MomentumAdvection(self.u)
                dudtR.fill("WENO-M")

                temp1 = self.u + self.dt*dudtR.tendency
                
                dudtR = MomentumAdvection(temp1)
                dudtR.fill("WENO-M")

                temp2 = .75*self.u + .25*temp1 +.25*self.dt*dudtR.tendency
            
                dudtR = MomentumAdvection(temp2)
                dudtR.fill("WENO-M")

                self.u = (1./3.)*self.u + (2./3.)*temp2 + (2./3.)*self.dt*dudtR.tendency
                self.time_elapsed += self.dt
                self.timesteps.append(self.time_elapsed)
                self.energy.append(self.Energy())

            elif self.flux_choice == "Z":
                dudtR = MomentumAdvection(self.u)
                dudtR.fill("WENO-Z")

                temp1 = self.u + self.dt*dudtR.tendency
                
                dudtR = MomentumAdvection(temp1)
                dudtR.fill("WENO-Z")

                temp2 = .75*self.u + .25*temp1 +.25*self.dt*dudtR.tendency
            
                dudtR = MomentumAdvection(temp2)
                dudtR.fill("WENO-Z")

                self.u = (1./3.)*self.u + (2./3.)*temp2 + (2./3.)*self.dt*dudtR.tendency
                self.time_elapsed += self.dt
                self.timesteps.append(self.time_elapsed)
                self.energy.append(self.Energy())

    def Energy(self):
        return .5 * np.sum((self.u)**2)






dudt = np.array([])
dudtm1 = np.array([])
dudtm2 = np.array([])

JSdudt = np.array([])
JSdudtm1 = np.array([])
JSdudtm2 = np.array([])

Mdudt = np.array([])
Mdudtm1 = np.array([])
Mdudtm2 = np.array([])

Zdudt = np.array([])
Zdudtm1 = np.array([])
Zdudtm2 = np.array([])
