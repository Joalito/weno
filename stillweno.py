import matplotlib.pyplot as plt
import numpy as np
from wenoInterpolation import JSWENO_interpolation, WENOZ_interpolation, WENOM_interpolation

nx = 60
dx = .5
dt = .1
delta = []
for _ in range(nx/3):
    delta.append(0)
for _ in range(nx/3):
    delta.append(1)
for _ in range(nx/3):
    delta.append(0)

delta = np.asarray(delta)
delta2 = delta
delta3 = delta


def positivestencil(delta):
    m2 = np.roll(delta, 2)
    m1 = np.roll(delta, 1)
    p1 = np.roll(delta, -1)
    p2 = np.roll(delta, -2)
    return (m2, m1, delta, p1, p2)

def negativestencil(delta):
    m1 = np.roll(delta, 1)
    p1 = np.roll(delta, -1)
    p2 = np.roll(delta, -2)
    p3 = np.roll(delta, -3)
    return (p3, p2, p1, delta, m1)

def flux_from_weno(phi1, phi2):
    global dx
    phiav = .5 * (phi1 + phi1)
    flux = -(phiav - np.roll(phiav, 1)) / dx
    return flux
    

for _ in range(50):

    phip = JSWENO_interpolation(*positivestencil(delta))
    phim = JSWENO_interpolation(*negativestencil(delta))

    phip2 = WENOZ_interpolation(*positivestencil(delta2))
    phim2 = WENOZ_interpolation(*negativestencil(delta2))

    phip3 = WENOM_interpolation(*positivestencil(delta3))
    phim3 = WENOM_interpolation(*negativestencil(delta3))


    deltatemp1 = delta + flux_from_weno(phip, phim)*dt
    delta2temp1 = delta2 + flux_from_weno(phip2, phim2)*dt
    delta3temp1 = delta3 + flux_from_weno(phip3, phim3)*dt

    phip = JSWENO_interpolation(*positivestencil(deltatemp1))
    phim = JSWENO_interpolation(*negativestencil(deltatemp1))

    phip2 = WENOZ_interpolation(*positivestencil(delta2temp1))
    phim2 = WENOZ_interpolation(*negativestencil(delta2temp1))

    phip3 = WENOM_interpolation(*positivestencil(delta3temp1))
    phim3 = WENOM_interpolation(*negativestencil(delta3temp1))


    deltatemp2 = .75*delta + .25* deltatemp1 + .25*dt*flux_from_weno(phip, phim)
    delta2temp2 = .75*delta2 + .25* delta2temp1 + .25*dt*flux_from_weno(phip2, phim2)
    delta3temp2 = .75*delta3 + .25* delta3temp1 + .25*dt*flux_from_weno(phip3, phim3)

    phip = JSWENO_interpolation(*positivestencil(deltatemp2))
    phim = JSWENO_interpolation(*negativestencil(deltatemp2))

    phip2 = WENOZ_interpolation(*positivestencil(delta2temp2))
    phim2 = WENOZ_interpolation(*negativestencil(delta2temp2))

    phip3 = WENOM_interpolation(*positivestencil(delta3temp2))
    phim3 = WENOM_interpolation(*negativestencil(delta3temp2))


    delta = (1./3.)*delta + (2./3.)*deltatemp2 + (2./3.)*dt*flux_from_weno(phip, phim)
    delta2 = (1./3.)*delta2 + (2./3.)*delta2temp2 + (2./3.)*dt*flux_from_weno(phip2, phim2)
    delta3 = (1./3.)*delta3 + (2./3.)*delta3temp2 + (2./3.)*dt*flux_from_weno(phip3, phim3)



fig = plt.figure()
ax = plt.axes(ylim=(-.25, 1.25))
plt.plot(delta, 'bo', delta2, 'ro', delta3, "go")

plt.show()
