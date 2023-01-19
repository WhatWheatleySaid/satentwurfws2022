
import numpy as np
import matplotlib.pyplot as plt

class Satellite(object):
    MASS = 3 #kg
    DRAG_COEFF = 2.5 # 1
    def __init__(self, a, e, i, om, OM, nu, mu, drag_area):
        self.a = a
        self.e = e
        self.i = i
        self.om = om
        self.OM = OM
        self._nu = nu
        self.mu = mu
        self.drag_area = drag_area
        self.rho = self.calcRho()
        self.getStateVectors()
        self.getDragVector()

    @property
    def nu(self):
        return self._nu

    @nu.setter
    def nu(self, val):
        self._nu = val
        self.getStateVectors()

    def getStateVectors(self):
        p = self.a * (1 - self.e**2)
        r = p / (1 + self.e*np.cos(self.nu))

        rr_pf = np.array([ np.cos(self.nu), np.sin(self.nu), 0 ])*r
        vv_pf = np.sqrt(self.mu/p)*np.array(
            [ -np.sin(self.nu), (self.e + np.cos(self.nu)), 0 ]
            )

        R_OM = np.array([ [np.cos(self.OM),np.sin(self.OM), 0], [-np.sin(self.OM),np.cos(self.OM),0], [0,0,1]])
        R_i = np.array([ [1,0,0],  [0,np.cos(self.i),np.sin(self.i)] , [0,-np.sin(self.i),np.cos(self.i)]])
        R_om = np.array([ [np.cos(self.om),np.sin(self.om),0] , [-np.sin(self.om),np.cos(self.om),0] , [0,0,1] ])

        R = np.matmul(np.matmul(R_om, R_i) , R_OM)

        self.rr = np.matmul(R, rr_pf)
        self.vv = np.matmul(R, vv_pf)

    def getDragVector(self):
        f = self.getAtmosDragForce()
        d_vec = - (f/self.MASS) * self.vv/np.linalg.norm(self.vv)
        self.drag_vec = d_vec
        return d_vec

    @property
    def rho(self):
        return self._rho

    @rho.getter
    def rho(self):
        self.calcRho()
        return self._rho
    @rho.setter
    def rho(self, val):
        self._rho = val

    def calcRho(self):
        # TODO: Write code for atmospheric model?
        return 2.51e-14

    def getAtmosDragForce(self):
        return .5*self.rho*self.DRAG_COEFF*self.drag_area*np.linalg.norm(self.vv)**2

# define our parameters
H = 540 # km
R_earth = 6371 #km
a = R_earth + H #km
c_d = 2.5
low_drag_area = 17025/1e6 #m**2
high_drag_area = 52640/1e6 #m**2

rho = 2.51e-14 # kg/m**3

mu_earth = 0.39860*1e6 # km**3/s**2

visviva = lambda r,a: mu_earth*((2/r) - (1/a))

# calculate force differential
v = visviva(a,a)

cubesat_mass = 3 # kg
# init
r_0_1 = np.array([540, 0, 0])

sat = Satellite(a, 0, 0, 0, 0, 0, mu_earth, low_drag_area)
print(sat.rr)
print(sat.vv)
print(sat.drag_vec)
