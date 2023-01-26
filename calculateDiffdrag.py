
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import RK45

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

    def calcNu(self):
        h = np.cross(self.rr, self.vv)
        e = np.cross(self.vv, h)/self.mu - self.rr/np.linalg.norm(self.rr)
        if np.dot(self.rr, self.vv) >= 0:
            self._nu = np.arccos(np.dot(e,self.rr)/(np.linalg.norm(self.e)*np.linalg.norm(self.rr) ))
        else:
            self._nu = 2*np.pi - np.arccos(np.dot(e,r)/(np.linalg.norm(self.e)*np.linalg.norm(self.rr) ))
        return self._nu

    def calcCircularAngle(self):
        return np.sin(self.rr[1] / self.rr[0])

# define our parameters
H = 540*1e3 # m
R_earth = 6371*1e3 #m
a = R_earth + H #m
c_d = 2.5
low_drag_area = 17025/1e6 #m**2 # XXX: unit!
high_drag_area = 52640/1e6 #m**2 # # XXX: unit

rho = 2.51e-14 # kg/m**3 # XXX wrong unit or convert other units

mu_earth = 3.986004418*1e14 # m**3/s**2

visviva = lambda r,a: mu_earth*((2/r) - (1/a))

# calculate force differential
v = visviva(a,a)

cubesat_mass = 3 # kg

sat1 = Satellite(a, 0, 0, 0, 0, 0, mu_earth, low_drag_area)
sat2 = Satellite(a, 0, 0, 0, 0, 0, mu_earth, high_drag_area)

initial_state = np.array([sat1.rr[0], sat1.rr[1], sat1.rr[2],
                         sat1.vv[0], sat1.vv[1], sat1.vv[2]])
def rhsfun(t, y, sat):
    r = y[0:3]
    v = y[3:]
    drag_f = sat.getAtmosDragForce()
    drag_accel_vec = (drag_f / sat.MASS) * (-v/np.linalg.norm(v))
    dr = v
    dv = -(sat.mu*r)/(np.linalg.norm(r,)**3) + drag_accel_vec
    dy = np.array([
              dr[0], dr[1], dr[2],
              dv[0], dv[1], dv[2]])
    return dy

end_time_s = 1*7*24*60*60
integrator1 = RK45(lambda t,y: rhsfun(t, y, sat1), 0, initial_state, end_time_s, max_step=60)
integrator2 = RK45(lambda t,y: rhsfun(t, y, sat2), 0, initial_state, end_time_s, max_step=60)
ts1 = []
ys1 = []
nu1 = []
while 1:
    integrator1.step()
    ts1.append(integrator1.t)
    ys1.append(integrator1.y)
    sat1.rr = np.array(integrator1.y[0:3])
    nu1.append(sat1.calcCircularAngle())
    assert integrator1.status != 'failed'
    if integrator1.status == 'finished':
        ys1 = np.array(ys1)
        ts1 = np.array(ts1)
        break
ts2 = []
ys2 = []
nu2 = []
while 1:
    integrator2.step()
    ts2.append(integrator2.t)
    ys2.append(integrator2.y)
    sat2.rr = np.array(integrator2.y[0:3])
    nu2.append(sat2.calcCircularAngle())
    assert integrator2.status != 'failed'
    if integrator2.status == 'finished':
        ys2 = np.array(ys2)
        ts2 = np.array(ts2)
        break

xs1 = ys1[:,0:3]
vs1 = ys1[:,3:]
xs2 = ys2[:,0:3]
vs2 = ys2[:,3:]
# ax = plt.figure().add_subplot(projection='3d')
# ax.plot(xs1[:,0], xs1[:,1], xs1[:,2])
plt.plot(ts1, nu1)
plt.show()
