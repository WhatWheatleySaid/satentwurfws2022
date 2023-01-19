


def atmosDragForce(rho, c_d, a, v):
    return .5*rho*c_d*a*v**2

# define our parameters
H = 400 # km
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

f_low = atmosDragForce(rho, c_d, low_drag_area, v)
f_high = atmosDragForce(rho, c_d, high_drag_area, v)

print(f_low)
print(f_high)
print(f_high-f_low)

cubesat_mass = 3 # kg
