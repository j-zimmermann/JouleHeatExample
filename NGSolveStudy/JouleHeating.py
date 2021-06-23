"""Simulate Joule (Resistive) Heating in NGSolve"""
import ngsolve as ng
from netgen.geom2d import Circle, CSG2d, EdgeInfo as EI, Solid2d

# radius
R = 0.1
# order
order = 4

geo = CSG2d()
rect = Solid2d([(0, 0),
                EI(bc="bottom"),  # set bc for segment (0,0)-(1,0)
                (1, 0),
                EI(bc="right"),  # set bc for segment (1,0)-(1,1)
                (1, 1),
                EI(bc="top"),  # set bc for segment (1,1)-(0,1)
                (0, 1),
                EI(bc="left"),  # set bc for segment (0,1)-(0,0)
                ], mat="rect")
circle = Circle(center=(0.5, 0.5), radius=R, mat="mat2", bc="circle")
geo.Add(rect - circle)
ngmesh = geo.GenerateMesh(maxh=0.4)
mesh = ng.Mesh(ngmesh)
mesh.Curve(order)
ng.Draw(mesh)

datafile = open('TemperatureDevelopment.dat', 'w')

"""Let's start with some functions and parameters needed for the whole model"""
T_0 = ng.CoefficientFunction(300.)
T_ref = ng.CoefficientFunction(293.)
V0 = ng.CoefficientFunction(.1)  # set boundary value of electrode to 10 V
f = ng.CoefficientFunction(0.0)
rho_0 = ng.CoefficientFunction(1.754e-8)  # reference resistivity
C_p = ng.CoefficientFunction(340.)  # heat capacity
rho = ng.CoefficientFunction(8930.)  # density  [kg/m^3]
k_iso = ng.CoefficientFunction(384.)  # thermal conductivity [W/(m*K)]
alpha = ng.CoefficientFunction(3.9e-3)  # coefficient for T-dependent resistivity [1/K]

t_max = 2000.  # time in s
dt = 50.  # time step
num_steps = int(t_max / dt)


# temperature-dependent resistivity
def sigma_T(Temp):
    return ng.CoefficientFunction(1. / (rho_0 * (1. + alpha * (Temp - T_0))))


# declare function space
V = ng.H1(mesh, order=order, dirichlet="left|right")
VT = ng.H1(mesh, order=order, dirichlet="left|right|circle")
print("Number of DOFs: {}".format(V.ndof))

# boundary conditions
# Dirichlet for assigned electric potential
bc_values_pot = {"top": V0,
                 "bottom": f,
                 "circle": 0,
                 "left": V0,
                 "right": f}
bc_cf_pot = ng.CoefficientFunction([bc_values_pot[bc] for bc in mesh.GetBoundaries()])

# Dirichlet for assigned temperature
bc_values_T = {"top": T_0,
               "bottom": T_0,
               "circle": T_0,
               "left": T_0,
               "right": T_0}
bc_cf_T = ng.CoefficientFunction([bc_values_T[bc] for bc in mesh.GetBoundaries()])

# Define initial value, which will later be the solution from the previous (n-th) time step
T_n = ng.GridFunction(VT, "T_n")
T_n.Set(T_ref)
ng.Draw(T_n)

# Define variational problem for temperature and potential
psi_T = V.TrialFunction()  # trial for potential
T_n1_T = V.TrialFunction()  # trial for temperature
v = V.TestFunction()
vT = VT.TestFunction()
# Time-stepping
t = 0

# write initial temperature
datafile.write("%f\t%f\n" % (t, T_n(mesh(0., 0.25))))
psi = ng.GridFunction(V, "pot")  # solution for potential
psi.Set(bc_cf_pot, ng.BND)
T_n1 = ng.GridFunction(VT, "T")  # solution for T
T_n1.Set(bc_cf_T, ng.BND)
ng.Draw(psi)
ng.Draw(T_n1)

# weak forms
# potential
a = ng.BilinearForm(V)
a += sigma_T(T_n) * ng.grad(psi_T) * ng.grad(v) * ng.dx
L = ng.LinearForm(V)
L += f * v * ng.dx  # keep in mind: f=0
# temperature
a_1 = ng.BilinearForm(VT)
a_1 += (rho * C_p * T_n1_T * vT * ng.dx
        + k_iso * dt * ng.grad(T_n1_T) * ng.grad(vT) * ng.dx)

L_1 = ng.LinearForm(VT)
L_1 += (rho * C_p * T_n + dt * sigma_T(T_n) * (ng.grad(psi) * ng.grad(psi))) * v * ng.dx

for n in range(num_steps):

    # assemble EM problem
    a.Assemble()
    L.Assemble()

    # solve EM problem
    res = L.vec.CreateVector()
    res.data = L.vec - a.mat * psi.vec
    psi.vec.data += a.mat.Inverse(V.FreeDofs()) * res

    # Update current time
    t += dt

    # assemble temperature problem
    L_1.Assemble()
    a_1.Assemble()

    # solve temperature problem
    res = L_1.vec.CreateVector()
    res.data = L_1.vec - a_1.mat * T_n1.vec
    T_n1.vec.data += a_1.mat.Inverse(VT.FreeDofs()) * res
    # Update previous solution
    T_n.Set(T_n1)

    # evaluate temperature at point of interest
    datafile.write("%f\t%f\n" % (t, T_n(mesh(0.5, 0.75))))
    ng.Redraw()

datafile.close()
hflux = -k_iso * ng.grad(T_n1)  # heat flux
# attemp to plot white arrows, does not work
ng.Draw(hflux, mesh, "flux")
