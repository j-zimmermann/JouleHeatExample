"""
Simulate Joule (Resistive) Heating in FEniCS
Copyright (C) 2019  Julius Zimmermann <julius.zimmermann@uni-rostock.de>
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
import dolfin as d
import matplotlib.pyplot as plt
import sys
import os.path

'''
GMSH boundaries (from .geo or .msh file)
1 stands for Physical Line,
2 stands for Physical Surface,
3 stands for Physical Volume
(nomenclature taken from Computational Reality book of Abali)
use the auto file!

1 1 "Insulation_1"
1 2 "Insulation_2"
1 3 "Boundary_1"
1 4 "Boundary_2"
1 5 "CircleCenter"
2 6 "Copper"
'''

d.set_log_active(False)
meshpath = '../Mesh/'
if len(sys.argv) > 1:
    meshname = sys.argv[1]
else:
    meshname = 'MeshDiskwithHoleAuto'
    # meshname = meshname + '_refined'
mesh_file = meshpath + meshname + '.h5'
if not os.path.isfile(mesh_file):
    print("Specify correct name for meshfile!")
    sys.exit(0)

mesh = d.Mesh()
hdf = d.HDF5File(mesh.mpi_comm(), mesh_file, "r")
hdf.read(mesh, "/mesh", False)
cells = d.MeshFunction("size_t", mesh, dim=2)
hdf.read(cells, "/subdomains")
facets = d.MeshFunction("size_t", mesh, dim=1)
hdf.read(facets, "/facets")

datafileHDF5 = d.HDF5File(mesh.mpi_comm(), "TemperatureDevelopment_test.h5", "w")
datafile = open('TemperatureDevelopment_test.dat', 'w')


d.plot(mesh)
plt.show()


"""Let's start with some functions and parameters needed for the whole model"""
T_0 = d.Constant(300.)
T_ref = d.Constant(293.)
V0 = d.Constant(.1)  # set boundary value of electrode to 10 V
f = d.Constant(0.0)
rho_0 = d.Constant(1.754e-8)  # reference resistivity
C_p = d.Constant(340.)  # heat capacity
rho = d.Constant(8930.)  # density  [kg/m^3]
k_iso = d.Constant(384.)  # thermal conductivity [W/(m*K)]
alpha = d.Constant(3.9e-3)  # coefficient for T-dependent resistivity [1/K]

t_max = 2000.  # time in s
dt = 50.  # time step
num_steps = int(t_max / dt)


# temperature-dependent resistivity
def sigma_T(Temp):
    return 1. / (rho_0 * (1. + alpha * (Temp - T_0)))


# declare function space
V = d.FunctionSpace(mesh, 'CG', 2)
print("Number of DOFs: {}".format(V.dim()))
dx = d.dx

# boundary conditions
# see declaration from GMSH above
# Dirichlet for assigned electric potential
bc0_e = d.DirichletBC(V, V0, facets, 3)
bc1_e = d.DirichletBC(V, f, facets, 4)

# Dirichlet for assigned temperature
bc0_t = d.DirichletBC(V, T_0, facets, 3)
bc1_t = d.DirichletBC(V, T_0, facets, 4)
bc2_t = d.DirichletBC(V, T_0, facets, 5)  # circle

bc_e = [bc0_e, bc1_e]
bc_t = [bc0_t, bc1_t, bc2_t]

# Define initial value, which will later be the solution from the previous (n-th) time step
T_n = d.Function(V)
T_n = d.interpolate(T_ref, V)

# Define variational problem for temperature and potential
psi = d.TrialFunction(V)  # trial for potential
T_n1 = d.TrialFunction(V)  # trial for temperature
v = d.TestFunction(V)
# potential
a = d.inner(sigma_T(T_n) * d.grad(psi), d.grad(v)) * dx
L = f * v * dx  # keep in mind: f=0

psi = d.Function(V)  # solution for potential
# temperature
a_1 = (rho * C_p * T_n1 * v * dx
       + k_iso * dt * d.dot(d.grad(T_n1), d.grad(v)) * dx)
L_1 = (rho * C_p * T_n
       + dt * sigma_T(T_n) * d.dot(d.grad(psi), d.grad(psi))) * v * dx
T_n1 = d.Function(V)  # solution for T

# Time-stepping
t = 0
# write initial temperature
datafile.write("%f\t%f\n" % (t, T_n(d.Point(0., 0.25))))
datafileHDF5.write(T_n, "/T_0")
# open ParaView file and write initial T
file = d.File('JouleHeatingT.pvd')
file << T_n

for n in range(num_steps):
    # solve EM problem
    d.solve(a == L, psi, bc_e)
    # p = d.plot(psi)
    # plt.colorbar(p)
    # plt.show()
    # Update current time
    t += dt
    # Compute solution
    d.solve(a_1 == L_1, T_n1, bc_t)
    # d.plot(T_n1)
    # plt.colorbar(p)
    # plt.show()
    # Update previous solution
    T_n.assign(T_n1)
    # evaluate temperature at point of interest
    datafile.write("%f\t%f\n" % (t, T_n(d.Point(0., 0.25))))
    datafileHDF5.write(T_n, "/T_{}".format(n + 1))
    file << T_n

datafile.close()
p = d.plot(T_n1)
plt.colorbar(p, format='%3.3f')
hflux = d.project(-k_iso * d.grad(T_n), d.VectorFunctionSpace(mesh, 'Lagrange', 1))
# attemp to plot white arrows, does not work
d.plot(hflux, mode="glyphs", color="white")
plt.savefig("TemperatureFlux_test.png")
plt.show()
