#!/usr/bin/env python3

import dolfin as d
import sys

if len(sys.argv) != 2:
    print("specify name of the mesh file (omit .xml)")
    sys.exit(0)

meshname = sys.argv[1]
mesh = d.Mesh(meshname + ".xml")
cd = d.MeshFunction('size_t', mesh, meshname + "_physical_region.xml")
fd = d.MeshFunction('size_t', mesh, meshname + "_facet_region.xml")
hdf = d.HDF5File(mesh.mpi_comm(), meshname + ".h5", "w")
hdf.write(mesh, "/mesh")
hdf.write(cd, "/subdomains")
hdf.write(fd, "/facets")
print("Converted mesh to .h5 format")
