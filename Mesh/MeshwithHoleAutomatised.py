# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v8.3.0 with dump python functionality
###

import sys
import os
import salome

salome.salome_init()
theStudy = salome.myStudy

import salome_notebook
notebook = salome_notebook.NoteBook(theStudy)
cur_path = os.getcwd()
print("Current path is {}".format(cur_path))
sys.path.insert(0, cur_path)

####################################################
##       Begin of NoteBook variables section      ##
####################################################

if len(sys.argv) != 5:
    print("Error - provide input as: meshname,height,width,radius")

# h = HEIGHT

meshname = sys.argv[1]
h = sys.argv[2]
w = sys.argv[3]
r = sys.argv[4]
# w = WIDTH
# r = RADIUS
notebook.set("h", h)
notebook.set("w", w)
notebook.set("r", r)
####################################################
##        End of NoteBook variables section       ##
####################################################
###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New(theStudy)

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Face_1 = geompy.MakeFaceHW("h", "w", 1)
Disk_1 = geompy.MakeDiskR("r", 1)
Cut_1 = geompy.MakeCutList(Face_1, [Disk_1], True)
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["VERTEX"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["VERTEX"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["EDGE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
listSubShapeIDs = geompy.SubShapeAllIDs(Cut_1, geompy.ShapeType["FACE"])
Insulation_1 = geompy.CreateGroup(Cut_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Insulation_1, [6])
Insulation_2 = geompy.CreateGroup(Cut_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Insulation_2, [8])
Boundary_1 = geompy.CreateGroup(Cut_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Boundary_1, [3])
Boundary_2 = geompy.CreateGroup(Cut_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Boundary_2, [10])
CircleCenter = geompy.CreateGroup(Cut_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(CircleCenter, [12])
Copper = geompy.CreateGroup(Cut_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(Copper, [1])
geompy.addToStudy(O, 'O')
geompy.addToStudy(OX, 'OX')
geompy.addToStudy(OY, 'OY')
geompy.addToStudy(OZ, 'OZ')
geompy.addToStudy(Face_1, 'Face_1')
geompy.addToStudy(Disk_1, 'Disk_1')
geompy.addToStudy(Cut_1, 'Cut_1')
geompy.addToStudyInFather(Cut_1, Insulation_1, 'Insulation_1')
geompy.addToStudyInFather(Cut_1, Insulation_2, 'Insulation_2')
geompy.addToStudyInFather(Cut_1, Boundary_1, 'Boundary_1')
geompy.addToStudyInFather(Cut_1, Boundary_2, 'Boundary_2')
geompy.addToStudyInFather(Cut_1, CircleCenter, 'CircleCenter')
geompy.addToStudyInFather(Cut_1, Copper, 'Copper')

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New(theStudy)
Mesh_1 = smesh.Mesh(Cut_1)
NETGEN_1D_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize( 0.031 )
NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
NETGEN_2D_Parameters_1.SetOptimize( 1 )
NETGEN_2D_Parameters_1.SetFineness( 4 )
NETGEN_2D_Parameters_1.SetMinSize( 0.0065181 )
NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )
Insulation_1_1 = Mesh_1.GroupOnGeom(Insulation_1,'Insulation_1',SMESH.EDGE)
Insulation_2_1 = Mesh_1.GroupOnGeom(Insulation_2,'Insulation_2',SMESH.EDGE)
Boundary_1_1 = Mesh_1.GroupOnGeom(Boundary_1,'Boundary_1',SMESH.EDGE)
Boundary_2_1 = Mesh_1.GroupOnGeom(Boundary_2,'Boundary_2',SMESH.EDGE)
CircleCenter_1 = Mesh_1.GroupOnGeom(CircleCenter,'CircleCenter',SMESH.EDGE)
Copper_1 = Mesh_1.GroupOnGeom(Copper,'Copper',SMESH.FACE)
isDone = Mesh_1.Compute()
[ Insulation_1_1, Insulation_2_1, Boundary_1_1, Boundary_2_1, CircleCenter_1, Copper_1 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
    Mesh_1.ExportMED(cur_path + '/' + meshname + '.med', auto_groups=0, minor=32, overwrite=1, meshPart=None, autoDimension=1)
    print('Mesh exported')
    pass
except:
    print('ExportToMEDX() failed. Invalid file name?')
smesh.SetName(Mesh_1, 'Mesh_1')

## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')
smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')
smesh.SetName(Copper_1, 'Copper')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Insulation_1_1, 'Insulation_1')
smesh.SetName(Boundary_1_1, 'Boundary_1')
smesh.SetName(Insulation_2_1, 'Insulation_2')
smesh.SetName(CircleCenter_1, 'CircleCenter')
smesh.SetName(Boundary_2_1, 'Boundary_2')


if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser(True)
