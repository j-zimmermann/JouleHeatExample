# !/bin/bash

FENICS='JouleHeating'
h=1
w=2
r='0.2'
mn='MeshDiskwithHole_test'

cd Mesh

SALOME='MeshwithHoleAutomatised'
# run SALOME with arguments: meshname, height, width, radius
salome -t "$SALOME".py args:$mn,$h,$w,$r
./MED2MSH.sh "$mn".med 2
./MSH2XML.sh "$mn".msh
./MeshtoHDF5_py3.py "$mn" 


cd ../FenicsStudy
python3 $FENICS.py $mn
pvpython ParaviewScript.py
convert -delay 5 -loop 0 movie.*.png movie.gif
