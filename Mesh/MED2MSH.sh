#!/bin/bash



medfile=$1
dim=$2
if [ "$#" -ne 2 ]; then
	if [ "$#" -eq 0 ]; then
		echo "Script to convert MED files to MSH files. Provide name of medfile and dimension of mesh."
	elif [ "$#" -eq 1 ]; then
	       	echo "Please provide also dimension of mesh"
	fi
fi

if [[ $medfile == *.med ]]; then
	mshfile=${medfile%.med}.msh2
	gmsh "$medfile" -"$dim" -v 0 -o "$mshfile"
	mshfilenew=${medfile%.med}.msh
	mv "$mshfile" "$mshfilenew"
else 
	echo "No .med file provided!"	
fi
