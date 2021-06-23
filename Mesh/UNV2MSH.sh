#!/bin/bash


unvfile=$1
dim=$2
if [ "$#" -ne 2 ]; then
	if [ "$#" -eq 0 ]; then
		echo "Script to convert UNV files to MSH files. Provide name of unvfile and dimension of mesh."
	elif [ "$#" -eq 1 ]; then
	       	echo "Please provide also dimension of mesh"
	fi
fi

if [[ $unvfile == *.unv ]]; then
	mshfile=${unvfile%.unv}.msh2
	gmsh "$unvfile" -"$dim" -v 0 -o "$mshfile"
	mshfilenew=${unvfile%.unv}.msh
	mv "$mshfile" "$mshfilenew"
else 
	echo "No .unv file provided!"	
fi
