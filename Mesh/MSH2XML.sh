#!/bin/bash

mshfile=$1
if [ "$#" -ne 1 ]; then
	echo "Script to convert MSH files to dolfin XML files. Provide name of mshfile."
fi

if [[ $mshfile == *.msh ]]; then
	xmlfile=${mshfile%.msh}.xml
	dolfin-convert "$mshfile" "$xmlfile"
else 
	echo "No .msh file provided!"	
fi
