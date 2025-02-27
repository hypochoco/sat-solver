#!/bin/bash

E_BADARGS=65
if [ $# -ne 1 ]
then
	echo "Usage: `basename $0` <input>"
	exit $E_BADARGS
fi
	
input=$1

./src/build/project1 $input