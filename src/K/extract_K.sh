#!/bin/bash -f

gawk '
{
	match($0, /Isothermal Compressibility Kappa\ += ([1-9].[0-9]+e[-+][0-9]+)/, arr)  # for scientific format
	K = arr[1]
	if(length(K) > 1)
		print 1/K
}
' < $1

#	match($0, /Isothermal Compressibility Kappa\ +=\ +([0-9]+.[0-9]+)/, arr)
