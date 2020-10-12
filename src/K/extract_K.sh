#!/bin/bash -f

gawk '
{
	match($0, /Isothermal Compressibility Kappa\ += ([1-9].[0-9]+e[-+][0-9]+)/, arr)
	K = arr[1]
	if(length(K) > 1)
		print K
}
' < $1
