#!/usr/bin/bash

##
## script for matrix modification.  Eliminate all weight-3 columns from G.
##  
##
##

vecdec=../src/vecdec 
fdem=surf_d3.dem
#fdem=../input/rep_code5.dem
nam=try # for matrices: try0H.mmx, try1H.mmx, etc.
cws=tmp # for codewords: tmp0.nz, tmp1.nz, ...

# read the DEM model and write the initial (level=0) matrices
$vecdec fdem=$fdem mode=3.28 fout=${nam}0

for (( level=0; level<1; level++ )) ; do
    echo "#######################" starting $level "##############################"
    $vecdec mode=2 finH=${nam}${level}H.mmx finL=${nam}${level}L.mmx finP=${nam}${level}P.mmx outC=${cws}${level}.nz steps=100000
done

