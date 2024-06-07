#!/bin/bash
# -*- shell-script -*-
#

vecdec=../src/vecdec
stim=../../Stim/out/stim
mtx_sub=../src/mtx_sub
d=./tmp # directory to use for tmp files 
mkdir $d
fH=${d}/tmpH.mmx

# generate the original DEM and OBS events (fixed error probability for now)
$vecdec fdem=surf_d3.dem steps=0 gdet=tmp/bigD.01 gobs=tmp/bigL.01 useP=0.002

# create original matrices H,L,P
$vecdec mode=3.28 fout=${d}/tmp fdem=surf_d3.dem

# establish block boundaries 
$mtx_sub debug=0  fin= tmp/tmpH.mmx rows=0,3 rows=4,11 rows=12,19 rows=20,23
# this gives the output (col ranges for each row block)
# 0 164
# 1 267
# 39 290
# 170 290
#
# cut two matrices: A0, A1 (diagonal blocks)
# 
# A0: rows=0,11 cols=0,267
$mtx_sub debug=0 fin=tmp/tmpH.mmx minR=0 maxR=11 minC=0 maxC=267 out=tmp/mA0
# B1: rows=12,23 cols=0,38 -- this matrix is empty 
$mtx_sub debug=0 fin=tmp/tmpH.mmx minR=12 maxC=38 out=tmp/mB1
# A1: rows=12,23 cols=39,290
$mtx_sub debug=0 fin=tmp/tmpH.mmx minR=12 minC=39 out=tmp/mA1

# prepare the first slice for decoding with A0 (cols from 1):
cut -c 1-12 tmp/bigD.01 > tmp/oneD.01
# do the actual decoding.  This creates the 268-col file perrA.01
$vecdec debug=0 steps=1000 useP=0.002 finH=tmp/mA0.mtx fdet=tmp/oneD.01 perr=tmp/perrA.01
# cut first 39 cols to final output 
cut -c 1-39 tmp/perrA.01 > tmp/poneA.01

# prepare the second slice for decoding with A1 
cut -c 13- tmp/bigD.01 > tmp/twoD.01
$vecdec debug=0 steps=1000 useP=0.002 finH=tmp/mA1.mtx fdet=tmp/twoD.01 perr=tmp/perrB.01
# combine for total predicted errors
paste --delimiters="" tmp/poneA.01 tmp/perrB.01 >tmp/perr.01

# now verify the decoding 
$vecdec debug=0 steps=0 useP=0.002 fdem=surf_d3.dem ferr=tmp/perr.01 gdet=tmp/outD.01 gobs=tmp/outL.01 
#diff tmp/bigD.01 tmp/outD.01 |grep -c -e "^>"
diff tmp/bigD.01 tmp/outD.01
