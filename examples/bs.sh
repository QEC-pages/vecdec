#!/bin/bash
# -*- shell-script -*-
#
# Bacon-Shor codes using RW's script 'prog.py'
#
# create file "bs2.in"
cat <<EOF > bs2.in
# bs code d=2, separate X and Z measurements 
# layout: 
# 0  1  
# 2  3
4 3 4 # n, number of repetitions, maximum operations on one ancilla qubit 
# separate X and Z measurements 
X0 X1 _  _ # g1 	 
X2 X3 _  _ # g2	
_ _ Z0 Z2  # g3     
_ _ Z1 Z3  # g4

# stabilizer generators from gauge generators
! 1 2 # g1g2  XXXX 
! 3 4 # g3g4: ZZZZ

: Z0 Z1 # logical operator 
EOF

# program variables to use 
vecdec=../src/vecdec
stim=../../Stim/out/stim
python=python3.11

# place to put the data in
outfile=bs.dat

# program parameters
dmin=2 # minimum code distance 
dmax=2 # maximum code distance
p0=0.032 # error probability for j=0 
jmin=0 # range for calculating p1 and p2 (using p0/2**j )
jmax=6
Ntotal=$((1024*4)) # total number of steps to use

echo "# running bs.sh" > $outfile
echo "# running bs.sh for Bacon-Shor code"
echo "# error probability p1" >> $outfile
echo "# columns: d p1 Nfails Ntotal"
index=0 # block index to use in gnuplot

for (( d0=$dmin; d0<=$dmax; d0+=1 )) do # distance loop
    for (( j1=$jmin; j1<=$jmax; j1++ )) ; do #   p1 loop (depolarizing)
      p1=`awk "BEGIN { print ${p0}*exp(-${j1}*log(2)) } "`
      echo "# columns: d p1 Nfails Ntotal" >> $outfile

      fnam=tmp # comment this line to keep all files 
      # now, generate the stim circuit (notice how 'p1' is substituted)
      $python  <<EOF >$fnam.stim
# on the fly python script 
from prog import *
circuit = do_code("bs${d0}.in", $p1, 0)
print(circuit)
quit()
EOF
      # analyze errors 
      $stim analyze_errors --in $fnam.stim > $fnam.dem
      
      $vecdec debug=0 steps=$((d0*d0*d0)) nvec=$Ntotal f=$fnam.dem > $fnam.out
      echo $d0 $p1 `cat $fnam.out` # show the output
      echo $d0 $p1 `cat $fnam.out` >> $outfile # save to big file        
    done # loop over p1
    echo ""
done


