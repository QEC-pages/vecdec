#!/bin/bash
# -*- shell-script -*-
#
# analyze how measurement errors affect threshold (p1-p2 model)
#
# program variables to use 
vecdec=../src/vecdec
stim=../../Stim/out/stim

# place to put the data in
outfile=toric.dat

# program parameters
dmin=3 # minimum code distance 
dmax=7 # maximum code distance
p0=0.032 # error probability for j=0 
jmin=0 # range for calculating p1 and p2 (using p0/2**j )
jmax=6
Ntotal=$((1024*100)) # total number of steps to use

echo "# running toric.sh" > $outfile
echo "# depolarizing probability p1, measurement error p2" >> $outfile
echo "# columns: d p1 p2 Nfails Ntotal"
index=0 # block index to use in gnuplot

for (( d0=$dmin; d0<=$dmax; d0+=2 )) do # distance loop
    for (( j1=$jmin; j1<=$jmax; j1++ )) ; do #   p1 loop (depolarizing)
      p1=`awk "BEGIN { print ${p0}*exp(-${j1}*log(2)) } "`
      if ((j1==jmax)) ; then p1=0 ; fi
      echo '#' d=$d0 p1=$p1 index=$index >> $outfile         
      echo '#' d=$d0 p1=$p1 index=$index
      echo "# columns: d p1 p2 Nfails Ntotal" >> $outfile

      for (( j2=$jmin; j2<=$jmax; j2++ )) ; do # p2 loop (measurement)
        p2=`awk "BEGIN { print ${p0}*exp(-${j2}*log(2)) } "`
        fnam=surf_d$d0_$j1_$j2 # filename to use
        fnam=tmp # comment this line to keep all files 
        # now, generate the stim circuit
        $stim gen --code toric_code --task rotated_memory_x \
          --distance $d0 --rounds $d0 \
          --after_clifford_depolarization $p1 \
          --after_reset_flip_probability $p1 \
          --before_round_data_depolarization $p1 \
          --before_measure_flip_probability $p2 \
          > $fnam.stim
        # analyze errors 
        $stim analyze_errors --in $fnam.stim > $fnam.dem
        
        $vecdec debug=0 steps=$((d0*d0*d0*10)) swait=$((d0*d0*d0)) ntot=$Ntotal nvec=8192 nfail=1000 f=$fnam.dem > $fnam.out
        echo $d0 $p1 $p2 `cat $fnam.out` # show the output
        echo $d0 $p1 $p2 `cat $fnam.out` >> $outfile # save to big file        
      done # loop over p2
      # use gnuplot to fit the last block of data as a function of 'p2' 
      gnuplot 2> "tmp.out" <<EOF 
set fit quiet
set fit errorvariables
fit a+b*x "$outfile" index $index us (log(\$3)):(log(\$4/\$5)):(1/sqrt(\$4)) yerror via a,b
print $d0 ,$p1 ,a ,b, a_err, b_err 
EOF
      echo "# fit results: d p1 a b a_err b_err (index=$index)"
      echo "#" `cat tmp.out` 
      echo "# fit results: d p1 a b a_err b_err" >> $outfile
      echo "#" `cat tmp.out` >> $outfile
      echo "" >> $outfile
      echo "" >> $outfile
      index=$((index+1))
    done # loop over p1 
done # loop over distance 

    echo $outfile