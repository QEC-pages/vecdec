#!/bin/bash
# -*- shell-script -*-
#
# use 01 simulation data from Stim 
#
# program variables to use 
vecdec=./vecdec
stim=../../Stim/out/stim

# place to put the data in
outfile=surf_new.dat
detf=dets.01
obsf=obs.01

# program parameters
dmin=3 # minimum code distance 
dmax=3 # maximum code distance
p0=0.004 # error probability for j=0 
jmin=0 # range for calculating p1 and p2 (using p0/2**j )
jmax=0
Ntotal=$((16)) # total number of steps to use

echo "# running surf3.sh" > $outfile
echo "# same depolarizing probability and measurement error p1" >> $outfile
echo "# columns: d p1 Nfails Ntotal"
echo "# columns: d p1 Nfails Ntotal" >> $outfile
index=0 # block index to use in gnuplot

for (( d0=$dmin; d0<=$dmax; d0+=2 )) do # distance loop
    fnam=surf_d$d0 # filename to use
    
    for (( j1=$jmin; j1<=$jmax; j1++ )) ; do #   p1 loop (depolarizing)
      p1=`awk "BEGIN { print ${p0}*exp(-${j1}*log(2)/2.0) } "`
      #      if ((j1==jmax)) ; then p1=0 ; fi
      # now, generate the stim circuit
      $stim gen --code surface_code --task rotated_memory_x \
        --distance $d0 --rounds $d0 \
        --after_clifford_depolarization $p1 \
        --after_reset_flip_probability $p1 \
        --before_round_data_depolarization $p1 \
        --before_measure_flip_probability $p1 \
        > $fnam.stim
      # analyze errors 
      $stim analyze_errors --in $fnam.stim > $fnam.dem
      # simulate
      $stim  sample_dem --shots $Ntotal --in $fnam.dem \
        --out $detf --out_format 01 \
        --obs_out $obsf --obs_out_format 01
      
      $vecdec debug=0 mode=1 fdet=$detf fobs=$obsf steps=$((d0*d0*d0)) lerr=1 swait=$((d0*d0)) \
        ntot=$Ntotal nvec=$((Ntotal)) nfail=5000 f=$fnam.dem > $fnam.out
      echo $d0 $p1 `cat $fnam.out` # show the output
      echo $d0 $p1 `cat $fnam.out` >> $outfile # save to big file
      
    done # loop over p1
    if (( 0 )) ; then
      gnuplot 2> "tmp.out" <<EOF 
set fit quiet
set fit errorvariables
fit [:log(0.007)] a+b*x "$outfile" index $index us (log(\$2)):(log(\$3/\$4)):(1/sqrt(\$3)) yerror via a,b
print $d0 ,a ,b, a_err, b_err 
EOF
      echo "# fit results: d a b a_err b_err (index=$index)"
      echo "#" `cat tmp.out` 
      echo "# fit results: d a b a_err b_err" >> $outfile
      echo "#" `cat tmp.out` >> $outfile
      echo "" >> $outfile
      echo "" >> $outfile
      index=$((index+1))
    fi
done # loop over distance 

    echo $outfile

# Binary variables, StDev(x)=sqrt((x2-x1^2/n)/n)=sqrt(p(1-p))
# delta(log(sum x/n))=delta(sum x/n)/(sum x/n)=
