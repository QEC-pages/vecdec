#!/bin/bash
# -*- shell-script -*-
#
# analyze how measurement errors affect threshold (p1-p2 model)
#
# program variables to use 
vecdec=../src/vecdec
stim=../../Stim/out/stim
pymatching=../../PyMatching/pymatching 
# place to put the data in
outfile=surf_est.dat

# program parameters
dmin=3 # minimum code distance 
dmax=9 # maximum code distance
p0=0.05 # error probability for j=0 
jmin=0 # range for calculating p1 and p2 (using p0/2**j )
jmax=15
Ntotal=$((1024*10)) # total number of steps to use
if (( $# > 0 )) ; then # we actually mean to do the calculation 
  echo "# running surf_est.sh" > $outfile
  echo "# depolarizing probability p1" >> $outfile
  echo "# using stim Ntotal=$Ntotal for comparison" >> $outfile 
  index=0 # block index to use in gnuplot
  
  for (( d0=$dmin; d0<=$dmax; d0+=2 )) ;  do # distance loop
      fnam=surf_d$d0 # filename to use
      detf=tmp_det.01
      obsf=tmp_obs.01
      
      echo "#" dist p1 estF wt LLR p_stim fails tot "idx=$index"    
      echo "#" dist p1 estF wt LLR p_stim fails tot "idx=$index" >> $outfile   
      for (( j1=$jmin; j1<=$jmax; j1++ )) ; do #   p1 loop (depolarizing)
        p1=`awk "BEGIN { print ${p0}*exp(-${j1}*log(2)/2.0) } "`
        #      if ((j1==jmax)) ; then p1=0 ; fi
        # now, generate the stim circuit
        if (( 1 )) ; then 
          $stim gen --code surface_code --task rotated_memory_x \
            --distance $d0 --rounds $d0 \
            --after_clifford_depolarization $p1 \
            --after_reset_flip_probability $p1 \
            --before_round_data_depolarization $p1 \
            --before_measure_flip_probability $p1 \
            > $fnam.stim
        else
          $stim gen --code surface_code --task rotated_memory_x \
            --distance $d0 --rounds 1 \
            --after_clifford_depolarization $p1 \
            > $fnam.stim
        fi
        # analyze errors 
        $stim analyze_errors --in $fnam.stim > $fnam.dem
        
        $stim  sample_dem --shots $Ntotal --in $fnam.dem \
          --out $detf --out_format 01 \
          --obs_out $obsf --obs_out_format 01
        
        $pymatching predict --dem $fnam.dem --in $detf --out p$obsf --in_format 01 --out_format 01
        pyOK=`paste -d " " p$obsf $obsf | grep "1 1\|0 0" | wc -l`
        pyF=`paste -d " " p$obsf $obsf | grep "0 1\|1 0" | wc -l`
        #      echo "# pymatching :" $pyOK fail $pyF
        $vecdec debug=0 steps=$((d0*d0*d0*d0)) mode=2 f=$fnam.dem > $fnam.out
        echo $d0 $p1 `cat $fnam.out` \
          `awk "BEGIN { print  $pyF/$Ntotal \" \" $pyF \" \" $pyF+$pyOK }"`
        
        echo $d0 $p1 `cat $fnam.out` \
          `awk "BEGIN { print  $pyF/$Ntotal \" \" $pyF \" \" $pyOK }"`\
          >> $outfile 
      done # loop over p1
      echo ""
      echo "" >> $outfile 
      echo "" >> $outfile 
      
  done # loop over distance 
else
  echo "# just the plot!"
fi
gnuplot <<EOF
set logs x
set logs y
plot "$outfile" index 0 us 2:3 tit "d=3 est" with lines lc "red"
replot "" index 0 us 2:6 tit "d=3 Stim est" with linespo lc "red"

replot "" index 1 us 2:3 tit "d=5 est" with lines lc "blue"
replot "" index 1 us 2:6 tit "d=5 Stim est" with linespo lc "blue"

replot "" index 2 us 2:3 tit "d=7 est" with lines lc "green"
replot "" index 2 us 2:6 tit "d=7 Stim" with linespo lc "green"

replot "" index 3 us 2:3 tit "d=9 est" with lines lc "black"
replot "" index 3 us 2:6 tit "d=9 Stim" with linespo lc "black"

set term push 
set term pdfcairo
set out "surf_est.pdf"
replot
set out 

EOF
    
# Binary variables, StDev(x)=sqrt((x2-x1^2/n)/n)=sqrt(p(1-p))
# delta(log(sum x/n))=delta(sum x/n)/(sum x/n)=
