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
dmax=7 # maximum code distance
p0=0.032 # error probability for j=0 
jmin=0 # range for calculating p1 and p2 (using p0/2**j )
jmax=15
Ntotal=$((1000*10)) # total number of steps to use
if (( $# > 0 )) ; then # we actually mean to do the calculation 
  echo "# running surf_est.sh" > $outfile
  echo "# depolarizing probability p1" >> $outfile
  echo "# using stim Ntotal=$Ntotal for comparison" >> $outfile 
  index=0 # block index to use in gnuplot
  
  for (( d0=$dmin; d0<=$dmax; d0+=2 )) ;  do # distance loop
      fnam=surf_d$d0 # filename to use
      detf=tmp_det.01
      obsf=tmp_obs.01

      p1=0.001
      $stim gen --code surface_code --task rotated_memory_x \
            --distance $d0 --rounds $d0 \
            --after_clifford_depolarization $p1 \
            --after_reset_flip_probability $p1 \
            --before_round_data_depolarization $p1 \
            --before_measure_flip_probability $p1 \
            > $fnam.stim

      $stim analyze_errors --in $fnam.stim > $fnam.dem

      n0=`grep -c -e "^[1-9]" $fnam.nz`;
      echo computing codewords for d=$d0
      for ((i=1; i<=25; i++)); do
	  $vecdec mode=2 steps=$((d0*d0*Ntotal)) fdem= $fnam.dem finC=${fnam}.nz outC=${fnam}.nz dW=2
	  nn=`grep -c -e "^[1-9]" ${fnam}.nz`;
	  if ((nn!=n0)); then
	      echo round $i of 25 d=$d0 n0=$n0 nn=$nn try again
	      n0=$nn;
	  else
	      echo done nn=$nn
	      break;
	  fi;
      done

      echo "#" dist p1 estF maxW=$d0 "idx=$index"    
      echo "#" dist p1 estF maxW=$d0 "idx=$index" >> $outfile   
      
      for (( j1=$jmin; j1<=$jmax+4; j1++ )) ; do #   p1 loop (depolarizing)
          p1=`awk "BEGIN { print ${p0}*exp(-${j1}*log(2)/2.0) } "`
	  
          $stim gen --code surface_code --task rotated_memory_x \
            --distance $d0 --rounds $d0 \
            --after_clifford_depolarization $p1 \
            --after_reset_flip_probability $p1 \
            --before_round_data_depolarization $p1 \
            --before_measure_flip_probability $p1 \
            > $fnam.stim
        # analyze errors 
          $stim analyze_errors --in $fnam.stim > $fnam.dem

	  out=`$vecdec debug=0 mode=2.2 steps=0 fdem= $fnam.dem finC=${fnam}.nz maxW=$d0`
	  echo $d0 $p1 $out >> $outfile
	  echo $d0 $p1 $out 	  	  
      done
      echo
      echo >> $outfile
      echo >> $outfile
      index=$((index+1))
      echo "#" dist p1 estF maxW=$((d0+1)) "idx=$index"    
      echo "#" dist p1 estF maxW=$((d0+1)) "idx=$index" >> $outfile   
      for (( j1=$jmin; j1<=$jmax+4; j1++ )) ; do #   p1 loop (depolarizing)
          p1=`awk "BEGIN { print ${p0}*exp(-${j1}*log(2)/2.0) } "`
	  
          $stim gen --code surface_code --task rotated_memory_x \
            --distance $d0 --rounds $d0 \
            --after_clifford_depolarization $p1 \
            --after_reset_flip_probability $p1 \
            --before_round_data_depolarization $p1 \
            --before_measure_flip_probability $p1 \
            > $fnam.stim
        # analyze errors 
          $stim analyze_errors --in $fnam.stim > $fnam.dem

	  out=`$vecdec debug=0 mode=2.2 steps=0 fdem= $fnam.dem finC=${fnam}.nz maxW=$((d0+1))`
	  echo $d0 $p1 $out >> $outfile
	  echo $d0 $p1 $out 	  	  
      done
      echo
      echo >> $outfile
      echo >> $outfile

      index=$((index+1))
      echo "#" dist p1 estF maxW=$((d0+2)) "idx=$index"    
      echo "#" dist p1 estF maxW=$((d0+2)) "idx=$index" >> $outfile   
      for (( j1=$jmin; j1 <= $jmax+4; j1++ )) ; do #   p1 loop (depolarizing)
          p1=`awk "BEGIN { print ${p0}*exp(-${j1}*log(2)/2.0) } "`
	  
          $stim gen --code surface_code --task rotated_memory_x \
            --distance $d0 --rounds $d0 \
            --after_clifford_depolarization $p1 \
            --after_reset_flip_probability $p1 \
            --before_round_data_depolarization $p1 \
            --before_measure_flip_probability $p1 \
            > $fnam.stim
        # analyze errors 
          $stim analyze_errors --in $fnam.stim > $fnam.dem

	  out=`$vecdec debug=0 mode=2.2 steps=0 fdem= $fnam.dem finC=${fnam}.nz maxW=$((d0+2))`
	  echo $d0 $p1 $out >> $outfile
	  echo $d0 $p1 $out 	  	  
      done
      echo
      echo >> $outfile
      echo >> $outfile

      if((d0==3)); then 
	  index=$((index+1))

	  echo "#" dist p1 fail_fract mode=0 "idx=$index"    
	  echo "#" dist p1 fail_fract mode=0 "idx=$index" >> $outfile   
	  for (( j1=$jmin; j1<=$jmax; j1++ )) ; do #   p1 loop (depolarizing)
              p1=`awk "BEGIN { print ${p0}*exp(-${j1}*log(2)/2.0) } "`
	      
              $stim gen --code surface_code --task rotated_memory_x \
		    --distance $d0 --rounds $d0 \
		    --after_clifford_depolarization $p1 \
		    --after_reset_flip_probability $p1 \
		    --before_round_data_depolarization $p1 \
		    --before_measure_flip_probability $p1 \
		    > $fnam.stim
              # analyze errors 
              $stim analyze_errors --in $fnam.stim > $fnam.dem
              shots=$((100*Ntotal))

	      out=`$vecdec debug=0 mode=0 ntot=$shots nfail=400 nvec=$((shots/10)) steps=100 lerr=1 fdem= $fnam.dem `
	      echo $d0 $p1 $out >> $outfile
	      echo $d0 $p1 $out 	  	  
	  done
	  echo
	  echo >> $outfile
	  echo >> $outfile
      fi
      
      index=$((index+1))            
      echo "#" dist p1 p_stim fails tot "idx=$index"    
      echo "#" dist p1 p_stim fails tot "idx=$index" >> $outfile   
      
      for (( j1=$jmin; j1<=$jmax; j1++ )) ; do #   p1 loop (depolarizing)
          p1=`awk "BEGIN { print ${p0}*exp(-${j1}*log(2)/2.0) } "`
	  
          $stim gen --code surface_code --task rotated_memory_x \
            --distance $d0 --rounds $d0 \
            --after_clifford_depolarization $p1 \
            --after_reset_flip_probability $p1 \
            --before_round_data_depolarization $p1 \
            --before_measure_flip_probability $p1 \
            > $fnam.stim
        # analyze errors 
          $stim analyze_errors --decompose_errors --in $fnam.stim > $fnam.dem
          shots=$((d0*d0*d0*d0*Ntotal))
          $stim  sample_dem --shots $shots --in $fnam.dem \
		 --out $detf --out_format 01 \
		 --obs_out $obsf --obs_out_format 01          

	  $pymatching predict --dem $fnam.dem --in $detf --out p$obsf --in_format 01 --out_format 01
          pyOK=`paste -d " " p$obsf $obsf | grep "1 1\|0 0" | wc -l`
          pyF=`paste -d " " p$obsf $obsf | grep "0 1\|1 0" | wc -l`
          #      echo "# pymatching :" $pyOK fail $pyF
	  #          $vecdec debug=0 steps=0 finC=$fnam.nz mode=2.2 f=$fnam.dem > $fnam.out
          echo $d0 $p1 `awk "BEGIN { print  $pyF/$shots \" \" $pyF \" \" $pyF+$pyOK }"`
        
          echo $d0 $p1 `awk "BEGIN { print  $pyF/$shots \" \" $pyF \" \" $pyF+$pyOK }"` >> $outfile
      done # loop over p1
      echo ""
      echo "" >> $outfile 
      echo "" >> $outfile
      index=$((index+1))
            
  done # loop over distance 
else
    echo "# just the plot!"
fi
gnuplot <<EOF
set logs x
set logs y
set key left top
plot [:][:1] "$outfile" index 0 us 2:3 tit "d=3 est maxW=3" with lines lc "red"
replot "" index 1  us 2:3 tit "d=3 est maxW=4" with lines lc "red" dashtype (3,2)
replot "" index 2  us 2:3 tit "d=3 est maxW=5" with lines lc "red" dashtype (4,2,1,2)
replot "" index 3  us 2:3 tit "d=3 mode=0" with lines lc "purple" 
replot "" index 4  us 2:3 tit "d=3 Stim" with linespo lc "red"
		   
replot "" index 5  us 2:3 tit "d=5 est maxW=5" with lines lc "blue"
replot "" index 6  us 2:3 tit "d=5 est maxW=6" with lines lc "blue" dashtype (3,2)
replot "" index 7  us 2:3 tit "d=5 est maxW=7" with lines lc "blue" dashtype (4,2,1,2)
replot "" index 8  us 2:3 tit "d=5 Stim est" with linespo lc "blue"
		   
replot "" index 9  us 2:3 tit "d=7 est maxW=7" with lines lc "green"
replot "" index 10 us 2:3 tit "d=7 est maxW=8" with lines lc "green" dashtype (3,2)
replot "" index 11 us 2:3 tit "d=7 est maxW=9" with lines lc "green" dashtype (4,2,1,2)
replot "" index 12 us 2:3 tit "d=7 Stim" with linespo lc "green"
		   
replot "" index 13 us 2:3 tit "d=9 est maxW=9" with lines lc "black"
replot "" index 14 us 2:3 tit "d=9 est maxW=10" with lines lc "black" dashtype (3,2)
replot "" index 15 us 2:3 tit "d=9 est maxW=11" with lines lc "black" dashtype (4,2,1,2)
replot "" index 16 us 2:3 tit "d=9 Stim" with linespo lc "black"

set term push 
set term pdfcairo
set out "surf_est.pdf"
replot
set out 
set term pop
EOF
    
# Binary variables, StDev(x)=sqrt((x2-x1^2/n)/n)=sqrt(p(1-p))
# delta(log(sum x/n))=delta(sum x/n)/(sum x/n)=
