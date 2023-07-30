#!/bin/bash
# -*- shell-script -*-
#
# try to estimate the value of Lambda 3/5 using `stim` and `vecdec` 
#

# For a run with `n` decoding rounds, fit logical error rate to exponential form $P_f=1-(1-2*eps)^n$.
# Calculate `eps(d)` for `d=3` and `d=5`.  Then we have `Lambda=eps(3)/eps(5)`.
#

# program variables to use 
vecdec=../src/vecdec
stim=../../Stim/out/stim
pymatching=../../PyMatching/pymatching 
# place to put the data in
outfile=lambda.dat

runstim=1
runvecd=1

# program parameters
d1=3 # minimum code distance for Lambda 
d2=5 # maximum code distance for Lambda 
p0=0.002 # error probability for j=0 
jmin=0 # range for calculating p1 and p2 (using p0/2**j )
jmax=0 # just one value
# nmin=dmin 
n_factor=2
Ntotal=$((64*1024*1024)) # total number of steps to use
echo "# running lambda.sh" > $outfile
echo "# depolarizing and measurement probability p1" >> $outfile
echo "# using stim Ntotal=$Ntotal for comparison" >> $outfile 

for (( j1=jmin; j1<=jmax; j1++)) ; do # p1 loop
  p1=`awk "BEGIN { print ${p0}*exp(-${j1}*log(2)/2.0) } "`
  fnam=surf_j$j1 # filename to use -- separate for each `j`
  if (( runstim )) ; then 
    echo "# output of 'lambda.sh' running stim/pymatching sim" > ${fnam}.out
  fi
  if (( runvecd )) ; then 
    echo "# output of 'lambda.sh' running vecdec mode=2 est" > v${fnam}.out
  fi
  for d0 in $d1 $d2 ; do 
    for (( n=$d2 ; n<=$((d2*n_factor)); n++ )) ; do # n-loop
      detf=tmp_det.01
      obsf=tmp_obs.01
        
      #        echo "#" dist p1 estF wt LLR p_stim fails tot "idx=$index"    
      #        echo "#" dist p1 estF wt LLR p_stim fails tot "idx=$index" >> $outfile
      if (( 0 )) ; then
        $stim gen --code surface_code --task rotated_memory_x \
          --distance $d0 --rounds $n \
          --after_clifford_depolarization $p1 \
          --after_reset_flip_probability $p1 \
          --before_round_data_depolarization $p1 \
          --before_measure_flip_probability $p1 \
          > $fnam.stim
      else
        $stim gen --code surface_code --task rotated_memory_x \
          --distance $d0 --rounds $n \
          --after_reset_flip_probability $p1 \
          --before_round_data_depolarization $p1 \
          --before_measure_flip_probability $p1 \
          > $fnam.stim
      fi

      if (( runstim )) ; then # run pymatching  

        # analyze errors
        $stim analyze_errors --in $fnam.stim > v$fnam.dem
        $stim  sample_dem --shots $((Ntotal)) --in v$fnam.dem \
          --out $detf --out_format 01 \
          --obs_out $obsf --obs_out_format 01
        $stim analyze_errors --decompose_errors --in $fnam.stim > $fnam.dem
                
        $pymatching predict --dem $fnam.dem --in $detf --out p$obsf \
          --in_format 01 --out_format 01 --in $fnam.stim > $fnam.dem_dec 
        
        pyOK=`paste -d " " p$obsf $obsf | grep "1 1\|0 0" | wc -l`
        pyF=`paste -d " " p$obsf $obsf | grep "0 1\|1 0" | wc -l`
        #      echo "# pymatching :" $pyOK fail $pyF
        #  $vecdec debug=0 steps=$((d0*d0*d0*d0)) mode=2 f=$fnam.dem > $fnam.out
        echo $n $d0 $p1 `awk "BEGIN { print  $pyF/$Ntotal \" \" $pyF \" \" $pyF+$pyOK }"`
        echo $n $d0 $p1 `awk "BEGIN { print  $pyF/$Ntotal \" \" $pyF \" \" $pyF+$pyOK }"` >> ${fnam}.out

      fi 
      if (( runvecd )) ; then # run vecdec 
        $stim analyze_errors --in $fnam.stim > v$fnam.dem
        $vecdec debug=0 steps=$((1000*d0*d0*n*n)) swait=$((10*d0*d0*n*n)) mode=2 f=v$fnam.dem > tmp.out 
        echo $n $d0 $p1 `cat tmp.out` >> v$fnam.out
        echo $n $d0 $p1 `cat tmp.out`
      fi
    done # n loop
    if (( runstim )); then 
      echo "" >> ${fnam}.out 
      echo "" >> ${fnam}.out
    fi
    if (( runvecd )) ; then 
      echo "" >> v${fnam}.out 
      echo "" >> v${fnam}.out
    fi
    echo ""
  done # d0 loop
  gnuplot <<EOF 
    set fit quiet
    set fit errorvariables
    fit a3+b3*x "${fnam}.out" index 0 us 1:4:(1/sqrt(\$6)) yerror via a3,b3
    fit a5+b5*x "" index 1 us 1:4:(1/sqrt(\$6)) yerror via a5,b5  
    fit va3+vb3*x "v${fnam}.out" index 0 us 1:4:(1/sqrt(\$6)) yerror via va3,vb3
    fit va5+vb5*x "" index 1 us 1:4:(1/sqrt(\$6)) yerror via va5,vb5  
    
    set logs x
    set logs y
    set style data linespo
    plot "${fnam}.out" index 0 us 1:4 tit "d=$d1" with linespoi lc "red"
    replot "" index 1 us 1:4 tit "d=$d2" with linespoi lc "blue"
    replot "v${fnam}.out" index 0 us 1:4 tit "d=$d1 vecdec" with linespoi lc "purple"
    replot "" index 1 us 1:4 tit "d=$d2 vecdec" with linespoi lc "green"
    replot a3+b3*x notit lc "red"
    replot a5+b5*x notit lc "blue"
    replot va3+vb3*x notit lc "purple"
    replot va5+vb5*x notit lc "green"
    print "d=3", $p1 ,a3 ,b3, a3_err, b3_err
    print "d=5", $p1 ,a5 ,b5, a5_err, b5_err
    print "d=3", $p1 ,va3 ,vb3, va3_err, vb3_err, " est"
    print "d=5", $p1 ,va5 ,vb5, va5_err, vb5_err, " est"
    print "lambda(stim)=",b3/b5, " lambda(est)=",vb3/vb5
    print "" 
    set term push
    set term pdfcairo
    set out "tmp.pdf
    replot
    set out
    set term pop
EOF
  if (( 0 )) ; then 
    gnuplot 2> "tmp.out" <<EOF 
    set fit quiet
    set fit errorvariables
    fit a3+b3*x "$fnam.out" index 0 us (log(\$1)):(log(\$4)):(1/sqrt(\$6)) yerror via a3,b3
    fit a5+b5*x "" index 1 us (log(\$1)):(log(\$4)):(1/sqrt(\$6)) yerror via a5,b5

    set logs x
    set logs y
    set style data linespo
    plot "" index 0 us 1:4 tit "d=$d1" lc "red"
    replot "" index 1 us 1:4 tit "d=$d2" lc "blue"
    replot exp(a3+b3*log(x)) notit lc "red"
    replot exp(a5+b5*log(x)) notit lc "blue"
    print $p1 ,a3 ,b3, a3_err, b3_err

    set term push
    set term pdfcairo
    set out "tmp.pdf
    replot
    set out
    set term pop
EOF
  echo "# fit results: d p1 a b a_err b_err (index=$index)"
  echo "#" `cat tmp.out` 
  echo "# fit results: d p1 a b a_err b_err" >> $outfile
  echo "#" `cat tmp.out` >> $outfile
  echo "" >> $outfile
  echo "" >> $outfile
  index=$((index+1))
  fi
  

done # p1 loop 


