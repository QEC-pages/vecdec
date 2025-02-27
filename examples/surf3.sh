#!/bin/bash
#
#! @brief surf3.sh - compare `vecdec` and pymatching decoders on surface codes.

# variables to use 
vecdec=../src/vecdec
stim=../../Stim/out/stim
pymatching=../../PyMatching/pymatching

# place to put the data in
outfile=surf_new.dat
detf=dets.01
obsf=obs.01

# program parameters
dmin=3 # minimum code distance 
dmax=7 # maximum code distance
p0=0.016 # error probability for j=0 
jmin=0 # range for calculating p1 and p2 (using p0/2**j )
jmax=12
#Ntotal=$((1024*128)) # total number of steps to use
Ntotal=$((1024*8)) 
nvec=$((1024*8))
factr=10 # NOTE: make it "1000" for an actual calculation ********************** !!!
#echo "# running surf3.sh" > $outfile
#echo "# same depolarizing probability and measurement error p1" >> $outfile
#echo "# columns: d p1 Nfails Ntotal"
#echo "# columns: d p1 Nfails Ntotal" >> $outfile
#index=0 # block index to use in gnuplot
echo "#" d p 'pF(vecdec) ' 
for (( ddd = $dmin; ddd <= $dmax; ddd += 2 )) ; do # distance loop
    fnam=surf_d$ddd # filename to use
    
    for (( j1=$jmin; j1<=$jmax; j1+=2 )) ; do #   p1 loop (depolarizing)
      p1=`awk "BEGIN { print ${p0}*exp(-${j1}*log(2)/2.0) } "`
      #      if ((j1==jmax)) ; then p1=0 ; fi
      # now, generate the stim circuit
      if (( 1 )) ; then 
        $stim gen --code surface_code --task rotated_memory_x \
          --distance $ddd --rounds $ddd \
          --after_clifford_depolarization $p1 \
          --after_reset_flip_probability $p1 \
          --before_round_data_depolarization $p1 \
          --before_measure_flip_probability $p1 \
          > $fnam.stim
      else # this is what is actually being run 
        $stim gen --code surface_code --task rotated_memory_x \
          --distance $ddd --rounds $ddd \
          --after_clifford_depolarization $p1 \
          > $fnam.stim
      fi
      # analyze errors 
      $stim analyze_errors --in $fnam.stim > $fnam.dem
      $stim analyze_errors --decompose_errors --in $fnam.stim > v$fnam.dem
      # simulate
      $stim  sample_dem --shots $Ntotal --in $fnam.dem \
        --out $detf --out_format 01 \
        --obs_out $obsf --obs_out_format 01
      # actually decode 
      $pymatching predict --dem $fnam.dem --in $detf --out p$obsf --in_format 01 --out_format 01
      pyF=`paste -d " " p$obsf $obsf | grep "0 1\|1 0" | wc -l`
      $pymatching predict --dem v$fnam.dem --in $detf --out p$obsf --in_format 01 --out_format 01
      pyFv=`paste -d " " p$obsf $obsf | grep "0 1\|1 0" | wc -l`            
      #      pyF=`paste -d " " p$obsf $obsf | grep "0 1\|1 0" | wc -l`
      echo "# pymatching :" $ddd $p1 `awk "BEGIN { print \"  \" $pyF \" \" $pyF/$Ntotal \"    \" $pyFv \" \" $pyFv/$Ntotal \" \" $Ntotal }"`
      # use the same sample in `vecdec`
      if ((ddd==3)) ; then steps_max=$((16*factr)) ;
      elif ((ddd==5)) ; then steps_max=$((64*factr)) ; # this is going to be very slow, but I want to know how bad `pymatching` actually is
      elif ((ddd==7)); then steps_max=$((256*factr)) ; fi
                           
      for (( steps=$((4*factr)) ; steps <= steps_max ; steps*= 2 )) ; do
        $vecdec debug=0 mode=0 fdet=$detf fobs=$obsf steps=$steps lerr=0 swait=0 \
          ntot=$Ntotal nvec=$((nvec)) nfail=0 f=$fnam.dem > $fnam.out
        echo "#" steps=$steps `awk '{print $1/$2 " " $1 " " $2 }' < $fnam.out`
      done                
      echo $ddd $p1 `awk '{print $1/$2 " " $1 }' < $fnam.out` \
        `awk "BEGIN { print  $pyF \" \" $pyF/$Ntotal \" \" $pyFv \" \" $pyFv/$Ntotal \" \" $Ntotal }"`      
    done # loop over p1
    echo ""
    echo ""
done # loop over distance 
echo "#" see results in $outfile
