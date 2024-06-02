#!/bin/bash
#
vecdec=../src/vecdec
Hx=../input/QX40.mtx
Hz=../input/QZ40.mtx
out=../input/try_mode2.dat
cws=../input/tmpQ40.nz
# program parameters
p0=0.032 # error probability for j=0 
jmin=0 # range for calculating p1 and p2 (using p0/2**j )
jmax=16
#Ntotal=$((1024*128)) # total number of steps to use
Ntotal=$((1024*1024*10)) 
nvec=$((1024*32))
m0steps=2000
maxW=12

if (( 0 )) ;  then 
    # generate the codewords
    cnt0=0
    echo "# generating codewords for finH=$Hx finG=$Hz steps=100000"
    for((i=1; i<=10; i++)) ; do
	$vecdec debug=0 mode=2 finH=$Hx finG=$Hz steps=100000 useP=0.01 maxW=$maxW finC=$cws outC=$cws dW=-1
	cnt=`grep -c -e "^[1-9]" $cws`
	echo i=$i $cnt0 $cnt 
	if((cnt==cnt0)) ; then
	    echo $cnt 
	    break;
	else
	    cnt0=cnt;
	fi
    done
    
    echo "#" estimated fail probability mode=2.2 Hx=$Hx Hz=$Hz cnt=$cnt maxW=$maxW > $out
    for (( j1=$jmin; j1<=$jmax; j1+=1 )) ; do 
	p1=`awk "BEGIN { print ${p0}*exp(-${j1}*log(2)/2.0) } "`
	str=`$vecdec debug=0 mode=2.2 finH=$Hx finG=$Hz useP=$p1 steps=0 maxW=$maxW finC=$cws`
	echo -e $p1 "\t" `$vecdec debug=0 mode=2.2 finH=$Hx finG=$Hz useP=$p1 steps=0 maxW=4  finC=$cws` maxW=4
	echo -e "   \t" `$vecdec debug=0 mode=2.2 finH=$Hx finG=$Hz useP=$p1 steps=0 maxW=6  finC=$cws` maxW=6
	echo -e "   \t" `$vecdec debug=0 mode=2.2 finH=$Hx finG=$Hz useP=$p1 steps=0 maxW=8  finC=$cws` maxW=8
	echo -e "   \t" `$vecdec debug=0 mode=2.2 finH=$Hx finG=$Hz useP=$p1 steps=0 maxW=10 finC=$cws` maxW=10 
	echo -e $p1 "\t" $str maxW=$maxW
	echo 
	echo -e $p1 "\t" $str maxW=$maxW >> $out
    done
    echo >> $out
    echo >> $out 
    echo 
    
    echo "#" computed fail probability mode=0 Hx=$Hx Hz=$Hz ntot=$Ntotal nvec=$nvec >> $out
    for (( j1=$jmin; j1<=$jmax; j1+=1 )) ; do 
	p1=`awk "BEGIN { print ${p0}*exp(-${j1}*log(2)/2.0) } "`
	str=`$vecdec debug=0 mode=0 finH=$Hx finG=$Hz useP=$p1 steps=$m0steps ntot=$Ntotal nvec=$nvec lerr=1 nfail=1000`
	grep $p1 $out 
	echo -e $p1 "\t" $str
	echo $p1 $str >> $out
    done
    echo
fi
echo >> $out
echo >> $out 
echo # now try plotting \"$out\"

for (( w=4 ; w<=6 ; w++ )); do 
    echo "#" estimated fail probability mode=2.2 Hx=$Hx Hz=$Hz cnt=$cnt maxW=$w >> $out
    echo "#" estimated fail probability mode=2.2 Hx=$Hx Hz=$Hz cnt=$cnt maxW=$w
    for (( j1=$jmin; j1<=$jmax; j1+=1 )) ; do 
	p1=`awk "BEGIN { print ${p0}*exp(-${j1}*log(2)/2.0) } "`
	str=`$vecdec debug=0 mode=2.2 finH=$Hx finG=$Hz useP=$p1 steps=0 maxW=$w finC=$cws`
	echo -e $p1 "\t" $str maxW=$w	
	echo -e $p1 "\t" $str maxW=$w >> $out
    done
    echo >> $out 
    echo >> $out
    echo
done
    
