prog=../src/vecdec
prog_dbl=../src/vecdec_dbl
code=PERMCODE1200.3.alist
outfile=tmp.out
p0=0.0$((3*1024))
debug=0;
echo "# running vecdec on $code " > $outfile
for((dd2=-1; dd2<=2; dd2++)); do
    if((dd2==-1)) ; then 
	program=$prog_dbl
	d2=0
    else
	program=$prog
	d2=$((300*dd2))
    fi
    echo "# use qllr2=$d2 $program" >> $outfile
    echo "# $dd2 $ii $pp d2=$d2 $program"
    for ((ii=0; ii <= 11; ii++)); do
	pp=`awk "BEGIN { print ${p0}*exp(-${ii}*log(2)*0.25) }"`
	$program debug=0 mode=1.1 finH=$code ntot=10000 nvec=$((100*10)) nfail=1000 useP=$pp qllr2=$d2 > tmp.tmp
	echo $pp `cat tmp.tmp`
	echo $pp `cat tmp.tmp` >> $outfile
	debug=0;
    done
    echo >> $outfile
    echo >> $outfile
done
