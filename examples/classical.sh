prog=../src/iter_dec
code=PERMCODE1200.3.alist
outfile=tmp.out
p0=0.0$((3*1024))
debug=0;
for ((ii=0; ii <= 12; ii++)); do
    pp=`awk "BEGIN { print ${p0}*exp(-${ii}*log(2)*0.25) }"`
    echo ii $pp
    $prog debug=$debug finH=$code nvec=$((1000*1000)) nfail=1000 useP=$pp qllr2=0 > tmp.tmp
    cat tmp.tmp
    echo $pp `cat tmp.tmp` >> $outfile
    debug=0;
done
echo >> $outfile
echo >> $outfile

