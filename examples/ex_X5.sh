num=5000
src/vecdec mode=1.12 seed=7 fdet=tmp_det.01 fobs=tmp_obs.01 \
	   fdem= ./examples/surf_d3.dem nvec=$num steps=10 lerr=-1 ntot=$num \
	   perr=tmp_E.01 pdet=tmp_D.01 pobs=tmp_L.01 
paste -d: tmp_det.01 tmp_obs.01 | head -$num > tmp1.01
paste -d: tmp_D.01 tmp_L.01 > tmp2.01
echo failed: `paste tmp1.01 tmp2.01 -d" "| grep -E -v -c '^([01]+:[01]+) \1$'`
