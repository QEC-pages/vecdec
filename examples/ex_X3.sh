src/vecdec mode=0 fdet=tmp_det.01 fobs=tmp_obs.01 \
	   fdem= ./examples/surf_d3.dem nvec=2500 steps=1000 ntot=5000 \
	   perr=tmp_E.01 pdet=tmp_D.01 pobs=tmp_L.01 
det_mism=`paste tmp_D.01 tmp_det.01 -d" "| grep -E -c -v '^([01]+) \1$'`
obs_mism=`paste tmp_L.01 tmp_obs.01 -d" "| grep -E -c -v '^([01]+) \1$'`
echo mismatched det: $det_mism obs: $obs_mism
