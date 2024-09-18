src/vecdec mode=0 fdet=tmp_det.01 fobs=tmp_obs.01 fdem= ./examples/surf_d3.dem nvec=2500 steps=1000 ntot=50000 \
  perr=tmp_E.01 pdet=tmp_D.01 pobs=tmp_L.01 
det_mism=`diff  tmp_D.01 tmp_det.01 | grep -c -e "^<"`
obs_mism=`diff  tmp_L.01 tmp_obs.01 | grep -c -e "^<"`
echo mismatch detector: $det_mism obs: $obs_mism
