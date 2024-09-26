
stim=../Stim/out/stim # location of `stim` command-line binary
$stim sample_dem --shots 5000 --in ./examples/surf_d3.dem \
	--out tmp_det.01 --out_format 01 \
	--obs_out tmp_obs.01 --obs_out_format 01
src/vecdec mode=0 fdet=tmp_det.01 fobs=tmp_obs.01 fdem= ./examples/surf_d3.dem nvec=2500 steps=1000 ntot=5000 
