# stim=../Stim/out/stim
dist=5
p1=0.003
# # construct `stim` file
# $stim gen --code surface_code --task rotated_memory_x \
#       --distance $dist --rounds $dist \
#       --after_clifford_depolarization $p1 \
#       --after_reset_flip_probability $p1 \
#       --before_round_data_depolarization $p1 \
#       --before_measure_flip_probability $p1 \
#       > examples/surf_d$dist.stim
# # generate the corresponding DEM file 
# $stim analyze_errors --in examples/surf_d${dist}.stim > examples/surf_d${dist}.dem
echo uW=2 dist=$dist p1=$p1
./src/vecdec debug=1 seed=7 mode=0 steps=5000 lerr=1 fdem=examples/surf_d${dist}.dem nvec=1000 ntot=1000 uW=2


