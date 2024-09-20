num=1000
# prepare the matrices end probability vectors (extension '.mtx' added)
src/mtx_sub fin=bigH.mmx maxR=83 maxC=1401 out=tryA1
src/mtx_sub fin=bigH.mmx minR=36 maxC=585  out=tryB2
src/mtx_sub fin=bigH.mmx minR=36 minC=586  out=tryA2
# split the probability vector (this is a hack)
# P1
grep %% bigP.mmx > tryP1.mmx # Matrix Market header
echo 1 1402 >> tryP1.mmx
tail +4 bigP.mmx | head -1402 >> tryP1.mmx
# P2
grep %% bigP.mmx > tryP2.mmx # Matrix Market header
echo 1 $((1679-586)) >> tryP2.mmx
tail -$((1679-586)) bigP.mmx >> tryP2.mmx
#
# generate the `det` and `obs` events (big file)
src/vecdec debug=0 gdet=big_det.01 gobs=big_obs.01 finH=bigH.mmx finP=bigP.mmx finL=bigL.mmx \
	   mode=0 steps=0 uW=-1 ntot=$num nvec=$num seed=7
#
# split `det` events
cut -c 1-84 big_det.01 > try_D1.01
cut -c 37-120 big_det.01 > try_D2.01
#
# first-step decoding
src/vecdec debug=0 mode=0 steps=1000 lerr=1 uW=2 \
	   fdet=try_D1.01 finH=tryA1.mtx finP=tryP1.mmx  \
	   ntot=$num nvec=$num seed=7 perr=try_E1.01
#
# cut out the first 586 cols of e1 
cut -c 1-586 try_E1.01 > try_E1x.01
#
# second step decoding
src/vecdec debug=0 mode=0 fdet=try_D2.01 finH=tryA2.mtx finP=tryP2.mmx  \
	   steps=1000 lerr=1 uW=2 ntot=$num nvec=$num seed=7 perr=try_E2.01 \
	   fer0=try_E1x.01 finA=tryB2.mtx
# create the combined error vector
paste try_E1x.01 try_E2.01 -d "" > out_E.01
# create the updated det file for syndrome verification (should be all-zero)
src/vecdec debug=0 mode=0 finH=bigH.mmx ntot=$num nvec=$num steps=0 uW=-1 useP=0.05 \
	   fdet=big_det.01 fer0=out_E.01 finA=bigH.mmx gdet=yyy.01
# create the obs vector
src/vecdec debug=0 mode=0 finH=bigH.mmx finL=bigL.mmx ntot=$num nvec=$num steps=0 uW=-1 useP=0.05 \
	   ferr=out_E.01 gobs=xxx.01
#
# output success count
paste -d: big_det.01 big_obs.01 | head -$num > tmp1.01
paste -d: yyy.01 xxx.01 > tmp2.01
#
echo two-step decoding failed: det `grep -c 1 yyy.01` \
     obs `paste big_obs.01 xxx.01 -d" "| grep -E -v -c '^([01]+) \1$'` out of $num 
#
# compare with full-block decoding (use more steps here since matrix is bigger)
src/vecdec mode=0 steps=10000 lerr=1 uW=2 \
	   fdet=big_det.01 fobs=big_obs.01 fdem=examples/surf_d5.dem \
	   ntot=$num nvec=$num seed=7 







