H=./examples/96.3.963.alist
# write out predicted errors and dets
src/vecdec mode=1.12 finH= $H seed=13 \
	   ntot=1000 nvec=1000 steps=40 useP=0.05 gdet=dets.01 perr=bperr.01 | tail -2 
# write out updated DET events
src/vecdec debug=1 mode=0 finH= $H ntot=1000 nvec=1000 steps=0 uW=-1 useP=0.05 \
	   fdet=dets.01 fer0=bperr.01 finA=$H gdet=yyy.01
echo "# number of convergence failures:" `grep -c 1 yyy.01`
