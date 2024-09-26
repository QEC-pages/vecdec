
for minor in 0 1 2  4 5 6 12 13 14 ; do
    echo "#"
    for steps in 10 50 ; do 
	echo -e minor=$minor steps=$steps " \t" \
            `./src/vecdec debug=0 mode=1.$minor maxosd=50 steps=$steps lerr=1 \
                   finH= ./examples/96.3.963.alist ntot=100000 nvec=100000 useP=0.05 seed=113`
    done
done
