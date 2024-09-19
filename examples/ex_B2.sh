
for lerr in 0 1 2 3 ; do
    echo "#"
    for steps in 5 10 25 50 ; do 
	echo -e lerr=$lerr steps=$steps " \t" \
            `./src/vecdec debug=0 mode=1 maxosd=50 steps=$steps lerr=$lerr  \
                   finH= ./examples/96.3.963.alist ntot=1000 nvec=1000 useP=0.05 seed=113`
    done
done
