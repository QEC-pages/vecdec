vecdec=./src/vecdec
for lerr in 0 1 ; do 
  for steps in 500 2000 5000 ; do 
    echo -e lerr=$lerr steps=$steps "\t" \
  	    `$vecdec debug=0 seed=113 mode=0 finH= ./examples/96.3.963.alist \
	             ntot=1000 nvec=1000 steps=$steps useP=0.05 lerr=$lerr` 
  done 
done
