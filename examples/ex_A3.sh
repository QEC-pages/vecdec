vecdec=./src/vecdec
for uW in 1 2 3 4; do 
    echo -e uW=$uW "\t" \
  	 `$vecdec debug=0 seed=113 mode=0 lerr=0 \
	     finH= ./examples/96.3.963.alist \
	     ntot=1000 nvec=1000 steps=5000 useP=0.05 uW=$uW ` 
done
