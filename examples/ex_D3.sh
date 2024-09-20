#!/usr/bin/bash -f
for ((w=6; w<=16; w++)) ; do 
  echo $w `src/vecdec debug=0 mode=2.1 steps=0 \
              finH= ./examples/96.3.963.alist useP=0.05 finC=tmp.nz maxW=$w` 
done
