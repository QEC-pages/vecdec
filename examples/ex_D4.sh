#!/usr/bin/bash -f
p1=0.01
for ((w=11; w<=15; w+=2)) ; do 
    echo maxW=$w mode=2.1 `src/vecdec debug=0 mode=2.1 steps=0 \
              finH= ./examples/96.3.963.alist useP=$p1 finC=tmp.nz maxW=$w` 
done
echo "#"
for ((w=11; w<=15; w+=2)) ; do 
    echo maxW=$w mode=2.2 `src/vecdec debug=0 mode=2.2 steps=0 \
              finH= ./examples/96.3.963.alist useP=$p1 finC=tmp.nz maxW=$w` 
done

