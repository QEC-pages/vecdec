#!/usr/bin/bash -f
./src/vecdec mode=2.1 steps=100000 seed=7 finH= ./examples/96.3.963.alist \
	     useP=0.05 outC=tmp.nz dW=10
