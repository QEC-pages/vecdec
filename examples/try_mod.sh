# example of available code transformations:
# removal of rows of weight 1, 2, and 3 from the G matrix.
# note: partition function verification is done with Mathematica 

# use these matrices 
mH=../examples/tryH.mmx
mL=../examples/tryL.mmx
pp=0.05
vecdec=../src/vecdec
mtx_sub=../src/mtx_sub

# generate codewords
$vecdec debug=0 steps=1000 mode=2 dW=-1 outC=tryC.nz finH=$mH useP=$pp
cat tryC.nz
echo 

# write matrices and calculate the partition function "prefix: tmp"
$vecdec debug=16 mode=3.31 finC=tryC.nz finH=$mH finL=$mL useP=$pp  fout=tmp dW=-1 
echo "original H:"
$mtx_sub debug=128 fin=tmpH.mmx 
echo "original L:"
$mtx_sub debug=128 fin=tmpL.mmx 
echo "original G:"
$mtx_sub debug=128 fin=tmpG.mmx 
echo "original K:"
$mtx_sub debug=128 fin=tmpK.mmx 

 # echo "export matrices for modified code"
$vecdec debug=32 mode=3.32 finH=$mH finL=$mL useP=$pp fout=tmp finC=tryC.nz dW=-1 maxW=3

 # generate codewords for the modified matrices
$vecdec debug=0 steps=1000 mode=2 dW=-1 outC=tmpC.nz finH=tmpH.mmx finP=tmpP.mmx

echo 
 # echo "# calculate new partition functions"
$vecdec debug=16 mode=3.31 finC=tmpC.nz finH=tmpH.mmx finP=tmpP.mmx finL=tmpL.mmx dW=-1 fout=try2 maxW=4
echo "modified H:"
$mtx_sub debug=128 fin=try2H.mmx 
echo "modified G:"
$mtx_sub debug=128 fin=try2G.mmx 
echo "modified L:"
$mtx_sub debug=128 fin=try2L.mmx 
echo "modified K:"
$mtx_sub debug=128 fin=try2K.mmx 
cat try2P.mmx 


