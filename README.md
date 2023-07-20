# vecdec - a simple vectorized decoder

## Installation 

The program uses `m4ri` library for binary linear algebra.
To install this library on a Ubuntu system, run
    `sudo apt-get install libm4ri-dev`
    
To use the example scripts, you will also need a command-line version of
[Stim](https://github.com/quantumlib/Stim)
    
For compilation *help*, change to the (vecdec/src/) directory and just run w/o
arguments  
    `make`
Since the program  is experimental, I recommend compiling with  
    `make vecdec EXTRA=""`
This will enable additional integrity checks.

## Error model

A detector error model (DEM) is a collection of independent `events`, each of
which is characterized by a probability `p`, a list of affected syndrome bits,
and a list of affected codewords.  Can be created by `stim`, see shell scripts
in the (vecdec/examples/) directory.  Notice that `stim` cycles are not
supported in the DEM file.  Only the lines starting with `error` are used; the
`detector` and `shift_detectors` are silently ignored, as well as any comments.
Any other entry will trigger an error.

```bash
# An example DEM file created by `stim`
error(0.125) D0
error(0.125) D0 D1
error(0.125) D0 D2
error(0.125) D1 D3
error(0.125) D1 L0
error(0.125) D2 D4
error(0.125) D3 D5
error(0.125) D4 D6
error(0.125) D5 D7
detector(1, 0) D0
detector(3, 0) D1
shift_detectors(0, 1) 0
detector(1, 0) D2
detector(3, 0) D3
shift_detectors(0, 1) 0
detector(1, 0) D4
detector(3, 0) D5
detector(1, 1) D6
detector(3, 1) D7
```

## How it works 

Given a binary vector `s` with detector events, the goal is to construct the
most likely binary vector `e` such that `H*e=s`.  The decoding is verified by
computing the vector of observable bits `L*e` and comparing it with the
corresponding result computed from the actual error.

The random errors and the corresponding detector/observable vectors can be
generated on the fly (for a given error model), or read from files in 01
format.

The program processes up to `nvec` syndrome vectors at a time.  The syndrome
columns are written as columns of a binary matrix `S`.  The original check
matrix `H` (extracted from the error model) is combined with `S`, to form
block-matrix `[H,S]`.  At each step, a column ordering `P` is randomly generated
(using values of event probabilities to help), the Gauss elimination is
performed on the rows of the combined matrix, creating the list of pivot columns
`[i0, i1, ...]`, with one entry per row.  Given the transformed syndrome column
`[s0, s1, ...]`, the output vector has values given by the list of pairs
[(i0,s0), (i1,s1), ...].  The energy is calculated as the sum of LLRs for
non-zero bits in `e` and recorded, along with the sparse form of `e` if the
energy is small enough.

The program stops after a sufficient number of attempts is made, and compares
the values `L*e` for each found error `e` with the similar ones generated from
the original errors; any mismatch is a logical error.

## How to run it

The `suggested` number of syndrome vectors to generate and process is given by
`ntot`.  In reality, the number of rounds is calculated as $\lceil$ `ntot/nvec`
$\rceil$, where `nvec` is another command line argument (default `nvec=16`).  In
each round, `nvec` different syndrome vectors will be generated and processed.
Ideally, `nvec` should be a factor `64`, since 64-bit integers are used to store
binary vectors internally.

The parameter `nfail`, when non-zero, will cause execution to stop after
accumulating a given number of logical errors.

The parameter `lerr`, when non-zero, specifies the maximum number of non-zero
error bits outside of the index set to try before generating a new permutation.
This is similar to OSD level.

Another important command-line parameter is `steps`.  It should be set to a
large number (experiment!) for decoding to be accurate, especially close to the
threshold.  The related parameter `swait` (if non-zero) specifies the number of
steps w/o any updates to any error vector to wait before terminating the cycle.

Use `debug=0` to suppress any output except for simulation results.  Use
`debug=1023` to output all possible debugging information (not all bits are used
at this time).

Use `fdem="filename"` (with or without quotes) or `fdem= "filename"` (with a
space) to specify the input file with the detector error model.  

With `mode=1` switch, you need to set the names of the input files, e.g., 
using the command-line arguments  
    `fdet=dets.01 fobs=obs_flips.01`

## All command-line arguments 

You can obtain these by running `vecdec --help`

```bash
src/vecdec: a simple vectorized random information set decoder.
  usage: src/vecdec param=value [[param=value] ... ]
	 Command line arguments are processed in the order given.
	 Supported parameters:
	 --help	: give this help (also '-h' or just 'help')
	 fdem=[string]	: name of the input file with detector error model
	 fdet=[string]	: input file with detector events (01 format)
	 fobs=[string]	: input file with observable flips (01 format)
	 steps=[integer]	: num of random window decoding steps (default: 1)
	 lerr=[integer]		: local search level after gauss (0, no search)
	 swait=[integer]	: steps w/o new errors to stop (0, do not stop)
	 nvec=[integer]		: max vector size for decoding (default: 16)
	 ntot=[integer]		: total syndromes to generate (default: 1)
	 nfail=[integer]	: total fails to terminate (0, do not terminate)
	 seed=[integer]		: RNG seed or use time(NULL) if 0 (default)
	 mode=[integer]		: bitmap for operation mode (default: 0)
		*   0: clear the entire mode bitmap to 0.
		*   1: read detector events/observables from 01 files 
		*   2: cycle global probabilities in error model from `pmin` to `pmax`
	 pmin = [double]	: min global probability with `mode&2` (-1)
	 pmax = [double]	: max global probability with `mode&2` (-1)
	 pstep =[double]	: step of probability with `mode&2` (1e-3)
	 debug=[integer]	: bitmap for aux information to output (default: 1)
		*   0: clear the entire debug bitmap to 0.
		*   1: output misc general info (on by default)
		*   2: output matrices for verification
	 See program documentation for input file syntax.
	 Multiple `debug` parameters are XOR combined except for 0.
	 Use debug=0 as the 1st argument to suppress all debug messages.
```

## Libraries 

The program uses `m4ri` library for binary linear algebra.

## Future

Eventually, the program will be able to process externally generated syndrome
data, or use other more efficient decoders, e.g., belief propagation.

## To do list 

- [ ] Finish local recursive search implementation 
- [ ] Ensure that local search is done in order of decreasing `p` (`?`)
- [ ] Test performance (statistics of updates): `A`: random window initial
      (sorted by `p`), `B`: random window secondary (unbiased permutations),
      `Aj`: local `j` after `A`, `Bj`: local `j` after `B`.
- [x] Remove non-DEM error model
- [ ] Add hashing table for triplets (s,e,energy) and implement ML (`mode=1`).
- [ ] Reserve `mode=4` for BP+OSD.  Since we are doing OSD which generates a
      local list, `mode=1` can also be used here
- [ ] Other minimum energy solvers? 
