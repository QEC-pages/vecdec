# vecdec - vectorized decoder and LER estimator

## Overview 

The program can do several kinds of calculations depending on the value of the
`mode` command-line argument:

- `mode=0` Run a vectorized Random Information Set decoder using the
  error data generated by the program itself or externally generated
  (e.g., using `Stim`).  See the [Vectorized decoder](#vectorized-decoder) section.
- `mode=1` Run belief propagation decoder.  See the [BP decoder](#bp-decoder) section
- `mode=2` Estimate logical error rate by enumerating most likely errors.  See
  the [LER estimator](#ler-estimator) section.
- `mode=3` Parse the DEM file and save the corresponding matrices to
  files.  See the [Export the matrices](#export-the-matrices) section.
  
In any mode, the program requires the parity check matrices and error
probabilities, e.g., as given in a Detector Error Model (DEM) file,
see the [Error model](#error-model) section.

For additional details, see the [section](#-all-command-line-arguments) on
command-line arguments, the source code in `vecdec/src` directory, and example
scripts in `vecdec/examples` and `vecdec/input` directories.


## Installation 

The program uses `m4ri` library for binary linear algebra.
To install this library on a Ubuntu system, run
    `sudo apt-get install libm4ri-dev`
    
To run scripts in `vecdec/examples` directory, you will need to install
command-line versions of [Stim](https://github.com/quantumlib/Stim) and
[PyMatching](https://github.com/oscarhiggott/PyMatching).
    
For compilation *help*, change to the (vecdec/src/) directory and just run w/o
arguments  
    `make`
Since the program  is experimental, I recommend compiling with  
    `make vecdec EXTRA=""` 
This will enable additional integrity
checks, and a lot of optional debugging information.

## Error model

A detector error model (DEM) is a collection of independent `events`, each of
which is characterized by a probability `p`, a list of affected syndrome bits,
and a list of affected codewords.  Can be created by `stim`, see shell scripts
in the (vecdec/examples/) directory.  Notice that `stim` cycles are not
supported in the DEM file.  Only the lines starting with `error` are used; the
`detector` and `shift_detectors` are silently ignored, as well as any comments.
Any other entry will trigger an error.

```
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

Alternatively, the error model can be specified in terms of CSS
matrices `H=Hx` and `L=Lx` (use command-line parameters `finH=` and
`finL=` to provide the file names).  For a quantum code, instead of
`L`, the dual CSS generator matrix `G=Hz` can be given.

All matrices with entries in `GF(2)` should have the same number of
columns, `n`, and obey the following orthogonality conditions: $$
H_XH_Z^T=0, \quad H_XL_Z^T=0, \quad L_XH_Z^T=0, \quad L_XL_Z^T=I,$$
where $I$ is an identity matrix.  Notice that matrix $L_Z$ is not used
by the program and it cannot be set.

The program can read matrices from files in sparse (`coordinate`)
MaTrix market exchange (`MTX`) format and David MacKay's `alist` format. 


## Vectorized decoder

This section describes operation with the command-line switch `mode=0`
(the default).

### How it works 

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
`[(i0,s0), (i1,s1), ...]`.  The energy is calculated as the sum of LLRs for
non-zero bits in `e` and recorded, along with the sparse form of `e` if the
energy is small enough.

The program stops after a sufficient number of attempts is made, and compares
the values `L*e` for each found error `e` with the similar ones generated from
the original errors; any mismatch is a logical error.

### How to run it

The `suggested` number of syndrome vectors to generate and process is given by
`ntot`.  In reality, the number of rounds is calculated as $\lceil$ `ntot/nvec`
$\rceil$, where `nvec` is another command line argument (currently, the default
`nvec=16`, although one should probably use at least `1023` for efficiency).  In
each round, `nvec` different syndrome vectors will be generated and processed.
Ideally, `nvec` should be a factor `64`, since 64-bit integers are used to store
binary vectors internally.

The parameter `nfail`, when non-zero, will cause execution to stop after
accumulating a given number of logical errors.

The parameter `lerr`, when non-zero, specifies the maximum number of
non-zero error bits outside of the index set to try before generating
a new permutation.  This is similar to OSD level.  For large codes,
recommended values are `lerr=0` and `lerr=1`; larger values make the
program prohibitively slow (run-time generally scales as code length
to the power of `lerr`).

Another important command-line parameter is `steps`.  It should be set
to a sufficiently large number (experiment!) for decoding to be
accurate, especially close to the threshold.  The related parameter
`swait` (if non-zero) specifies the number of steps w/o any updates to
any error vector to wait before terminating the cycle.

Use `debug=0` to suppress any output except for simulation results.
Use `debug=1023` to output all possible debugging information (not all
bits are used at this time).  The default value `debug=1` causes the
program to print information about basic steps of the calculation.

Use `fdem=filename` (with or without quotes) or `fdem= filename` (with
a space) to specify the input file with the detector error model.

You may also want to set the names of the input files with externally
generated errors, e.g., using the command-line arguments 
`fdet=dets.01 fobs=obs_flips.01`

## BP decoder 

Activate it with the command-line parameters `mode=1`.  Several
variants of the BP decoder are implemented; these are controlled by
the `submode` bitmap, and also by parameters `bpalpha`, `bpretry`,
`lerr`, and `maxosd`.

### submode bitmap values 

The `submode` bitmap is specified after a decimal point in the mode
parameter, e.g., `mode=1.12` has `submode=12=2+4+8`, which corresponds
to set bits `1`, `2`, and `3`.  As detailed below, this gives `serial-V`
BP schedule based on `average LLR`.

- bit `0` in the bitmap controls the use of regular (instantaneous)
  LLR to check the convergence.  Namely, when bit `0` is set,
  instantaneous LLR will be used.
  
- bit `1` controls the use of average LLR to check the convergence.

  When both bit `0` and bit `1` are sets, both sets of LLR values will
  be computed and used sequentially after every round of BP to check
  the convergence.  For convenience, when `neither` bit `0` nor bit
  `1` is set, both sets of LLR values will be computed and used.
  
- when bit `2` is set, `serial` BP update schedule is used, otherwise
  `parallel` schedule.  Additional options for `serial` schedule are
  set by bits `3` and `4` (these only have effect with serial BP
  schedule, i.e., when bit `2` is set):
  - when bit `3` is set, `serial-V` order is selected, otherwise `serial-C`
  - when bit `4` is set, node randomization is done at every round, as
    opposed to once per run with a given syndrome vector.
	
### additional parameters affecting BP 

- parameter `bpalpha` (allowed range from 0.0 to 1.0; the default
  value is `0.5`) sets the coefficient for calculating average LLR
  values.  Smaller values correspond to shorter memory; setting
  `bpalpha=0` gives the same result at instantaneous LLR.

- parameter `bpretry` (integer, default value `1`, do not retry)
  allows to repeat BP several times using different node permutation.
  The prescribed number of attempts is made until a solution which
  satisfies the syndrome is found.  This is useless with the
  `parallel` schedule (TODO: randomization of prior probabilities).
  
- the use of optional ordered-statistics decoding (OSD) after BP
  failed to converge is controlled by the parameters `lerr` (default
  `-1`) and `maxosd` (default `100`).  With OSD, the variable nodes
  are sorted by decreasing error probability (according to BP), and
  Gauss elimination is used to find valid error vectors.  With
  `lerr=0`, the unique solution corresponds to all-zero information
  set, with `lerr=1`, in addition to all-zero information set,
  information sets with a single non-zero bit are also checked, and
  the minimum-weight vector is selected; with `lerr=2` up to two
  non-zero bits are selected, etc.  When `maxosd` parameter is
  non-zero, the search with multiple non-zero information bits is
  restricted to this range of columns.

### Quantized LLR and `min-sum` vs. `sum-product` BP updates

When the program is compiled with `USE_QLLR` (the default), parameters
of integer-based Quantized LLR module can be set using `qllr#=ubt`,
where `#` can be 1, 2, or 3.  By default, the recommended values
`qllr1=12 qllr2=300 qllr3=7` are used.

The power `1<<qllr1` determines how integral LLRs relate to real LLRs
(`to_double=(1<<qllr1)*int_llr`).  Table resolution is
`2^(-(qllr1-qllr3))`.  The parameter `qllr2` is the number of entries
in a table for LLR operations.  With `qllr2=0` (no table), the
`Sum-Product` algorithm used in the program effectively reduces to a
`Min-Sum` algorithm.

## LER and code distance estimator

This section describes operation with the command-line switch `mode=2`.

Given the error model, i.e., the matrices $H$, $L$, and the vector of column
probability values $p_i$, the program tries to enumerate the likely binary
`codewords` $c$ such that $He=0$, $Le\neq0$, 
while the associate log-likelihood probability ratio (LLR) 
$\sum_i c_i \ln(1/p_i-1)$
is not too large. 

In fact, right now it just tries to find the single most likely codeword $c$ and
outputs an estimate of the corresponding contribution to the logical error
probability, $$\prod_{i\in \mathop{\rm supp}c} 2[ p_i(1-p_i)]^{1/2}.$$

The third column in the output is the (upper bound on the) code
distance.  This is particularly useful if you want to compare the
circuit distance with that of the original code.

To speed up the distance calculation, you can use the parameter `dmin`
(by default, `dmin=0`).  When non-zero, if a code word of weight `w`
$\le$ `dmin` is found, the distance calculation is terminated
immediately, and the result `-w` with a negative sign is printed.
This is useful if, e.g., we are trying to construct a fault-tolerant
measurement circuit whose circuit distance should coincide with the
distance `d` of the original code.  In this case, if we specify
`dmin=(d-1)`, distance calculation will be terminated immediately so
that a different circuit can be studied.

Additional parameters relevant for this mode: 

- `finC` optional filename to read the list of codewords from (default: empty, do not read)
- `outC` optional filename t write the full list of codewords to (default: empty, do not write)
- `maxC` maximum number of codewords to read/write/create (default: 0, no limit)
- `nfail` maximum weight of codewords above the minimum weight
  (default: 0, just keep the minimum weight codewords).  **TODO:
  rename this parameter -- right now it is shared with mode=0 and
  mode=1 but has different meaning!**

## Export the matrices 

This section describes operation with the command-line switch
`mode=3`.  In this case the program does not try to run anything and
just parses the DEM file and saves the corresponding parity-check `H`,
observables `L` matrices and the error-probabilities `P` vector.  In
addition, the program constructs the `G` matrix (whose rows are
orthogonal to rows of `H` and `L` matrices and the total rank equals
to `n`, the number of columns in these matrices).

The matrices are written in the Matrix Market format to files with names
`${fout}H.mmx`, `${fout}G.mmx`, `${fout}L.mmx`, and `${fout}P.mmx`,
where the header string is defined by the value of the `fout=`
command-line argument.  An exceptional value is `fout=stdout`, in
which case the contents of the files is directed (surprise!) to
`stdout`.

**Note:** For some obscure reasons, the constructed `G` matrix will
 only have weight-3 rows.  This is done by finding triplets of columns
 in matrices `H` and `L` which sum to zero.  This works for a DEM file
 created from a circuit which contains one- or two-qubit depolarizing
 noise.  With insufficient rank in weight-3 rows, the program will
 currently fail.


## All command-line arguments 

You can generate the list of supported command line arguments by running 
`vecdec --help`.

```bash
./vecdec:  vecdec - vectorized decoder and LER estimator
  usage: ./vecdec param=value [[param=value] ... ]
	 Command line arguments are processed in the order given.
	 Supported parameters:
	 --help		: give this help (also '-h' or just 'help')
	 --morehelp	: give more help
	 fdem=[string]	: name of the input file with detector error model
	 finH=[string]	: file with parity check matrix Hx (mm or alist)
	 finG=[string]	: file with dual check matrix Hz (mm or alist)
	 finL=[string]	: file with logical dual check matrix Lx (mm or alist)
	 finK=[string]	: file with logical check matrix Lz (mm or alist)
	 finP=[string]	: input file for probabilities (mm or a column of doubles)
	 finC=[string]	: input file name for codewords in `nzlist` format
	 outC=[string]	: output file name for codewords in `nzlist` format
			 (if same as finC, the file will be updated)
	 useP=[double]	: fixed probability value (override values in DEM file)
		 for a quantum code specify 'fdem' OR 'finH' and ( 'finL' OR 'finG' );
		 for classical just 'finH' (and optionally the dual matrix 'finL')
	 ferr=[string]	: input file with error vectors (01 format)
	 fdet=[string]	: input file with detector events (01 format)
	 fobs=[string]	: input file with observables (01 matching lines in fdet)
		 specify either 'ferr' OR a pair of 'ferr' and 'fdet' (or none for internal)
	 fout=[string]	: header for output file names ('tmp', see 'mode=3')
		 (space is OK in front of file names to enable shell completion)
	 steps=[integer]	: num of RIS or BP decoding steps (default: 50)
	 lerr =[integer]	: OSD search level (-1, only implemented with `mode=0`, `1`)
	 maxosd=[integer]	: max column for OSD2 and above (100)
	 bpalpha=[float]	: average LLR scaling coefficient for BP (default 0.5)
	 bpretry=[integer]	: retry BP up to this many times per syndrome (1)
	 swait=[integer]	: Gauss steps w/o new errors to stop (0, do not stop)
	 nvec =[integer]	: max vector size for decoding (default: 1024)
			 (list size for distance or energy calculations)
	 ntot =[long long int]	: total syndromes to generate (default: 1)
	 nfail=[long long int]	: total fails to terminate (0, do not terminate)
	 dmin=[integer]	: terminate distance calculation immediately when
			 a vector of weight w<=dmin is found, return '-w' (default: 0)
	 seed= [long long int]	: RNG seed or use time(NULL) if 0 (default)
	 qllr1=[integer]	: if 'USE_QLLR' is set, parameter 'd1' (12)
	 qllr2=[integer]	: if 'USE_QLLR' is set, parameter 'd2' (300)
	 qllr3=[integer]	: if 'USE_QLLR' is set, parameter 'd3' (7)
		 These are used to speed-up LLR calculations, see 'qllr.h'
		 Use 'qllr2=0' for min-sum.
	 mode=int[.int]	: operation mode[.submode] (default: 0.0)
		* 0: use basic vectorized (random information set) decoder
			 read detector events from file 'fdet' if given, otherwise
			 generate 'ntot' detector events and matching observable flips;
			 decode in chunks of size 'nvec'. 
			 Read observable flips from file 'fobs' if given
		* 1: Belief Propagation decoder.  By default (no submode specified)
			   parallel BP using both regular and average LLRs.
			   Other options are selected by 'submode' bitmap: 
			 .1 (bit 0) use regular LLR
			 .2 (bit 1) use average LLR - these take precendence
			 .4 (bit 2) use serial BP schedule (not parallel)
			 .8 (bit 3) with serial, use V-based order (not C-based)
			 .16 (bit 4) with serial, randomize node order in each round 
			     (by default randomize only once per run)
		* 2: generate most likely fault vectors, estimate Prob(Fail).
			 Use up to 'steps' random information set (RIS) decoding steps
			 unless no new fault vectors have been found for 'swait' steps.
			 Keep vectors of weight up to 'nfail' above min weight found.
			 Generate up to 'ntot' unique min-energy fault vectors.
			 (if the corresponding parameters are set to non-zero values)
			 If 'outC' is set, write full list of CWs to this file.
			 If 'finC' is set, read initial set of CWs from this file.
		* 3: Read in the DEM file and output the corresponding 
			 H, G, L, and K matrices and the probability vector P.
			 Use 'fout=' command line argument to generate file names
			 ${fout}H.mmx, ${fout}G.mmx, ${fout}L.mmx, ${fout}K.mmx, and ${fout}P.mmx
			 with 'fout=stdout' all output is sent to 'stdout'
	 debug=[integer]	: bitmap for aux information to output (default: 1)
		*   0: clear the entire debug bitmap to 0.
		*   1: output misc general info (on by default)
		*   2: output matrices for verification
			 see the source code for more options
	 See program documentation for input file syntax.
	 Multiple 'debug' parameters are XOR combined except for 0.
	 Use debug=0 as the 1st argument to suppress all debug messages.
```

## Libraries 

The program uses `m4ri` library for binary linear algebra.  To install
under Ubuntu, run
```
sudo apt-get update -y
sudo apt-get install -y libm4ri-dev
```

`Tiny Mersenne Twister` written by Mutsuo Saito and Makoto Matsumoto
is used for random number generation (header file `src/tinymt64.h` is
included with the distribution).

`uthash` by Troy D. Hanson and Arthur O'Dwyer is used for hashing
storage (header file `src/uthash.h` is included with the
distribution).

Shell scripts in the `examples/` directory assume command-line
versions of `PyMatching` and `Stim` packages compiled and located in
`../PyMatching` and `../Stim` (with respect to the location of this
file).  These excellent packages written by Oscar Higgott and Craig
Gidney, respectively, are available from `GitHub`. Please refer to the
documentation of these packages for the installation instructions.

## Future

Eventually, more sophisticated algorithms will be added to the
program.  Currently in the works are serial BP decoder and stat-mech
maximum-likelihood decoder, as well as some degeneracy-reducing matrix
transformations for `mode=3`.

