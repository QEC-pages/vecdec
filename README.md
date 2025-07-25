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

In addition to regular decoders, a fast cluster-based pre-decoder
similar to a single-step Union Find decoder is available for both RIS
and BP decoders.  See the [Pre-decoder](#pre-decoder) section.

Some examples of using `vecdec` in various modes are given in the
[Common tasks](#common-tasks) section.

For additional details, see the [section](#-all-command-line-arguments) on
command-line arguments, the source code in `vecdec/src` directory, and example
scripts in `vecdec/examples` and `vecdec/input` directories.


## Installation 

### Libraries and related programs
The program uses `m4ri` library for binary linear algebra.
To install this library on a Ubuntu system, run
    `sudo apt-get install libm4ri-dev`
	
Alternatively, download `m4ri` library from
[github](https://github.com/malb/m4ri/) using, e.g., 
```
git clone https://github.com/malb/m4ri/
```
and follow the instructions to compile and install it.  You may need
to edit your `C_INCLUDE_PATH`, `CPLUS_INCLUDE_PATH`, or `LIBRARY_PATH`
if you install outside normal directory structure.
Alternatively, you can also set `INC` variable in the `Makefile`, or set it during 
compilation as an argument to `make` (this assumes `m4ri` was compiled
but not installed in `../../m4ri` directory relative the
`vecdec/src`):
```
make -j all DINC="-I../../m4ri -L../../m4ri/.libs"
```

To run scripts in `vecdec/examples` directory, you will also need to
install command-line versions of  
[Stim](https://github.com/quantumlib/Stim) and  
[PyMatching](https://github.com/oscarhiggott/PyMatching).  

### Compilation

Get the source code from
[https://github.com/QEC-pages/vecdec](https://github.com/QEC-pages/vecdec),
e.g., using `git`, and compile by running 
```
git clone https://github.com/QEC-pages/vecdec
cd vecdec/src
make -j all
```

The program should compile without warnings on a reasonably recent
Linux system.

For compilation *help*, change to the (vecdec/src/) directory and just run w/o
arguments 
    `make`

Normal compilation defines the variable `NDEBUG`.
If you run into trouble, you way want to recompile without it, by running
    `make clean && make vecdec EXTRA=""`  
This will enable additional integrity checks, and a lot of optional
debugging information.  Some of the additional checks may be
expensive; the program runs slower in this mode.

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
columns, `n`, and obey the following orthogonality conditions:
$$H_XH_Z^T=0,\quad H_XL_Z^T=0,\quad L_XH_Z^T=0,\quad L_XL_Z^T=I,$$
where $I$ is an identity matrix.  Notice that the latter identity is
not required; it is sufficient that `Lx` and `Lz` matrices have the
same full row rank `=k`, the dimension of the code, each row of `Lx`
has a non-zero scalar product with a row of `Lz`, and vice versa.

The program can read matrices from files in sparse (`coordinate`)
MaTrix market exchange (`MTX`) format and David MacKay's `alist`
format (the format is recognized automatically).

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
the `submode` bitmap, and also by parameters `bpgamma`, `bpretry`,
`lerr`, and `maxosd`.

### submode bitmap values 

The `submode` bitmap is specified after a decimal point in the mode
parameter, e.g., `mode=1.14` has `submode=14` ($=2+4+8=2^1+2^2+2^3$), which corresponds
to set bits `1`, `2`, and `3`.  As detailed below, this gives `serial-V`
BP schedule based on `average LLR`.

- bit `0` in the bitmap controls the use of regular (instantaneous)
  LLR to check the convergence.  Namely, when bit `0` is set,
  instantaneous LLR will be used.
  
- bit `1` controls the use of average LLR to check the convergence.

  When both bit `0` and bit `1` are set, both sets of LLR values will
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

- parameter `bpgamma` (allowed range from 0.0 to 1.0; the default
  value is `0.5`) sets the coefficient for calculating average LLR
  values.  Smaller values correspond to shorter memory; setting
  `bpgamma=0` gives the same result at instantaneous LLR.

- parameter `bpretry` (integer, default value `1`, do not retry)
  allows to repeat BP several times using different node permutations.
  The prescribed number of attempts is made until a solution which
  satisfies the syndrome is found.  This is useless with the
  `parallel` schedule (**TODO:** randomization of prior probabilities).
  
- the use of optional ordered-statistics decoding (OSD) after BP
  failed to converge is controlled by the parameters `lerr` (default
  `-1`) and `maxosd` (default `100`).  With OSD, the variable nodes
  are sorted by decreasing error probability (according to BP), and
  Gauss elimination is used to find valid error vectors.  With
  `lerr=0`, the unique solution corresponds to all-zero information
  set (variable nodes corresponding to pivot columns), with `lerr=1`,
  in addition to the all-zero information set, information sets with a
  single non-zero bit are also checked, and the minimum-weight vector
  is selected; with `lerr=2` up to two non-zero bits are selected,
  etc.  When `maxosd` parameter is non-zero, the search with
  *multiple* non-zero information bits is restricted to this range of
  columns (default value `maxosd=100`).

### Quantized LLR and `min-sum` vs. `sum-product` BP updates

When the program is compiled with `USE_QLLR` (the default), parameters
of integer-based Quantized LLR module can be set using `qllr#=int`,
where `#` can be 1, 2, or 3.  By default, the recommended values
`qllr1=12 qllr2=300 qllr3=7` are used.

The power `1<<qllr1` determines how integral LLRs relate to real LLRs
(`to_double=(1<<qllr1)*int_llr`).  Table resolution is
`2^(-(qllr1-qllr3))`.  The parameter `qllr2` is the number of entries
in a table for LLR operations.  With `qllr2=0` (no table), the
`Sum-Product` algorithm used in the program effectively reduces to a
`Min-Sum` algorithm.

## Pre-decoder

To speed-up the RIS and BP decoders, especially at small error
probabilities, a fast pre-decoder based on syndrome cluster
decomposition (similar to a single-step Union Find decoder) is now
available.  Namely, some error vectors and associated syndromes are
stored in a hash.  The initial syndrome vector is decomposed into
connected clusters based on syndrome connectivity graph; the look up
decoding is attempted for each cluster separately.  

The errors in hash are controlled by parameters `uW` (maximum error weight; use
`0` to just skip zero-syndrome vectors and `-1` to skip the pre-decoder
altogether; default value is `2`), `uR` (maximum v-v graph distance between
non-zero bits in stored error vectors; use `0` for no limit; default value is
`1`), and `maxU` (maximum number of vectors in the hash; the default value is
`0` for no limit).  

The bitmap parameter `uX` controls additional options which may or may not
improve the decoding; the default value is `0`.
 - With non-zero bit-0 (e.g., `uX=1`), the program will attempt to match the
   entire syndrome as a single `low-weight` error; this may give a better
   decoding accuracy for some codes and, in some cases, a marginal speed-up (the
   cluster-matching pre-decoder is usually fast).
 - With non-zero bit-1 (e.g., `uX=2`), if some of the clusters failed to match,
   the program will remove the matched clusters and pass to the decoder only the
   unmatched clusters.  This may improve the decoding probability accuracy,
   e.g., with a BP decoder which tends to fail often due to degeneracy.  Also,
   this option is not recommended for surface codes, but may work well for codes
   with larger syndrome distances where mismatched clusters are
   not as likely.

## LER and code distance estimator

**This information it currently out of date.  Please see the output of `vecdec
--help` and `vecdec mode=2 --help` for reference.**

This section describes operation with the command-line switch `mode=2`.
Summary of additional mode options:
* Without any submode specified, the program just tries to find a
  minimum-weight codeword.  In this regime, it is recommended to add
  `useP=-1` to the command line, to skip energy calculations. 
* With submode 1 (i.e., `mode=2.1`), the program will use the computed
  codewords to estimate the fail probability.  Two quantities will be
  computed: the sum of estimated fail probabilities over the codewords
  and the maximum single-codeword fail probability.  The results are
  not expected to be accurate since estimated prefactors will be used.
* The submode 2 (i.e., `mode=2.2`) is similar, except that exact
  prefactors will be used (for a codeword of weight `w` a prefactor is
  computed by examining all $2^w$ binary error patterns).

Given the error model, i.e., the matrices $H$, $L$, and the vector of column
probability values $p_i$, the program tries to enumerate the likely binary
`codewords` $c$ such that $He=0$, $Le\neq0$, 
while the associate log-likelihood probability ratio (LLR) 
$\sum_i c_i \ln(1/p_i-1)$
is not too large. 

The following gives a more detailed description for `mode=2.1`:
The program stores the list of codewords found in a hash, and
outputs four numbers: 
- The sum of estimated contributions to the logical error probability
  from all codewords found, $\displaystyle\sum_c\prod_{i\in
  \mathop{\rm supp}c} 2[ p_i(1-p_i)]^{1/2}$.  If the list of codewords
  is large enough, this gives an upper bound on the fail probability.
- Maximum fail probability from a single codeword (maximum term
  contributing to the sum above).  Up to a prefactor, this gives a
  lower bound on the fail probability.
- Minimum weight of the codeword found.  This gives the code distance
  (`Z`-distance for a CSS code with `H=Hx`, `G=Hz`).  This is
  particularly useful if you want to compare the circuit distance with
  that of the original code.
- Number of codewords contributing to the total.  With `dW=0`, this
  gives the total number of distinct minimum-weight codewords found.

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

- `finC` optional filename to read the list of codewords from
  (default: empty, do not read).  Notice that the energy of the
  codewords is not stored in the file, it will be recomputed.  The
  codewords read are verified after reading; error will result if any
  codeword does not satisfy the orthogonality conditions, `H*c=0`,
  `L*c!=0` (this may happen, e.g., if an incorrect file was specified).
- `outC` optional filename to write the full list of codewords to
  (default: empty, do not write).
- `maxC` maximum number of codewords to read/write/create (default: `0`, no limit).
- `dW` maximum weight of codewords above the minimum weight (default:
  `0`, just keep the minimum weight codewords).  Setting `dW=-1`
  suppresses the upper limit on the weight of codewords stored.
- `dE` maximum energy of codewords above the minimum energy (default:
  `-1`, no limit on the energy of codewords found).

## Export the matrices 

This section describes operation with the command-line switch
`mode=3`.  In this case the program does not try to run anything and
just parses the DEM file and saves the corresponding parity-check
`H=Hx`, observables `L=Lx` matrices and the error-probabilities `P`
vector.  In addition, the program constructs the `G=Hz` and `K=Lz`
matrices.  As a reminder, rows of `G` are orthogonal to rows of `H`
and `L`, while rows of `K` are orthogonal to the rows of `H`, are
linearly independent from rows of `G`, and each has a non-zero product
with a row of `L`.  Also, `rank H + rank G = n-k`, `rank L = rank K =
k`, where `n` is the number of variable nodes (matrix columns) and `k`
is the dimension of the code.  For a classical code, matrix `G` is
trivial (has zero rank), while the `L` matrix can be selected to have
all rows of weight `1`.

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

You can generate the list of supported command line arguments by
running `vecdec --help`.  Additional information is available by
running `vecdec mode=[int] help` or `vecdec --morehelp`.  Here is the
corresponding output:

### run `./vecdec --help`
```bash
./src/vecdec:  vecdec - vectorized decoder and LER estimator
  usage: ./src/vecdec param=value [[param=value] ... ]
	 Command line arguments are processed in the order given except
	 for 'mode' and 'debug' (these are processed first).
	 Supported parameters:
	 --help		: give this help (also '-h' or just 'help')
	 mode=[integer] help		: help for specific mode
	 --morehelp	: give more help on program conventions
	 fdem=[string]	: name of the input file with detector error model
	 finH=[string]	: file with parity check matrix Hx (mm or alist)
	 finHT=[string]	: file with transposed parity check matrix Hx (mm only)
	 finG=[string]	: file with dual check matrix Hz (mm or alist)
	 finL=[string]	: file with logical dual check matrix Lx (mm or alist)
	 finK=[string]	: file with logical check matrix Lz (mm or alist)
	 finP=[string]	: input file for probabilities (mm or a column of doubles)
	 finQ=[string]	: input file for alt probabilities (for use with mode 2.16)
	 finA=[string]	: additional matrix to correct syndromes (mm or alist)
	 finC=[string]	: input file name for codewords in `nzlist` format
		 (space is OK in front of file names to enable shell completion)
	 outC=[string]	: output file name for codewords in `nzlist` format
			 (if same as finC, the file will be updated)
	 maxC=[long long int]	: max number of codewords to read/write/store
	 uW=[integer]	: max weight of an error cluster in hash (default: 2)
		 ('0': no hash but skip zero-weight syndrome vectors; '-1': do not skip)
	 uR=[integer]	: max range of v-v neighbors for errors in syndrome hash
		 (use '0' for no limit; recommended default: 1)
	 uX=[integer]	: bitmap for cluster-based predecoder options (default: 0)
		 1 (bit 0) try low-weight error w/o cluster decomp (recommend with 'uR=0'); 
		 2 (bit 1) use partial cluster matches (experimental)
	 maxU=[long long integer]	: max number of syndrome vectors in hash
		 for pre-decoding (default: '0', no limit)
	 epsilon=[double]	: small probability cutoff (default: 1e-8)
	 useP=[double]	: fixed probability value (override values in DEM file)
		 default: 0, do not override; negative value = weight-only mode
	 mulP=[double]	: scale probability values from DEM file
		 for a quantum code specify 'fdem' OR 'finH' and ( 'finL' OR 'finG' );
		 for classical just 'finH' or 'finHT' (and optionally the dual matrix 'finL')
	 useQ=[double]	: alt fixed probability value for use with mode 2.16
	 refQ=[double]	: reference error probability for use with mode 2.16
	 ferr=[string]	: input file with error vectors (01 format)
	 fer0=[string]	: add'l error to correct det events 's+A*e0' (01 format)
		 where matrix 'A' is given via 'finA', 's' via 'fdet', and 'e0'
		 are the additional error vectors.
	 fobs=[string]	: input file with observables (01 matching lines in fdet)
	 fdet=[string]	: input file with detector events (01 format)
		 specify either 'ferr' OR a pair of 'ferr' and 'fdet' (or none for internal)
	 gobs, gdet=[string]	: out file for generated vectors (01 format)
	 perr, pobs, pdet=[string]	: out file for predicted vectors (01 format)
	 fout=[string]	: header for output file names ('tmp', see 'mode=3')
	 steps=[integer]	: num of RIS or BP decoding steps (default: 50)
	 lerr =[integer]	: OSD search level (-1, only implemented with `mode=0`, `1`)
	 maxosd=[integer]	: max column for OSD2 and above (100)
	 bpgamma=[float]	: average LLR scaling coefficient for BP (default 0.5)
	 bpretry=[integer]	: retry BP up to this many times per syndrome (1)
	 swait=[integer]	: Gauss steps w/o new errors to stop (0, do not stop)
	 nvec =[integer]	: max vector size for decoding (default: 1024)
			 (list size for distance or energy calculations)
	 ntot =[long long int]	: total syndromes to generate (default: 1)
	 nfail=[long long int]	: total fails to terminate (0, do not terminate)
	 dW=[integer]	: if 'dW>=0', may keep vectors of weight up to 'minW+dW' (0)
	 maxW=[integer]	: if non-zero, skip any vectors above this weight (0)
	 dE=[double]	: if 'dE>=0', may keep vectors of energy up to 'minE+dE'
			 (default value: -1, no upper limit on energy)
	 dmin=[integer]	: terminate distance calculation immediately when
			 a vector of weight 'W<=dmin' is found, return '-w' (default: 0)
	 seed= [long long int]	: RNG seed or automatic if <=0 (default: 0)
	 qllr1=[integer]	: if 'USE_QLLR' is set, parameter 'd1' (12)
	 qllr2=[integer]	: if 'USE_QLLR' is set, parameter 'd2' (300)
	 qllr3=[integer]	: if 'USE_QLLR' is set, parameter 'd3' (7)
		 These are used to speed-up LLR calculations, see 'qllr.h'
		 Use 'qllr2=0' for min-sum.
	 mode=int[.int]	: operation mode[.submode] (default: 0.0)
		* 0: use basic vectorized (random information set) decoder
		* 1: Belief Propagation decoder.
		* 2: Generate most likely fault vectors, estimate Prob(Fail).
		* 3: Read in the DEM file and optionally write the corresponding 
			 G, K, H, and L matrices and the probability vector P.
	 debug=[integer]	: bitmap for aux information to output (default: 1)
		*   0: clear the entire debug bitmap to 0.
		*   1: output misc general info (on by default)
		*   2: additional info; calculate code confinement
			 see the source code for more options
	 See program documentation for input file syntax.
	 Multiple 'debug' parameters are XOR combined except for 0.
```

### run `./vecdec mode=0 --help`
```sh
 mode=0 : use basic vectorized (random information set) decoder
	 No 'submode' can be used with this mode. 
	 This decoder is exponentially slow but it is not specific to 
	 quantum LDPC codes.  Accuracy and performance are 
	 determined by parameters 'steps' (number of RIS rounds) and 'lerr'.
	 Long codes may require exponentially large number of 'steps'.
	 Values 'lerr>1' can be slow for long codes.
	 Specify a single DEM file 'fdem', or 'finH', 'finL', and 'finP'
	 separately (either 'finL' or 'finG' is needed for a quantum code).
	 Use 'useP' to override error probability values in DEM file.   
	 Use 'mulP' to scale error probability values from DEM file.   
	 Errors can be generated internally or read from 01 file 'ferr'.
	 Alternatively, files with detector events and observables 
	 can be specified via 'fdet' and 'fobs'. 
	 Set 'nfail' and/or 'swait' for early termination.
	 Total of 'ntot' errors will be read or generated in chunks of 'nvec'.
	                                                       

	 With 'uW' non-negative, use hash storage to store likely syndrome
		 vectors to speed up the decoding.  Parameter 'maxU>0' sets the limit on the
		 number of syndrome vectors in the hash; no limit if 'maxU=0'.  
		 Parameter 'uR>0' sets the limit on the v-v graph distance between non-zero positions
		 in an error vector stored in the hash; no limit if 'uR=0'
		 Bitmap 'uX', when non-zero, enables experimental options for cluster-based
		 predecoder: try to match syndrome as a whole (bit 0), and using partially
		 matched syndrome clusters, in which case only residual error is sent to the
		  main decoder (default: 0) (use with caution, fail rates may increase!)
		 With debug&2 non-zero and uW>0, print out the confinement function
```

### run `./vecdec mode=1 --help`
```sh
 mode=1.[submode] : use one of several iterative decoder versions
	Submode bitmap values:
			 .1 (bit 0) use regular LLR
			 .2 (bit 1) use average LLR - these take precendence
			 .4 (bit 2) use serial BP schedule (not parallel)
			 .8 (bit 3) with serial, use V-based order (not C-based)
			 .16 (bit 4) with serial, randomize node order in each round 
			     (by default randomize only once per run)
			 .32 (bit 5) with serial, suppress all node order randomization
			     (must have bpretry=1 and bit 4 not set in submode)
	 For convenience, 'submode=0' is equivalent to 'submode=3'. 
	 This decoder may experience convergence issues.
	 Accuracy and performance are determined by parameters 
	 'steps' (number of BP rounds), 'lerr' (OSD level, defaul=-1, on OSD).
	 and 'maxosd', the number of columns for OSD in levels 2 and above.
	 Using 'steps' not higher than 50 is recommended.
	 Use 'bpgamma' to specify how averaging is done (default: 0.5).
	 Use 'bpretry' to specify how many times to retry BP (default: 1, do not retry) 
	 Use 'qllr' parameters to set LLR quantization
	   or compile with 'VER=""' option for double LLR values. 
	 With 'qllr2=0' the Sum-Product algorithm reduces to a 'Min-Sum'.
	 Specify a single DEM file 'fdem', or 'finH', 'finL', and 'finP'
	 separately (either 'finL' or 'finG' is needed for a quantum code).
	 Use 'useP' to override error probability values in DEM file.   
	 Use 'mulP' to scale error probability values from DEM file.   
	 Errors can be generated internally or read from 01 file 'ferr'.
	 Alternatively, files with detector events and observables 
	 can be specified via 'fdet' and 'fobs'. 
	 Set 'nfail' and/or 'swait' for early termination.
	 Total of 'ntot' errors will be read or generated in chunks of 'nvec'.
	                                                       

	 With 'uW' non-negative, use hash storage to store likely syndrome
		 vectors to speed up the decoding.  Parameter 'maxU>0' sets the limit on the
		 number of syndrome vectors in the hash; no limit if 'maxU=0'.  
		 Parameter 'uR>0' sets the limit on the v-v graph distance between non-zero positions
		 in an error vector stored in the hash; no limit if 'uR=0'
		 Bitmap 'uX', when non-zero, enables experimental options for cluster-based
		 predecoder: try to match syndrome as a whole (bit 0), and using partially
		 matched syndrome clusters, in which case only residual error is sent to the
		  main decoder (default: 0) (use with caution, fail rates may increase!)
		 With debug&2 non-zero and uW>0, print out the confinement function
```

### run `./vecdec mode=2 --help`
```sh
 mode=2 : Generate most likely fault vectors, estimate Prob(Fail).
	Submode bitmap values:
			 .1  (bit 0) calculate original fail probability estimate
			 .2  (bit 1) fail prob estimate using average LLRs and cw count
			 .4  (bit 2) calculate exact greedy probability estimate
			 .8  (bit 3) approx greedy prob estimate with prefactor
			 .16 (bit 4) use reference `refQ/finQ/useQ` to calculate fail
				 probability estimates (as opposed to direct summations)
	 Use up to 'steps' random information set (RIS) steps
	 unless no new codewords (fault vectors) have been found for 'swait' steps.
	 Use 'steps=0' to just use the codewords from the file 
	 With `dW>=0`, keep vectors of weight up to 'dW' above min weight found.
	 With `dE>=0`, keep vectors of energy up to 'dE' above minimum E found (sum of LLRs).
	 When 'maxC' is non-zero, generate up to 'maxC' unique codewords.
	 If 'outC' is set, write full list of CWs to this file.
	 If 'finC' is set, read initial set of CWs from this file.
	 Accuracy and performance are determined by parameters 
	 'steps' (number of RIS rounds), 'lerr' (OSD level, defaul=-1, no OSD).
	 Specify a single DEM file 'fdem', or 'finH', 'finL', and 'finP'
	 separately (either 'finL' or 'finG' is needed for a quantum code).
	 Use 'useP' to override error probability values in DEM file.   
	 Use 'mulP' to scale error probability values from DEM file.   
	 Similarly, use 'finQ' or 'useQ' arguments to specify alternative
	 probability vectors with mode 2.16
```

### run `./vecdec mode=3 --help`
```sh
 mode=3 : Export matrices associated with the code.
	 Read in the DEM file and optionally write the corresponding 
	 G, K, H, and L matrices and the probability vector P.
	 By default (submode&31=0) output everything, otherwise
	 .1 (bit 0) write G=Hz matrix
	 .2 (bit 1) write K=Lz matrix
	 .4 (bit 2) write H=Hx matrix
	 .8 (bit 3) write L=Lx matrix
	 .16 (bit 4) write P vector
	 Codewords file 'finC', if given, will be used to create 'G' and 'K'
	   matrices with rows of smallest possible weight.
	 In addition, mode=3.32 (just bit 5 set) in combination with
	  codewords file 'finC' forces code transformation mode.
	 Similarly, use mode=3.64 (just bit 6 set) to create DEM file '${fout}D.dem'
	 Use 'fout=' command line argument to generate file names
	 ${fout}H.mmx, ${fout}G.mmx, ${fout}L.mmx, ${fout}K.mmx, and ${fout}P.mmx
	 with 'fout=stdout' all output is sent to 'stdout'
```

### run `./vecdec --morehelp`
```sh
   Matrices used by ./src/vecdec:
	 We have a CSS code with binary generator matrices Hx=H, Hz=G,
	 and logical-operator generating matrices Lx=L and Lz=K.  These
	 matrices have 'n' columns each and satisfy orthogonality properties
		 Hx*Hz^T=0, Lx*Hz^T=0, Lx*Lz^T=0, and Lx*Lz^T=I (identity matrix).
	 We are trying to correct binary Z errors, whose probabilities
	 are given by the n-component double vector P
	  (or can be overridden with the command-line parameter 'useP').
		   A detector error model (DEM) file, in the format produced by Stim,
	 contains matrices Hx and Lx, and the error probability vector P.
	 The same code can be obtained by specifying the matrices 
	 independently, via separate files.

   The dual CSS matrix Hz can be specified instead of Lx.
	 In such a case, the internal error generator must be used
	 (an attempt to specify 'fdet' and 'fobs' files will result in an error).

   For a classical code, just give the parity check matrix Hx=H.
	 In this case G matrix is trivial (has zero rank), and
	 Lx has all rows of weight '1'.  
	 Only the internal error generator can be used for classical codes
	                                                       
   Note: detection events (syndrome bits) are given by the product 'H*e'
	       observables are given by 'L*e'
	                                                          
	 Parameter(s) used by all modes:                         
	 seed=[integer] : when negative or zero, combine provided value
		 with 'time(null)' and 'pid()' for more randomness.
	                                                       
```		 


### Additional details

- The command-line argument `seed` is used when a positive value is
  specified.  When a negative value `-x` is specified, a combination
  of `time(NULL)`, `getpid()` and `x` is used to provide maximum
  randomness.  Namely, 
  ```C
  seed = x + time(NULL) + 1000000ul * getpid();
  ```
- The command-line argument `maxW` can be used to ensure that all
  codewords of higher weight are omitted.  This is useful with
  `mode=3.32` (matrix transformation) as the complexity grows
  exponentially with the weight.  [Currently only transformations
  corresponding to codewords of weight `1`, `2`, and `3` are
  implemented.]

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

