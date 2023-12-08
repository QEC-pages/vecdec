---
title: Notes on `vecdec` 
author: Leonid P. Pryadko
date: 2023/07/11
bibliography: refs.bib
---

# Conventional decoding algorithms 

## ML list decoding 

Given a sufficient number of decoding steps (different column permutations $P$
and different information sets of weight up to $t$), the RW-$t$ algorithm
generates a large list of vectors matching the syndromes.  The idea of the
*maximum-likelihood* (ML) *list decoding* is to use a small-weight subset of the
generated error vectors to approximately calculate the entropy of each outcome,
i.e., compare the free energies $\ln Z(e)$, where $Z(e)$ is the sum of
probabilities of all errors degenerate with $e$.  Secondly, if we also count the
number of `hits` for every error vector, we can estimate the probability of
missing the "correct" vector of weight smaller than $e_{\rm min}$ found so far.

In practice, we want to implement hashing storage to keep $[s,e,x,p(e)]$, where
$p(e)$ is the error probability and $x\equiv Le$ is the logical operator
outcome.  Then, for a given syndrome $s$, we should compare $Z_x(s)\equiv
\sum_{e:s,x}p(e)$ for different values of $x$ and choose the largest.  After
decoding, some of the entries may be moved to the immediate decoding list
containing the pairs $[s,e_{\rm ML}]$.  In practice, it is not clear how much
memory would be required.  With sufficiently small median error probability
$\bar p$, it may do to store all errors of weight up to $3n\bar p$, say.

## Information set enumeration (up to level `t`)

For a given column permutation $P$ the program constructs the reduced row
echelon form (RREF) of the parity check matrix $H$ with the help of the Gauss
algorithm.  The syndrome vectors are transposed simultaneously as columns of the
matrix $S$.  The set of pivot columns in $H$ give the index set $J$; its
complement $I$ in classical coding theory is known as the *information set*.  Up
to a column permutation, the RREF of $H$ has the form $(E,A)$, where $E$ is the
identity matrix on the set $J$ (of size $r\equiv\mathop{\rm rank}H$) and $A$ is
a matrix of size $r$ by $n-r$.  The dual matrix has the form $(A^T,E')$.

For a given (transformed) syndrome vector $s$, the 0th level random-window
decoding (RW-0) amounts to choosing positions from the index set $J$ to match
the non-zero bits in a transformed syndrome vector $s$.  The intuition is that
the vector weight is restricted by $r$, i.e., we are preferentially drawing
small-weight vectors.  Furthermore, if the actual error $e$ that happened is
supported on $J$, it is easy to verify that the decoding is accurate.  The
corresponding probability is high for errors of small weight.

With level-$t$ random-window decoding (RW-$t$), up to $t$ non-zero positions
from the information set $I$ are chosen, the corresponding columns are added to
the syndrome.  This is supposed to address the situations where all but $t'\le
t$ bits of the error are supported in $J$.  The complexity of enumerating
information set vectors of weight $w$ is ${n-r\choose w}$, while the benefit
(improved success probability, i.e., decreased complexity) does not scale as
rapidly; typically one uses $t\le 2$ in order to preserve the overall complexity
scaling of each RW step.

Unfortunately, especially with LDPC codes, the probabilities of encountering
different information sets (the numbers of permutations $P$ that give $I$) may
differ a lot.  Second difficulty is that of non-uniform probabilities $p_j$ for
errors in different columns of $H$.  Ideally, we should come up with a random
process which generates the column according to a given set of probabilities in
the error model.

## Belief propagation with OSD

Belief propagation (BP) is notoriously difficult for highly-degenerate quantum
codes, as it gets stuck often at degenerate configurations.  On the other hand,
it tends to clear small-weight errors efficiently, after ${\cal O}(1)$ steps with
overall complexity linear, so that the remaining errors can be cleaned up with a
secondary decoder.  Ordered-statistics decoder (OSD) has been suggested for this
purpose, which is essentially an RW decoder using the list of *aposteriori*
probabilities returned by BP.

# Logical fail rate predictors 

## Asymptotic minimum-weight fail rate 

With a given detector error model (specified by detector-fault incidence matrix
$J$, observable-fault incidence matrix $L$, and a set of independent
probabilities $p_i$ for each column of $J$), we would like to know the logical
fail rate of the minimum-weight decoder in the regime of small probabilities.
For simplicity, let us first assume all error probabilities be the same and
equal to $p$.  Then, according to Fowler [@Fowler-2013], the asymptotic logical
error rate at small $p$ has the form $P_L=Bp^{\lceil d/2\rceil}$, where the
prefactor $B$ can be calculated by enumerating the fault paths ("codewords") of
weight $d$.  Moreover, while different fault paths are not uncorrelated, in the
case of a surface code, ignoring the correlations gives $B$ with around 12%
accuracy.  With $p_i$ different, the weight of a single codeword $c$ can be
approximated as $1/2\prod_{i:c_i\neq 0}2[p_i(1-p_i)]^{1/2}$.

## Weighted Hashimoto matrix approach 

In the case of graph errors, summation over paths can be done approximately with
the help of weighted Hashimoto matrix, defined in terms of directed edges
(arcs), $$ M_{a=i\to j,b=i'\to j'}= W_a \delta_{j,i'}(1-\delta_{i,j'}),\quad
W_a\equiv 2[p_a(1-p_a)]^{1/2}, $$ where $W_a$ is the weight of a single edge.
Namely, to enumerate the sum of all non-backtracking paths of length $m$
starting at arc $a$ and ending at arc $b$, we just write $[M^m]_{a,b}W_b$ (no
summation).  Notice that at sufficiently small error probabilities $p_a$, the
logical fail probability is going to be dominated by the leading-order terms with $m=d$.

In practice, when $H$ corresponds to a graph (e.g., in the case of a surface
code), we just have to construct a vector $x$ of in-arcs and a vector $y$ of
weighted out-arcs, and calculate $P_L=\alpha x M^d y$, where the prefactor
coefficient $\alpha$ can be used to estimate the corrections, e.g., due to the
distribution of error probabilities $p_a$.

## Estimate the correction due to overlaps

The issue with the expansion over fault-lines (codewords) is non-locality.
Indeed, $p_0=(1-p)^n\le e^{-np}$ can be small unless $pn\to0$, and $q_{t+1}$ is
also strongly suppressed by similar factors, to the point of not being useful.
To construct a more useful expansion, consider a local notion of failure.
Namely, given an index set $J$ and a given error vector $e$, consider the
probability $P_J(e)$ of a minimum-energy fault strictly on $J$.  That is, we
assume that there exists a codeword $c$ (irreducible or not) supported on $J$
such that $P(e+c)>P(e)$, while for every irreducible codeword $c'$ not entirely
supported on $J$, $p(e+c')\le p(e)$.

The function $P_A(e)$ has nice properties: 
- It is increasing, meaning that if $A\subset B$, $P_A(e)\le P_B(e)$
- When sets $A$ and $B$ have partial overlaps, $A\neq A\cap B\neq\emptyset$ and
  $B\neq A \cap B$, as long as $P_{A\cap B}(e)=0$, $P_{A\cup
  B}(e)=P_{A}(e)+P_{B}(e)$ **(???)**

Furthermore, consider a set of non-trivial irreducible binary codewords $c$ such
that $Hc=0$ and $Lc\neq0$; necessarily, $\mathop{\rm wgt} c\ge d$.  Irreducible
means that $c$ cannot be separated into a pair of non-zero binary vectors with
disjoint supports and zero syndromes.  In the case of a surface code,
irreducible $c$ is a homologically non-trivial chain without self-intersections.

## Actual to-do list 
- [x] Implement recursive enumeration of vectors from the information set.
- [x] Implement reading of externally generated simulation data (e.g., from
      Stim) using separate files with detection events and the corresponding
      observables.  Here is a sample command-line to produce binary vectors in
      `01` format: 

```bash
stim sample_dem \
  --shots 5 \
  --in example.dem \
  --out dets.01 \
  --out_format 01 \
  --obs_out obs_flips.01 \
  --obs_out_format 01
```
- [x] Implement a mode for finding the list of most likely zero-syndrome error
      vectors.
- [ ] Implement hashing storage for small-weight errors (say, up to weight 2)
      and the corresponding syndromes for fast lookup decoding.  Subsequently,
      use decoding results to add most likely vectors to this list.
- [ ] Implement hashing storage for near-ML decoding using lists of small-weight
      errors.  Perhaps use `vecdec` in *scalar* mode for this, working with just
      one (or a few) syndrome vectors at a time, and keeping the corresponding
      hashing storage separate, so that only the error vectors (in sparse form?)
      and the corresponding probabilities need to be store.  E.g., try using the
      [`uthash` library](https://github.com/troydhanson/uthash).
- [ ] Implement BP decoding with OSD
  - [x] Actual BP steps 
  - [x] Add `BoxPlus()` from `it++` library
  - [ ] Add error estimation
  - [ ] BP with randomization
  - [ ] BP with Freezing / Stabilizer inactivation (Savin et al)
  - [ ] BP version by Kung, Kuo, Lai (http://arXiv.org/abs/2305.03321) and/or Kuo, Lai (http://arXiv.org/abs/2104.13659)
  - [ ] Add OSD
  - [ ] Stopping sets analysis?
  - [ ] Initial BP acceleration ? 
  - [ ] Other tricks from papers by Narayanan; Kuo+Lai; Roffe; Valentin Savin
- [ ] Come up with alternative simulation mechanisms, e.g., in the regime of
      small error probabilities $p$.
- [ ] Hashimoto matrix implementation steps
  - [ ] Given a DEM in a *graph form,* construct the sparse weighted Hashimoto
        matrix (of size $2n\times 2n$).
  - [ ] Construct an in-vector $x$ and a weighted out-vector $y$
  - [ ] Calculate the corresponding logical fault-rates
  - [ ] Estimate the fudge-factor $\alpha$, given the statistics of error probabilities
  - [ ] Estimate the effect of correlations between the trajectories.
- [ ] Code transformations reducing the degeneracy for `mode=3`
  - [ ] Remove `w=3` degeneracies (rows of weight 3 in `G`)
  - [ ] Remove `w=4` degeneracies (rows of weight 4 in `G`)
  - [ ] Code for arbitrary row weight (exponentially large matrices may result)

