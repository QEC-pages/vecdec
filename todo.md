---
title: Notes on `vecdec` 
author: Leonid P. Pryadko
date: 2023/07/11
---

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
it tends to clear small-weight error efficiently, after ${\cal O}(1)$ steps with
overall complexity linear, so that the remaining errors can be cleaned up with a
secondary decoder.  Ordered-statistics decoder (OSD) has been suggested for this
purpose, which is essentially an RW decoder using the list of *aposteriori*
probabilities returned by BP.

## Actual to-do list 
- [ ] Implement recursive enumeration of vectors from the information set.
- [ ] Implement hashing storage for small-weight errors (say, up to weight 2)
      and the corresponding syndromes for fast lookup decoding.
- [ ] Implement hashing storage for errors corresponding to different syndrome
      vectors.  Perhaps use `vecdec` in *scalar* mode for this, working with
      just one (or a few) syndrome vectors at a time, and keeping the
      corresponding hashing storage separate, so that only the error vectors
      (in sparse form?) and the corresponding probabilities need to be stored.
- [ ] Implement BP decoding with OSD
- [ ] Come up with alternative simulation mechanisms, e.g., in the regime of
      small error probabilities $p$.
