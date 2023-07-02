# vecdec - a simple vectorized decoder

## Error model

An error model is a collection of rows, or independent `events`, each of which
is characterized by a probability `p`, a list of affected syndrome bits, and a
list of affected codewords.  Can be read from a file similar to the following
example: 
```bash
# end-of line comments can be used anywhere 
# [rows `r` in the check matrix] [number of codewords, `k`] [number of entries `n`]
3 2 4
# each row is the probability followed by the list of syndrome bits, 
# followed by the list of affected codewords separated by a semicolon.
0.001 0 1 ; 0
0.01 1 2 ; 1
0.005 0 1 2 ; 0 1 
0.0001 2 0; 1 
# second semicolon can be used to start new block of data if wanted 
# use as many rows as necessary
# rows with identical entries will be automatically combined
```

Given a syndrome vector `s`, the goal is to construct the most likely binary
vector `e` such that `H*e=s`.  For each input row, the program outputs a binary
vector `L*e` with the found `e`. Most common operation as a filter: syndrome
vectors on `stdin`, output vectors on `stdout`.

## Syndrome vectors

Each syndrome vector is a row of of exactly `r` zerows and ones.  Rows starting
with `#` are considered comments and ignored. An empty input row produces an
empty row on output.

## How it works

The program processes a up to `n` syndrome vectors at a time.  The syndrome
columns are written after all columns of the original check matrix `H`.  At each
step, a column ordering `P` is randomly generated (using values of event
probabilities to help), the Gauss elimination is performed on the rows of the
combined matrix, creating the list of pivot columns `[i0, i1, ...]`, with one
entry per row.  Given the transformed syndrome column `[s0, s1, ...]`, the
output vector has values given by the list of pairs [(i0,s0), (i1,s1), ...].
The energy is calculated as the sum of LLRs for non-zero bits in `e` and
recorded, along with the sparse form of `e` if the energy is small enough.

The program stops after a sufficient number of attempts is made, and outputs the
values `L*e`, one per syndrome row, in this order.

## Libraries 

The program uses `m4ri` library for binary linear algebra.


