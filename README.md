# vecdec - a simple vectorized decoder

## Error model

An error model is a collection of rows, or independent `events`, each of which is characterized by a probability `p`, a list of affected syndrome bits, and a list of affected codewords.  Can be read from a file similar to the following example:
```
# end-of line comment 
# [rows `r` in the check matrix] [number of independent codewords, `m`]
3 2
# each row is the probability followed by the list of syndrome bits, followed by the list of codewords.
# use a semicolon to separate blocks of data.
0.001 0 1 ; 0
0.01 1 2 ; 1
0.005 0 1 2 ; 0 1
# use as many rows as necessary
# rows with identical entries will be automatically combined
```

Given a syndrome vector `s`, the goal is to construct the most likely binary vector `e` such that `H*e=s`.  For each input row, the program outputs a binary vector `L*e` with the found `e`. Most common operation as a filter: syndrome vectors on `stdin`, output vectors on `stdou`.

## Syndrome vectors

Each syndrome vector is a row of of exactly `r` zerows and ones.
Rows starting with `#` are considered comments and ignored. An empty row produces an empty row on output.
