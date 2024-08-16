
## 2024/08/15

Implemented a cluster-based pre-decoder using hash look-up tables.
With this, one can dramatically reduce the number of errors to be
processed with RIS or BP decoder.  Also, in the case of BP, the net
error rate may go down.

The new command-line parameters are `uW` (max error cluster weight,
default `2`), `uR` (max graph distance between set bits in a cluster,
default `4` -- set it to `0` to have no limit -- but it would be too
slow on large codes), and `maxU` (limit on the number of error vectors
in hash), default `0` (no limit).

If you set `uW=-1`, you get more or less the old behavior.

with `uW=0`, no clusters will be generated, but the program would skip
0-weight error vectors (it was already done with BP decoders mode=1).

with `uW>1`, the stored syndromes will be used to decompose the
syndrome.  The graph-based algorithms are pretty efficient (they were
invented for percolation theory), so with small error rates, the
number of syndrome vectors that actually need decoding goes down
dramatically.
