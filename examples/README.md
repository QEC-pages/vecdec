b# examples 

Edit the shell scripts `surf12.sh` and/or `bs.sh` to ensure that `stim`,
`python`, and `gnuplot` variables point to correct locations, then run, e.g.,

```sh
sh bs.sh
```
it should automatically create d=2 Bacon-Shor circuit, and, for every value of
`p1`, create an error model, and run `vecdec` to generate decoding data. 


# classical LDPC codes 
files `1920.1280.3.303.alist` and `PERMCODE1200.3.alist` are classical
LDPC codes taken from David MacKay's website
(https://www.inference.org.uk/mackay/CodesFiles.html)
