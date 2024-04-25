tsp
==============

Ising model - by "Heat Bath" and Wolff clustering [1][2]. 

References:
-----------
[1] MCMC from Scratch, Masanori Hanada and So Matsuura <br/>
[2] [MCMC](https://github.com/masanorihanada/MCMC-Sample-Codes) <br/>

Run
-----

Greedy nearest neighbour ordering:
```
% cargo run --release -- --greedy assets/cities734.dat 

```

Output is saved to "route.dat", which has the same format as the input and can be used to
continue reordering any number of times, e.g.:

```
% cargo run --release 

#iteration  #avg spin   #avg energy
  10000     -0.0702     -1.3466
  20000     -0.0043     -1.3503
  30000      0.0437     -1.3444
  40000      0.0281     -1.3340

```


![GIF](https://raw.githubusercontent.com/jesper-olsen/mcmc/master/ising/animated_Wolf_h0_t2_3.gif) 

