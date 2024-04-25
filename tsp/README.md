tsp
==============

Travelling Salesman Problem - by replica exchange MCMC [1][2]. 

City examples (lattitude,longitude) in assets are from [3]:
* cities16862: Italy,
* cities10639: Finland, 
* cities734: Uruguay,
* cities194: Qatar

References:
-----------
[1] MCMC from Scratch, Masanori Hanada and So Matsuura <br/>
[2] [MCMC](https://github.com/masanorihanada/MCMC-Sample-Codes) <br/>
[3] [som-tsp](https://github.com/diego-vicente/som-tsp) <br/>

Run
-----

Greedy nearest neighbour ordering:
```
% cargo run --release -- --greedy assets/cities734.dat 

```

Output is saved to "route.dat", which has the same format as the input and can be used to
continue reordering any number of times, e.g.:

```
% cargo run --release -- --niter 10000 --dbeta 0.5 --num_replicas 200 route.dat

```

![PNG](https://raw.githubusercontent.com/jesper-olsen/mcmc/master/tsp/cities734.png) 

