mcmc
==============

Ising model - by "Heat Bath" and Wolff clustering [1][2]. 

References:
-----------
[1] MCMC from Scratch, Masanori Hanada and So Matsuura <br/>
[2] [MCMC](https://github.com/masanorihanada/MCMC-Sample-Codes) <br/>

Run
-----

Output is saved to "route.dat", which has the same format as the input and can be used to
continue reordering any number of times, e.g.:

```
% cargo run -- -h

Usage: ising [OPTIONS]

Options:
  -d, --dim <D>         lattice dimensions (width,height) [default: 750,600]
  -t, --temp <T>        temperature [default: 2.3]
      --coupling_h <H>  external magnetic field [default: 0.1]
  -n, --niter <N>       iterations [default: 1000000]
  -u, --u <U>           update rule - 0: heat bath, 1: Wolff clustering [default: 1]
  -h, --help            Print help
  -V, --version         Print version
```

```
% cargo run --release 

#iteration  #avg spin   #avg energy
  10000     -0.0702     -1.3466
  20000     -0.0043     -1.3503
  30000      0.0437     -1.3444
  40000      0.0281     -1.3340
```

```
% sh ising.sh
```


![GIF](https://github.com/jesper-olsen/mcmc/blob/main/ising/Assets/animated_Wolf_h0_t2_3.gif) 

![GIF](https://github.com/jesper-olsen/mcmc/blob/main/ising/Assets/animated_Wolf_h0_t2_2.gif) 
