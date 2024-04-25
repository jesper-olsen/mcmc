// 2d Ising model - heat bath and  Wolff algorithm
/// U.Wolff, Collective Monte Carlo updating for spin systems. Phys.Rev.Lett.62(4),361(1989)
use clap::Parser;
use core::str::FromStr;
use std::fmt;
use std::fs::File;
use std::io::{BufWriter, Write};

pub mod marsaglia;

const COUPLING_J: f64 = 1.0; // constant

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(short, long="dim", default_value_t = String::from("750,600"))]
    ///lattice dimensions (width,height)
    d: String,

    #[arg(short, long = "temp", default_value_t = 2.30)]
    ///temperature
    t: f64,

    #[arg(long = "coupling_h", default_value_t = 0.1)]
    ///external magnetic field
    h: f64,

    #[arg(short, long = "niter", default_value_t = 1_000_000)]
    ///iterations
    n: usize,

    #[arg(short, long, default_value_t = 1)]
    ///update rule - 0: heat bath, 1: Wolff clustering
    u: usize,
}

struct Lattice {
    spin: Vec<bool>,
    xdim: usize,
    ydim: usize,
}

impl Lattice {
    pub fn new(xdim: usize, ydim: usize) -> Self {
        let spin = vec![true; xdim * ydim];
        Lattice { spin, xdim, ydim }
    }

    pub fn get_spin(&self, x: usize, y: usize) -> i32 {
        let i = y * self.xdim + x;
        if self.spin[i] {
            1
        } else {
            -1
        }
    }

    pub fn set_spin(&mut self, x: usize, y: usize, v: bool) {
        let i = y * self.xdim + x;
        self.spin[i] = v;
    }

    pub fn flip_spin(&mut self, x: usize, y: usize) {
        let i = y * self.xdim + x;
        self.spin[i] = !self.spin[i];
    }

    pub fn total_spin(&self) -> i32 {
        self.spin.iter().map(|&b| if b { 1 } else { -1 }).sum()
    }

    pub fn action(&self, temperature: f64, coupling_h: f64) -> f64 {
        let mut sum = 0;
        for ix in 0..self.xdim {
            let ixp1 = (ix + 1) % self.xdim; // boundary condition
            for iy in 0..self.ydim {
                let iyp1 = (iy + 1) % self.ydim; // boundary condition
                sum += self.get_spin(ix, iy) * (self.get_spin(ixp1, iy) + self.get_spin(ix, iyp1));
            }
        }
        (sum as f64 * COUPLING_J + self.total_spin() as f64 * coupling_h) / -temperature
    }
}

impl fmt::Display for Lattice {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for iy in 0..self.ydim {
            for ix in 0..self.xdim {
                write!(f, "{} ", self.get_spin(ix, iy))?
            }
            writeln!(f)?
        }
        writeln!(f)
    }
}

/// Calculation of the probability when the spin at (ix,iy) is updated
fn heat_bath_prob(
    lat: &Lattice,
    (ix, iy): (usize, usize),
    temperature: f64,
    coupling_h: f64,
) -> f64 {
    let ixp1 = (ix + 1) % lat.xdim; // boundary condition!
    let iyp1 = (iy + 1) % lat.ydim;
    let ixm1 = (ix + lat.xdim - 1) % lat.xdim;
    let iym1 = (iy + lat.ydim - 1) % lat.ydim;
    let v: i32 = [(ixp1, iy), (ix, iyp1), (ixm1, iy), (ix, iym1)]
        .iter()
        .map(|(x, y)| lat.get_spin(*x, *y))
        .sum();
    let ep = -(coupling_h + v as f64 * COUPLING_J) / temperature;
    let em = -ep;
    (-ep).exp() / ((-ep).exp() + (-em).exp())
}

fn make_cluster(
    lat: &Lattice,
    rnd: &mut marsaglia::Marsaglia,
    temperature: f64,
) -> (i32, Vec<(usize, usize)>) {
    let mut in_out = vec![true; lat.xdim * lat.ydim];
    let rand_site = rnd.uni() * lat.xdim as f64 * lat.ydim as f64;
    let (ix, iy) = (rand_site as usize / lat.ydim, rand_site as usize % lat.ydim);
    in_out[iy * lat.xdim + ix] = false;
    let mut i_cluster = Vec::new();
    i_cluster.push((ix, iy));
    let prob = 1.0 - (-2.0 * COUPLING_J / temperature).exp();
    let mut k = 0;
    while k < i_cluster.len() {
        let (ix, iy) = i_cluster[k];
        let ixp1 = (ix + 1) % lat.xdim; // boundary condition!
        let iyp1 = (iy + 1) % lat.ydim;
        let ixm1 = (ix + lat.xdim - 1) % lat.xdim;
        let iym1 = (iy + lat.ydim - 1) % lat.ydim;
        for (x, y) in [(ixp1, iy), (ix, iyp1), (ixm1, iy), (ix, iym1)] {
            if lat.get_spin(x, y) == lat.get_spin(ix, iy)
                && in_out[y * lat.xdim + x]
                && rnd.uni() < prob
            {
                i_cluster.push((x, y));
                in_out[y * lat.xdim + x] = false;
            }
        }
        k += 1;
    }
    (lat.get_spin(ix, iy), i_cluster)
}

fn wolf_update(
    lat: &mut Lattice,
    rng: &mut marsaglia::Marsaglia,
    temperature: f64,
    coupling_h: f64,
) {
    let (cspin, i_cluster) = make_cluster(&lat, rng, temperature);
    let metropolis = rng.uni();
    if (-2.0 / temperature * coupling_h * cspin as f64 * i_cluster.len() as f64).exp() > metropolis
    {
        for (ix, iy) in i_cluster.iter() {
            lat.flip_spin(*ix, *iy);
        }
    }
}

fn heat_bath_update(lat: &mut Lattice, rng: &mut marsaglia::Marsaglia, temp: f64, coupling_h: f64) {
    let rand_site = rng.uni() * lat.xdim as f64 * lat.ydim as f64;
    let (ix, iy) = (rand_site as usize / lat.ydim, rand_site as usize % lat.ydim);
    lat.set_spin(
        ix,
        iy,
        heat_bath_prob(&lat, (ix, iy), temp, coupling_h) > rng.uni(),
    ); // metropolis
}

fn parse_number_pair<T: FromStr>(s: &str, separator: char) -> Result<(T, T), &str> {
    let parts: Vec<&str> = s.split(separator).collect();

    if parts.len() != 2 {
        return Err("Invalid format. Use NUMBER1,NUMBER2");
    }

    Ok((
        T::from_str(parts[0]).map_err(|_| "Invalid number")?,
        T::from_str(parts[1]).map_err(|_| "Invalid number")?,
    ))
}

fn parse_pair<T: FromStr>(s: &str, label: &str) -> (T, T) {
    match parse_number_pair::<T>(s, ',') {
        Ok((x, y)) => (x, y),
        Err(msg) => {
            println!("failed to parse {label}: {msg}");
            std::process::exit(1)
        }
    }
}

fn main() {
    let args = Args::parse();
    let (width, height) = parse_pair::<usize>(&args.d, "dimensions");
    let temperature = args.t;
    let coupling_h = args.h;
    let niter = args.n;
    let log_interval = niter / 10;
    let algorithm = args.u;

    let mut rng = marsaglia::Marsaglia::new(12, 34, 56, 78);
    let mut lat = Lattice::new(width, height);
    println!("#iteration  #avg spin   #avg energy");
    for iter in 1..=niter {
        match algorithm {
            0 => heat_bath_update(&mut lat, &mut rng, temperature, coupling_h),
            _ => wolf_update(&mut lat, &mut rng, temperature, coupling_h),
        }

        if iter % (niter / 1000) == 0 {
            let energy = lat.action(temperature, coupling_h) * temperature;
            let n = (lat.xdim * lat.ydim) as f64;
            println!(
                "{iter:>7}  {:>10.4}  {:>10.4}",
                lat.total_spin() as f64 / n,
                energy / n
            );
        }
        if log_interval > 0 && iter % log_interval == 0 {
            let fname = format!("Frames/frame{}.txt", iter / log_interval);
            let file = File::create(fname).unwrap();
            let mut writer = BufWriter::new(file);
            writeln!(writer, "{lat}").unwrap();
        }
    }
}
