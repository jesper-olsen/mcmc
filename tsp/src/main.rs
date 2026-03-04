use clap::Parser;
use gnuplot::{AxesCommon, Caption, Color, Figure};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Error, Write};
use stmc_rs::marsaglia::Marsaglia;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(long = "niter", default_value_t = 5000)]
    ///number of iterations
    niter: usize,

    #[arg(short, long = "num_replicas", default_value_t = 200)]
    ///number of replicas
    num_replicas: usize, //b: usize,

    #[arg(short, long = "dbeta", default_value_t = 0.5)]
    ///step size: dbeta
    dbeta: f64,

    #[arg(short, long = "greedy", default_value_t = false)]
    ///greedy reordering
    greedy: bool,

    /// filename - x y coordinates, one pair per line
    fname: String,
}

fn main() -> Result<(), Error> {
    let args = Args::parse();
    println!(
        "fname: {} niter: {} dbeta: {} num_replicas: {} greedy {}",
        args.fname, args.niter, args.dbeta, args.num_replicas, args.greedy
    );
    let cities = read_cities(&args.fname)?;
    let mut route = if args.greedy {
        greedy_ordering(&mut cities.clone())
    } else {
        tsp(&cities, args.niter, args.dbeta, args.num_replicas)
    };

    println!("Route length: {}", calc_distance(&cities, &route));
    let fname = "route.dat";
    println!("Saving ordering to {fname}");

    let file = File::create(fname)?;
    let mut writer = BufWriter::new(file);
    for city_idx in &route {
        let (x, y) = cities[*city_idx];
        writeln!(writer, "{x} {y}")?;
    }

    // add link from end to start
    route.push(0);
    let x2: Vec<(f64, f64)> = route.iter().map(|i| cities[*i]).collect();
    plot(&x2);
    Ok(())
}

// format: one line per city, (x,y) coordinates
fn read_cities(fname: &str) -> Result<Vec<(f64, f64)>, Error> {
    let file = File::open(fname)?;
    Ok(BufReader::new(file)
        .lines()
        .filter_map(|line| {
            let l = line.ok()?;
            let mut parts = l.split_whitespace();
            let x = parts.next()?.parse().ok()?;
            let y = parts.next()?.parse().ok()?;
            Some((x, y))
        })
        .collect())
}

fn euclid_sq(a: &(f64, f64), b: &(f64, f64)) -> f64 {
    let dx = a.0 - b.0;
    let dy = a.1 - b.1;
    dx * dx + dy * dy
}

/// return index of element closest to x[0]
fn nearest(x: &[(f64, f64)]) -> usize {
    let mut imin = 0;
    let mut dmin = f64::MAX;
    for (i, e) in x[1..].iter().enumerate() {
        let d = euclid_sq(e, &x[0]);
        if d < dmin {
            dmin = d;
            imin = i + 1;
        }
    }
    imin
}

fn greedy_ordering(x: &mut [(f64, f64)]) -> Vec<usize> {
    let mut ordering: Vec<usize> = (0..x.len()).collect();
    for i in 0..x.len() - 1 {
        let j = nearest(&x[i..]);
        if i + 1 != j + i {
            x.swap(i + 1, j + i);
            ordering.swap(i + 1, j + i);
        }
    }
    ordering
}

fn calc_distance(x: &[(f64, f64)], ordering: &[usize]) -> f64 {
    let n = ordering.len();
    if n < 2 {
        return 0.0;
    }

    (0..n)
        .map(|i| {
            let city_a = x[ordering[i]];
            let city_b = x[ordering[(i + 1) % n]]; // Wraps index n back to 0

            let dx = city_a.0 - city_b.0;
            let dy = city_a.1 - city_b.1;
            (dx * dx + dy * dy).sqrt()
        })
        .sum()
}

fn tsp(x: &[(f64, f64)], niter: usize, dbeta: f64, num_replicas: usize) -> Vec<usize> {
    let ncity = x.len();
    let mut minimum_distance = f64::MAX;
    let mut minimum_ordering = vec![0; ncity];

    let mut ordering = vec![vec![0; ncity]; num_replicas];
    for e in ordering.iter_mut() {
        (0..e.len()).for_each(|icity| e[icity] = icity);
    }

    let mut rng = Marsaglia::default();
    for iter in 0..niter {
        // MCMC - 1 or more steps
        for _ in 0..1 {
            for (i, replica) in ordering.iter_mut().enumerate() {
                let (mut k, mut l) = (0, 0);
                while k == l {
                    k = (rng.uni() * ((ncity - 1) as f64)) as usize + 1;
                    l = (rng.uni() * ((ncity - 1) as f64)) as usize + 1;
                }

                // Metropolis for each replica
                let beta = (i + 1) as f64 * dbeta;
                let action_init = calc_distance(x, replica) * beta;
                replica.swap(k, l);
                let action_fin = calc_distance(x, replica) * beta; // TODO: optimise - only four edges can change when two cities are swapped...

                // Metropolis test
                if (action_init - action_fin).exp() <= rng.uni() {
                    // reject
                    replica.swap(k, l);
                }
            }
        }

        // Exchange replicas
        for i in 1..num_replicas {
            let beta1 = i as f64 * dbeta;
            let beta2 = (i + 1) as f64 * dbeta;
            let action_init =
                calc_distance(x, &ordering[i - 1]) * beta1 + calc_distance(x, &ordering[i]) * beta2;
            let action_fin =
                calc_distance(x, &ordering[i - 1]) * beta2 + calc_distance(x, &ordering[i]) * beta1;
            // Metropolis test
            if (action_init - action_fin).exp() > rng.uni() {
                ordering.swap(i - 1, i);
            }
        }

        let distance = calc_distance(x, &ordering[num_replicas - 1]);
        if distance < minimum_distance {
            minimum_distance = distance;
            minimum_ordering[..].copy_from_slice(&ordering[num_replicas - 1][..])
        }

        if iter % 500 == 0 {
            println!("iter {iter:4}    distance: {minimum_distance:.4}");
        }
    }

    println!("MO: {minimum_ordering:?}");
    minimum_ordering
}

fn plot(route: &[(f64, f64)]) {
    let mut fg = Figure::new();

    //fg.set_terminal("qt", "");

    fg.set_title("Traveling Salesman Problem");
    fg.axes2d().set_x_label("X", &[]).set_y_label("Y", &[]);

    fg.axes2d().points(
        route.iter().map(|&(x, _)| x),
        route.iter().map(|&(_, y)| y),
        &[Caption("Route"), Color(gnuplot::RGBString("blue"))],
    );

    fg.axes2d().lines(
        route.iter().map(|&(x, _)| x),
        route.iter().map(|&(_, y)| y),
        &[Caption("Route"), Color(gnuplot::RGBString("red"))],
    );

    fg.show().unwrap();
}

#[cfg(test)]
mod tests {
    #[test]
    fn tsp_test() {
        let mut x = crate::read_cities(&"assets/cities10.dat").unwrap();
        let r = crate::tsp(&x, 5000, 0.5, 200);
        let expect = vec![0, 5, 4, 3, 2, 1, 9, 6, 7, 8];
        assert_eq!(r, expect);
    }

    #[test]
    fn test_greedy() {
        let mut x = vec![(0.0, 0.0), (3.0, 3.0), (2.0, 2.0), (1.0, 1.0)];
        let ordering = crate::greedy_ordering(&mut x);
        assert_eq!(ordering, vec![0, 3, 2, 1]);
    }
}
