use clap::Parser;
use gnuplot::{Figure, Caption, Color, AxesCommon};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

pub mod marsaglia;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    //#[arg(short, long = "niter", default_value_t = 5000)]
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

fn main() {
    let args = Args::parse();
    println!(
        "fname: {} niter: {} dbeta: {} num_replicas: {} greedy {}",
        args.fname, args.niter, args.dbeta, args.num_replicas, args.greedy
    );

    let mut x = read_cities(&args.fname);
    let route = if args.greedy {
        let mut route = greedy_ordering(&mut x.clone());
        x.push(x[0]);
        route.push(0);
        route
    } else {
        x.push(x[0]);
        tsp(&x, args.niter, args.dbeta, args.num_replicas)
    };     

    println!("Route length: {}", calc_distance(&x, &route));
    let x2: Vec<(f64,f64)> = route.iter().map(|i| x[*i]).collect();
    plot(&x, &x2);
    let fname = "route.dat";
    println!("Saving ordering to {fname}");
    let file = File::create(fname).unwrap();
    let mut writer = BufWriter::new(file);
    for icity in 0..x.len()-1 {
        let p = &x[route[icity]];
        writeln!(writer, "{} {}", p.0, p.1).unwrap();
    }
}

fn read_cities(fname: &str) -> Vec<(f64, f64)> {
    let file = File::open(fname).unwrap();
    let reader = BufReader::new(file);
    let mut x: Vec<(f64, f64)> = Vec::new(); // city coordinates - 2d plane
    for r in reader.lines() {
        let line = r.expect("Failed to read line");
        let mut iter = line.split_whitespace();
        let (f0, f1) = (
            iter.next().unwrap().parse().unwrap(),
            iter.next().unwrap().parse().unwrap(),
        );
        x.push((f0, f1));
    }
    x
}

fn euclid_sq(a: &(f64, f64), b: &(f64, f64)) -> f64 {
    let dx = a.0 - b.0;
    let dy = a.1 - b.1;
    dx * dx + dy * dy
}

/// return index of element closest to x[0]
fn nearest(x: &[(f64,f64)]) -> usize {
    let mut imin = 0;
    let mut dmin = f64::MAX;
    for (i,e) in x[1..].iter().enumerate() {
        let d = euclid_sq(e,&x[0]);
        if d<dmin {
            dmin = d;
            imin = i+1;
        }
    }
    imin
}

fn greedy_ordering(x: &mut [(f64,f64)]) -> Vec<usize> {
    let mut ordering: Vec<usize> = (0..x.len()).map(|i| i).collect();
    for i in 0..x.len()-1 {
        let j = nearest(&x[i..]);
        if i+1!=j+i {
            x.swap(i+1,j+i);
            ordering.swap(i+1,j+i);
        }
    }
    ordering
}

fn calc_distance(x: &[(f64, f64)], ordering: &[usize]) -> f64 {
    ordering
        .windows(2)
        .map(|pair| {
            let (i, j) = (pair[0], pair[1]);
            let (r1, r2) = (x[i].0 - x[j].0, x[i].1 - x[j].1);
            (r1 * r1 + r2 * r2).sqrt()
        })
        .sum()
}

fn tsp(x: &Vec<(f64, f64)>, niter: usize, dbeta: f64, num_replicas: usize) -> Vec<usize> {
    let ncity = x.len() - 1;
    let mut minimum_distance = 0.0;
    let mut minimum_ordering = vec![0; ncity + 1];

    let mut ordering = vec![vec![0; ncity + 1]; num_replicas];
    for e in ordering.iter_mut() {
        (0..e.len()).for_each(|icity| e[icity] = icity);
    }

    let mut rng = marsaglia::Marsaglia::new(12, 34, 56, 78);
    for iter in 0..niter {
        // MCMC - 1 or more steps
        for _ in 0..1 {
            for (i,replica) in ordering.iter_mut().enumerate() {
                let (mut k, mut l) = (0, 0);
                while k == l {
                    k = (rng.uni() * ((ncity - 1) as f64)) as usize + 1;
                    l = (rng.uni() * ((ncity - 1) as f64)) as usize + 1;
                }

                // Metropolis for each replica
                let beta = (i + 1) as f64 * dbeta;
                let action_init = calc_distance(&x, replica) * beta;
                replica.swap(k, l);
                let action_fin = calc_distance(&x, replica) * beta; // TODO: optimise

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
            let action_init = calc_distance(&x, &ordering[i - 1]) * beta1
                + calc_distance(&x, &ordering[i]) * beta2;
            let action_fin = calc_distance(&x, &ordering[i - 1]) * beta2
                + calc_distance(&x, &ordering[i]) * beta1;
            // Metropolis test
            if (action_init - action_fin).exp() > rng.uni() {
                ordering.swap(i - 1, i);
            }
        }

        let distance = calc_distance(&x, &ordering[num_replicas - 1]);
        if iter == 0 || distance < minimum_distance {
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

fn plot(cities: &Vec<(f64,f64)>, route: &Vec<(f64,f64)>) {
    let mut fg = Figure::new();

    fg.set_terminal("qt", "");

    fg.set_title("Traveling Salesman Problem");
    fg.axes2d()
      .set_x_label("X",&[])
      .set_y_label("Y",&[]);

    fg.axes2d()
      .points(
        cities.iter().map(|&(_, y)| y),
        cities.iter().map(|&(x, _)| x),
        &[Caption(""), Color("blue"),],
    );

    fg.axes2d()
      .lines(
        route.iter().map(|&(_, y)| y),
        route.iter().map(|&(x, _)| x),
        &[Caption("Route"), Color("red"),],
    );

    fg.show().unwrap();
}


#[cfg(test)]
mod tests {
    #[test]
    fn tsp_test() {
        let mut x = crate::read_cities(&"assets/cities10.dat");
        x.push(x[0]); 
        let r = crate::tsp(&x, 5000, 0.5, 200);
        let expect = vec![0, 5, 4, 3, 2, 1, 9, 6, 7, 8, 10];
        assert_eq!(r, expect);
    }

    #[test]
    fn test_greedy() {
        let mut x = vec![(0.0, 0.0), (3.0, 3.0), (2.0, 2.0), (1.0, 1.0), (0.0,0.0)];
        let ordering=crate::greedy_ordering(&mut x);
        assert_eq!(ordering, vec![0,4,3,2,1]);
    }
}
