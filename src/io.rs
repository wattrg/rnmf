use crate::mesh::cartesian::*;
use std::fs::File;
use std::io::prelude::*;
use vtkio::model::*;

pub use indicatif::{ProgressBar, ProgressStyle};

pub trait RnmfProgressBar {
    fn create (increments: usize) -> ProgressBar;
    fn increment (&mut self, size: usize);
    fn finish (&mut self);
}


impl RnmfProgressBar for ProgressBar {
    fn create (increments: usize) -> ProgressBar {
        let sty = ProgressStyle::default_bar()
            .template("{spinner:.yellow} [{elapsed_precise}] [{bar:40.green/yellow}] {pos:>7}/{len:7} (eta: {eta_precise})")
            .progress_chars("==> ");
        let pb = ProgressBar::new(increments as u64);
        pb.set_style(sty);
        pb
        
    }
    fn increment (&mut self, size: usize) {
        self.inc(size as u64);
    }
    fn finish(&mut self){
        self.finish_with_message("done");
    }
}

pub fn write_csv(name: &str, mut data: CartesianDataFrame){
    let mut file = File::create(format!("{}.csv",name)).unwrap();
    let n = data.underlying_mesh.n[0] as isize;
    print!("  writing csv...");
    for (indx,val) in data.into_iter().enumerate_index(){
        file.write_all(format!("{}", val).as_bytes()).unwrap();
        if indx.0 == n-1 {
            file.write_all(b"\n").unwrap();
        }
        else {
            file.write_all(b",").unwrap();
        }
    }
    println!(" done");
}

// pub fn write_vtk(iter: usize, data: CartesianDataFrame){

// }