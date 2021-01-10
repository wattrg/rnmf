use crate::mesh::cartesian2d::*;
use std::fs::File;
use std::io::prelude::*;
//use vtkio::model::*;
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

pub fn write_csv_2d(name: &str, mut data: CartesianDataFrame2D){
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

// pub fn write_vtk(name: &str, iter: usize, data: CartesianDataFrame2D){

//     let file = File::create(format!("{}_{}.csv", name, iter)).unwrap();
//     let mesh = data.underlying_mesh;

//     // compute the extent of the mesh
//     let mut hi_2d: i32 = 0;
//     let mut hi_3d: i32 = 0;
//     if mesh.dim >= 2 {
//         hi_2d = mesh.n[1] as i32;
//         if mesh.dim == 3{
//             hi_3d = mesh.n[2] as i32;
//         }
//     }
//     let extent = Extent::Ranges([0..=mesh.n[0] as i32, 0..=hi_2d as i32, 0..=hi_3d]);
//     let vtk_data = Attribute::scalars("psi", 1).add_field_data(FieldArray::new("A",1));
//     let vtk = Attributes::new();
// }