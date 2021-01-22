use crate::mesh::cartesian2d::*;
use std::fs::File;
use std::io::prelude::*;
use vtkio::model::*;
pub use indicatif::{ProgressBar, ProgressStyle};
use colored::*;

/// Extension to `indicatif::ProgressBar` to simplify its use
pub trait RnmfProgressBar {
    fn create (increments: usize) -> ProgressBar;
    fn increment (&mut self, size: usize);
    fn finish (&mut self);
}


impl RnmfProgressBar for ProgressBar {
    fn create (increments: usize) -> ProgressBar {
        let sty = ProgressStyle::default_bar()
            .template("Progress: [{elapsed_precise}] [{bar:40.green/yellow}] {pos:>7}/{len:7} (eta: {eta_precise})")
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



pub fn write_vtk(loc: &str, 
                 name: &str, 
                 iter: usize, 
                 data: CartesianDataFrame2D, 
                 pb: &ProgressBar){
    let file_name = format!("{}_{:0>8}.vti", name, iter);
    pb.println(format!("    Writing vtk file '{}'", file_name));
    let file = format!("{}_{}", loc, file_name);

    let mesh = &data.underlying_mesh;

    // generate information required for vtk file format
    let extent = Extent::Ranges([0..=mesh.n[0] as i32, 0..=mesh.n[1] as i32 as i32, 0..=1]);
    let version = Version::new((0,1));
    let byte_order = ByteOrder::LittleEndian;

    // generate iobuffer for the data
    let psi_iobuffer = IOBuffer::F64(data.get_valid_cells());

    // create the DataArrayBase
    let psi_data_array_base = DataArrayBase {
        name: String::from("psi"),
        elem: ElementType::Scalars{num_comp: 1, lookup_table: Option::None},
        data: psi_iobuffer, 
    };

    // create the attribute/s
    let psi_attribute = Attribute::DataArray(DataArray::from(psi_data_array_base));
    let vtk_attributes = Attributes{
        cell: vec![psi_attribute.clone()],
        point: vec![],
    };

    // piece of image data
    let vtk_image_data_piece = Piece::Inline(Box::new(
        ImageDataPiece {
            extent: extent.clone(),
            data: vtk_attributes,
        })
    );

    // data set
    let vtk_data_set = DataSet::ImageData{
        extent,
        origin: [mesh.lo[0] as f32, mesh.lo[1] as f32, 0.0],
        spacing: [mesh.dx[0] as f32, mesh.dx[1] as f32, 1.0],
        meta: Option::None,
        pieces: vec![vtk_image_data_piece]
    };

    // overall vtk model
    let vtk_model = Vtk {
        version,
        title: String::from(""),
        byte_order,
        data: vtk_data_set,
    };



    // write file
    match vtkio::export(vtk_model, file){
        Ok(_) => {}
        Err(vtkio::Error::IO(err))  => {
            pb.println(format!("      {} failed writing vtk file: {}", "Warning:".yellow(), err));
        }
        Err(vtkio::Error::Write(err)) => {
            pb.println(format!("      {} failed writing vtk file: {:?}", "Warning:".yellow(), err));
        }
        Err(vtkio::Error::XML(err)) => {
            pb.println(format!("      {} failed writing vtk file: {}", "Warning:".yellow(), err));
        }
        Err(vtkio::Error::UnknownFileExtension(err)) => {
            pb.println(format!("      {} failed writing vtk file, invalid extension: {:?}", "Warning:".yellow(), err));
        }
        _ => {
            pb.println(format!("      {} un-expected error writing vtk file", "Warning:".yellow()));
        }
    }



}