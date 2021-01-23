use crate::mesh::cartesian2d::*;
//use std::fs::File;
//use std::io::prelude::*;
use vtkio::model::*;
pub use indicatif::{ProgressBar, ProgressStyle};
use colored::*;
use crate::Real;
use std::convert::TryFrom;

pub type OutputCallBack = fn(&CartesianDataFrame2D) -> Result<VtkData, String>;

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




pub struct VtkData {
    data: Vec<Real>,
    data_type: ElementType,
}

impl TryFrom<CartesianDataFrame2D> for VtkData {
    type Error = String;
    fn try_from(df: CartesianDataFrame2D) -> Result<Self, Self::Error>{
        match df.n_comp {
            1 => {
                Ok(VtkData{
                    data: df.get_valid_cells(),
                    data_type: ElementType::Scalars{num_comp: 1, lookup_table: Option::None},
                })
            }
            2 => {
                let mut reshaped_data: Vec<Real> = Vec::new();
                for i in 0..df.underlying_mesh.n[0] as isize{
                    for j in 0..df.underlying_mesh.n[1] as isize {
                        reshaped_data.push(df[(i,j,0)]);
                        reshaped_data.push(df[(i,j,1)]);
                        reshaped_data.push(0.0);
                    }
                }
                Ok(VtkData{
                    data: reshaped_data,
                    data_type: ElementType::Vectors,
                })
            }
            3 => {
                let mut reshaped_data: Vec<Real> = Vec::new();
                for i in 0..df.underlying_mesh.n[0] as isize{
                    for j in 0..df.underlying_mesh.n[1] as isize {
                        reshaped_data.push(df[(i,j,0)]);
                        reshaped_data.push(df[(i,j,1)]);
                        reshaped_data.push(df[(i,j,2)]);
                    }
                }
                Ok(VtkData{
                    data: reshaped_data,
                    data_type: ElementType::Vectors,
                })
            }
            _ => {Err(format!("Incorrect number of components, should be 1, 2, or 3, not {}", df.n_comp))}
        }
    }
}

/// Write a vti file
/// Required arguments
///     loc:  folder to place file in
///     name: name of the file to write (iteration number will be appended by this function) 
///     df:   the data frame containing the data to use to generate the output (this will eventually become a data frame container)
///     data: the variables to write to the file
///     pb:   the progress bar, so that messages can be printed properly
pub fn write_vtk(loc: &str, 
                 name: &str, 
                 iter: usize, 
                 df: &CartesianDataFrame2D,
                 data: &std::collections::HashMap<String, OutputCallBack>,
                 pb: &ProgressBar){

    // generating file name information, and print message letting user know we are writing a vtk file
    let file_name = format!("{}_{:0>8}.vti", name, iter);
    pb.println(format!("    Writing vtk file '{}'", file_name));
    let file = format!("{}_{}", loc, file_name);


    // generate information that is constant across all the Attributes
    let mesh = &df.underlying_mesh;
    let extent = Extent::Ranges([0..=mesh.n[0] as i32, 0..=mesh.n[1] as i32 as i32, 0..=1]);
    let version = Version::new((0,1));
    let byte_order = ByteOrder::LittleEndian;

    // loop through all items in the data hash map, and generate the cell data Attributes
    let mut cell_data: Vec<Attribute> = Vec::with_capacity(data.len());
    for (name, get_data) in data.iter() {
        let vec_data = get_data(df);

        match vec_data {
            Ok(dat) => {
                let iobuffer = IOBuffer::F64(dat.data);
                
                let dab = DataArrayBase{
                    name: String::from(name),
                    elem: dat.data_type,
                    data: iobuffer,
                };

                cell_data.push(Attribute::DataArray(DataArray::from(dab)));

            }
            Err(msg) => {
                pb.println(format!("        {} Couldn't convert {} to vtk format, so not including it: '{}'"
                                , "Warning:".yellow(), name, msg))
            }
        }
        
    }

    // create the attribute/s
    let vtk_attributes = Attributes{
        cell: cell_data,
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