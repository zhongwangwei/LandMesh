# LandMesh

This repository contains the grid generation tool for the Common Land Model version 2024 (CoLM2024). It allows for the creation of unstructured meshes with adaptive refinement based on various surface characteristics.

## Key Features

- Generates initial icosahedral or hexagonal grids
- Performs adaptive mesh refinement based on configurable thresholds
- Supports refinement criteria like land type heterogeneity, topography, LAI, soil properties, etc.
- Outputs grid files compatible with CoLM2024


## Usage

1. Edit `mkgrd.nml` to configure grid generation options
2. Compile the Fortran code 
3. Run the executable
4. Output files will be generated in the specified directory

## Dependencies

- NetCDF
- A Fortran compiler (tested with gfortran)

## Contributing

Contributions to improve the grid generation tool are welcome. Please submit issues and pull requests on GitHub.

## License

This project is licensed under the GNU General Public License v2.0. See the LICENSE file for details.

## Citation

If you use this tool in your research, please cite:

[Fan, H., Xu, Q., Bai, F., Wei, Z., Zhang, Y., Lu, X., Wei, N., Zhang, S., Yuan, H., Liu, S. and Li, X., 2024. An unstructured mesh generation tool for efficient high‐resolution representation of spatial heterogeneity in land surface models. Geophysical Research Letters, 51(6), p.e2023GL107059.]

## Contact

For questions or support, please contact zhongwang wei (zhongwang007@gmail.com).

## Author
Hanwen Fan and Zhongwang Wei @ SYSU


## REVISION HISTORY

! 2024.07.19 Zhongwang Wei

! 2023.10.28  Hanwen Fan and Zhongwang Wei @ SYSU

! 2023.02.21  Zhongwang Wei @ SYSU

! 2021.12.02  Zhongwang Wei @ SYSU 

! 2020.10.01  Zhongwang Wei @ SYSU

###
## Directory Structure

//////////////////////////////////////////////////////////////////////////////////

[makegrid subfolder]    [Description]
source_data     Stores input data for makegrid
refine          Code related to constructing unstructured grids
mkgrd.nml       Execution file for mkgrd folder

/////////////////////////////////////////////////////////////////////////////////

[How to execute mkgrd.nml]
./mkgrd/mkgrd.x mkgrd.nml

/////////////////////////////////////////////////////////////////////////////////

[mkgrd folder]
- `mkgrd.F90`: Main program 
- `icosahedron.F90`: Generates initial icosahedral grid
- `MOD_refine_sjx.F90`: Module for refining triangular grids
- `MOD_refine_lbx.F90`: Module for refining hexagonal grids
- `MOD_refine_lbx_step2.F90`: 2nd, 3rd... refinement of polygonal unstructured grid
- `MOD_Get_Contain_Patch.F90`: Calculates grid overlaps and patch types
- `precision.F90   `: Variable format abbreviation settings
- `define.h `:  Macro definitions
- `GetContain.F90 `:  Calculate containment relationship between unstructured and structured grids
- `GetThreshold.F90 `:  Calculate threshold file for unstructured grid refinement
- `Get_Contain_PatchId.F90 `:  Calculate containment relationship and patchtype needed for MPI

- Other supporting modules and utilities


///////////////////////////////////////////////////////////////////////////////

[output folder]
tmp                              Intermediate output files for polygonal grid refinement
threshold                        Threshold files
sfcfiles                         Initial unrefined grid files
result                           Files from final iteration of polygonal grid refinement
patchtypes                       Files needed for MPI execution
contain                          Files for containment relationship between unstructured and structured grids

//////////////////////////////////////////////////////////////////////////////

[Execution logic]: Please modify the namelist 

1. If no refinement is needed, the program will directly run Get_Contain_PatchId.F90 to calculate the containment relationship and patchtype file for MPI

2. If refinement is needed,  the program will run GetContain.F90 to calculate the containment relationship

3. Then will run GetThreshold.F90 to calculate triangular grid thresholds (since refinement is based on triangular grids)

4. Run refine_lbx.F90 for one level of refinement (or directly execute refine_sjx.F90 for triangular grids)

5. If only one level of refinement is needed, run Get_Contain_PatchId.F90 to calculate the containment relationship and patchtype file for MPI

6. If multiple levels of refinement are needed, run refine_lbx_step2.F90 to create multi-level refined grid files

6. Run Get_Contain_PatchId.F90 to calculate containment relationship and patchtype file for MPI

7. Main final output files: unstructured grid file, containment relationship file, patchtype file

///////////////////////////////////////////////////////////////////////////////
