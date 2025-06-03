# EarthMesh

EarthMesh is a grid generation tool for land surface, ocean, and atmospheric models. This is an update from v1 which was primarily for land surface models. The `mesh_type` variable in `v2/mkgrd.F90` (which can be `landmesh`, `oceanmesh`, or `earthmesh`) determines the type of grid generated.

This repository contains the grid generation tool for the Common Land Model version 2024 (CoLM2024). It allows for creating unstructured meshes with adaptive refinement based on various surface characteristics.

## Key Features

- Generates initial icosahedral or hexagonal grids
- Performs adaptive mesh refinement based on configurable thresholds
- Supports refinement criteria like land type heterogeneity, topography, LAI, soil properties, etc. (for `landmesh`)
- Supports generation of ocean and global earth meshes.
- Outputs grid files compatible with CoLM2024 and other models

## Dependencies

- NetCDF
- A Fortran compiler (tested with gfortran)

## Author
Rui Zhang (V2), Hanwen Fan (V1) and Zhongwang Wei @ SYSU

## Directory Structure for V2

The primary code for v2 is in the `v2/` directory. Key files/components within `v2/` include:

- `mkgrd.F90`: The main program for grid generation.
- `*.nml` (e.g., `landmeshv7.nml`, `nxp36.nml`): Namelist files for configuring grid generation options.
- `Makefile`: Used for compiling the v2 code.
- `run.slurm`, `make.slurm`: Example SLURM submission scripts for job scheduling.

The `postprocessing/` directory contains scripts for visualizing or analyzing generated grids.

## Usage (V2)

### Compilation
1. Navigate to the `v2/` directory: `cd v2`
2. Run `make` to compile the Fortran code: `make`
   (Ensure a Fortran compiler and NetCDF libraries are installed and configured in your environment).

### Configuration
1. Edit the relevant `*.nml` file in the `v2/` directory (e.g., `landmeshv7.nml`) to set parameters for grid generation.
2. Key parameters include `mesh_type` which can be set to `landmesh`, `oceanmesh`, or `earthmesh` to determine the model domain.

### Execution
1. From within the `v2/` directory, run the executable with a namelist file as an argument:
   `./mkgrd.x <your_namelist_file.nml>`
   For example: `./mkgrd.x landmeshv7.nml`
2. Alternatively, if using a SLURM environment, you can submit the job using the provided scripts:
   `sbatch run.slurm` (after configuring `run.slurm` appropriately).

## Output Folders (V2)

Output directories are created within the `file_dir` specified in the namelist file (typically `base_dir/expnme/`). These include:

- `contain/`: Stores containment relationship files, mapping grid cells between different resolutions or unstructured and structured grids.
- `gridfile/`: Stores the main grid definition files (e.g., cell vertices, connectivity).
- `patchtype/`: Stores patch type files, which are often needed for domain decomposition in MPI-parallel models.
- `result/`: Stores final mesh files from the last iteration of refinement. These are typically the primary output grids.
- `tmpfile/`: Stores intermediate output files generated during polygonal grid refinement steps.
- `threshold/`: Stores threshold files used for adaptive refinement if refinement is enabled.

## Execution Logic (V2)

The core execution logic for v2, including the handling of different `mesh_type` values (`landmesh`, `oceanmesh`, `earthmesh`) and the adaptive refinement processes, is implemented in the `v2/mkgrd.F90` file. The main program block within this file orchestrates the grid generation steps.

## Contributing

Contributions to improve the grid generation tool are welcome. Please submit issues and pull requests on GitHub.

## License

This project is licensed under the GNU General Public License v2.0. See the LICENSE file for details.

## Citation

If you use this tool in your research, please cite:

[Fan, H., Xu, Q., Bai, F., Wei, Z., Zhang, Y., Lu, X., Wei, N., Zhang, S., Yuan, H., Liu, S. and Li, X., 2024. An unstructured mesh generation tool for efficient high‐resolution representation of spatial heterogeneity in land surface models. Geophysical Research Letters, 51(6), p.e2023GL107059.]

## Contact

For questions or support, please contact Zhongwang Wei (zhongwang007@gmail.com).

## REVISION HISTORY
! 2025.01.09 The Initial version of V2 is being developed by Rui Zhang
! 2025.01.09 V1 has been archived by Zhongwang Wei
! 2024.07.19 Zhongwang Wei
! 2023.10.28  Hanwen Fan and Zhongwang Wei @ SYSU
! 2023.02.21  Zhongwang Wei @ SYSU
! 2021.12.02  Zhongwang Wei @ SYSU
! 2020.10.01  Zhongwang Wei @ SYSU
