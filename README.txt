VPRM Shapeshifter
Michal Galkowski
MPI-BGC Jena, October 2018
michal.galkowski@bgc-jena.mpg.de

Tool to preprocess data of MODIS indices for VPRM in WRF-GHG v3.9.1.1.
Designed to run on Mistral DKRZ cluster, but should be easily ported
to any linux machine.



Changelog =========================================================
v1.5.1 - New option: no_colons, FALSE by default; setting the output filename format
v1.5 - Using netcdf4 compression as default now; some small bugfixes in description.
v1.4 - Added possibility of renaming Kaplan input files (additional function options); added .docx presentation about installation
v1.3 - Added netcdf parameters necessary for appropriate Land-Use setup in WRF-Simulations.



Requirements ======================================================
R                - version > 3.0 (tested on 3.3.3 and 3.5.0)
R ncdf4 library  - available on R startup on Mistral
git              - available on Mistral by default



Instructions of usage ==============================================

Valid on Mistral supercomputer (DKRZ), on other machines these
need to be adopted (user responsibility, limited assistance possible).

Note:
R ncdf4 library is available on R startup on Mistral, git is available from the command line by default

1. Prepare environment
$ module load r

2. Create vprm_shapeshifter directory in a chosen location
$ cd {YOUR TOOLS DIRECTORY}
$ mkdir vprm_shapeshifter
$ cd vprm_shapeshifter

3. Get the code from git repository by executing the following command:
$ git clone /work/mj0143/b301033/Projects/REPOSITORY/vprm_shapeshifter.rep/.

4. Change the directory paths in the file:
run_vprm_shapeshifter.R
Note: Kaplan input for CH4 biogenic emissions is optional

5. To run the script, either:
a) execute, from command line:
$ Rscript run_vprm_shapeshifter.R
b) OR, execute, from inside R:
> source("run_vprm_shapeshifter.R")

6. Output will be written into the output dir in the form of daily WRF inputs
named similar to:
vprm_input_d01_2000_01_01.nc
These can be linked to the WRF run directory and included in the namelist.input file.
Run with chem_opt = 16 or 17 (tested with 17)

Please inform me (michal.galkowski@bgc-jena.mpg.de) about any bugs in this software
