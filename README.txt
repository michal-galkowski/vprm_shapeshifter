VPRM Shapeshifter v1.2
Michal Galkowski
MPI-BGC Jena, October 2018
michal.galkowski@bgc-jena.mpg.de

Tool to preprocess data of MODIS indices for VPRM in WRF-GHG v3.9.1.1.
Designed to run on Mistral DKRZ cluster, but should be easily ported
to any linux machine. For help with running on other machines, please
contact the author.

Requirements:
R                - version > 3.0 (only tested on 3.3.3 and 3.5.0)
R ncdf4 library  - available on R startup on Mistral
git              - available on Mistral by default

Instructions of usage @Mistral

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

5. Run the script from:
a) command line:
$ Rscript run_vprm_shapeshifter.R
b) inside R:
> source("run_vprm_shapeshifter.R")

6. Output will be written into the output dir in the form of daily WRF inputs
named similar to:
vprm_input_d01_2000_01_01.nc
These can be linked to the WRF run directory and included in the namelist.input file.
Run with chem_opt = 16 or 17 (tested with 17)

Please inform me about any bugs in this software
