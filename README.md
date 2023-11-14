# Benchmarks for Membrane Protein Energy Functions
##### This is under construction ###################
This is a set of scientific benchmark tests for evaluating membrane protein modeling energy functions. The test probe an energy function's ability to capture membrane protein orientation, stability, sequence, and structure. The methods are described in detail in the citation below. 

 - Alford RF, Samanta R & Gray JJ; "Diverse scientific benchmarks for implicit membrane energy functions", JCTC, 2021
 - Samanta R & Gray JJ; "Implicit model to capture electrostatic features of membrane environment", bioRxiv, 2023

## Manifest

 - `tests/` - Scientific benchmark scripts for generation, post-processing, and analysis
 - `targets/` - Experimental and model datasets for testing
 - `example/` - Example for analyzing output data from franklin2019
 - `LICENSE` - MIT license for benchmark code
 - `config.txt` - Path and platform information for Rosetta

## Prerequisites

#### System software

The test framework requires Python version 3.7+, PyRosetta version 3.7+ . In addition, the data generation stage takes advantange of high-performance computing resources. the default setup is configured for a slurm cluster. We also support conda clusters. For other setups, please email the author for assistance. 

#### Molecular modeling software

The tests use both the command-line and python interfaces to the Rosetta macromolecular modeling suite. Rosetta is available to academic users for free and to comercial users for a fee. 

To get Rosetta, obtain a license and download the package at <https://www.rosettacommons.org/software/license-and-download>. To compile the code, navigate to the `Rosetta/main/source/` directory and run the command below. 

```
./scons.py -jX bin mode=release 
```

Here, "X" is the number of processors to use during compilation. For compilation on a laptop, the recommended number of processors is 1. If you are working on a larger workstation or high performance computing cluster, we recommend scaling up to 8-24 processors. More information can be found in the [Rosetta Build Documentation](https://www.rosettacommons.org/docs/wiki/build_documentation/Build-Documentation#setting-up-rosetta-3_basic-setup). 

To get PyRosetta, follow instructions at: 

```
https://www.pyrosetta.org/downloads
```

## Setup

To begin, modify the `config.txt` file to include information about your Rosetta installation, the location of this repository, and the platform. Each variable is described in the example below. 

```
benchmark = /path/to/benchmark  	# Path to this repo
rosettadir = /path/to/rosetta   	# Path to Rosetta bin
platform = linux 			# can be linux or mac
buildenv = release			# can be release or debug
compiler = gcc				# can be gcc or clang
```

## Documentation

The implicit membrane energy function benchmarks involves three steps: (1) data generation, (2) post-processing, and (3) analysis. We will walk through each step below. 

#### Step 1: Generate benchmark data

The first step performs all initial PyRosetta and Rosetta modeling calculations via a computing cluster. The generation script takes about 30min to run. Afterward, job completion requires roughly 1K CPU hours for all 12 tests. To run the generation step, use the command line below. 

	To test franklin2019: './generate_test_data.py --energy_fxn franklin2019 --which_tests all'
	To test franklin2023: './generate_test_data.py --energy_fxn franklin2023 --which_tests all'


The `--energy_fxn` flag sets the energy function to test, referred to by the name of the weights file in the Rosetta database (must be present in both Rosetta & PyRosetta). The `--which_tests` flag indicates which tests to be run. This can be `all` or a comma-separated list of the tests as given below. 

| Test                        | #  | Description 													   |
|-----------------------------|----|-------------------------------------------------------------------|
| tm-peptide-tilt-angle       | 1  | Tilt angle of single-span transmembrane peptides           	   |
| adsorbed-peptide-tilt-angle | 2  | Tilt angle of surface-adsorbed peptides    				   |
| protein-tilt-angle          | 3  | Tilt angle of multi-pass membrane proteins 					   |
| hydrophobic-length          | 4  | Hydrophobic thickness of multi-pass membrane proteins             |
| ddG-of-insertion            | 5  | Energetic cost of transfering a peptide from water to bilayer     |
| ddG-of-pH-insertion         | 6  | Energetic cost of pH-dependent water-to-biolayer peptide transfer |
| ddG-of-mutation             | 7  | Energetic cost of single-point mutations in the membrane          |
| sequence-recovery           | 8  | Recovery of sequence features after full fixed-backbone redesign  |
| sc-distribution             | 9  | Depth-dependent distribution of side chains relative to native    |
| decoy-discrimination        | 10 | Native structure discrimination          						   |
| helix-kink                  | 11 | Helix kink angle prediction            						   |
| protein-protein-docking     | 12 | Membrane protein-protein docking           					   |

The outputs are then organized in a `data/` directory created by the script. The first subdirectory is the name of the energy function in use (e.g., franklin2019/franklin2023). Then, this folder contains 12 subdirectories for the data output from each test. 

#### For test 1-3 #######
We divide the run into 6 small executions which calculates the enrgy landscape over a small distance. However we then run one more step to recompile the energy landscape into a single file with the following code.
Test1: ./combiningfiles_peptide_tilt_angle.py --energy_fxn franklin2023 --which_tests tm-peptide-tilt-angle
Test2: ./combiningfiles_peptide_tilt_angle.py --energy_fxn franklin2023 --which_tests adsorbed-peptide-tilt-angle
Test3: ./combiningfiles_protein_tilt_angle.py --energy_fxn franklin2023 --which_tests adsorbed-peptide-tilt-angle

#### Step 2: Post-Process benchmark data

The sequence recovery and structure prediction benchmark tests require an imtermediate post-processing step before final data analysis. To run the post-processing step, run the command line below. The flags are described in step #1. 

	./process_test_data.py --energy_fxn franklin2023 --which_tests all

#### Step 3: Analyze benchmark data 

The final step is to visualize and analyze the results of each benchmark tests. To do so, we provide a package of R scripts for generating the appropriate plots. While R can be run for the command line, we recommend downloading the R studio IDE from (https://rstudio.com/). You can run the example analysis script for `franklin2019` data from the `example/` directory through R studio or with the command line below. 

	Rscript analyze_f19_tests.R 

If you want to avoid using R package you can also keep using the python. 
The following code is particularly written for tests 1,2,5,7 and 9 (Results for Franklin2023 paper). 
	'./plot_benchmark_dataset --energy_fxn franklin2023 --which_tests ddG-of-mutation'

#### Determining the weights of score function #####
./calculate_franklin2023_weights.py