## Fast mDLAG (multi-population DLAG)

[![][license-img]][license-url]

[license-img]: https://img.shields.io/github/license/mashape/apistatus.svg
[license-url]: https://github.com/egokcen/fast-mdlag/blob/main/LICENSE

Fast mDLAG contains methods for dramatically accelerating the fitting of
Delayed Latents Across (multiple) Groups (mDLAG) models. mDLAG is a dimensionality
reduction framework for characterizing the multi-dimensional, concurrent flow of signals
across multiple groups of time series (for example, the recorded activity of multiple
neuronal populations). As preliminaries to this work, see 
[this paper](https://nips.cc/virtual/2023/poster/70171) and its
accompanying [repository](https://github.com/egokcen/mDLAG).

### Table of Contents

- [Citing this work](#citing-this-work)
- [Overview](#overview)
- [System requirements](#system-requirements)
- [Installation guide](#installation-guide)
- [Instructions for use](#instructions-for-use)
- [Contact](#contact)
- [License](#license)

## Citing this work

This repository accompanies the following paper:
- Gokcen, E., Jasper, A. I., Kohn, A., Machens, C. K., & Yu, B. M.
  Fast multi-group Gaussian process factor models. _Neural Computation_ (2025).

Please read it carefully before using the code, as it describes all of the
terminology and usage modes. For proper attribution, please cite the above reference
and [this codepack](CITATION.cff) if using any portion of this code for your own 
purposes.

## Overview

This codepack includes implementations of the primary contributions of the paper: a
method to fit mDLAG via inducing variables (mDLAG-inducing), and a method to fit mDLAG
via the frequency domain (mDLAG-frequency). To keep the codepack self-contained, the
baseline approach to fitting mDLAG (mDLAG-time) is also included (original repository
[here](https://github.com/egokcen/mDLAG)). See also 
[this repository](https://github.com/egokcen/DLAG) for a frequency
domain approach to fitting [DLAG](https://www.nature.com/articles/s43588-022-00282-5.epdf?sharing_token=hyFsoNFZfDyh-EVI76hoddRgN0jAjWel9jnR3ZoTv0Nm1Wps5WZ1Nog-dORHPLUG97YnGS0JZBkvhpO7c5pblBICIHRXMKZ04hmmro2Tn12HIbx2e2LrperSJc6bwzqptnPIaVOrqvl8DcloXzDaOBhLlAqzUvwM4uMyl96KvTE%3D).

The codepack directory structure is as follows:
- `mDLAG`

  Self-contained implementations of mDLAG-time and mDLAG-frequency.
  `mDLAG/demo/demo_mdlag.m` is a script that provides a generic demonstration of
  mDLAG-time. `mDLAG/demo/demo_mdlag_freq.m` provides a generic demonstration of
  mDLAG-frequency and comparison to mDLAG-time.

- `smDLAG`

  A self-contained implementation of mDLAG-inducing. `smDLAG/demo/demo_smdlag.m` is a
  script that provides a generic demonstration of mDLAG-inducing.
  
  There will be some conflicting function definitions across `mDLAG` and `smDLAG`. We 
  take care of these namespace issues throughout the codepack. Just beware of these
  potential conflicts if you plan to use both `mDLAG` and `smDLAG` in your own project.
  
  Note that mDLAG-inducing is referred to as smDLAG throughout the codepack, for
  "sparse mDLAG". In the Gaussian process (GP) literature, inducing variables fall under
  the umbrella of sparse GP methods. In the paper, we chose to avoid "sparse mDLAG" in
  favor of mDLAG-inducing to distinguish between different flavors of sparsity, for
  example the sparsity imposed on the loading matrices by automatic relevance
  determination.

- `simulation`

  The code and results to reproduce analyses and figures related to simulated
  experiments. The original datasets are themselves too large to store in this
  repository. However, the scripts and random seeds needed to reproduce those datasets
  locally on your own machine are included. The relevant subdirectories and the figures
  to which they relate are as follows:
    - `basic_demo`: Figs. 1 and 2
    - `scaling_seqlength`: Figs. 3 and 7, Supplementary Figs. 1, 2, and 8
    - `scaling_numgroups`: Fig. 4, Supplementary Fig. 7
    - `circulant_demo`: Figs. 1 and 6, Supplementary Fig. 9
    - `bias_modelselection`: Fig. 8
    - `scaling_samplingrate`: Supplementary Fig. 7
    - `scaling_neurons`: Supplementary Fig. 7

- `npx`

  Part of the data, code, and results to reproduce analyses and figures related to the
  Neuropixels recordings. The related
  figures are as follows:
    - Fig. 5, Supplementary Figs. 3–6 and 10

## System requirements

This codepack was written in Matlab (The MathWorks, Inc.), and must be run in 
Matlab.

It has been tested on Matlab 2023b and 2024a, on Linux 
(Pop!_OS 22.04 LTS, Red Hat Enterprise Linux release 7.9) operating systems.

Note that all runtimes provided below are estimates based on these tests. Your
runtimes may vary.

Some functions rely on the C/MEX Matlab interface for speedup, whereby C code 
can be compiled and used from within Matlab. A native C compiler is necessary
to take advantage of this functionality. Windows users may require extra 
installation of, for example, Microsoft Visual C++ or MinGW. The code is written
to default to a native Matlab code implementation if mex files cannot be 
properly executed, and will still work correctly.

## Installation guide

Simply download and extract the
[latest release](https://github.com/egokcen/fast-mdlag/releases) of this codepack to
your desired local directory.

Install [Matlab](https://www.mathworks.com/products/matlab.html).

For C/MEX Compilation:
1. Enter `mDLAG` directory. Run `startup.m`.
2. Enter `smDLAG` directory. Run `startup.m`.

Assuming Matlab is already installed on your machine, setup should not take 
more than a few minutes.

## Instructions for use

### mDLAG-time generic demo

1. Enter `mDLAG` directory.
2. Run `startup.m` to add all necessary dependencies to the Matlab path.
3. Open `mDLAG/demo/demo_mdlag.m` (remain in `mDLAG` directory).
4. Run `demo_mdlag.m` cell-by-cell. The script contains directions and 
   descriptions of relevant user-defined parameters.
   
With the current data and settings, `demo_mdlag.m` should take ~5 min to run.

### mDLAG-frequency generic demo

1. Enter `mDLAG` directory.
2. Run `startup.m` to add all necessary dependencies to the Matlab path.
3. Open `mDLAG/demo/demo_mdlag_freq.m` (remain in `mDLAG` directory).
4. Run `demo_mdlag_freq.m` cell-by-cell. The script contains directions and 
   descriptions of relevant user-defined parameters.
   
With the current data and settings, `demo_mdlag_freq.m` should take seconds to run.

### mDLAG-inducing generic demo

1. Enter `smDLAG` directory.
2. Run `startup.m` to add all necessary dependencies to the Matlab path.
3. Open `smDLAG/demo/demo_smdlag.m` (remain in `smDLAG` directory).
4. Run `demo_smdlag.m` cell-by-cell. The script contains directions and 
   descriptions of relevant user-defined parameters.
   
With the current data and settings, `demo_smdlag.m` should take ~10 min to run.

### Reproduce the elements of Fig. 1B and 1D

1. Enter `simulation/basic_demo` directory.
2. Open and run `smdlag_matrix_notebook.m`.

### Reproduce the elements of Fig. 2

1. Enter `simulation/basic_demo` directory.
2. Run `startup.m` to set up common directories and other constants.
3. To see pre-saved results, open `results_summary_basic_demo.m` and run it
   cell-by-cell.
4. To fit each method from scratch:
   1. mDLAG-time: Run `fit_mdlag_basic_demo_time.m`. The script should take ~1.2 hours
      to run.
   2. mDLAG-inducing: Run `fit_smdlag_basic_demo.m`. The script should take ~10 min to
      run.
   3. mDLAG-frequency: Run `fit_mdlag_basic_demo_freq.m`. The script should take ~30
      sec to run.

### Reproduce the elements of Fig. 3

1. Enter the `simulation/scaling_seqlength` directory.
2. Run `startup.m` to set up common directories and other constants.
3. (*Only needed for steps 4.ii and 6*) Run `generate_data_scaling_seqlength.m` to
   generate the simulated datasets (with T = 500 time points per trial). The script
   should take seconds to run.
4. To reproduce Fig. 3A, open `performance_evaluation/prediction_seqlength.m`.
   Remain in the `simulation/scaling_seqlength` directory.
   1. To see pre-saved performance calculations, skip to Line 80 and run the script
      cell-by-cell.
   2. To re-calculate performance from existing model fits, run the script in full.
      The script should take ~10 min to run.
5. To reproduce Fig. 3B-D, open `performance_evaluation/runtime_seqlength.m`. Remain
   in the `simulation/scaling_seqlength` directory. Run the script.
6. To perform the base experiments (with up to T = 500 time points per trial) from
   scratch:
   1. mDLAG-time: Run `fit_mdlag_scaling_seqlength_time.m`. With all runs parallelized,
      the script should take ~10 hours to run.
   2. mDLAG-inducing: Run `fit_smdlag_scaling_seqlength.m`. With all runs parallelized,
      the script should take ~30 min to run.
   3. mDLAG-frequency: Run `fit_mdlag_scaling_seqlength_freq.m`. With all runs
      parallelized, the script should take ~5 min to run.
7. To perform the extended experiments (with up to T = 5,000 time points per trial) from
   scratch:
   1. Run `generate_data_longtrials.m` to generate the simulated datasets 
      (with T = 5,000 time points per trial). The script should take seconds to run.
   2. mDLAG-frequency: Run `fit_mdlag_longtrials_freq.m`. With all runs parallelized,
      the script should take ~10 min to run.

### Reproduce the elements of Fig. 4

1. Enter the `simulation/scaling_numgroups` directory.
2. Run `startup.m` to set up common directories and other constants.
3. (*Only needed for steps 4.ii and 6*) Run `generate_data_scaling_numgroups.m` to
   generate the simulated datasets. The script should take seconds to run.
4. To reproduce Fig. 4A, open `performance_evaluation/prediction_numgroups.m`.
   Remain in the `simulation/scaling_numgroups` directory.
   1. To see pre-saved performance calculations, skip to Line 80 and run the script
      cell-by-cell.
   2. To re-calculate performance from existing model fits, run the script in full.
      The script should take ~2 min to run.
5. To reproduce Fig. 4B-D, open `performance_evaluation/runtime_numgroups.m`. Remain
   in the `simulation/scaling_numgroups` directory. Run the script.
6. To perform these experiments from scratch:
   1. mDLAG-time: Run `fit_mdlag_scaling_numgroups_time.m`. With all runs parallelized,
      the script should take ~41 hours to run.
   2. mDLAG-inducing: Run `fit_smdlag_scaling_numgroups.m`. With all runs parallelized,
      the script should take ~2 min to run.
   3. mDLAG-frequency: Run `fit_mdlag_scaling_numgroups_freq.m`. With all runs
      parallelized, the script should take ~1 min to run.

### Reproduce the elements of Fig. 5 and Supplementary Figs. 3–6 and 10

1. Enter `npx` directory.
2. Run `startup.m` to set up common directories and other constants.
3. To reproduce Fig. 5B and Supplementary Fig. 4A, open and run `pred_summary.m`.
4. To reproduce Fig. 5C-E and Supplementary Fig. 4B, open `runtime_summary.m`:
   1. Fig. 5C-E: Set `NUM_INDUCE = 20`. Run the script.
   2. Supplementary Fig. 4B: Set `NUM_INDUCE = 32`. Run the script.
5. To reproduce Supplementary Fig. 3, open and run `gp_timescale_summary.m`.
6. To reproduce Supplementary Fig. 5, open and run `runtime_summary_seqlength.m`.
7. To reproduce Supplementary Fig. 6, open and run `runtime_summary_numgroups.m`.
8. To reproduce Supplementary Fig. 10, open and run `finetune_example.m`.

### Reproduce the elements of Fig. 6, Fig. 1F, and Supplementary Fig. 9

1. Enter the `simulation/circulant_demo` directory.
2. Run `startup.m` to set up common directories and other constants.
3. Open and run `circulant_bias_notebook.m` cell-by-cell:
   - Sections 1b and 3a reproduce Fig. 1F
   - Sections 1a and 1c reproduce Fig. 6A and 6B
   - Section 2a reproduces Fig. 6C
   - Sections 3, 3c, and 3d reproduce Supplementary Fig. 9

### Reproduce the elements of Fig. 7 and Supplementary Fig. 8

1. Enter the `simulation/scaling_seqlength` directory.
2. Run `startup.m` to set up common directories and other constants.
3. (*Only needed for step 7, if not already done*) Run
   `generate_data_scaling_seqlength.m` to generate the simulated datasets. The script
   should take seconds to run.
4. Open `performance_evaluation/bias_seqlength.m`. Remain in the
   `simulation/scaling_seqlength` directory. Note three constants toward the top of the
   script: `FIX_TYPE`, `PLOT_SPARSE`, and `PLOT_TAPER`.
5. To reproduce Fig. 7, set `FIX_TYPE = 0`, `PLOT_SPARSE = 0`, and `PLOT_TAPER = 1`.
   Run the script.
6. To reproduce Supplementary Fig. 8:
   1. Fix timescales: Set `FIX_TYPE = 1`. Run the script.
   2. Fix delays: Set `FIX_TYPE = 2`. Run the script.
7. To perform these experiments from scratch:
   1. mDLAG-frequency with tapering: Run `fit_mdlag_scaling_seqlength_freq_taper.m`.
      With all runs parallelized, the script should take ~5 min to run.
   2. mDLAG-frequency, fixed timescales: Run
      `fit_mdlag_scaling_seqlength_freq_fixtau.m`. With all runs parallelized, the
      script should take ~5 min to run.
   3. mDLAG-frequency, fixed delays: Run
      `fit_mdlag_scaling_seqlength_freq_fixD.m`. With all runs parallelized, the
      script should take ~5 min to run.

### Reproduce the elements of Fig. 8

1. Enter the `simulation/bias_modelselection` directory.
2. Run `startup.m` to set up common directories and other constants.
3. (*Only needed for step 6*) Run `generate_data_modelselection.m` to
   generate the simulated datasets. The script should take seconds to run.
4. To reproduce Fig. 8A and 8B, open `performance_evaluation/bias_modelselection.m`.
   Remain in the `simulation/bias_modelselection` directory. Run the script.
5. To reproduce Fig. 8C, open `performance_evaluation/bias_modelselection_taper.m`.
   Remain in the `simulation/bias_modelselection` directory. Run the script.
6. To perform these experiments from scratch:
   1. mDLAG-time: Run `fit_mdlag_modelselection_time.m`. With all runs parallelized, the
      script should take ~7 hours to run.
   2. mDLAG-inducing: Run `fit_smdlag_modelselection.m`. With all runs parallelized, the
      script should take ~5 hours to run.
   3. mDLAG-frequency: Run `fit_mdlag_modelselection_freq.m`. With all runs
      parallelized, the script should take ~25 min to run.
   4. mDLAG-frequency with tapering: Run `fit_mdlag_modelselection_freq_taper.m`.
      With all runs parallelized, the script should take ~25 min to run.

### Reproduce the elements of Supplementary Fig. 1

1. Enter the `simulation/scaling_seqlength` directory.
2. Run `startup.m` to set up common directories and other constants.
3. Open `performance_evaluation/single_dataset_seqlength.m`. Remain in the
   `simulation/scaling_seqlength` directory. Run the script.

### Reproduce the elements of Supplementary Fig. 2

1. Enter the `simulation/scaling_seqlength` directory.
2. Run `startup.m` to set up common directories and other constants.
3. (*Only needed for step 5, if not already done*) Run
   `generate_data_scaling_seqlength.m` to generate the simulated datasets. The script
   should take seconds to run.
4. Open `performance_evaluation/bias_inducingpoints.m`. Remain in the
   `simulation/scaling_seqlength` directory. Run the script.
5. To perform these experiments from scratch, run `fit_smdlag_bias_inducingpoints.m`.
   With all runs parallelized, the script should take ~30 min to run.

### Reproduce the elements of Supplementary Fig. 7A

1. Enter the `simulation/scaling_numgroups` directory.
2. Run `startup.m` to set up common directories and other constants.
3. Run `generate_data_scaling_numgroups.m` to generate the simulated datasets. The
   script should take seconds to run.
4. To reproduce Supplementary Fig. 7A, open `performance_evaluation/bias_numgroups.m`.
   Remain in the `simulation/scaling_numgroups` directory. Run the script.
5. To perform the tapering experiments from scratch, run
   `fit_mdlag_scaling_numgroups_freq_taper.m`. With all runs parallelized, the script
   should take ~1 min to run.

### Reproduce the elements of Supplementary Fig. 7B

1. Enter the `simulation/bias_samplingrate` directory.
2. Run `startup.m` to set up common directories and other constants.
3. (*Only needed for step 5*) Run `generate_data_bias_samplingrate.m` to generate the
   simulated datasets. The script should take seconds to run.
4. To reproduce Supplementary Fig. 7B, open 
   `performance_evaluation/bias_samplingrate.m`. Remain in the 
   `simulation/bias_samplingrate` directory. Run the script.
5. To perform these experiments from scratch:
   1. mDLAG-time: Run `fit_mdlag_bias_samplingrate_time.m`. With all runs parallelized,
      the script should take ~10 min to run.
   2. mDLAG-frequency: Run `fit_mdlag_bias_samplingrate_freq.m`. With all runs
      parallelized, the script should take ~1 min to run.
   3. mDLAG-frequency with tapering: Run `fit_mdlag_bias_samplingrate_freq_taper.m`.
      With all runs parallelized, the script should take ~1 min to run.

### Reproduce the elements of Supplementary Fig. 7C

1. Enter the `simulation/scaling_neurons` directory.
2. Run `startup.m` to set up common directories and other constants.
3. (*Only needed for step 5*) Run `generate_data_scaling_neurons.m` to generate the
   simulated datasets. The script should take seconds to run.
4. To reproduce Supplementary Fig. 7C, open `performance_evaluation/bias_neurons.m`.
   Remain in the `simulation/scaling_neurons` directory. Run the script.
5. To perform these experiments from scratch:
   1. mDLAG-time: Run `fit_mdlag_scaling_neurons_time.m`. With all runs parallelized,
      the script should take ~4 hours to run.
   2. mDLAG-inducing: Run `fit_smdlag_scaling_neurons.m`. With all runs parallelized,
      the script should take ~1 hour to run.
   3. mDLAG-frequency: Run `fit_mdlag_scaling_neurons_freq.m`. With all runs
      parallelized, the script should take ~1.5 hours to run.
   4. mDLAG-frequency with tapering: Run `fit_mdlag_scaling_neurons_freq_taper.m`.
      With all runs parallelized, the script should take ~1.5 hours to run.
   
## Contact
For questions, please contact Evren Gokcen at egokcen@cmu.edu.

## License
[MIT](LICENSE)