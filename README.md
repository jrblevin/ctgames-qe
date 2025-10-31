# Replication Package for "Identification and Estimation of Continuous Time Dynamic Discrete Choice Games"

[Identification and Estimation of Continuous Time Dynamic Discrete Choice Games][ctgames]  
[Jason R. Blevins][jblevins]

## Overview

The code in this replication package implements structural estimation and
Monte Carlo experiments for two continuous-time dynamic discrete choice
models analyzed in the paper.

The first model (`mc1p` directory) is a continuous-time single-agent renewal
model inspired by Rust (1987) that analyzes optimal bus engine replacement
decisions.  There are two main programs for this model: `rustct` carries out
the structural estimation of the model with real data while `mc1p` carries
out the Monte Carlo experiments.

The second model (`mcnp` directory) is a continuous-time quality ladder model
based on Ericson and Pakes (1995) that analyzes oligopoly dynamics with entry,
investment, and exit decisions.  A single program, `mcnp`, carries out the
Monte Carlo experiments.

## Data Availability

The Monte Carlo experiments in this paper do not involve analysis of external
data.  The data used are generated via simulation in the code itself.

The NFXP data used in the empirical example is from Rust (1987).  It is
included in this replication package and is freely available under the
MIT License.

### NFXP Data

The data for the Nested Fixed Point Algorithm (NFXP) (i.e., bus engine
replacement data) used in the single-agent renewal model (`mc1p`) consists of
odometer readings and dates of bus engine replacements for 162 buses in the
fleet of the Madison Metropolitan Bus Company that were in operation during
the period December 1974 to May 1985.  This data was originally collected and
used by Rust (1987).  The data files represent different bus models/vintages
and are provided in ASCII format.  Each file is named by bus model (e.g.,
`d309.asc` for Davidson model 309 buses, `g870.asc` for Grumman model 870
buses).

The original data and detailed documentation can be downloaded from John
Rust's website at <https://editorialexpress.com/jrust/nfxp.html> or from the
Zenodo repository at <https://doi.org/10.5281/zenodo.3374587>.  The data files
are also included in this replication package for convenience.

Data files: `mc1p/data/*.asc`

## Computational Requirements

### Software Requirements

  - Fortran 2008 or later compiler:
      - GNU Fortran (gfortran)
      - Intel Fortran Compiler
  - LAPACK and BLAS libraries (or Intel MKL as alternative)
  - OpenMP support (included with above compilers)
  - GNU Make
  - Bash shell
  - R (for statistical calculations in Table 2 generation)

The code has been tested with both GNU Fortran and Intel Fortran.
Substantial parts of the code use bash scripting, which may require Linux or
macOS.

### Controlled Randomness

For reproducibility, random number generator seeds for data generation are set
within the Monte Carlo experiment control programs: `mc1p/mc1p.f90` (line 130)
and `mcnp/mcnp.f90` (line 128).  However, differences in compilers, library
versions, compiler optimization settings, hardware, etc. may still lead to
slightly different results.

### Hardware Requirements

#### Summary

| Table   | Description                         | Hardware Used                    | Runtime |
|---------|-------------------------------------|----------------------------------|---------|
| Table 1 | Rust (1987) Sample Characteristics  | MacBook Pro (12 core M2 Max)     | 30 sec  |
| Table 2 | Estimates from Rust (1987) Data     | MacBook Pro (12 core M2 Max)     | 30 sec  |
| Table 3 | Single Agent Renewal Model          | MacBook Pro (12 core M2 Max)     | 10 min  |
| Table 4 | Quality Ladder Model Specifications | Mac Pro (28 core Intel Xeon W)   | 20 hr   |
| Table 5 | Quality Ladder Model (Partial)      | HPC Cluster (48 core Intel Xeon) | 10 hr   |
| Table 5 | Quality Ladder Model (Full)         | HPC Cluster (48 core Intel Xeon) | 5 days  |

Notes:

- Runtimes in the table above are for replicating complete tables, taking
  advantage of parallelism on the hardware shown.  Actual runtime may vary
  significantly by Fortran compiler used, hardware capabilities, and the
  level of multiprocessing used.

- The high performance computing (HPC) cluster nodes we used have 48 Intel
  Xeon Platinum 8260 2.40GHz CPU cores and 192 GB memory.

- Table 5 is very computationally intensive to fully replicate.  Please
  review the replication instructions below for details and recommendations.

    - For example, on the HPC cluster node used the final table row with N=8
      and Δ=1.0 required approximately 4 hours per Monte Carlo trial and a
      full replication requires 100 trials.  With 48 cores at our disposal,
      we executed trials concurrently using 6 parallel processes each with
      8 cores used for parallel processing within a trial.

#### Single Agent Model (Tables 1-3)

It is possible to run the single-agent renewal model estimation and Monte
Carlo experiments on a standard 2025 laptop or desktop computer.

- Tables 1 and 2: The `rustct` structural estimation code for the single agent
  model has been tested and easily completes in under one minute on a 2023
  MacBook Pro (12-core M2 Max processor).

- Table 3: The `mc1p` single agent Monte Carlo experiments and timing
  exercises have been tested on a 2019 Mac Pro workstation (2.5 GHz 28-Core
  Intel Xeon W processor), where the full set of experiments completed in
  around 5 minutes. It completes in around 10 minutes on the same 2023
  MacBook Pro mentioned above.

#### Quality Ladder Model (Tables 4-5)

The quality ladder model Monte Carlo experiments are very computationally
intensive and were carried out in parallel on an HPC cluster over several
days.  The `mcnp` quality ladder model Monte Carlo experiments involve
repeatedly solving for the equilibrium of a complex dynamic game with many
firms and computing the matrix exponential for a large matrix to calculate
the log likelihood function for each trial parameter value during the
optimization.  This becomes feasible due to the computational benefits of
continuous time games as well as the use of high performance, vectorized
Fortran code.

However, with more time it is also possible to run these experiments on a
high-end workstation.  For our largest experiment with 8 firms, a single
Monte Carlo replication takes about 5 hours on a 2019 Mac Pro.  Completing
all 100 replications reported in the paper would therefore take around 21
days.  On the other hand, for the smaller 4 firm model all 100 replications
can be completed in around 8 hours.

Due to the increased computational requirements for Table 5 as reported in
the paper, we also provide partial replication instructions and results below
which are based on the first 10 out of 100 Monte Carlo trials.

## Description of Code

### Single-Agent Renewal Model (`mc1p/`)

Main Programs:

  - `rustct.f90`: Structural estimation using bus engine replacement data.
  - `mc1p.f90`: Monte Carlo experiment control program.

Core Modules:

  - `rust_model.f90`: Single-agent continuous-time dynamic discrete choice
    model.
  - `rust_data.f90`: Data loading and processing for bus fleet data.
  - `dataset.f90`: Underlying data structure.

Libraries:

  - `lbfgsb/`: [L-BFGS-B][] optimization routines.
  - `expokit/`: [Expokit][] for dense Matrix exponential calculations.

Control Files:

  - `control/mc-<delta>-<nm>.ctl`: Monte Carlo experiment control files.

Automated replication scripts:

  - `table_1_2.sh`: Replicates Tables 1 and 2.
  - `table_3.sh`: Replicates Table 3.

### Quality Ladder Model (`mcnp/`)

Main Program:

  - `mcnp.f90`: Monte Carlo experiment control program.

Core Modules:

  - `model.f90`: Continuous-time dynamic oligopoly (quality ladder) model
    with entry, exit, and investment.
  - `encoding.f90`: Efficient state space encoding/decoding.
  - `dataset.f90`: Data structure for storing simulated observations.
  - `sparse.f90`: Sparse matrix operations for computational efficiency.

Libraries:

  - `lbfgsb/`: [L-BFGS-B][] optimization routines.

Control Files:

  - `control/example.ctl`: A small-scale example control file for testing on
    a standard laptop or desktop computer.
  - `control/mc-<nn>-<delta>.ctl`: Monte Carlo simulation control files for
    different numbers of firms `nn` and data sampling interval `delta`.
  - `control/time-<nn>.ctl`: Timing control files for different numbers of
    firms `nn`.

Automated replication scripts:

  - `table_4.sh`: Replicates Table 4.
  - `table_5.sh`: Replicates Table 5.
  - `run_parallel.sh`: Helper script for Table 5 that handles automatic
    distribution of Monte Carlo trials across parallel processes.
  - `table_5_save.sh`: Helper script for Table 5 that collects results
    from Monte Carlo replications and creates LaTeX table.

## Instructions to Replicators

### Setup

1. Ensure you have one of the supported Fortran compilers installed.
2. Ensure BLAS and LAPACK libraries are available (or Intel MKL).
3. Ensure that R is installed for computing statistics reported in tables.
4. Clone or extract this repository to a local directory.

### Quick Start

For users with a working Fortran development system looking to quickly test
the code:

```bash
# Build and estimate single-agent model (heterogeneous specification)
cd mc1p && make && ./rustct-heterogeneous

# Run a single-agent Monte Carlo experiment
./mc1p control/mc-1.00-200.ctl

# Run the simple example quality ladder model Monte Carlo experiment
cd ../mcnp && make
OMP_NUM_THREADS=1 ./mcnp control/example.ctl
```

### Installing Dependencies

On Debian/Ubuntu Linux:

```bash
sudo apt-get install gfortran liblapack-dev libblas-dev r-base
```

On CentOS/RHEL Linux:

```bash
sudo yum install gcc-gfortran lapack-devel blas-devel R
```

On macOS with Homebrew:

```bash
brew install gcc r
```

macOS provides BLAS and LAPACK through the Accelerate framework.

### Building the Code

For the single-agent model (`mc1p`):

```bash
cd mc1p
make clean
make
```

For the quality ladder model (`mcnp`):

```bash
cd mcnp
make clean
make
```

The Makefiles will use GNU Fortran by default.  To use the Intel Fortran
compiler, prefix the usual `make` command with `SYSTEM=intel`:

```bash
SYSTEM=intel make
```

**Note**: The table-specific replication scripts (`table_1_2.sh`,
`table_3.sh`, `table_4.sh`, `table_5.sh`) will automatically build the
required binaries if they are missing, so you may skip this step if you plan
to use those scripts directly.

### Single Agent Renewal Model Empirical Results (Tables 1 and 2)

**Estimated time**: 30 seconds on a modern laptop.

To automatically generate both Tables 1 and 2:

```bash
cd mc1p
./table_1_2.sh
```

This script:

- Builds binaries for three model variants using conditional compilation:
    - `rustct-abbe`: Fixed λ specification
    - `rustct-homogeneous`: Homogeneous λ specification
    - `rustct-heterogeneous`: Heterogeneous λ specification
- Runs each model variant using the Rust (1987) bus data and saves log files:
    - `./rustct-abbe > logs/rustct-abbe.log`
    - `./rustct-homogeneous > logs/rustct-homogeneous.log`
    - `./rustct-heterogeneous > logs/rustct-heterogeneous.log`
- Extracts appropriate results from the log files
- Computes the LR statistics and p-values.
    - This step requires `Rscript`. If it is not available, p-values will be
      omitted from the table.
- Generates `mc1p/results/table_1.tex` (sample characteristics)
- Generates `mc1p/results/table_2.tex` (parameter estimates and
  heterogeneity tests)

You may wish to manually compile the binaries (`rustct-abbe`,
`rustct-homogeneous`, and `rustct-heterogeneous`) used to produce these
results (for example, to change the Fortran compiler used). You can use
`make` to build the programs directly:

```bash
cd mc1p

# Default settings (GNU Fortran)
make rustct-abbe rustct-homogeneous rustct-heterogeneous

# Use Intel Fortran compiler
SYSTEM=intel make rustct-abbe rustct-homogeneous rustct-heterogeneous
```

### Single-Agent Monte Carlo Experiments (Table 3)

**Estimated time**: 10 minutes on a modern laptop.

To execute all the single-agent model Monte Carlo experiments reported in
Table 3, you can use the `table_3.sh` script:

```bash
cd mc1p
./table_3.sh
```

This script:

- Compiles the `mc1p` binary used for the Monte Carlo experiments.
- Runs each of the specifications reported in the paper and saves individual
  log files:
    - `./mc1p control/mc-<delta>-<nm>.ctl > logs/mc-<delta>-<nm>.log`
    - For each specification, by default all cores on your machine will be
      used to carry out Monte Carlo replications in parallel.
- Collects results from each specification and produces Table 3.
- Saves the complete LaTeX table to `mc1p/results/table_3.tex`
- Displays the results to the console.

Alternatively, you can manually compile the program and run individual
specifications by using the appropriate control file.  For example:

```bash
cd mc1p
make mc1p
./mc1p control/mc-1.00-3200.ctl
```

Individual control files are named as `mc-<delta>-<nm>.ctl` where `delta` is
the observation time interval (`0.00` for continuous time) and `nm` is the
number of markets simulated.

The `mc1p` program supports OpenMP parallelization, executing individual Monte
Carlo replications in parallel across threads.  To explicitly set the number
of threads, set the `OMP_NUM_THREADS` environment variable:

```bash
export OMP_NUM_THREADS=8   # Adjust number of cores for your system
```

### Quality Ladder Model: Table 4 (Timing Results)

**Estimated time**: 20 hours on a modern, multicore workstation

The "Obtain V" column reports wall clock times for solving the model using
value iteration and is not intended to be replicated exactly. The timing will
depend on your system characteristics.  To produce the complete table
automatically:

```bash
cd mcnp
./table_4.sh
```

This script runs timing experiments for N=2 through N=30 firms using the
corresponding control files (e.g., `./mcnp control/time-02.ctl`), saves log
files for each specification, extracts model size and timing information from
the output, and produces Table 4 in LaTeX format at
`mcnp/results/table_4.tex`.

_Note:_ Experiments with N≥26 firms require substantial RAM (~32 GB).
If experiments fail with an error message such as "Killed", this likely
indicates that the process ran out of memory.

### Quality Ladder Model: Table 5 (Monte Carlo Results)

Table 5 is the most computationally intensive to replicate.  We provide three
distinct replication approaches depending on your computational resources and
time constraints:

1. **Full Replication**: Complete reproduction with 100 Monte Carlo trials
   (~5 days in parallel on HPC cluster)
2. **Incremental Replication**: Run selected specifications individually
   (flexible timing)
3. **Partial Replication**: Quick validation with 10 trials instead of 100
   (~10 hours on HPC cluster node)

The table below lists the recommended number of threads/cores used per Monte
Carlo trial and the approximate time per trial.

| N   | Sampling   | Threads Per Trial | Time Per Trial |
|-----|------------|-------------------|----------------|
| 2   | Continuous | 2                 | 6 sec          |
|     | Δ = 1.0    | 2                 | 12 sec         |
| 4   | Continuous | 4                 | 2 min          |
|     | Δ = 1.0    | 4                 | 4 min          |
| 6   | Continuous | 4                 | 16 min         |
|     | Δ = 1.0    | 4                 | 32 min         |
| 8   | Continuous | 8                 | 1 hr           |
|     | Δ = 1.0    | 8                 | 4 hr           |

To calculate the approximate time required to run all 100 Monte Carlo trials
_sequentially_ on similar hardware, multiply the Time per Trial by 100.
To reduce the computational time required, we recommend carrying out these
experiments using nested parallelism:

- Run multiple parallel Monte Carlo experiments (_processes_)
- Each experiment uses multi-threading (_threads_)
- Total number of cores used concurrently: _total = processes × threads_

More specifically, the `mcnp` program supports OpenMP parallelization, using
multiple _threads_ within each simulation to compute the vectorized value
function, log likelihood function, etc.  Setting the `OMP_NUM_THREADS`
environment variable appropriately for each specification is handled
internally by the replication script `table_5.sh` and `run_parallel.sh`.

If you have multiple cores available, the `table_5.sh` script will
automatically distribute Monte Carlo trials across multiple _processes_,
based on the recommended number of Threads Per Trial in the table above.
In that case, divide the expected time by the number of parallel processes.

_Example:_ On the HPC cluster node used, the final table row with N=8 and
Δ=1.0 required approximately 4 hours per Monte Carlo trial, and a full
replication requires 100 trials.  We executed trials concurrently using 6
parallel _processes_ (each with 8 _threads_ for parallel processing within
a trial).

Finally, we note that the scripts described below will automatically compile
the `mcnp` binary using GNU Fortran, but if you wish to manually build it you
can use `make` like so:

```bash
cd mcnp
make clean
make
```

#### Approach 1: Full Replication (100 Monte Carlo Trials Per Specification)

This approach reproduces all results as reported in the paper.

**Requirements**: HPC cluster  
**Estimated time**: ~5 days on 48-core system

```bash
cd mcnp

# Optional: Set number of physical cores (auto-detected if not set)
export OMP_NUM_PROCS=48

# Run all specifications with 100 Monte Carlo trials each
./table_5.sh --full
```

**Note:** The `OMP_NUM_PROCS` variable overrides automatic CPU detection.
Use your physical core count, not hyperthreaded logical cores. To check:
`lscpu` on Linux or `sysctl -n hw.physicalcpu` on macOS.

The script automatically:

- Builds the `mcnp` binary used for the Monte Carlo specifications.
- Runs all 8 specifications (2, 4, 6, 8 firms with both continuous and
  discrete sampling).
    - Carries out each specification by running `mcnp` with the appropriate
      control files and saves log files.
    - Carries out 100 Monte Carlo trials per specification.
    - Uses parallel execution of trials across multiple processes with
      recommended thread allocation per process (handled internally by
      `run_parallel.sh`).
    - Each parallel process will use its own log file.
    - Skips already-completed specifications (resumes after interruption).
- Collects all results and generates the LaTeX table at
  `results/table_5.tex`.

Check status of completed specifications:

```bash
./table_5.sh --status --full
```

Resume interrupted runs:

```bash
# Resumes automatically, skipping completed specifications
./table_5.sh --full
```

#### Approach 2: Incremental Replication (Selected Specifications)

**Requirements**: HPC cluster or high-end, multicore workstation
with ≥32 GB RAM  
**Estimated time**: Variable, depending on hardware and selection

```bash
cd mcnp

# Optional: Set number of physical cores (auto-detected if not set)
export OMP_NUM_PROCS=48

# Choose which specific specifications to run individually.
# For convenience, all 8 available specifications are listed here:
./table_5.sh --full --experiment mc-02-0.0
./table_5.sh --full --experiment mc-02-1.0
./table_5.sh --full --experiment mc-04-0.0
./table_5.sh --full --experiment mc-04-1.0
./table_5.sh --full --experiment mc-06-0.0
./table_5.sh --full --experiment mc-06-1.0
./table_5.sh --full --experiment mc-08-0.0
./table_5.sh --full --experiment mc-08-1.0

# Generate LaTeX table from completed specifications only
./table_5.sh --full --save
```

This approach allows you to:

- Validate that the code works correctly using smaller specifications.
- Distribute computation across multiple machines or sessions.

#### Approach 3: Partial Replication (10 Monte Carlo Trials Per Specification)

**Requirements**: HPC cluster or high-end, multicore workstation
with ≥32 GB RAM  
**Estimated time**: 10 hours on 48-core HPC cluster node

```bash
cd mcnp

# Optional: Set number of physical cores (auto-detected if not set)
export OMP_NUM_PROCS=48

# Run all specifications with only 10 trials each
./table_5.sh --partial
```

This uses special control files (`control/mc-*-partial.ctl`) with `nmc = 10`
instead of 100.

As above, specific specifications can also run in partial replication mode
(first 10 Monte Carlo trials only):

```bash
./table_5.sh --partial --experiment mc-08-0.0
```

The following table displays the expected partial replication results
corresponding to Table 5 in the paper.  It shows the means and standard
deviations of the parameter estimates for only the first 10 of the full
100 Monte Carlo trials reported in the paper.

**Table: Quality Ladder Model Monte Carlo Results (Partial Replication)**

| N   | K      | Sampling       |      | λ_L   | λ_H   | γ     | κ     | η     | μ     |
|-----|--------|----------------|------|-------|-------|-------|-------|-------|-------|
|     |        | DGP            | True | 1.000 | 1.200 | 0.400 | 0.800 | 4.000 | 0.900 |
| 2   | 56     | Continuous     | Mean | 1.003 | 1.210 | 0.401 | 0.809 | 4.035 | 0.899 |
|     |        |                | S.D. | 0.015 | 0.015 | 0.009 | 0.029 | 0.093 | 0.021 |
|     |        | Δ = 1.0        | Mean | 1.120 | 1.320 | 0.399 | 0.939 | 4.398 | 0.950 |
|     |        |                | S.D. | 0.198 | 0.205 | 0.007 | 0.276 | 0.723 | 0.073 |
| 4   | 840    | Continuous     | Mean | 1.004 | 1.205 | 0.405 | 0.802 | 4.067 | 0.900 |
|     |        |                | S.D. | 0.017 | 0.022 | 0.014 | 0.052 | 0.258 | 0.033 |
|     |        | Δ = 1.0        | Mean | 0.952 | 1.147 | 0.400 | 0.692 | 3.796 | 0.891 |
|     |        |                | S.D. | 0.095 | 0.092 | 0.006 | 0.151 | 0.372 | 0.034 |
| 6   | 5,544  | Continuous     | Mean | 1.003 | 1.198 | 0.409 | 0.791 | 4.017 | 0.900 |
|     |        |                | S.D. | 0.015 | 0.023 | 0.024 | 0.034 | 0.200 | 0.022 |
|     |        | Δ = 1.0        | Mean | 0.990 | 1.191 | 0.398 | 0.780 | 3.934 | 0.898 |
|     |        |                | S.D. | 0.084 | 0.086 | 0.006 | 0.143 | 0.335 | 0.025 |
| 8   | 24,024 | Continuous     | Mean | 1.001 | 1.207 | 0.412 | 0.805 | 4.013 | 0.898 |
|     |        |                | S.D. | 0.015 | 0.021 | 0.019 | 0.024 | 0.179 | 0.021 |
|     |        | Δ = 1.0        | Mean | 1.003 | 1.204 | 0.399 | 0.792 | 3.937 | 0.905 |
|     |        |                | S.D. | 0.131 | 0.131 | 0.004 | 0.172 | 0.377 | 0.052 |

### Automated Full Replication (All Tables)

For users with access to HPC resources, a `main.sh` script is provided that
automates the complete replication of all tables:

```bash
./main.sh
```

This script runs each table in order (Tables 1-5), with parallelization
within each table as appropriate.

Important: Given the computational requirements, particularly for Table 5,
this script should only be used on HPC systems with substantial resources
(≥32 GB RAM, ≥40 cores).

For most users, we recommend running the individual table scripts separately
as described above, which provides more flexibility and control over the
replication process.

## List of Tables and Programs

The provided code reproduces all tables in the paper.  Figures in the paper do
not require code.  The following table summarizes the complete list of tables
and programs.  Please see the detailed replication instructions above for each.

| Output  | Program                     | Description                                      | Output File                |
|---------|-----------------------------|--------------------------------------------------|----------------------------|
| Table 1 | `mc1p/table_1_2.sh`         | Rust (1987) Sample Characteristics               | `mc1p/results/table_1.tex` |
| Table 2 | `mc1p/table_1_2.sh`         | Estimates Based on Data from Rust (1987)         | `mc1p/results/table_2.tex` |
| Table 3 | `mc1p/table_3.sh`           | Single Agent Renewal Monte Carlo Results         | `mc1p/results/table_3.tex` |
| Table 4 | `mcnp/table_4.sh`           | Quality Ladder Monte Carlo Specifications        | `mcnp/results/table_4.tex` |
| Table 5 | `mcnp/table_5.sh --partial` | Quality Ladder Monte Carlo (Partial Replication) | `mcnp/results/table_5.tex` |
| Table 5 | `mcnp/table_5.sh --full`    | Quality Ladder Monte Carlo (Full Replication)    | `mcnp/results/table_5.tex` |

**Notes:**

- Runtimes given are for replicating the complete tables, taking advantage of
  parallelism on the hardware shown.  Actual runtime may vary significantly by
  Fortran compiler used, hardware capabilities, and the level of
  multiprocessing used.

- The HPC cluster nodes we used have 48 Intel Xeon Platinum 8260 2.40GHz CPU
  cores and 192 GB memory.

- Table 5 is the most computationally intensive to replicate.  We executed
  Monte Carlo trials in parallel across processors on an HPC cluster node.
  See the full replication instructions above for details and recommendations.

- Original log files and results obtained by the author are provided in the
  respective `mc1p/original_logs/` and `mcnp/original_logs/` directories.

## References

- Byrd, R. H., P. Lu, and J. Nocedal (1995).
  [A limited memory algorithm for bound constrained optimization](https://doi.org/10.1137/0916069).
  _SIAM Journal on Scientific and Statistical Computing_ 16, 1190–1208.

- Ericson, R. and A. Pakes (1995).
  [Markov-Perfect Industry Dynamics: A Framework for Empirical Work](https://doi.org/10.2307/2297841).
  _Review of Economic Studies_ 62, 53–82.

- Rust, J. (1987).
  [Optimal Replacement of GMC Bus Engines: An Empirical Model of Harold Zurcher](https://doi.org/10.2307/1911259).
  _Econometrica_ 55, 999–1033.

- Sidje, R. B. (1998).
  [Expokit: A software package for computing matrix exponentials](https://doi.org/10.1145/285861.285868).
  _ACM Transactions on Mathematical Software_ 24, 130–156.

- Zhu, C., R. H. Byrd, P. Lu, and J. Nocedal (1997).
  [Algorithm 778: L-BFGS-B, FORTRAN routines for large scale bound constrained optimization](https://doi.org/10.1145/279232.279236).
  _ACM Transactions on Mathematical Software_ 23, 550-560.

## License

The replication code is licensed under the BSD License. See [LICENSE](LICENSE)
for details.

The [L-BFGS-B][] license is provided in the `mc1p/lbfgsb/` directory:
[License.txt](lbfgsb/License.txt).

The [Expokit][] license is provided in the `mc1p/expokit/` directory:
[LICENSE](expokit/LICENSE).


[jblevins]: https://jblevins.org/
[ctgames]: https://jblevins.org/research/ctgames
[L-BFGS-B]: https://users.iems.northwestern.edu/~nocedal/lbfgsb.html
[Expokit]: https://www.maths.uq.edu.au/expokit/
