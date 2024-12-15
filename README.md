
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` 
up-to-date. 
`devtools::build_readme()` is handy for this. -->

# The Multimodal Marginal Epistasis test <img src="man/figures/logo.png" align="right" height="200" alt="" />

<!-- badges: start -->
<!-- badges: end -->

The `mmer` package implements a computationally and statistically
efficient method for detecting marginal epistasis in genome-wide
association studies (GWAS). Find the full package documentation
including examples and articles here: [Multimodal Marginal Epistasis
test Documentation](https://lcrawlab.github.io/mmer/).

## Key Features

- Hutchinson’s stochastic trace estimator: efficient and scalable
  computation
- Mailman algorithm: fast vector-by-matrix operation
- Linear mixed model: controls for population structure
- Multimodal Input: incorporates additional data from HDF5 files to
  improve power in detecting gene-by-gene interactions.
- Optimize for Memory Constraints: Highly configurable block wise
  processing of the data allows to make the most of available resources.
  See also [How To Optimize the Memory Requirements of
  MME](articles/tutorial-memory-optimization.html).
- Parallelization: Utilizes OpenMP for multi-threaded processing.

## Installation

You can install the development version of mmer from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("lcrawlab/mmer")
```

## Dependencies

System requirements of the package:

- GNU make
- R (\>= 4.4)
- Rhdf5lib (from BioConductor)
- OpenMP (optional)

To install `Rhdf5lib`, first install the tool `BiocManager` from CRAN,
then install the library using this tool.

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rhdf5lib")
```

The full list of R dependencies can be found in the [DESCRIPTION
file](https://github.com/lcrawlab/mmer/blob/main/DESCRIPTION).

### OpenMP

For OS X and Linux, the OpenMP library can be installed via one of the
(shell) commands specified below:

| System | Command |
|:---|:---|
| **OS X (using Homebrew)** | `brew install libomp` |
| **Debian-based systems (including Ubuntu)** | `sudo apt-get install libomp-dev` |

To enable openMP, it may be necessary to configure the compiler flags
`SHLIB_OPENMP_CXXFLAGS` and `LDFLAGS` in the `~/.R/Makevars` file.

| System | Required Flags           |
|--------|--------------------------|
| OS X   | `-Xclang -fopenmp -lomp` |
| Linux  | `-fopenmp -lomp`         |

## Known Issues

Compiling the package requires the compiler to find the libraries for
the dependencies. For unix systems, the libraries are typically
installed at `/usr/local/lib` and `/usr/local/include`. For users using
OS X and homebrew, the libraries are typically installed at
`/opt/homebrew/lib` and `/opt/homebrew/include`.

Non-standard library paths need to be configured. The `src/Makevars`
file configures the compiler flags and considers the `LDFLAGS` and
`CPPFLAGS` from the `~/.R/Makevars` file.

## References

- Stamp J, Crawford L (2024). mvMAPIT: Multivariate Genome Wide Marginal
  Epistasis Test. R package version 2.0.3,
  <https://lcrawlab.github.io/mvMAPIT/>,
  <https://github.com/lcrawlab/mvMAPIT>.
- Stamp et al. (2023): Leveraging genetic correlation between traits for
  epistasis detection in GWAS. G3: Genes, Genomes, Genetics.
- Fu, B., Pazokitoroudi, A., Xue, A., Anand, A., Anand, P., Zaitlen, N.,
  & Sankararaman, S. (2023). A biobank-scale test of marginal epistasis
  reveals genome-wide signals of polygenic epistasis. bioRxiv.
- Crawford et al. (2017): Detecting epistasis with the marginal
  epistasis test. PLoS Genetics.
- Devresse et al. (2024): HighFive - Header-only C++ HDF5 interface.
  <https://zenodo.org/records/13120799>
