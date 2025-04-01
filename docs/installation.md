
## Installation


We recommend running deconvATAC within virtual environments, such as Conda, to prevent conflicts.


### Create conda environment

```bash
conda create -n deconvATAC python=3.9 r-base=4.3.0
conda activate deconvATAC
```

### Installing deconvATAC

First, clone the directory: 
```bash
git clone https://github.com/theislab/deconvATAC.git
```

Install the package: 
```bash
cd deconvATAC
pip install .
```
.. note:: 
    If you encounter issues with `glibc` during the installation you can try to install it using conda:
    `conda create -n deconvATAC python=3.9 r-base=4.3.0 gcc_linux-64 gxx_linux-64`
### Installing optional dependencies

deconvATAC is only installed with the packages needed for the simulation, highly variable peak selection, and metrics. For running the deconvolution methods, we recommend to work with a different environment for each method to prevent dependency conflicts, with deconvATAC installed in each. 
You can install the dependencies needed for the python-based deconvolution methods with: 

```bash
pip install .[cell2location] # note: for zsh shell, please use brackets: '.[cell2location]'
pip install .[tangram]
pip install .[destvi]
```

#### RCTD

For installing RCTD, please use the following 
```bash
conda install bioconda::r-spacexr
```
In your R terminal, install
```bash
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("S4Vectors")
BiocManager::install("SingleCellExperiment")
```

#### SpatialDWLS

For SpatialDWLS, the Giotto package needs to be installed. Please follow the installation guidelines in the [Giotto documentation](https://drieslab.github.io/Giotto_website/articles/installation.html) for installation of the package. 


