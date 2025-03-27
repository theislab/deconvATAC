# deconvATAC

The deconvATAC package provides code used in our benchmarking study for deconvoluting spatialATAC data via deconvolution tools designed for spatial transcriptomics. In our study, we benchmark five top-performing spatial transcriptomics deconvolution methods. deconvATAC additionally provides a framework for simulating spatial multi-modal data from dissociated single-cell data, as well as metrics for evaluating the performance of deconvolution. 


Please refer to the [documentation][link-docs].

Data used in this study is available on [Zenodo]()

<p align="left">
<img src="https://github.com/theislab/deconvATAC/blob/docs_sarah/docs/figure1.png/?raw=true" alt="Study overview" width="700"/>



## Installation


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

### Installing optional dependencies

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



## Citation

> t.b.a


[issue-tracker]: https://github.com/theislab/deconvATAC/issues
[link-docs]: https://deconvATAC.readthedocs.io
[link-api]: https://deconvatac.readthedocs.io/en/latest/autoapi/index.html
