# deconvATAC

The deconvATAC package provides code used in our benchmarking study for deconvoluting spatialATAC data via deconvolution tools designed for spatial transcriptomics. In our study, we benchmark five top-performing spatial transcriptomics deconvolution methods. deconvATAC additionally provides a framework for simulating spatial multi-modal data from dissociated single-cell data, as well as metrics for evaluating the performance of deconvolution. 


Please refer to the [documentation][link-docs].

Data used in this study is available on [Zenodo](https://zenodo.org/records/15089738)

<p align="left">
<img src="https://github.com/theislab/deconvATAC/blob/main/docs/figure1.png/?raw=true" alt="Study overview" width="700"/>



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
> [!NOTE]  
> If you encounter issues with `glibc` during the installation you can try to install it using conda:
> ```conda create -n deconvATAC python=3.9 r-base=4.3.0 gcc_linux-64 gxx_linux-64```
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


> **Spatial transcriptomics deconvolution methods generalize well to spatial chromatin accessibility data**
>
> Sarah Ouologuem, Laura D Martens, Anna C Schaar, Maiia Shulman, Julien Gagneur, Fabian J Theis
>
> _Bioinformatics_, Volume 41, Issue Supplement_1, July 2025, Pages i314–i322 doi: [10.1093/bioinformatics/btaf268](https://doi.org/10.1093/bioinformatics/btaf268).


[issue-tracker]: https://github.com/theislab/deconvATAC/issues
[link-docs]: https://deconvATAC.readthedocs.io
[link-api]: https://deconvatac.readthedocs.io/en/latest/autoapi/index.html
