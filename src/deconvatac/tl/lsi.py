
import muon as mu 
import numpy as np 


def lsi(adata, num_components=50):
    """
    Runs LSI on adata.X and uses HVFs if adata.var["highly_variable"] is present. Otherwise, runs LSI on all features.

    :param adata: anndata to run LSI on
    :param num_components: number of components to compute
    """

    if "highly_variable" in adata.var.columns:
        adata_hvfs = adata[:, adata.var["highly_variable"]]
        mu.atac.tl.lsi(adata_hvfs, n_comps=num_components)
        # save output to original anndata
        adata.obsm["X_lsi"] = adata_hvfs.obsm['X_lsi']
        adata.uns["lsi"] = adata_hvfs.uns['lsi']
        adata.varm["LSI"] = np.zeros(shape=(adata.n_vars, adata_hvfs.varm["LSI"].shape[1]))
        adata.varm["LSI"][adata.var["highly_variable"]] = adata_hvfs.varm['LSI']
    else:
        mu.atac.tl.lsi(adata, n_comps=num_components)    
