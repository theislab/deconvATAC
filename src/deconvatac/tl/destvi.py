import os
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt


def destvi(
    adata_spatial,
    adata_ref,
    labels_key=None,
    layer_spatial=None,
    layer_ref=None,
    use_gpu=True,
    max_epochs_spatial=2000,
    max_epochs_ref=300,
    return_adatas=False,
    plots=True,
    results_path="./destvi_results",
    model_ref_kwargs={},
    train_ref_kwargs={},
    model_spatial_kwargs={}, 
    train_spatial_kwargs={}
):
    '''
    Run DestVI

    Parameters 
    -----------

    adata_spatial : AnnData
        AnnData of the spatial data, filtered by highly variable features. Feature space needs to be the same as the one of adata_ref. 
    adata_ref : AnnData 
        AnnData of the reference data, filtered by highly variable features. Feature space needs to be the same as the one of adata_spatial.
    labels_key : str
        Cell type key in adata_ref.obs for label information
    layer_spatial : str
        Layer of adata_spatial to use for deconvolution. If None, uses adata_spatial.X.
    layer_ref : str
        Layer of adata_ref to use for deconvolution. If None, uses adata_ref.X.
    use_gpu : bool
        Whether to use the GPU.
    max_epochs_spatial: int
        Number of epochs for the stLVM.
    max_epochs_ref: int 
        Number of epochs for the scLVM. 
    return_adatas: bool 
        Whether to return AnnDatas with deconvolution results. Returns tupel: (adata_spatial, adata_ref).
    plots: bool 
        Whether to plot ELBO plots and UMAP of scLVM latent space.
    results_path: str
        Path to save estimated cell type abundances to. 
    model_ref_kwargs: dict
        Parameters for scvi.model.CondSCVI()
    train_ref_kwargs: dict
        Parameters for scvi.model.CondSCVI.train()
    model_spatial_kwargs: dict
        Parameters for scvi.model.DestVI.from_rna_model()
    train_spatial_kwargs: dict
        Parameters for scvi.model.DestVI.train()
    
        
    Returns
    --------

    - Saves estimated proportions as csv-file to results_path.
    - If return_adatas=True, returns tupel (adata_spatial, adata_ref) with saved deconvolution results. 
    '''
    import scvi
    
   # 1. Fit scLVM
    adata_ref_copy = adata_ref.copy()
    scvi.model.CondSCVI.setup_anndata(adata=adata_ref_copy, layer=layer_ref, labels_key=labels_key)
    model_ref = scvi.model.CondSCVI(adata_ref_copy,**model_ref_kwargs)
    if plots:
        model_ref.view_anndata_setup()

    model_ref.train(max_epochs=max_epochs_ref,use_gpu=use_gpu,**train_ref_kwargs)
    if plots:
        model_ref.history["elbo_train"].iloc[5:].plot()
        plt.show()
        # plot latent space
        adata_ref_copy.obsm["scLVM_latent"] = model_ref.get_latent_representation()
        sc.pp.neighbors(adata_ref_copy, use_rep="scLVM_latent")
        sc.tl.umap(adata_ref_copy)
        sc.pl.umap(adata_ref_copy, color=labels_key, frameon=False, title="UMAP of scLVM latent space")

    # 2. Fit stLVM
    adata_spatial_copy = adata_spatial.copy()
    scvi.model.DestVI.setup_anndata(adata_spatial_copy, layer=layer_spatial)
    model_st = scvi.model.DestVI.from_rna_model(st_adata=adata_spatial_copy, sc_model=model_ref,**model_spatial_kwargs)
    model_st.train(max_epochs=max_epochs_spatial,use_gpu=use_gpu,**train_spatial_kwargs)
    if plots:
        model_st.history["elbo_train"].iloc[10:].plot()
        plt.show()

    # 3. Extract results
    proportions = model_st.get_proportions()

    if not os.path.exists(results_path):
        os.mkdir(results_path)
    proportions.to_csv(results_path + "/predicted_proportions.csv")

    if return_adatas:
        adata_spatial_copy.obsm["proportions"] = proportions
        ct_list = proportions.columns.values
        for ct in ct_list:
            data = adata_spatial_copy.obsm["proportions"][ct].values
            adata_spatial_copy.obs[ct] = np.clip(data, 0, np.quantile(data, 0.99))
        return adata_spatial_copy, adata_ref_copy