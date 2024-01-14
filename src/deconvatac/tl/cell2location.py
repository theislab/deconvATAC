import os
import cell2location as c2l


def cell2location(
    adata_spatial,
    adata_ref,
    N_cells_per_location,
    detection_alpha,
    labels_key=None,
    layer_spatial=None,
    layer_ref=None,
    use_gpu=True,
    max_epochs_spatial=30000,
    max_epochs_ref=None,
    return_adatas=False,
    plots=True,
    results_path="./cell2location_results",
    setup_ref_kwargs={},
    train_ref_kwargs={},
    setup_spatial_kwargs={},
    train_spatial_kwargs={},
):
    """
    Run Cell2Location

    Parameters
    -----------

    adata_spatial : AnnData
        AnnData of the spatial data, filtered by highly variable features. Feature space needs to be the same as the one of adata_ref.
    adata_ref : AnnData
        AnnData of the reference data, filtered by highly variable features. Feature space needs to be the same as the one of adata_spatial.
    N_cells_per_location : float
        Expected cell number per location.
    detection_alpha : float
        Regularisation of per-location normalisation.
    labels_key : str
        Cell type key in adata_ref.obs for label information
    layer_spatial : str
        Layer of adata_spatial to use for deconvolution. If None, uses adata_spatial.X.
    layer_ref : str
        Layer of adata_ref to use for deconvolution. If None, uses adata_ref.X.
    use_gpu : bool
        Whether to use the GPU.
    max_epochs_spatial: int
        Number of epochs for the spatial mapping model. If None, defaults to np.min([round((20000 / n_cells) * 400), 400]).
    max_epochs_ref: int
        Number of epochs for the reference model. If None, defaults to np.min([round((20000 / n_cells) * 400), 400]).
    return_adatas: bool
        Whether to return AnnDatas with deconvolution results. Returns tupel: (adata_spatial, adata_ref).
    plots: bool
        Whether to plot QC and ELBO plots.
    results_path: str
        Path to save estimated cell type abundances to.
    setup_ref_kwargs: dict
        Parameters for cell2location.models.RegressionModel.setup_anndata()
    train_ref_kwargs: dict
        Parameters for cell2location.models.RegressionModel.train()
    setup_spatial_kwargs: dict
        Parameters for cell2location.models.Cell2location.setup_anndata()
    train_spatial_kwargs: dict
        Parameters for cell2location.models.Cell2location.train()


    Returns
    --------

    - Saves 'q05_cell_abundance_w_sf' and 'means_cell_abundance_w_sf' as csv-files to results_path.
    - If return_adatas=True, returns tupel (adata_spatial, adata_ref) with saved deconvolution results.
    """

    adata_ref_copy = adata_ref.copy()
    adata_spatial_copy = adata_spatial.copy()

    # 1. Fit sc model
    c2l.models.RegressionModel.setup_anndata(
        adata=adata_ref_copy, layer=layer_ref, labels_key=labels_key, **setup_ref_kwargs
    )
    model_ref = c2l.models.RegressionModel(adata_ref_copy)
    if plots:
        model_ref.view_anndata_setup()

    model_ref.train(max_epochs=max_epochs_ref, use_gpu=use_gpu, **train_ref_kwargs)
    if plots:
        model_ref.plot_history(20)

    adata_ref_copy = model_ref.export_posterior(adata_ref_copy)
    if plots:
        model_ref.plot_QC()

    if "means_per_cluster_mu_fg" in adata_ref_copy.varm.keys():
        inf_aver = adata_ref_copy.varm["means_per_cluster_mu_fg"][
            [f"means_per_cluster_mu_fg_{i}" for i in adata_ref_copy.uns["mod"]["factor_names"]]
        ].copy()
    else:
        inf_aver = adata_ref_copy.var[
            [f"means_per_cluster_mu_fg_{i}" for i in adata_ref_copy.uns["mod"]["factor_names"]]
        ].copy()
    inf_aver.columns = adata_ref_copy.uns["mod"]["factor_names"]

    # 2. Fit spatial model
    c2l.models.Cell2location.setup_anndata(adata=adata_spatial_copy, layer=layer_spatial, **setup_spatial_kwargs)
    model_st = c2l.models.Cell2location(
        adata_spatial_copy,
        cell_state_df=inf_aver,
        N_cells_per_location=N_cells_per_location,
        detection_alpha=detection_alpha,
    )
    if plots:
        model_st.view_anndata_setup()

    model_st.train(max_epochs=max_epochs_spatial, use_gpu=use_gpu, **train_spatial_kwargs)
    if plots:
        model_st.plot_history(1000)

    adata_spatial_copy = model_st.export_posterior(adata_spatial_copy)
    if plots:
        model_st.plot_QC()

    # 3. Save results
    if not os.path.exists(results_path):
        os.mkdir(results_path)
    adata_spatial_copy.obsm["q05_cell_abundance_w_sf"].to_csv(results_path + "/q05_cell_abundance_w_sf.csv")
    adata_spatial_copy.obsm["means_cell_abundance_w_sf"].to_csv(results_path + "/means_cell_abundance_w_sf.csv")

    if return_adatas:
        return adata_spatial_copy, adata_ref_copy
