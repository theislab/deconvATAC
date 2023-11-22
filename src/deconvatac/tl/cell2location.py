import os
import cell2location as c2l

def cell2location(
    adata_st,
    adata_ref,
    N_cells_per_location,
    detection_alpha,
    return_adatas=False,
    use_gpu=True,
    plots=True,
    result_path="./cell2location_results",
    layer_st=None,
    layer_ref=None,
    batch_key_st=None,
    batch_key_ref=None,
    labels_key_st=None,
    labels_key_ref=None,
    categorical_covariate_keys_st=None,
    categorical_covariate_keys_ref=None,
    continuous_covariate_keys_st=None,
    continuous_covariate_keys_ref=None,
    max_epochs_st=30000,
    max_epochs_ref=None,
    batch_size_st=None,
    batch_size_ref=2500,
    train_size_st=1,
    train_size_ref=1,
    lr_st=0.002,
    lr_ref=0.002,
):
    ## TODO: add parameter descriptions
    # expects adata_st and adata_ref to be already filtered by highly variable features

    # 1. Fit sc model
    c2l.models.RegressionModel.setup_anndata(
        adata=adata_ref,
        layer=layer_ref,
        batch_key=batch_key_ref,
        labels_key=labels_key_ref,
        categorical_covariate_keys=categorical_covariate_keys_ref,
        continuous_covariate_keys=continuous_covariate_keys_ref,
    )
    model_ref = c2l.models.RegressionModel(adata_ref)  # also add parameters of this function?
    if plots:
        model_ref.view_anndata_setup()
    model_ref.train(
        max_epochs=max_epochs_ref, use_gpu=use_gpu, batch_size=batch_size_ref, train_size=train_size_ref, lr=lr_ref
    )
    if plots:
        model_ref.plot_history(20)
    adata_ref = model_ref.export_posterior(
        adata_ref,
        use_quantiles=True,
        add_to_obsm=["q05", "q50", "q95", "q0001"],
        sample_kwargs={"batch_size": batch_size_ref, "use_gpu": use_gpu},
    )

    if plots:
        model_ref.plot_QC()
    if "means_per_cluster_mu_fg" in adata_ref.varm.keys():
        inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][
            [f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]
        ].copy()
    else:
        inf_aver = adata_ref.var[[f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]].copy()
    inf_aver.columns = adata_ref.uns["mod"]["factor_names"]  # also save this table?

    # 2. Fit spatial model
    c2l.models.Cell2location.setup_anndata(
        adata=adata_st,
        layer=layer_st,
        batch_key=batch_key_st,
        labels_key=labels_key_st,
        categorical_covariate_keys=categorical_covariate_keys_st,
        continuous_covariate_keys=continuous_covariate_keys_st,
    )

    model_st = c2l.models.Cell2location(
        adata_st, cell_state_df=inf_aver, N_cells_per_location=N_cells_per_location, detection_alpha=detection_alpha
    )
    if plots:
        model_st.view_anndata_setup()
    model_st.train(
        max_epochs=max_epochs_st, use_gpu=use_gpu, batch_size=batch_size_st, train_size=train_size_st, lr=lr_st
    )
    if plots:
        model_st.plot_history(1000)
    adata_st = model_st.export_posterior(
        adata_st, sample_kwargs={"num_samples": 1000, "batch_size": model_st.adata.n_obs, "use_gpu": use_gpu}
    )
    if plots:
        model_st.plot_QC()

    # 3. Save results
    if not os.path.exists(result_path):
        os.mkdir(result_path)
    adata_st.obsm["q05_cell_abundance_w_sf"].to_csv(result_path + "/q05_cell_abundance_w_sf.csv")
    adata_st.obsm["means_cell_abundance_w_sf"].to_csv(result_path + "/means_cell_abundance_w_sf.csv")

    if return_adatas:
        return adata_st, adata_ref

