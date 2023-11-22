import os
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import scvi


def destvi(
    adata_st,
    adata_ref,
    return_adatas=False,
    use_gpu=True,
    plots=True,
    result_path="./destvi_results",
    # parameters for setup
    layer_ref=None,
    layer_st=None,
    labels_key_ref=None,
    n_hidden_ref=128,
    n_latent_ref=5,
    n_layers_ref=2,
    weight_obs_ref=False,
    dropout_rate_ref=0.05,
    vamp_prior_p_st=15,
    l1_reg_st=0.0,
    # training parameters:
    max_epochs_ref=300,
    max_epochs_st=2000,
    lr_ref=0.001,
    lr_st=0.003,
    train_size_ref=1,
    train_size_st=1.0,
    validation_size_ref=None,
    validation_size_st=None,
    shuffle_set_split_ref=True,
    shuffle_set_split_st=True,
    batch_size_ref=128,
    batch_size_st=128,
    plan_kwargs_ref=None,
    plan_kwargs_st=None,
    n_epochs_kl_warmup_st=200,
):
    ## TODO: add parameter descriptions
    # expects adata_st and adata_ref to be already filtered by highly variable features

    # 1. Fit scLVM
    adata_ref_copy = adata_ref.copy()
    scvi.model.CondSCVI.setup_anndata(adata=adata_ref_copy, layer=layer_ref, labels_key=labels_key_ref)
    model_ref = scvi.model.CondSCVI(
        adata_ref_copy,
        n_hidden=n_hidden_ref,
        n_latent=n_latent_ref,
        n_layers=n_layers_ref,
        weight_obs=weight_obs_ref,
        dropout_rate=dropout_rate_ref,
    )
    if plots:
        model_ref.view_anndata_setup()

    model_ref.train(
        max_epochs=max_epochs_ref,
        lr=lr_ref,
        use_gpu=use_gpu,
        train_size=train_size_ref,
        validation_size=validation_size_ref,
        shuffle_set_split=shuffle_set_split_ref,
        batch_size=batch_size_ref,
        plan_kwargs=plan_kwargs_ref,
    )
    if plots:
        model_ref.history["elbo_train"].iloc[5:].plot()
        plt.show()
        # plot latent space
        adata_ref_copy.obsm["scLVM_latent"] = model_ref.get_latent_representation()
        sc.pp.neighbors(adata_ref_copy, use_rep="scLVM_latent")
        sc.tl.umap(adata_ref_copy)
        sc.pl.umap(adata_ref_copy, color=labels_key_ref, frameon=False, title="UMAP of scLVM latent space")

    # 2. Fit stLVM
    adata_st_copy = adata_st.copy()
    scvi.model.DestVI.setup_anndata(adata_st_copy, layer=layer_st)
    model_st = scvi.model.DestVI.from_rna_model(
        st_adata=adata_st_copy, sc_model=model_ref, vamp_prior_p=vamp_prior_p_st, l1_reg=l1_reg_st
    )
    model_st.train(
        max_epochs=max_epochs_st,
        lr=lr_st,
        use_gpu=use_gpu,
        train_size=train_size_st,
        validation_size=validation_size_st,
        shuffle_set_split=shuffle_set_split_st,
        batch_size=batch_size_st,
        n_epochs_kl_warmup=n_epochs_kl_warmup_st,
        plan_kwargs=plan_kwargs_st,
    )
    if plots:
        model_st.history["elbo_train"].iloc[10:].plot()
        plt.show()

    # 3. Extract results
    proportions = model_st.get_proportions()

    if not os.path.exists(result_path):
        os.mkdir(result_path)
    proportions.to_csv(result_path + "/predicted_proportions.csv")

    if return_adatas:
        adata_st_copy.obsm["proportions"] = proportions
        ct_list = proportions.columns.values
        for ct in ct_list:
            data = adata_st_copy.obsm["proportions"][ct].values
            adata_st_copy.obs[ct] = np.clip(data, 0, np.quantile(data, 0.99))
        return adata_st_copy, adata_ref_copy