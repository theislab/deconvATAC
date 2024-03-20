import os
import numpy as np
import pandas as pd
import scanpy as sc


def tangram(
    adata_spatial,
    adata_ref,
    labels_key,
    run_rank_genes=False,
    layer_rank_genes=None,
    num_epochs=1000,
    device="cpu",
    return_adatas=False,
    result_path="./tangram_results",
    **kwargs,
):
    """
    Run Tangram

    Parameters
    -----------

    adata_spatial : AnnData
        AnnData of the spatial data.
    adata_ref : AnnData
        AnnData of the reference data.
    labels_key : str
        Cell type key in adata_ref.obs for label information
    run_rank_genes: bool
        If true, will run sc.tl.rank_genes_groups on reference anndata followed by tg.pp_adatas.
        If false, will only run tg.pp_adatas with all peaks of the input anndatas. Thus, expects anndatas to be already filtered by HVFs.
    layer_rank_genes: str
        Only if run_rank_genes is true. Layer to use for sc.tl.rank_genes_groups.
    num_epochs: int
        Number of epochs.
    device : string or torch.device
        Which device to use.
    return_adatas: bool
        Whether to return AnnDatas with deconvolution results. Returns tupel: (adata_spatial, adata_ref).
    results_path: str
        Path to save estimated cell type abundances to.
    **kwargs:
        Parameters for tangram.mapping_utils.map_cells_to_space()

    Returns
    --------

    - Saves estimated proportions as csv-file to results_path.
    - If return_adatas=True, returns tupel (adata_spatial, adata_ref) with saved deconvolution results.
    """
    import tangram as tg

    adata_ref_copy = adata_ref.copy()
    adata_spatial_copy = adata_spatial.copy()

    # 1. Preprocess anndatas
    if run_rank_genes:
        sc.tl.rank_genes_groups(adata_ref_copy, groupby=labels_key, layer=layer_rank_genes, method="wilcoxon")
        markers_df = pd.DataFrame(adata_ref_copy.uns["rank_genes_groups"]["names"]).iloc[0:1000, :]
        markers = list(np.unique(markers_df.melt().value.values))
        tg.pp_adatas(adata_sc=adata_ref_copy, adata_sp=adata_spatial_copy, genes=markers)
    else:
        tg.pp_adatas(adata_sc=adata_ref_copy, adata_sp=adata_spatial_copy, genes=list(adata_ref_copy.var.index))

    # 2. Run tangram
    adata_results = tg.mapping_utils.map_cells_to_space(
        adata_sc=adata_ref_copy, adata_sp=adata_spatial_copy, num_epochs=num_epochs, device=device, **kwargs
    )

    # 3. Extract and save results
    tg.project_cell_annotations(adata_results, adata_spatial_copy, annotation=labels_key)

    if not os.path.exists(result_path):
        os.mkdir(result_path)
    adata_spatial_copy.obsm["tangram_ct_pred"].to_csv(result_path + "/tangram_ct_pred.csv")

    if return_adatas:
        return adata_spatial_copy, adata_ref_copy
