import os
import numpy as np
import pandas as pd
import scanpy as sc
import tangram as tg





def tangram(
    adata_st,
    adata_ref,
    labels_key,
    return_adatas=False,
    run_rank_genes=False,
    layer_rank_genes=None,
    result_path="./tangram_results",
    num_epochs=1000,
    device="cpu",
    **kwargs,
):
    """
    run_rank_genes: If true, will run sc.tl.rank_genes_groups on reference anndata followed by tg.pp_adatas
                    If false, will only run tg.pp_adatas with all peaks of the input anndatas. Thus, expects anndatas to be already filtered by HVFs
    layer_rank_genes: Only if run_rank_genes is true. Layer to use for sc.tl.rank_genes_groups
    """

    # 1. Preprocess anndatas
    if run_rank_genes:
        sc.tl.rank_genes_groups(adata_ref, groupby=labels_key, layer=layer_rank_genes, method="wilcoxon")
        markers_df = pd.DataFrame(adata_ref.uns["rank_genes_groups"]["names"]).iloc[0:1000, :]
        markers = list(np.unique(markers_df.melt().value.values))
        tg.pp_adatas(adata_sc=adata_ref, adata_sp=adata_st, genes=markers)
    else:
        tg.pp_adatas(adata_sc=adata_ref, adata_sp=adata_st, genes=list(adata_ref.var.index))

    # 2. Run tangram
    adata_results = tg.mapping_utils.map_cells_to_space(
        adata_sc=adata_ref, adata_sp=adata_st, num_epochs=num_epochs, device=device, **kwargs
    )

    # 3. Extract and save results
    tg.project_cell_annotations(adata_results, adata_st, annotation=labels_key)

    if not os.path.exists(result_path):
        os.mkdir(result_path)
    adata_st.obsm["tangram_ct_pred"].to_csv(result_path + "/tangram_ct_pred.csv")

    if return_adatas:
        return adata_st, adata_ref