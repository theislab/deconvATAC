from rpy2 import robjects
import anndata2ri
import pandas as pd

def card(
    adata_spatial, 
    adata_ref, 
    labels_key,
    spatial_count_layer="X",
    sc_count_layer="X",
    r_lib_path=None, 
    results_path="./card_results", 
    card_kwargs={}
): 
    '''
    Run CARD (Cell type Assignment by Reference-based Deconvolution) method

    Parameters 
    -----------
    adata_spatial : AnnData
        AnnData object for spatial transcriptomics data.
    adata_ref : AnnData 
        AnnData object for reference single-cell RNA-seq data.
    labels_key : str
        Key in adata_ref.obs for cell type labels.
    spatial_count_layer : str
        Layer in adata_spatial for count data. Defaults to 'X'.
    sc_count_layer : str
        Layer in adata_ref for count data. Defaults to 'X'.
    r_lib_path : str
        Path to R library if necessary.
    results_path : str
        Directory to save the results.
    card_kwargs : dict
        Additional arguments for CARD functions.

    Returns
    --------
    - Saves estimated cell type proportions to results_path.
    '''
    # Setup R environment
    if r_lib_path is not None:
        robjects.r.assign("lib_path", r_lib_path)
        robjects.r(".libPaths(lib_path)")
    robjects.r("library(CARD)")

    # Convert Python objects to R
    anndata2ri.activate()
    scexp_ref = anndata2ri.py2rpy(adata_ref)
    scexp_spatial = anndata2ri.py2rpy(adata_spatial)
    
    # Prepare and run CARD
    robjects.r.assign("scexp_ref", scexp_ref)
    robjects.r.assign("scexp_spatial", scexp_spatial)
    robjects.r.assign("labels_key", labels_key)
    robjects.r.assign("results_path", results_path)
    robjects.r.assign("spatial_count_layer", spatial_count_layer)
    robjects.r.assign("sc_count_layer", sc_count_layer)

    card_r_command = f"""
        # Prepare data
        sc_data <- assay(scexp_ref, "{sc_count_layer}")
        spatial_data <- assay(scexp_spatial, "{spatial_count_layer}")

        # Create CARD object
        card_obj <- createCARDObject(
            sc_count = sc_data,
            sc_meta = colData(scexp_ref),
            spatial_count = spatial_data,
            spatial_location = scexp_spatial@int_colData@listData$reducedDims$spatial,
            ct.varname = "{labels_key}",
            {', '.join(f'{k}={v}' for k, v in card_kwargs.items())}
        )

        # Run CARD deconvolution
        card_obj <- CARD_deconvolution(CARD_object = card_obj)
        
        # Save results
        write.csv(card_obj@Proportion_CARD, file = paste0(results_path, "/estimated_proportions.csv"))
    """
    robjects.r(card_r_command)
    anndata2ri.deactivate()


# Example usage
# card_deconvolution(adata_spatial=my_spatial_data, adata_ref=my_sc_data, labels_key="cell_type")
