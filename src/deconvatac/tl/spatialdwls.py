from rpy2 import robjects
import anndata2ri
import scanpy as sc
import muon as mu 
import numpy as np 


def spatialdwls(
    adata_spatial, 
    adata_ref, 
    labels_key,
    n_cell=50,
    highly_variable_key_ref="highly_variable",
    highly_variable_key_spatial="highly_variable",
    layer_ref="X",
    layer_spatial="X",
    dimred = "LSI",
    r_lib_path=None, 
    results_path="./spatialdwls_results", 
): 
    '''
    Run SpatialDWLS

    Parameters 
    -----------

    adata_spatial : AnnData
        AnnData of the spatial data.
    adata_ref : AnnData 
        AnnData of the reference data.
    labels_key : str
        Cell type key in adata_ref.obs for label information.
    n_cell : int
        Number of cells per spot.
    highly_variable_key_ref: str
        Key in adata_ref.var of highly variable features. 
    highly_variable_key_spatial: str
        Key in adata_spatial.var of highly variable features. If None, uses all features for PCA/LSI and leiden clustering
    layer_spatial : str
        Layer of the normalized counts in adata_spatial to use for deconvolution. 
    layer_ref : str
        Layer of the normalized counts in adata_ref to use for deconvolution. 
    dimred : str ["PCA", "LSI"]
        Dimensionality reduction to run and use for leiden clustering. 
    r_lib_path : str
        Path to R library.   
    results_path : str
        Path to save estimated cell type abundances to. 
    
        
    Returns
    --------

    - Saves estimated proportions as csv-file to results_path.
    '''

    adata_spatial_copy = adata_spatial.copy()
    adata_ref_copy = adata_ref.copy()

    if r_lib_path is not None:
        robjects.r.assign("lib_path", r_lib_path)
        robjects.r(".libPaths(lib_path)")

    # Load spatialdwls package
    robjects.r("library(Giotto)")

    # Load reference into R and create signature matrix: 
    scexp_sc = anndata2ri._py2r.py2rpy_anndata(adata_ref_copy)
    robjects.r.assign("scexp_sc", scexp_sc)
    robjects.r.assign("labels_key", labels_key)
    robjects.r.assign("layer_ref", layer_ref) # layer of normalized counts 
    robjects.r.assign("highly_variable_key_ref", highly_variable_key_ref) # key in adata_ref.var with HVFs 
    robjects.r("""
                    counts_sc = assay(scexp_sc,layer_ref)
                    meta_sc = scexp_sc@colData
                    giotto_obj_ref = createGiottoObject(raw_exprs=counts_sc, cell_metadata=meta_sc)
                    giotto_obj_ref@norm_expr = giotto_obj_ref@raw_exprs
                    
                    var_ref = rowData(scexp_sc)
                    highly_variable = rownames(var_ref[var_ref[[highly_variable_key_ref]],])
                    signature_matrix = makeSignMatrixDWLS(giotto_obj_ref,
                                        expression_values = "normalized",
                                        reverse_log = FALSE,
                                        sign_gene = highly_variable,
                                        cell_type_vector = giotto_obj_ref@cell_metadata[[labels_key]])
                """
            ) 
      
    # Preprocess spatial data & compute leiden clustering 
    if highly_variable_key_spatial is not None and highly_variable_key_spatial != "highly_variable":
            adata_spatial_copy.var["highly_variable"] = adata_spatial_copy.var[highly_variable_key_spatial]
    if layer_spatial != "X":
         adata_spatial_copy.X = adata_spatial_copy.layers[layer_spatial]
    if dimred == "PCA":
        sc.pp.pca(adata_spatial_copy)
        sc.pp.neighbors(adata_spatial_copy)
        sc.tl.leiden(adata_spatial_copy)
    elif dimred == "LSI": 
        lsi(adata_spatial_copy)
        sc.pp.neighbors(adata_spatial_copy, use_rep="X_lsi")
        sc.tl.leiden(adata_spatial_copy)

    # Load spatial data into R
    scexp_st = anndata2ri._py2r.py2rpy_anndata(adata_spatial_copy)
    robjects.r.assign("scexp_st", scexp_st)
    robjects.r.assign("layer_spatial", layer_spatial) 
    robjects.r("""
                    counts_st = assay(scexp_st,layer_spatial) 
                    meta_st = scexp_st@colData
                    gene_meta_st = rowData(scexp_st)
                    coords = as.data.frame(scexp_st@int_colData@listData$reducedDims$spatial)
                    colnames(coords) = c("x", "y")
                    rownames(coords) = colnames(counts_st)
            
                    giotto_obj_spatial = createGiottoObject(raw_exprs=counts_st, cell_metadata=meta_st, gene_metadata = gene_meta_st, spatial_locs=coords)
                    giotto_obj_spatial@norm_expr = giotto_obj_spatial@raw_exprs
                    
                """
            )

    # Run deconvolution and save results 
    robjects.r.assign("n_cell", n_cell)
    robjects.r.assign("results_path", results_path)
    robjects.r("""
            giotto_obj_spatial = runDWLSDeconv(giotto_obj_spatial, cluster_column="leiden", sign_matrix=signature_matrix, n_cell=n_cell)
            
            dir.create(results_path)
            write.csv(giotto_obj_spatial@spatial_enrichment$DWLS, paste0(results_path, "/estimated_proportions.csv"))
                    
                """
            )


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

