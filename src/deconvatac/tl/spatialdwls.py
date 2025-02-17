

def spatialdwls(
    adata_spatial, 
    adata_ref, 
    labels_key,
    cluster_key = "leiden",
    n_cell=50,
    tfidf=True,
    r_lib_path=None, 
    results_path="./spatialdwls_results", 
): 
    
    '''
    Run SpatialDWLS

    Parameters 
    -----------

    adata_spatial : AnnData
        AnnData of the spatial data , filtered by highly variable features. Feature space needs to be the same as the one of adata_ref.
        Normalized counts are expected to be saved in .X. Names of observations and features are expected as row names of adata_spatial.obs and adata_spatial.var. 
    adata_ref : AnnData 
        AnnData of the reference data, filtered by highly variable features. Feature space needs to be the same as the one of adata_spatial.
         Normalized counts are expected to be saved in .X. Names of observations and features are expected as row names of adata_ref.obs and adata_ref.var. 
    labels_key : str
        Cell type key in adata_ref.obs for label information.
    cluster_key : str
        Cluster key in adata_spatial.obs for cluster information.
    n_cell : int
        Number of cells per spot.
    tfidf : bool
        Whether the normalized counts are TFIDF-normalized.
    r_lib_path : str
        Path to R library.  
    results_path : str
        Path to save estimated cell type abundances to. 
    
        
    Returns
    --------

    - Saves estimated proportions as csv-file to results_path.
    '''

    from rpy2 import robjects
    import anndata2ri
    import scanpy as sc
    from deconvatac.tl import lsi
    
    adata_spatial_copy = adata_spatial.copy()
    adata_ref_copy = adata_ref.copy()

    if r_lib_path is not None:
        robjects.r.assign("r_path", r_lib_path)
        robjects.r(".libPaths(r_path)")

    # Load spatialdwls package
    robjects.r("""
               library(Giotto)
               library(data.table)
               """)

    # Load reference into R and create signature matrix: 
    print("Loading reference into R and creating signature matrix")
    scexp_sc = anndata2ri._py2r.py2rpy_anndata(adata_ref_copy)
    robjects.r.assign("scexp_sc", scexp_sc)
    robjects.r.assign("labels_key", labels_key)
    robjects.r.assign("tfidf", tfidf)
    robjects.r("""
                counts_sc = assay(scexp_sc,"X")
                meta_sc = as.data.table(scexp_sc@colData)
                meta_sc$cell_ID = rownames(scexp_sc@colData)
                var_sc = as.data.table(rowData(scexp_sc))
                var_sc$feat_ID = rownames(rowData(scexp_sc))
               
                giotto_obj_ref = createGiottoObject(expression=counts_sc,raw_exprs= counts_sc,cell_metadata=meta_sc, feat_metadata = var_sc)
                norm_expr = create_expr_obj(name = 'normalized',exprMat = counts_sc)
                giotto_obj_ref = set_expression_values(gobject = giotto_obj_ref,values = norm_expr)
                highly_variable = rownames(rowData(scexp_sc))
               
                signature_matrix = makeSignMatrixDWLS(giotto_obj_ref,
                                    expression_values = "normalized",
                                    sign_gene = highly_variable,reverse_log = !tfidf,
                                    cell_type_vector = meta_sc[, get(labels_key)])
                """)
        
    # Load spatial data into R
    print("Loading spatial data into R")
    scexp_st = anndata2ri._py2r.py2rpy_anndata(adata_spatial_copy)
    robjects.r.assign("scexp_st", scexp_st)
    robjects.r("""  
                counts_st = assay(scexp_st,"X") 
                meta_st = as.data.table(scexp_st@colData)
                meta_st$cell_ID = rownames(scexp_st@colData)
                gene_meta_st = as.data.table(rowData(scexp_st))
                gene_meta_st$feat_ID = rownames(rowData(scexp_st))
                coords = as.data.frame(scexp_st@int_colData@listData$reducedDims$spatial)
                colnames(coords) = c("x", "y")
                rownames(coords) = colnames(counts_st)
               
                giotto_obj_spatial = createGiottoObject(expression=counts_st, cell_metadata=meta_st, feat_metadata = gene_meta_st, spatial_locs=coords)
                norm_expr = create_expr_obj(name = 'normalized', exprMat = counts_st)
                giotto_obj_spatial = set_expression_values(gobject = giotto_obj_spatial,values = norm_expr)
                """)

    # Run deconvolution 
    print("Running deconvolution")
    robjects.r.assign("n_cell", n_cell)
    robjects.r.assign("cluster_col", cluster_key)
    robjects.r("giotto_obj_spatial = runDWLSDeconv(giotto_obj_spatial, cluster_column=cluster_col, sign_matrix=signature_matrix, n_cell=n_cell)")

    # Save results 
    print("Saving results")
    robjects.r.assign("results_path", results_path)
    robjects.r("""
                dir.create(results_path)
                results_dt = getSpatialEnrichment(giotto_obj_spatial, output='data.table')
                write.csv(results_dt, paste0(results_path, '/proportions.csv'))
                """)


