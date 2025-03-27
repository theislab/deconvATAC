

def rctd(
    adata_spatial, 
    adata_ref, 
    labels_key,
    doublet_mode = 'full',
    r_lib_path=None, 
    results_path="./rctd_results", 
    create_rctd_kwargs = {}
): 
    '''
    Run RCTD

    Parameters 
    -----------

    adata_spatial : AnnData
        AnnData of the spatial data.
    adata_ref : AnnData 
        AnnData of the reference data.
    labels_key : str
        Cell type key in adata_ref.obs for label information
    doublet_mode: str ["doublet", "full"]
        On which mode to run RCTD:  'doublet' (at most 1-2 cell types per pixel),
                                    'full' (no restrictions on number of cell types)
    r_lib_path : str
        Path to R library in which RCTD is installed.   
    results_path : str
        Path to save estimated cell type abundances to. 
    create_rctd_kwargs : dict 
        Parameters for create.RCTD().
        
    Returns
    --------

    - Saves estimated proportions as csv-file to results_path.
    '''
    from rpy2 import robjects
    import anndata2ri

    if r_lib_path is not None:
        robjects.r.assign("lib_path", r_lib_path)
        robjects.r(".libPaths(lib_path)")

    # load rctd package
    robjects.r("library(spacexr)")

    # load single-cell data into R & create 'Reference' object
    scexp_sc = anndata2ri._py2r.py2rpy_anndata(adata_ref)
    robjects.r.assign("scexp_sc", scexp_sc)
    robjects.r.assign("labels_key", labels_key)
    robjects.r(
        """
                counts_sc = assay(scexp_sc,"X")
                meta_sc = scexp_sc@colData[[labels_key]]
                reference = Reference(counts_sc, meta_sc)
                """
    )

    # load spatial data into R & create 'SpatialRNA' object
    scexp_st = anndata2ri._py2r.py2rpy_anndata(adata_spatial)
    robjects.r.assign("scexp_st", scexp_st)
    robjects.r(
        """
                counts_st = assay(scexp_st,"X") 
                coords = as.data.frame(scexp_st@int_colData@listData$reducedDims$spatial)
                colnames(coords) = c("x", "y")
                rownames(coords) = colnames(counts_st)
                puck = SpatialRNA(coords, counts_st)
                """
    )

    # load create_rctd_kwargs into R 
    robjects.r("""
           create_rctd_keys = list()
           create_rctd_values = list()
           """)
    for key, value in create_rctd_kwargs.items(): 
        robjects.r.assign("temp_key", key)
        robjects.r.assign("temp_value", value)
        robjects.r("""
            create_rctd_keys = append(create_rctd_keys, temp_key)
            create_rctd_values = append(create_rctd_values, temp_value)
            """)
    robjects.r("names(create_rctd_values) = create_rctd_keys")

    # run deconvolution and save results
    robjects.r.assign("results_path", results_path)
    robjects.r.assign("doublet_mode", doublet_mode)
    robjects.r(
        """
                myRCTD = do.call(create.RCTD, c(list(spatialRNA=puck, reference=reference, max_cores=1),create_rctd_values))
                myRCTD = run.RCTD(myRCTD, doublet_mode=doublet_mode)
                """
    )
    if doublet_mode == "full":
        robjects.r("""
                    weights = myRCTD@results$weights
                    norm_weights <- normalize_weights(weights)
                    dir.create(results_path)
                    write.csv(as.data.frame(as.matrix(norm_weights)), paste0(results_path, "/estimated_proportions.csv"))
                    """)
    elif doublet_mode == "doublet":
        robjects.r("""
                    results = myRCTD@results
                    dir.create(results_path)
                    write.csv(results$results_df, paste0(results_path, "/estimated_proportions.csv"))
                    """)
