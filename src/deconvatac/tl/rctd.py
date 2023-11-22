from rpy2 import robjects
import anndata2ri


def rctd(adata_st, adata_sc, label_key, r_lib_path=None, result_path="./rctd_results"):
    ### TODO: add parameters of the functions create.RCTD and run.RCTD
    ### for now expects raw counts to be saved in adata.X

    if r_lib_path is not None:
        robjects.r.assign("lib_path", r_lib_path)
        robjects.r(".libPaths(lib_path)")

    # load rctd package
    robjects.r("library(spacexr)")

    # load single-cell data into R & create 'Reference' object
    scexp_sc = anndata2ri._py2r.py2rpy_anndata(adata_sc)
    robjects.r.assign("scexp_sc", scexp_sc)
    robjects.r.assign("label_key", label_key)
    robjects.r(
        """
                counts_sc = assay(scexp_sc,"X")
                meta_sc = scexp_sc@colData[[label_key]]
                reference = Reference(counts_sc, meta_sc)
                print("Done with reference")
                """
    )

    # load spatial data into R & create 'SpatialRNA' object
    scexp_st = anndata2ri._py2r.py2rpy_anndata(adata_st)
    robjects.r.assign("scexp_st", scexp_st)
    robjects.r(
        """
                counts_st = assay(scexp_st,"X") 
                coords = as.data.frame(scexp_st@int_colData@listData$reducedDims$spatial)
                colnames(coords) = c("x", "y")
                rownames(coords) = colnames(counts_st)
                puck = SpatialRNA(coords, counts_st)
                print("Done with spatial")
                """
    )

    # run deconvolution and save results
    robjects.r.assign("result_path", result_path)
    robjects.r(
        """
                print("Running RCTD")
                myRCTD = create.RCTD(puck, reference, CELL_MIN_INSTANCE = 0, max_cores = 1)
                myRCTD = run.RCTD(myRCTD, doublet_mode = 'doublet')
                results = myRCTD@results
                print("Saving results")
                dir.create(result_path)
                write.csv(results$results_df, result_path)
                """
    )

    # return adata_st with results saved in adata.obs ?