from anndata import AnnData
import numpy as np



def highly_variable_peaks(
          adata: AnnData, 
          cluster_key: str, 
          layer: str = None,
          scale: float = 1, 
          n_top_features: int = 20000
          ):
    '''
    Selects highly variable features the "var" way
    Adapted from: https://github.com/GreenleafLab/ArchR/blob/c61b0645d1482f80dcc24e25fbd915128c1b2500/R/IterativeLSI.R#L1015

        scatac: anndata object of the (reference) scATAC data
        cluster: name of column in adata.obs containing the clusters 
        scale: scale factor, i.e. log2((sums/feature_sums) *scale + 1)
        n_top_features: how many features to select 
    '''
    # step 1: get matrix of shape (clusters x features) (i.e. sum up features per cluster)
    clusters = adata.obs[cluster_key].unique()
    if layer is not None:
        matrix = adata.layers[layer]
    else:
        matrix = adata.X

    group_matrix = [matrix[adata.obs[cluster_key] == clus].sum(axis=0) for clus in clusters]
    group_matrix = np.asarray(np.stack(group_matrix, axis=0))
    
    # step 2: log normalize 
    # divide each count by sum(cluster) 
    group_matrix = group_matrix/group_matrix.sum(axis=1, keepdims=True)
    # scale & log normalize 
    group_matrix = np.log2(group_matrix * scale + 1)
    
    # step 3: 
    # compute variance of each feature 
    var = np.var(group_matrix, axis=0)
    
    # step 4
    # get indices of the top n features with the highest variance 
    idx = np.argpartition(var,-n_top_features)[-n_top_features:]
    
    # step 5: save HVF to anndata
    hv_features = [adata.var.index[i] for i in idx]
    adata.var["highly_variable"] = adata.var.index.isin(hv_features)
    
    return 