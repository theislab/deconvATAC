import os
import mudata
import moscot as mt
from moscot.problems.space import MappingProblem
import seaborn as sns
import scanpy as sc
import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt

def moscot_analysis(
    ref_path,
    data_path,
    output_path="./moscot_results",
    rna_key="rna",
    atac_key="atac",
    highly_variable_key="highly_variable",
    highly_accessible_key="highly_accessible",
    obsm_key="X_pca",
    cell_type_key="cell_type",
    return_adatas=False,
    plots=True,
    # Solve parameters
    alpha=0.5,
    epsilon=0.01,
    tau_a=1.0,
    tau_b=1.0,
    rank=-1,
    scale_cost='mean',
    batch_size=None,
    stage=('prepared', 'solved'),
    initializer=None,
    initializer_kwargs={},
    jit=True,
    min_iterations=None,
    max_iterations=None,
    threshold=0.001,
    linear_solver_kwargs={},
    device=None
):
    # Load the data
    ref = mudata.read(ref_path)
    data = mudata.read(data_path)
    
    adata_ref_rna = ref[rna_key]
    adata_ref_atac = ref[atac_key]
    adata_rna = data[rna_key]
    adata_atac = data[atac_key]

    adata_ref_rna_hv = adata_ref_rna[:, adata_ref_rna.var[highly_variable_key]]
    adata_ref_atac_hv = adata_ref_atac[:, adata_ref_atac.var[highly_variable_key]]
    adata_ref_atac_ha = adata_ref_atac[:, adata_ref_atac.var[highly_accessible_key]]

    # Define the problem and solve it  
    mp = MappingProblem(adata_sc=adata_ref_rna_hv, adata_sp=adata_rna)
    mp = mp.prepare(sc_attr={"attr": "obsm", "key": obsm_key})
    mp = mp.solve(
        alpha=alpha,
        epsilon=epsilon,
        tau_a=tau_a,
        tau_b=tau_b,
        rank=rank,
        scale_cost=scale_cost,
        batch_size=batch_size,
        stage=stage,
        initializer=initializer,
        initializer_kwargs=initializer_kwargs,
        jit=jit,
        min_iterations=min_iterations,
        max_iterations=max_iterations,
        threshold=threshold,
        linear_solver_kwargs=linear_solver_kwargs,
        device=device
    )

    transport_matrix = mp.problems[("src", "tgt")].solution.transport_matrix
    df_transport = pd.DataFrame(transport_matrix)

    if not os.path.exists(output_path):
        os.mkdir(output_path)
        
    df_transport.to_csv(os.path.join(output_path + 'predicted_proportions.csv'), index=False)

    deconv_results = df_transport
    cell_types = adata_ref_atac.obs[cell_type_key].values
    deconv_results.columns = cell_types
    deconv_results = deconv_results.div(deconv_results.sum(axis=1), axis=0)
    deconv_results_aggregated = deconv_results.groupby(by=deconv_results.columns, axis=1).sum()
    deconv_results_aggregated.to_csv(os.path.join(output_path + 'normalized_cell_type_predicted_proportions.csv'), index=False)

    if return_adatas:
        return adata_ref, adata

    if plots:
        # Add any plotting code here if needed
        sns.heatmap(transport_matrix)
        plt.show()

# Example usage:
# moscot_analysis(ref_path="/vol/storage/data/simulations/human_developing_cerebral_cortex.h5mu",
#                 data_path="/vol/storage/data/simulations/Brain_1.h5mu")
