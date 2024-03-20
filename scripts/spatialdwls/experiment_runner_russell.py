from sacred import Experiment
import seml
import mudata as mu
import scanpy as sc
from deconvatac.tl import spatialdwls


ex = Experiment()
seml.setup_logger(ex)


@ex.post_run_hook
def collect_stats(_run):
    seml.collect_exp_stats(_run)


@ex.config
def config():
    overwrite = None
    db_collection = None
    if db_collection is not None:
        ex.observers.append(seml.create_mongodb_observer(db_collection, overwrite=overwrite))



class ExperimentWrapper:
    """
    A simple wrapper around a sacred experiment, making use of sacred's captured functions with prefixes.
    This allows a modular design of the configuration, where certain sub-dictionaries (e.g., "data") are parsed by
    specific method. This avoids having one large "main" function which takes all parameters as input.
    """

    def __init__(self, init_all=True):
        if init_all:
            self.init_all()

    @ex.capture(prefix="data")
    def init_dataset(self, mdata_spatial_path, mdata_reference_path, var_HVF_column, labels_key, modality):

        self.spatial_path =  mdata_spatial_path
        self.adata_spatial = mu.read_h5mu(mdata_spatial_path).mod[modality]
        self.adata_reference = mu.read_h5mu(mdata_reference_path).mod[modality]
        # subset on HVFs
        self.adata_spatial = self.adata_spatial[:, self.adata_reference.var[var_HVF_column]]
        self.adata_reference = self.adata_reference[:, self.adata_reference.var[var_HVF_column]]

        self.modality =modality
        if self.modality == "rna":
            self.tfidf = False
            self.adata_spatial.X = self.adata_spatial.layers["log_norm"]
            self.adata_reference.X = self.adata_reference.layers["log_norm"]
            self.cluster_key = "leiden_pca"
        elif self.modality == "atac": 
            self.tfidf = True
            self.adata_spatial.X = self.adata_spatial.layers["tfidf_normalized"]
            self.adata_reference.X = self.adata_reference.layers["tfidf_normalized"]
            self.cluster_key = "leiden_lsi"

        self.labels_key = labels_key
        
    @ex.capture(prefix="method")
    def init_method(self, method_id):
        self.method_id =  method_id
    
    def init_all(self):
        self.init_dataset()
        self.init_method()


    @ex.capture(prefix="model")
    def run(self,output_path, n_cell, r_lib_path):
        
        dataset = self.spatial_path.split("/")[-1].split(".")[0]
        output_path = output_path + self.modality + '/' + dataset
        spatialdwls(adata_spatial=self.adata_spatial, 
                    adata_ref=self.adata_reference, 
                    labels_key = self.labels_key,
                    n_cell=n_cell,
                    tfidf=self.tfidf,
                    cluster_key=self.cluster_key,
                    r_lib_path=r_lib_path,
                    results_path=output_path)

        results = {
            "result_path": output_path + "/proportions.csv", 
            "dataset": dataset, 
            "modality": self.modality
        }
        return results


@ex.command(unobserved=True)
def get_experiment(init_all=False):
    print("get_experiment")
    experiment = ExperimentWrapper(init_all=init_all)
    return experiment


@ex.automain
def train(experiment=None):
    if experiment is None:
        experiment = ExperimentWrapper()
    return experiment.run()
