from sacred import Experiment
import seml
import scanpy as sc
from deconvatac.tl import cell2location


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
    def init_dataset(self, adata_spatial_path, adata_reference_path, var_HVF_column, labels_key):

        self.spatial_path = adata_spatial_path
        self.adata_spatial = sc.read_h5ad(adata_spatial_path)
        self.adata_reference = sc.read_h5ad(adata_reference_path)
        # subset on HVFs
        self.adata_spatial = self.adata_spatial[:, self.adata_reference.var[var_HVF_column]]
        self.adata_reference = self.adata_reference[:, self.adata_reference.var[var_HVF_column]]

        self.labels_key = labels_key
        
    
    def init_all(self):
        self.init_dataset()

    @ex.capture(prefix="model")
    def run(self, detection_alpha, N_cells_per_location, output_path, use_gpu):
        
        dataset = self.spatial_path.split("/")[-1].split(".")[0]
        output_path = output_path + dataset
        cell2location(adata_spatial=self.adata_spatial,
                        adata_ref=self.adata_reference,
                        N_cells_per_location=N_cells_per_location,
                        detection_alpha=detection_alpha,
                        labels_key=self.labels_key,
                        use_gpu=use_gpu,
                        results_path=output_path, 
                        plots = False)

        results = {
            "save_path": output_path
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
