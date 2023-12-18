from typing import Optional

import anndata
import numpy as np
from scipy.sparse import issparse


def reads_to_fragments(
    adata: anndata.AnnData,
    read_layer: Optional[str] = None,
    fragment_layer: str = "fragments",
) -> None:
    """Convert scATAC-seq read counts to appoximate fragment counts.

    Parameters
    ----------
    adata
        AnnData object that contains read counts.
    read_layer
        Key in`.layer` that the read counts are stored in.
    fragment_layer
        Key in`.layer` that the fragment counts will be stored in.

    Returns
    -------
    Adds layer with fragment counts in `.layers[fragment_layer]`.
    """
    adata.layers[fragment_layer] = adata.layers[read_layer].copy() if read_layer else adata.X.copy()
    if issparse(adata.layers[fragment_layer]):
        adata.layers[fragment_layer].data = np.ceil(adata.layers[fragment_layer].data / 2)
    else:
        adata.layers[fragment_layer] = np.ceil(adata.layers[fragment_layer] / 2)
