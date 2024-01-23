from typing import Union

import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon


def rmse(true: Union[pd.DataFrame, np.ndarray], predicted: Union[pd.DataFrame, np.ndarray]) -> float:
    """Compute RMSE on true and predicted cell type proportions.

    Parameters
    ----------
    true
        True cell type proportions.
    predicted
        Predicted cell type proportions.

    Returns
    -------
    Root mean squared error.
    """
    rmse = np.sqrt(np.mean((true - predicted) ** 2))
    return rmse


def jsd(true: Union[pd.DataFrame, np.ndarray], predicted: Union[pd.DataFrame, np.ndarray]) -> float:
    """Compute Jensen-Shannon divergence on true and predicted cell type proportions.

    Parameters
    ----------
    true
        True cell type proportions.
    predicted
        Predicted cell type proportions.

    Returns
    -------
    Jensen-Shannon divergence.
    """
    jsd = jensenshannon(true, predicted, axis=1, base=2)
    return jsd.mean()
