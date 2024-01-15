import math
import os
from typing import Optional

import anndata as ad
import muon as mu
import numba as nb
import numpy as np
import pandas as pd
import scipy


def conway_maxwell_poisson(lambda_, nu):
    """Sample from the Conway-Maxwell-Poisson distribution."""
    # Calculate normalizing constant
    lambda_ = int(lambda_)
    nu = float(nu)
    C = np.sum([(pow(lambda_, k) / math.factorial(k)) ** nu for k in range(1000)])

    # Sample from the distribution
    u, sum_p, k = np.random.rand(), 0, 0
    while sum_p < u:
        sum_p += (pow(lambda_, k) / math.factorial(k)) ** nu / C
        k += 1
    return k - 1


class Sampler:
    """Class to sample cells and clusters from a given dataset."""

    def __init__(
        self,
        reference: [mu.MuData, ad.AnnData],
        cell_type_key: str,
        num_spots: int,
        n_regions: int,
        region_type: str = "stripes",
        cell_number_mean: [int, list] = 6,
        cell_number_nu: [float, list] = 20.0,
        cell_type_number: [int, list] = 4,
        balance: Optional[str] = "balanced",
    ):
        self.reference = reference
        self.cell_type_key = cell_type_key
        self.num_spots = num_spots
        self.obs = reference.obs if isinstance(reference, ad.AnnData) else reference[list(reference.mod.keys())[0]].obs

        self.n_regions = n_regions
        self.region_type = region_type

        # Check that cell_number_mean, cell_number_nu and cell_type_number are lists of length n_regions
        if isinstance(cell_number_mean, list):
            if self.region_type == "gradient_number":
                assert len(cell_number_mean) == 2
            else:
                assert len(cell_number_mean) == n_regions
        if isinstance(cell_number_nu, list):
            assert len(cell_number_nu) == n_regions
        if isinstance(cell_type_number, list):
            assert len(cell_type_number) == n_regions

        if isinstance(cell_number_mean, int):
            self.cell_number_mean = np.ones(n_regions, dtype=int) * cell_number_mean
        else:
            self.cell_number_mean = np.array(cell_number_mean)
        if isinstance(cell_number_nu, float):
            self.cell_number_nu = np.ones(n_regions) * cell_number_nu
        else:
            self.cell_number_nu = np.array(cell_number_nu)
        if isinstance(cell_type_number, int):
            self.cell_type_number = np.ones(n_regions, dtype=int) * cell_type_number
        else:
            self.cell_type_number = np.array(cell_type_number)

        if balance not in ["balanced", "unbalanced"]:
            raise ValueError('balance must be one of ["balanced", "unbalanced"].')
        self.init_sample_prob(balance=balance)

    def init_sample_prob(self, balance="unbalanced"):
        """
        Initialize the sample probabilities based on cell type counts.

        Args:
            cell_type_key (str): Key for the cell type column in the dataset.

        Returns
        -------
            None
        """
        cell_counts = self.obs[self.cell_type_key].value_counts(normalize=True)

        if balance == "unbalanced":
            self.cluster_p = cell_counts
            self.cell_p = self.obs[self.cell_type_key].map(self.cluster_p).astype(float)
        elif balance == "balanced":
            self.cluster_p = pd.Series(1 / len(cell_counts), index=cell_counts.index)
            self.cell_p = 1 / self.obs[self.cell_type_key].map(cell_counts).astype(float) / len(cell_counts)

        self.clusters = self.obs[self.cell_type_key].cat.categories
        self.cluster_p = self.cluster_p[self.clusters]

    def define_regions(self, used_clusters):
        """Define the regions to sample from."""
        if self.region_type == "circles":
            self.circle_regions()
        elif self.region_type == "stripes":
            self.stripe_regions()
        elif self.region_type == "gradient_number":
            self.gradient_number_regions()
        elif self.region_type == "gradient_celltype":
            self.gradient_celltype_regions(used_clusters)
        else:
            raise ValueError(
                "region_type must be one of ['stripes', 'circles', 'gradient_number', 'gradient_celltype']"
            )

    def stripe_regions(self):
        """Define regions as stripes."""
        X, Y = self.get_coords()
        # Calculate the height of each stripe
        stripe_height = X.shape[0] / self.n_regions
        # Assign points to stripe regions based on their y-coordinates
        self.regions = np.floor_divide(Y, stripe_height).ravel().astype(int)

    def circle_regions(self):
        """Define regions in circles."""
        X, Y = self.get_coords()
        grid_size = X.shape[0]
        # Define the center of the circles
        circle_center = (grid_size - 1) / 2
        # Calculate the radius of the largest circle to include all points
        largest_circle_radius = min(circle_center, grid_size - 1 - circle_center)

        # Create concentric circles
        circle_radii = np.linspace(0, largest_circle_radius, self.n_regions, endpoint=False)
        regions = np.zeros(X.shape, dtype=int)
        for i in range(self.n_regions - 1):
            mask = np.logical_and(
                np.sqrt((X - circle_center) ** 2 + (Y - circle_center) ** 2) > circle_radii[i],
                np.sqrt((X - circle_center) ** 2 + (Y - circle_center) ** 2) <= circle_radii[i + 1],
            )
            regions[mask] = i

        regions[np.sqrt((X - circle_center) ** 2 + (Y - circle_center) ** 2) > circle_radii[-1]] = self.n_regions - 1
        self.regions = regions.ravel()

    def gradient_number_regions(self):
        """Define regions as a gradient."""
        # one region but we change the cell type probabilities
        self.n_regions = 10
        self.cell_number_mean = np.linspace(self.cell_number_mean[0], self.cell_number_mean[1], self.n_regions).astype(
            int
        )
        self.cell_number_nu = np.ones(self.cell_number_mean.shape) * self.cell_number_nu[0]
        self.stripe_regions()

    def gradient_celltype_regions(self, used_clusters):
        """Define regions as a gradient."""
        # one region but we change the cell type probabilities
        self.n_regions = 10
        self.cell_number_nu = np.ones(self.n_regions) * self.cell_number_nu[0]
        self.cell_number_mean = np.ones(self.n_regions) * self.cell_number_mean[0]
        self.stripe_regions()
        # change the cell type probabilities

        cluster_mask = self.obs[self.cell_type_key].isin([used_clusters[0][1]]).values
        self.cell_p = {region: self.cell_p.copy() for region in range(self.n_regions)}
        for region in range(self.n_regions):
            self.cell_p[region][cluster_mask] = self.cell_p[region][cluster_mask] + (region)

    def sample_data(self):
        """Sample data from the given dataset."""
        used_clusters = {
            region: np.random.choice(self.cluster_p.index, size=(cell_type_number), p=self.cluster_p, replace=False)
            for region, cell_type_number in enumerate(self.cell_type_number)
        }

        self.define_regions(used_clusters)

        # Define cell density (in different regions)
        cell_count = np.array(
            [
                conway_maxwell_poisson(self.cell_number_mean[region], self.cell_number_nu[region])
                for region in self.regions
            ]
        )
        cell_count[cell_count == 0] = 1  # avoid zero cells

        # print which clusters are used in each region and the mean number of cells per cluster
        for i, _ in enumerate(self.cell_type_number):
            print(f"Region {i}: {used_clusters[i]} (mean cells per cluster: {cell_count[self.regions == i].mean()})")

        if "gradient" in self.region_type:
            used_clusters = [used_clusters[0] for _ in self.regions]
        else:
            used_clusters = [used_clusters[region] for region in self.regions]

        params = list(zip(cell_count, used_clusters, self.regions))

        return self.sample_spots(
            params,
        )

    def sample_spots(self, params):
        """
        Sample cells based on given parameters.

        Args:
            params (np.ndarray): Array of cell and cluster counts.

        Returns
        -------
            tuple: Tuple of expression and density arrays.
        """
        sample_exp = {"tmp": self.reference} if isinstance(self.reference, ad.AnnData) else self.reference.mod
        exp = {key: np.zeros((len(params), adata.shape[1])) for key, adata in sample_exp.items()}
        density = np.zeros((len(params), len(self.clusters)))

        for i, (num_cell, used_clusters, region) in enumerate(params):
            cluster_mask = self.obs[self.cell_type_key].isin(used_clusters).values

            if self.region_type == "gradient_celltype":
                p = self.cell_p[region][cluster_mask] / self.cell_p[region][cluster_mask].sum()
            else:
                p = self.cell_p[cluster_mask] / self.cell_p[cluster_mask].sum()
            sampled_cells = np.random.choice(
                self.obs.index[cluster_mask],
                size=num_cell,
                p=p,
            )
            for key, adata in sample_exp.items():
                exp[key][i, :] = adata[sampled_cells, :].X.sum(axis=0)
            density[i, :] = self.obs.loc[sampled_cells, self.cell_type_key].value_counts().loc[self.clusters].values

        return exp, density

    def get_coords(self):
        """Get the coordinates for the spots."""
        # Create a grid
        grid_size = int(np.sqrt(self.num_spots))
        x = np.arange(0, grid_size)
        y = np.arange(0, grid_size)

        # Create a meshgrid for the coordinates
        X, Y = np.meshgrid(x, y)
        return X, Y


def generate_spatial_data(
    reference: [mu.MuData, ad.AnnData],
    cell_type_key: str,
    num_spots: int = 1024,
    n_regions: int = 5,
    balance: Optional[str] = None,
    cell_number_mean: [int, list] = 6,
    cell_number_nu: [float, list] = 20.0,
    cell_type_number: [int, list] = 4,
    **kwargs,
) -> [ad.AnnData, mu.MuData]:
    """
    Generate spatial data.

    Parameters
    ----------
    reference : Union[mu.MuData, ad.AnnData]
        The reference dataset used for generating spatial data.
    cell_type_key : str
        The key in the reference dataset that specifies the cell type information.
    num_spots : int, optional
        The number of spots (locations) to generate, by default 1024.
    n_regions : int, optional
        The number of spatial regions to generate, by default 5.
    balance : str, optional
        The balancing method to use for generating spatial data, by default None.
    cell_number_mean : Union[int, list], optional
        The mean number of cells per spot, by default 6.
    cell_number_nu : Union[float, list], optional
        The dispersion parameter for the negative binomial distribution used to model cell numbers, by default 20.0.
    cell_type_number : Union[int, list], optional
        The number of cell types to generate, by default 4.
    **kwargs
        Additional keyword arguments.

    Returns
    -------
    Union[ad.AnnData, mu.MuData]
        The generated spatial data.

    Notes
    -----
    This function generates spatial data based on a reference dataset. It uses a sampling approach to generate
    synthetic spatial data with specified characteristics such as cell type composition, cell numbers, and spatial
    organization.

    The generated spatial data is returned as an `AnnData` object if the reference dataset is an `AnnData` object,
    or as a `MuData` object if the reference dataset is a `MuData` object.
    """
    np.random.seed(0)
    sampler = Sampler(
        reference=reference,
        cell_type_key=cell_type_key,
        num_spots=num_spots,
        n_regions=n_regions,
        cell_number_mean=cell_number_mean,
        cell_number_nu=cell_number_nu,
        cell_type_number=cell_type_number,
        balance=balance,
        **kwargs,
    )

    exp, density = sampler.sample_data()
    spatial_mod = {}

    X, Y = sampler.get_coords()
    coords = np.vstack((X.flatten(), Y.flatten())).T * 10
    for key, data in exp.items():
        spatial_ann = ad.AnnData(scipy.sparse.csr_matrix(data))
        spatial_ann.var.index = reference.var_names if isinstance(reference, ad.AnnData) else reference[key].var_names
        spatial_ann.obs["cell_count"] = density.sum(axis=1)
        spatial_ann.uns["density"] = pd.DataFrame(density, columns=sampler.clusters, index=spatial_ann.obs_names)
        spatial_ann.obsm["proportions"] = density / density.sum(axis=1)[:, None]
        spatial_ann.uns["proportion_names"] = sampler.clusters.values
        spatial_ann.obsm["spatial"] = coords.astype("int")

        spatial_mod[key] = spatial_ann

    if isinstance(reference, ad.AnnData):
        return spatial_mod["tmp"]
    elif isinstance(reference, mu.MuData):
        return mu.MuData(spatial_mod)
