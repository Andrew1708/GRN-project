import scanpy as sc
from anndata import AnnData

def normalize(adata):
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    return adata

def select_hvgs(adata, n_top_genes=None):
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    adata = adata[:, adata.var.highly_variable].copy()
    return adata

def preprocess(adata, n_top_genes=None):
    adata = normalize(adata)
    adata = select_hvgs(adata, n_top_genes=n_top_genes)
    sc.pp.scale(adata)
    return adata