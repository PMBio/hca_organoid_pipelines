#!env python
# Author: Jeongbin Park

import anndata
import scanpy as sc
import zarr
import os
import gzip

def run_harmony(ad):
    sc.pp.pca(ad, use_highly_variable=True, n_comps=50)
    sc.external.pp.harmony_integrate(ad, "sample_id", max_iter_harmony=100)
    return ad

def clr_normalize_each_cell(adata, inplace=True):
    """Normalize count vector for each cell, i.e. for each row of .X"""

    import numpy as np
    import scipy

    def seurat_clr(x):
        # TODO: support sparseness
        s = np.sum(np.log1p(x[x > 0]))
        exp = np.exp(s / len(x))
        return np.log1p(x / exp)

    if not inplace:
        adata = adata.copy()

    # apply to dense or sparse matrix, along axis. returns dense matrix
    adata.X = np.apply_along_axis(
        seurat_clr, 1, (adata.X.A if scipy.sparse.issparse(adata.X) else adata.X)
    )
    return adata

def main(argv):
    if len(argv) < 3:
        print("Usage: ./firstline_downstream_analysis.py {sample_name} {analysis_path} {output_directory}")
    sample = argv[1]
    analysis_path = argv[2]
    outdir = argv[3]

    print("Running first-line data analysis...")
    data = sc.read_10x_h5(analysis_path, gex_only=False)
    data.var_names_make_unique()
    data.layers["counts"] = data.X.copy()
    
    protein = data[:, data.var["feature_types"] == "Antibody Capture"].copy()
    rna = data[:, data.var["feature_types"] == "Gene Expression"].copy()

    sc.pp.filter_genes(protein, min_counts=1)
    sc.pp.filter_genes(rna, min_counts=1)
    
    protein.layers["counts"] = protein.X.copy()
    clr_normalize_each_cell(protein)
    sc.pp.log1p(protein)
    sc.pp.neighbors(protein, n_neighbors=30)
    sc.tl.leiden(protein, key_added="protein_leiden", resolution=0.2)

    protein.obsp["protein_connectivities"] = protein.obsp["connectivities"].copy()
    sc.tl.umap(protein)
    
    rna.var["mito"] = rna.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(rna, qc_vars=["mito"], inplace=True)
    rna = rna[rna.obs['pct_counts_mito'] < 20, :]
    
    rna.layers["counts"] = rna.X.copy()
    sc.pp.normalize_total(rna)
    sc.pp.log1p(rna)
    sc.pp.pca(rna)
    #rna = run_harmony(rna)
    #sc.pp.neighbors(rna, n_neighbors=30, use_rep="X_pca_harmony")
    sc.pp.neighbors(rna, n_neighbors=30)
    sc.tl.umap(rna)
    sc.tl.leiden(rna, key_added="rna_leiden")
    
    rna.obsm["protein"] = protein.to_df()
    rna.obsm["protein_umap"] = protein.obsm["X_umap"]
    rna.obs["protein_leiden"] = protein.obs["protein_leiden"]
    rna.obsp["rna_connectivities"] = rna.obsp["connectivities"].copy()
    rna.obsp["protein_connectivities"] = protein.obsp["protein_connectivities"]
    
    sc.tl.umap(rna)
    #sc.pl.umap(rna, color=["rna_leiden", "protein_leiden"], size=100, title=[f"RNA UMAP colored by RNA Leiden clusters ({sample})", f"RNA UMAP colored by Protein Leiden clusters ({sample})"])
    
    protein.write(f"{outdir}/{sample}_protein.h5ad"
    rna.write(f"{outdir}/{sample}_rna.h5ad"

    if os.path.exists(f"{outdir}/samples.json.gz"):
        with gzip.open(f"../public/samples.json.gz", "rt") as f:
            samples = json.loads(f.read())
    else:
        samples = []
    samples.append(sample)
    with gzip.open(f"{outdir}/samples.json.gz", "wt") as f:
        json.dump(samples, f)

    if not os.path.exists(f'{outdir}/{sample}'):
        os.mkdir(f'{outdir}/{sample}')

    ad = rna
    with gzip.open(f"{outdir}/{sample}/metadata.json.gz", "wt") as f:
        obj = {
            'obsm': list(ad.obsm.keys()),
            'obs': list(ad.obs.columns),
            'var': list(ad.var.columns),
        }
        json.dump(obj, f)
    # Zarr.js doesn't support strings
    for name, col in ad.obs.iteritems():
        print(f"Saving obs metadata {name}...")
        with gzip.open(f"{outdir}/{sample}/obs_{name}.json.gz", "wt") as f:
            dump_column(col, f, True)
    for name, col in ad.var.reset_index(level=0).iteritems():
        if name[:9] == 'gene_ids-':
            if name == 'gene_ids-0':
                name = 'gene_ids'
            else:
                continue
        print(f"Saving var metadata {name}...")
        with gzip.open(f"{outdir}/{sample}/var_{name}.json.gz", "wt") as f:
            dump_column(col, f, False)
    print("Creating zarr group...")
    zg = zarr.group(f'{outdir}/{sample}/data.zarr', 'w')
    print("Saving zarr array 'X'...")
    zg.array('X', ad.X, dtype='>f4')
    for obsm_key in ad.obsm.keys():
        print("Saving zarr array '{obsm_key}'...")
        zg.array(obsm_key, ad.obsm[obsm_key], dtype='>f4', chunks=(int(np.ceil(ad.obsm[obsm_key].shape[0]/10)), 1))
    print("Done!")

if __name__ == "__main__":
    import sys
    main(sys.argv)