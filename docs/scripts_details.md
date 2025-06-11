# Scripts Details
 
## `create_cistopic_object.py`

### Parameters

| Name | Description |
|------|-------------|
| `--path_to_fragments` | Path to a `.tsv` file containing fragment information (Chromosome, Start, End, Name, Strand). |
| `--path_to_regions` | Path to the region file defining accessible regions (used by pycisTopic). |
| `--path_to_barcodes` | Path to a file listing valid barcodes (cell identifiers). |
| `--path_to_blacklist` | Path to a blacklist file of genomic regions to exclude. |
| `--metadata_path` | Path to a `.tsv` metadata file containing cell-level annotations. Must include `Cell_Barcode` and `Library_ID` columns. |
| `--project` | Name of the project — used in naming the resulting object file. |
| `--temp_dir` | Temporary working directory (can be used for intermediate steps). |
| `--out_dir` | Output directory where the final `.pkl` file will be saved. |

---

### Outputs

| Name / Destination | Type | Description |
|--------------------|------|-------------|
| `<out_dir>/<project>_pycistopic_obj.pkl` | Pickled Python object | Serialized pycisTopic object enriched with metadata, ready for downstream use. |

---

## `pycistopic-obj-with-model-creation.py`

### Parameters

| Name | Description |
|------|-------------|
| `--cistopic_path` | Path to the existing `.pkl` file containing a `pycisTopic` object (without an LDA model). |
| `--mallet_path` | Path to the `mallet` binary (e.g., `/path/to/mallet/bin/mallet`). |
| `--temp_dir` | Temporary directory for intermediate files produced by MALLET. |
| `--out_dir` | Directory to save the updated `.pkl` file with the best LDA model. |

---

### Outputs

| Name / Destination | Type | Description |
|--------------------|------|-------------|
| `<out_dir>/<project>_model_pycistopic.pkl` | Pickled Python object | The original `pycisTopic` object now enriched with a selected LDA model. |
| `stdout` | Console log | Prints the filename of the saved model-enhanced object, e.g., `Saved: myproject_model_pycistopic.pkl`. |

---

## `scenicplus-preprocessing.py`
### Parameters

| Name | Description |
|------|-------------|
| `--cistopic_object` | Path to a `.pkl` file containing the pycisTopic object with an LDA model. |
| `--out_dir` | Output directory for results (will contain region sets and preprocessed RNA data). |
| `--temp_dir` | Temporary directory for intermediate files during DAR computation. |
| `--rna_adata` | Path to the `.h5ad` file containing RNA-seq data in AnnData format. |
| `--sample_name` | Sample name used to tag output files. |

---

### Outputs

| Name / Destination | Type | Description |
|--------------------|------|-------------|
| `<out_dir>/regions.bed` | `.bed` | BED-format file listing all regions from the input `cistopic_object`. |
| `<out_dir>/region_sets/Topics_otsu/*.bed` | `.bed` | Topic regions binarized using Otsu's method. |
| `<out_dir>/region_sets/Topics_top_3k/*.bed` | `.bed` | Top 3000 regions for each topic using ntop binarization. |
| `<out_dir>/region_sets/DARs_cell_type/*.bed` | `.bed` | Differentially accessible regions for each cell type (`Classified_Celltype`). Empty files are automatically deleted. |
| `<out_dir>/rna_adata.h5ad` | `.h5ad` | Preprocessed RNA AnnData file with PCA, neighbors, and UMAP computed. |


# `scdart-preprocessing.py`

---

## Parameters

| Name | Description |
|------|-------------|
| `--cistopic_path` | Path to the `.pkl` file containing a `pycisTopic` object with fragment matrix and region names. |
| `--h5ad_path` | Path to the `.h5ad` file with the RNA data (AnnData format). |
| `--chromsizes_path` | Path to the chromosome sizes `.tsv` file (e.g., UCSC style). |
| `--tss_path` | Path to the gene annotation `.bed` file containing TSS positions. |
| `--out_dir` | Directory where the output matrices will be saved. |

---

## Outputs

| Name / Destination | Type | Description |
|--------------------|------|-------------|
| `<out_dir>/<project>_region2gene_dense.csv` | CSV | Binary matrix linking accessible regions to genes using filtered gene activity weights. |
| `<out_dir>/<project>_RNA_counts_dense.csv` | CSV | Dense expression matrix for selected genes across cells. |
| `<out_dir>/<project>_ATAC_counts_dense.csv` | CSV | Dense accessibility matrix for selected regions across cells. |
| `stdout` | Console log | Prints confirmation when matrices are saved successfully. |
---

### `create_scenic_config.py`

## Parameters

| Argument | Description |
|----------|-------------|
| `--cisTopic_obj` | Path to the `.pkl` file containing the cisTopic object. |
| `--scrna_data` | Path to the `.h5ad` file containing preprocessed scRNA-seq data. |
| `--region_set_folder` | Directory containing genomic region sets (used by pycisTarget). |
| `--ctx_db` | Path to the context-based motif enrichment database `.feather` file. |
| `--dem_db` | Path to the differential enrichment motif database `.feather` file. |
| `--motif_annotations` | Path to the motif annotation file (e.g., `.tbl` from AertsLab collection). |
| `--scp_out_dir` | Directory where SCENIC+ outputs will be written. |
| `--out_dir` | Output path to write the final `config.yaml` file. |
| `--temp_dir` | Directory for temporary intermediate files (default: `/tmp/scenic_temp`). |
| `--is_multiome` | Whether the dataset is multiomic (contains both scRNA and scATAC). Accepts `True` or `False`. |
| `--nr_cells_per_metacells` | Number of cells per metacell. If 0, no aggregation is applied (used in multiome mode). |
| `--metacell_key` | AnnData key under which metacell groupings are stored. Used for aggregation and downstream analysis. |

---

## Output

| File | Type | Description |
|------|------|-------------|
| `<out_dir>` | YAML file | SCENIC+ configuration file, defining all input paths, output destinations, and key parameter blocks: general settings, data preparation, motif enrichment, and inference configuration. |
