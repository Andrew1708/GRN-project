# 📁 Project Name: *Your Project Title Here*

> A Snakemake-based reproducible pipeline for *[brief purpose here, e.g. scRNA + scATAC analysis]*.

---

## 📦 Repository Structure

This project follows the [Snakemake Best Practices](https://snakemake.readthedocs.io/en/stable/project_info/tutorial.html). Below is an overview of the folder structure:

```
├── .gitignore                # Git ignored files
├── README.md                 # This file
├── LICENSE.md                # License file
├── workflow/                 # Main Snakemake pipeline directory
│   ├── rules/                # Modular Snakemake rule files
│   │   ├── module1.smk
│   │   └── module2.smk
│   ├── envs/                 # Conda environment YAMLs
│   │   ├── tool1.yaml
│   │   └── tool2.yaml
│   ├── scripts/              # Python/R scripts used in the workflow
│   │   ├── script1.py
│   │   └── script2.R
│   ├── notebooks/            # Analysis/visualization notebooks
│   │   ├── notebook1.py.ipynb
│   │   └── notebook2.r.ipynb
│   ├── report/               # Snakemake report configuration
│   │   ├── plot1.rst
│   │   └── plot2.rst
│   └── Snakefile             # Entry point to the Snakemake workflow
├── config/                   # Configuration and sample sheets
│   ├── config.yaml
│   └── some-sheet.tsv
├── results/                  # Output directory for final results
└── resources/                # Reference data, indices, etc.
```

---

## 🧬 Snakemake Rules

| Rule Name       | File             | Description                            |
|----------------|------------------|----------------------------------------|
| `create_pycistopic_obj`   | `pycistopic.smk`    | *Create a pycistopic object with cell metadata.*         |
| `create_pycistopic_obj_with_model`   | `pycistopic.smk`    | *Train pycistopic model with different number of topics (n=50) and create a pycistopic with the better model (optimal number of topics).*       |
| `create_pycistopic_obj` | `module1.smk` | Generates pycistopic object from input fragments and metadata. |

## 📌 *How to use each rule?*
### create_pycistopic_obj
   ```bash
   snakemake --directory workflow --cores 1 --use-conda {pycistopic_out_dir}/{sample_wildcard}_pycistopic_obj.pkl
   ```

### create_pycistopic_obj_with_model
   ```bash
   snakemake --directory workflow --cores 1 --use-conda {pycistopic_model_out_dir}/{sample_wildcard}_model_pycistopic.pkl
   ```

### preprocessing
   ```bash
   snakemake --directory workflow --cores 1 --use-conda {scenicplus_preprocessing_out_dir}/{sample_wildcard}
   ```

### prepare_fasta
   ```bash
   snakemake --directory workflow --cores 1 --use-conda {scenicplus_preprocessing_out_dir}/{sample}/consensus_fasta.fa
   ```

### create_cistarget_db
   ```bash
   snakemake --directory workflow --cores 1 --use-conda {scenicplus_preprocessing_out_dir}/{{sample}}/cisTarget_db
   ```

---

## 🧾 Scripts

| Script Name                                | Description                              |
|--------------------------------------------|------------------------------------------|
| `pycistopic-obj-creation.py`               | *This script generates a `pycisTopic` object from fragments, region, and barcode data, adds cell metadata, and serializes the final object using `pickle`.*         |
| `pycistopic-obj-with-model-creation.py`     | *This script loads a serialized `pycisTopic` object, fits multiple LDA topic models using MALLET, selects the best model based on evaluation, attaches it to the object, and saves the result.* |
| `scenicplus-preprocessing.py`     | *Performs necessary preprocessing steps for SCENIC+ pipeline. For scRNA: chooses high variable genes, normalizes and scales data. For scATAC: performs topic binarization and computes differential accessibility regions* |
| `scdart-preprocessing.py`     | *> This script performs preprocessing steps required for diagonal integration of multiome data (scRNA-seq and scATAC-seq). It processes gene expression and chromatin accessibility data to generate **dense matrices** for RNA, ATAC, and a **region-to-gene binary association matrix** using pycisTopic’s gene activity weights. These matrices can be used in downstream integration or modeling tasks such as scDART.* |
| `scdart-preprocessing.py`     | *> This script performs preprocessing steps required for diagonal integration of multiome data (scRNA-seq and scATAC-seq). It processes gene expression and chromatin accessibility data to generate **dense matrices** for RNA, ATAC, and a **region-to-gene binary association matrix** using pycisTopic’s gene activity weights. These matrices can be used in downstream integration or modeling tasks such as scDART.* |


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

## 📓 Notebooks

| Notebook             | Language | Purpose                                          |
|----------------------|----------|--------------------------------------------------|
| `notebook1.py.ipynb` | Python   | *E.g., Data QC and visualization.*               |
| `notebook2.r.ipynb`  | R        | *E.g., Final plots for differential analysis.*   |

> 📊 *Mention what data is loaded, expected input, and what it outputs (figures, tables, etc).*

---

## 🚀 Running the Workflow

1. **Edit `config/config.yaml`** with your paths and parameters.
2. **Create Conda environments** (automatically handled by Snakemake):

   ```bash
   snakemake --directory workflow --cores 1 --use-conda {file_you_wish_to_create}
   ```

---

## 📚 Dependencies & Credits

These scripts rely on the following libraries:

- [pycisTopic](https://pycistopic.readthedocs.io/en/latest/index.html)
- [SCENIC+](https://scenicplus.readthedocs.io/en/latest/index.html)
- 

> **Disclaimer:** The code was primarily developed while working through official tutorials and community guides for the libraries above.

---

## 📄 License

[MIT License](LICENSE.md)
