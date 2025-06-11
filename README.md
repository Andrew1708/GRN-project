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
   snakemake --directory workflow --cores 1 --use-conda {scenicplus_preprocessing_out_dir}/{sample}/cisTarget_db
   ```

### generate_scenic_config
   ```bash
   snakemake --directory workflow --cores 1 --use-conda {scenicplus_preprocessing_out_dir}/{sample}/{nr_cells_per_metacells}
   ```

### scenicplus
   ```bash
   snakemake --directory workflow --cores 1 --use-conda {scenicplus_out_dir}/{sample}_metacell_{nr_cells_per_metacells}/scplusmdata.h5mu
   ```

### benchmark
   ```bash
   snakemake --directory workflow --cores 1 --use-conda {benchmark_out_dir}/{grn_tool}/{sample_detailed}/scplusmdata.h5mu
   ```
   Note: if you run scenicplus {sample_detailed} = {sample}\_metacell\_{nr_cells_per_metacells} (e.g.)

### ⚠️ Note:
   The usage of snakemake allows to use wildcards - variables that are dynamically defined and allocated after user input. For this pipeline the following wildcard were used:
   - {sample} identifier for each example (e.g. AML12_DX, AML12_DX_MO_scDART, etc.)
   - {nr_cells_per_metacells} used in SCENIC+. Defines the number of cells to merge and form a metacell. If "0" is provided it is assumed that the dataset is a multiome, otherwise cells are merged based on a provided metacell_key (usage: snakemake --config metacell_key="XXXX"; Default: "Metacell_Key")

---

## 🧾 Scripts

| Script Name                                | Description                              |
|--------------------------------------------|------------------------------------------|
| `pycistopic_obj_creation.py`               | *This script generates a `pycisTopic` object from fragments, region, and barcode data, adds cell metadata, and serializes the final object using `pickle`.*         |
| `pycistopic_obj_with_model_creation.py`     | *This script loads a serialized `pycisTopic` object, fits multiple LDA topic models using MALLET, selects the best model based on evaluation, attaches it to the object, and saves the result.* |
| `scenicplus_preprocessing.py`     | *Performs necessary preprocessing steps for SCENIC+ pipeline. For scRNA: chooses high variable genes, normalizes and scales data. For scATAC: performs topic binarization and computes differential accessibility regions* |
| `scdart_preprocessing.py`     | *> This script performs preprocessing steps required for diagonal integration of multiome data (scRNA-seq and scATAC-seq). It processes gene expression and chromatin accessibility data to generate **dense matrices** for RNA, ATAC, and a **region-to-gene binary association matrix** using pycisTopic’s gene activity weights. These matrices can be used in downstream integration or modeling tasks such as scDART.* |
| `scdart_-_preprocessing.py`     | *> This script performs preprocessing steps required for diagonal integration of multiome data (scRNA-seq and scATAC-seq). It processes gene expression and chromatin accessibility data to generate **dense matrices** for RNA, ATAC, and a **region-to-gene binary association matrix** using pycisTopic’s gene activity weights. These matrices can be used in downstream integration or modeling tasks such as scDART.* |
| `create_scenic_config.py`     | *> This script generates a YAML configuration file for the SCENIC+ pipeline, defining all input data paths, parameter settings, and output locations for downstream gene regulatory network inference.* |


See more details at [scripts details](docs/scripts_details.md).

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
