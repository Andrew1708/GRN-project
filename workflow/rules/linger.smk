rule prepare_linger_inputs:
    input:
        rna_h5ad= lambda wildcards: f"{config['scrna_source_dir']}/{wildcards.sample}_filtered.h5ad",
        atac_mtx= lambda wildcards: f"{config['scatac_source_dir']}/{wildcards.sample}_filtered_peak_matrix.mtx",
        atac_barcodes= lambda wildcards: f"{config['scatac_source_dir']}/{wildcards.sample}_filtered_barcodes.tsv",
        atac_peaks= lambda wildcards: f"{config['scatac_source_dir']}/{wildcards.sample}_peaks.bed"
    output:
        rna=f"{config['linger_preprocessing_out_dir']}/{{sample}}/RNA.txt",
        atac=f"{config['linger_preprocessing_out_dir']}/{{sample}}/ATAC.txt",
        label=f"{config['linger_preprocessing_out_dir']}/{{sample}}/label.txt"
    params:
        label_col= config['celltype_col'],  # change if your label column differs
        outdir= lambda wildcards: f"{config['linger_preprocessing_out_dir']}/{wildcards.sample}",
    conda:
        "../envs/linger.yaml"
    shell:
        """
        python scripts/linger_prepare_input.py \
            --rna_h5ad {input.rna_h5ad} \
            --label_column {params.label_col} \
            --atac_mtx {input.atac_mtx} \
            --atac_barcodes {input.atac_barcodes} \
            --atac_peaks {input.atac_peaks} \
            --output_dir {params.outdir}
        """

rule run_linger:
    input:
        rna= lambda wildcards: f"{config['linger_preprocessing_out_dir']}/{wildcards.sample}/RNA.txt",
        atac= lambda wildcards: f"{config['linger_preprocessing_out_dir']}/{wildcards.sample}/ATAC.txt",
        label= lambda wildcards: f"{config['linger_preprocessing_out_dir']}/{wildcards.sample}/label.txt"
    output:
        tf_re= f"{config['linger_out_dir']}/{{sample}}/tmp/cell_population_TF_RE_binding.txt",
        re_tg= f"{config['linger_out_dir']}/{{sample}}/tmp/cell_population_cis_regulatory.txt",
        tf_tg= f"{config['linger_out_dir']}/{{sample}}/tmp/cell_population_trans_regulatory.txt",
        rna = f"{config['linger_out_dir']}/{{sample}}/tmp/linger_wrk/data/adata_RNA.h5ad",
        atac = f"{config['linger_out_dir']}/{{sample}}/tmp/linger_wrk/data/adata_ATAC.h5ad"
    log:
    params:
        data_dir= lambda wildcards: f"{config['linger_preprocessing_out_dir']}/{wildcards.sample}",  # Directory with RNA.txt, ATAC.txt, label.txt
        datadir=f"{config['linger_data_bulk']}/" ,                      # Placeholder for unused argument
        outdir= lambda wildcards: f"{config['linger_out_dir']}/{wildcards.sample}/tmp/",  # Output directory for results
    conda:
        "../envs/linger.yaml"
    shell:
        """
        python scripts/linger/linger_grn_inference.py \
            --data_dir {params.data_dir} \
            --datadir {params.datadir} \
            --outdir {params.outdir} \
        """

rule create_linger_grn:
    input:
        tf_re=lambda wildcards: f"{config['linger_out_dir']}/{wildcards.sample}/tmp/cell_population_TF_RE_binding.txt",
        re_tg=lambda wildcards: f"{config['linger_out_dir']}/{wildcards.sample}/tmp/cell_population_cis_regulatory.txt",
        tf_tg=lambda wildcards: f"{config['linger_out_dir']}/{wildcards.sample}/tmp/cell_population_trans_regulatory.txt",
        metadata= lambda wildcards: f"{config['scrna_source_dir']}/{wildcards.sample}_filtered.h5ad",
        rna=lambda wildcards: f"{config['linger_out_dir']}/{wildcards.sample}/tmp/linger_wrk/data/adata_RNA.h5ad",
        atac=lambda wildcards: f"{config['linger_out_dir']}/{wildcards.sample}/tmp/linger_wrk/data/adata_ATAC.h5ad"
    output:
        mudata= f"{config['linger_out_dir']}/{{sample}}/mdata.h5mu",
    params:
        outdir= lambda wildcards: f"{config['linger_out_dir']}/{wildcards.sample}",
        batch_size=1500,
        celltype_col=config['celltype_col'],  # Column name for cell types in RNA/ATAC obs
    conda:
        "../envs/linger_merge.yaml"
    shell:
        """
        python scripts/linger/linger_merge_grn.py \
            --tf_re {input.tf_re} \
            --re_tg {input.re_tg} \
            --tf_tg {input.tf_tg} \
            --rna {input.rna} \
            --atac {input.atac} \
            --outdir {params.outdir} \
            --batch_size {params.batch_size} \
            --meta {input.metadata} \
            --celltype_col {params.celltype_col}
        """
