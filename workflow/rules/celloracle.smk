rule run_cicero:
    input:
        matrix_file = lambda wildcards: f"{config['scatac_source_dir']}/{wildcards.sample}_filtered_peak_matrix.mtx",
        barcodes_file = lambda wildcards: f"{config['scatac_source_dir']}/{wildcards.sample}_filtered_barcodes.tsv",
        peaks_file = lambda wildcards: f"{config['scatac_source_dir']}/{wildcards.sample}_peaks.bed",
        chrom_sizes_file = lambda wildcards: config['chromsizes'],
    output:
        all_peaks = f"{config['co_preprocessing_out_dir']}/{{sample}}/all_peaks.csv",
        connections = f"{config['co_preprocessing_out_dir']}/{{sample}}/cicero_connections.csv"
    params:
        output_folder = lambda wildcards: f"{config['co_preprocessing_out_dir']}/{wildcards.sample}",
        env = config['cicero_env']
    shell:
        """
        conda run --no-capture-output -n {params.env} Rscript rscripts/run_cicero.R \
            {input.matrix_file} \
            {input.barcodes_file} \
            {input.peaks_file} \
            {input.chrom_sizes_file} \
            {params.output_folder}
        """

rule process_peaks:
    input:
        peaks = lambda wildcards: f"{config['co_preprocessing_out_dir']}/{wildcards.sample}/all_peaks.csv",
        cicero = lambda wildcards: f"{config['co_preprocessing_out_dir']}/{wildcards.sample}/cicero_connections.csv"
    output:
        processed_peaks = f"{config['co_preprocessing_out_dir']}/{{sample}}/processed_peak_file.csv"
    params:
        ref_genome = "hg38"
    conda:
        "../envs/celloracle.yaml"
    shell:
        """
        python scripts/celloracle_process_peak_data.py \
            --peaks {input.peaks} \
            --cicero {input.cicero} \
            --output {output.processed_peaks} \
            --ref-genome {params.ref_genome}
        """


rule scan_motifs:
    input:
        peaks = lambda wildcards: f"{config['co_preprocessing_out_dir']}/{wildcards.sample}/processed_peak_file.csv"
    output:
        grn = f"{config['co_preprocessing_out_dir']}/{{sample}}/base_GRN_dataframe.parquet"
    params:
        ref_genome="hg38",
        fpr=0.02,
        score_threshold=10.0,
        out_dir = lambda wildcards: f"{config['co_preprocessing_out_dir']}/{wildcards.sample}",
    conda:
        "../envs/celloracle.yaml"
    shell:
        """
        python scripts/celloracle_scan_motifs.py \
            --peaks {input.peaks} \
            --output-prefix {params.out_dir} \
            --ref-genome {params.ref_genome} \
            --fpr {params.fpr} \
            --score-threshold {params.score_threshold}
        """

rule co_rna_preprocessing:
    input:
        h5ad= lambda wildcards: f"{config['scrna_source_dir']}/{wildcards.sample}_filtered.h5ad"
    output:
        out= f"{config['co_preprocessing_out_dir']}/{{sample}}/rna_preprocessed.h5ad"
    params:
        n_top_genes=2000
    conda:
        "../envs/celloracle.yaml"
    shell:
        """
        python scripts/celloracle_scrna_preprocessing.py \
            --input {input.h5ad} \
            --output {output.out} \
            --n_top_genes {params.n_top_genes}
        """

# rule run_celloracle:
#     input:
#         h5ad = lambda wildcards: f"{config['co_preprocessing_out_dir']}/{wildcards.sample.replace('_without_atac', '').replace('_with_atac', '')}/rna_preprocessed.h5ad"
#     output:
#         h5mu = f"{config['co_out_dir']}/{{sample}}/mdata.h5mu"
#     params:
#         base_grn_flag = lambda wildcards: (
#             f"--use_atac --base_grn {config['co_preprocessing_out_dir']}/{wildcards.sample.replace('_without_atac', '').replace('_with_atac', '')}/base_GRN_dataframe.parquet"
#             if "without_atac" not in wildcards.sample else ""
#         )
#     conda:
#         "../envs/celloracle.yaml"
#     shell:
#        """
#         python scripts/celloracle_grn_inference.py \
#             --input {input.h5ad} \
#             {params.base_grn_flag} \
#             --output {output.h5mu}
#         """

rule run_celloracle_with_atac:
    wildcard_constraints:
        sample=".*_with_atac"
    input:
        h5ad = lambda wildcards: f"{config['co_preprocessing_out_dir']}/{wildcards.sample.replace('_with_atac', '')}/rna_preprocessed.h5ad",
        base_grn = lambda wildcards: f"{config['co_preprocessing_out_dir']}/{wildcards.sample.replace('_with_atac', '')}/base_GRN_dataframe.parquet",
        atac_peak = lambda wildcards: f"{config['scatac_source_dir']}/{wildcards.sample.replace('_with_atac', '')}_peaks.bed",
        atac_barcodes = lambda wildcards: f"{config['scatac_source_dir']}/{wildcards.sample.replace('_with_atac', '')}_filtered_barcodes.tsv",
        atac_matrix = lambda wildcards: f"{config['scatac_source_dir']}/{wildcards.sample.replace('_with_atac', '')}_filtered_peak_matrix.mtx"
    output:
        h5mu = f"{config['co_out_dir']}/{{sample}}/mdata.h5mu"
    conda:
        "../envs/celloracle.yaml"
    shell:
        """
        python scripts/celloracle_grn_inference.py \
            --input {input.h5ad} \
            --use_atac --base_grn {input.base_grn} \
            --atac_peak {input.atac_peak} \
            --atac_barcodes {input.atac_barcodes} \
            --atac_matrix {input.atac_matrix} \
            --output {output.h5mu}
        """

rule run_celloracle_without_atac:
    wildcard_constraints:
        sample=".*_without_atac"
    input:
        h5ad = lambda wildcards: f"{config['co_preprocessing_out_dir']}/{wildcards.sample.replace('_without_atac', '')}/rna_preprocessed.h5ad"
    output:
        h5mu = f"{config['co_out_dir']}/{{sample}}/mdata.h5mu"
    conda:
        "../envs/celloracle.yaml"
    shell:
        """
        python scripts/celloracle_grn_inference.py \
            --input {input.h5ad} \
            --output {output.h5mu}
        """
