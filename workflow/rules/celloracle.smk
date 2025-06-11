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
    conda:
        "../envs/cicero.yaml"
    shell:
        """
        Rscript -e "options(repos = c(CRAN='https://cloud.r-project.org'))" \
                -e "source('rscripts/run_cicero.R')" \
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
        "../envs/celloracle_env.yaml"
    shell:
        """
        python scripts/process_peaks.py \
            --peaks {input.peaks} \
            --cicero {input.cicero} \
            --output {output.processed_peaks} \
            --ref-genome {params.ref_genome}
        """



