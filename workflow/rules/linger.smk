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
        tfb_potential= f"{config['linger_out_dir']}/{{sample}}/cell_population_TF_RE_binding.txt",
        cis_network= f"{config['linger_out_dir']}/{{sample}}/cell_population_cis_regulatory.txt",
        trans_network= f"{config['linger_out_dir']}/{{sample}}/cell_population_trans_regulatory.txt",
        tfb_potential_all =  f"{config['linger_out_dir']}/{{sample}}/cell_population_TF_RE_binding_all.txt",
        cis_network_all = f"{config['linger_out_dir']}/{{sample}}/cell_type_specific_cis_regulatory_all.txt",
        trans_network_all = f"{config['linger_out_dir']}/{{sample}}/cell_type_specific_trans_regulatory_all.txt"
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
