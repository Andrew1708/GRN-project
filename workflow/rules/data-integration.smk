rule scdart_preprocessing:
    input:
        cistopic = lambda wildcards: f"{pycistopic_model_out_dir}/{wildcards.sample}_model_pycistopic.pkl",
        h5ad = lambda wildcards: f"{source_rna_dir}/{wildcards.sample}_filtered.h5ad",
    
    output:
        region2gene_mat = f"{scdart_temp_dir}/{{sample}}/{{sample}}_region2gene_dense.csv",
        rna_counts = f"{scdart_temp_dir}/{{sample}}/{{sample}}_RNA_counts_dense.csv",
        atac_counts = f"{scdart_temp_dir}/{{sample}}/{{sample}}_ATAC_counts_dense.csv",

    params:
        out_dir = lambda wildcards: f"{scdart_temp_dir}/{wildcards.sample}",
        chromsizes = config["chromsizes"],
        tss = config["tss"],
    
    conda:
        "../envs/scenicplus.yaml",

    shell:
        """
        python scripts/scdart-preprocessing.py \
            --cistopic_path {input.cistopic} \
            --h5ad_path {input.h5ad} \
            --chromsizes_path {params.chromsizes} \
            --tss_path {params.tss} \
            --temp_dir {params.out_dir} \
        """

rule scdart:
    input:
        rna = lambda wildcards: f"{scdart_temp_dir}/{wildcards.sample}/{wildcards.sample}_RNA_counts_dense.csv",
        atac = lambda wildcards: f"{scdart_temp_dir}/{wildcards.sample}/{wildcards.sample}_ATAC_counts_dense.csv",
        reg = lambda wildcards: f"{scdart_temp_dir}/{wildcards.sample}/{wildcards.sample}_region2gene_dense.csv",

    output:
        f"{scdart_temp_dir}/{{sample}}/{{sample}}_matches.csv"
    
    params:
        out_dir = lambda wildcards: f"{scdart_temp_dir}/{wildcards.sample}",
        scdart_dir= config["scdart_dir"],
    
    conda:
        "../envs/scdart.yaml",

    shell:
        """
        python scripts/scdart-integration.py \
            --project_name {wildcards.sample} \
            --rna {input.rna} \
            --atac {input.atac} \
            --reg {input.reg} \
            --out_dir {params.out_dir} \
            --scdart_dir {params.scdart_dir} \
            """