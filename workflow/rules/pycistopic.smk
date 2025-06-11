rule create_pycistopic_obj:
    input:
        fragments= lambda wildcards: f"{source_atac_dir}/{wildcards.sample}_filtered_fragments.tsv",
        regions= lambda wildcards: f"{source_atac_dir}/{wildcards.sample}_peaks.bed",
        barcodes= lambda wildcards: f"{source_atac_dir}/{wildcards.sample}_filtered_barcodes.tsv",
        metadata= lambda wildcards: f"{source_atac_dir}/{wildcards.sample}_filtered_metadata.tsv"
    output:
        f"{pycistopic_out_dir}/{{sample}}_pycistopic_obj.pkl"
    params:
        temp_dir=pycistopic_temp,
        blacklist=blacklist,
        out_dir=pycistopic_out_dir,
        project=lambda wildcards: wildcards.sample

    conda:
        "../envs/scenicplus.yaml"
    shell:
        """
        python scripts/pycistopic_obj_creation.py \
            --temp_dir {params.temp_dir} \
            --out_dir {params.out_dir} \
            --path_to_fragments {input.fragments} \
            --path_to_regions {input.regions} \
            --path_to_barcodes {input.barcodes} \
            --path_to_blacklist {params.blacklist} \
            --project {params.project} \
            --metadata_path {input.metadata}
        """

rule create_pycistopic_obj_with_model:
    input:
        pkl_obj= lambda wildcards: f"{pycistopic_out_dir}/{wildcards.sample}_pycistopic_obj.pkl"

    output:
        f"{pycistopic_model_out_dir}/{{sample}}_model_pycistopic.pkl"
    
    params:
        temp_dir= lambda wildcards: f"{mallet_temp}/{wildcards.sample}_model",
        out_dir=pycistopic_model_out_dir,
        mallet_path=pycistopic_mallet_path
    
    conda:
        "../envs/scenicplus.yaml"

    shell:
        """
        python scripts/pycistopic_obj_with_model_creation.py \
            --cistopic_path {input.pkl_obj} \
            --temp_dir {params.temp_dir} \
            --out_dir {params.out_dir} \
            --mallet_path {params.mallet_path} \
        """
