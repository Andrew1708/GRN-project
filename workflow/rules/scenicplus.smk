rule preprocessing:
    input:
        scrna= lambda wildcards: f"{scp_scRNA_input_dir}/{wildcards.sample}_filtered.h5ad",
        scatac= lambda wildcards: f"{scp_scATAC_input_dir}/{wildcards.sample}_model_pycistopic.pkl"

    output:
        region_bed = f"{scp_preprocessing_out_dir}/{{sample}}/regions.bed",
        adata = f"{scp_preprocessing_out_dir}/{{sample}}/rna_adata.h5ad",
        #add
    
    params:
        temp_dir= lambda wildcards: f"{scp_preprocessing_out_dir}/{wildcards.sample}_temp",
        out_dir = lambda wildcards: f"{scp_preprocessing_out_dir}/{wildcards.sample}"

    conda:
        "../envs/scenicplus.yaml"

    shell:
        """
        python scripts/scenicplus-preprocessing.py \
            --cistopic_object {input.scatac} \
            --out_dir {params.out_dir} \
            --temp_dir {params.temp_dir} \
            --rna_adata {input.scrna} \
            --sample_name {wildcards.sample} 
        """

rule prepare_fasta:
    input:
        script_dir= create_db_script_dir,
        genome_fasta= f"{pycistarget_dir}/hg38.fa",
        chromsizes= f"{pycistarget_dir}/hg38.chrom.sizes",
        region_bed= lambda wildcards: f"{scp_preprocessing_out_dir}/{wildcards.sample}/regions.bed",
    
    output:
        fasta = f"{scp_preprocessing_out_dir}/{{sample}}/consensus_fasta.fa",
    
    params:
        padding      = 1000,
        shuffle_bg   = "yes",
    
    conda:
        "../envs/scenicplus.yaml"

    shell:
        r"""
        set -euo pipefail
        bash {input.script_dir}/create_fasta_with_padded_bg_from_bed.sh \
             {input.genome_fasta} \
             {input.chromsizes} \
             {input.region_bed} \
             {output.fasta} \
             {params.padding} \
             {params.shuffle_bg}
        """

rule create_cistarget_db:
    input:
        script_dir= create_db_script_dir,
        fasta= lambda wildcards: f"{scp_preprocessing_out_dir}/{wildcards.sample}/consensus_fasta.fa",
        motif_list= f"{pycistarget_dir}/motifs.txt",
        cbuster_dir= f"{pycistarget_dir}/aertslab_motif_colleciton/v10nr_clust_public/singletons",
        cbuster_exec = f"{pycistarget_dir}/cbust"
    
    output:
        outdir = directory(f"{scp_preprocessing_out_dir}/{{sample}}/cisTarget_db/"),


    conda:
        "../envs/scenicplus.yaml"

    shell:
        r"""
        set -euo pipefail
        python {input.script_dir}/create_cistarget_motif_databases.py \
                -f {input.fasta} \
                -M {input.cbuster_dir} \
                -m {input.motif_list} \
                -c {input.cbuster_exec} \
                -o {output.outdir} \
                --bgpadding 1000 \
                -t 20
        """