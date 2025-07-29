rule preprocessing:
    input:
        scrna= lambda wildcards: f"{scp_scRNA_input_dir}/{wildcards.sample}_filtered.h5ad",
        scatac= lambda wildcards: f"{scp_scATAC_input_dir}/{wildcards.sample}_model_pycistopic.pkl"

    output:
        region_bed = f"{scp_preprocessing_out_dir}/{{sample}}/regions.bed",
        adata = f"{scp_preprocessing_out_dir}/{{sample}}/rna_adata.h5ad",
        cistopic_object = f"{scp_preprocessing_out_dir}/{{sample}}/cistopic_object.pkl",
        region_set_folder = directory(f"{scp_preprocessing_out_dir}/{{sample}}/region_sets/"),
    
    params:
        temp_dir= lambda wildcards: f"{scp_preprocessing_out_dir}/{wildcards.sample}_temp",
        out_dir = lambda wildcards: f"{scp_preprocessing_out_dir}/{wildcards.sample}"

    conda:
        "../envs/scenicplus.yaml"

    shell:
        """
        python scripts/scenicplus_preprocessing.py \
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
        ctx_db = f"{scp_preprocessing_out_dir}/{{sample}}/cisTarget_db.regions_vs_motifs.rankings.feather",
        dem_db = f"{scp_preprocessing_out_dir}/{{sample}}/cisTarget_db.regions_vs_motifs.scores.feather",

    params:
        outdir = lambda wildcards: directory(f"{scp_preprocessing_out_dir}/{wildcards.sample}/cisTarget_db"),

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
                -o {params.outdir} \
                --bgpadding 1000 \
                -t 20
        """

rule generate_scenic_config:
    input: 
        cisTopic_obj = lambda wildcards: f"{scp_preprocessing_out_dir}/{wildcards.sample}/cistopic_object.pkl",
        scrna_data = lambda wildcards: f"{scp_preprocessing_out_dir}/{wildcards.sample}/rna_adata.h5ad",
        region_set_folder = lambda wildcards: f"{scp_preprocessing_out_dir}/{wildcards.sample}/region_sets/",
        ctx_db = lambda wildcards: f"{scp_preprocessing_out_dir}/{wildcards.sample}/cisTarget_db.regions_vs_motifs.rankings.feather",
        dem_db = lambda wildcards: f"{scp_preprocessing_out_dir}/{wildcards.sample}/cisTarget_db.regions_vs_motifs.scores.feather",

    params:
        motif_annotations = config["motif_collection"],
        scp_out_dir = lambda wildcards: f"{scp_out_dir}/{wildcards.sample}_metacell_{wildcards.nr_cells_per_metacells}",
        temp_dir = lambda wildcards: f"{scp_temp_dir}/{wildcards.sample}",
        nr_cells_per_metacells = lambda wildcards: wildcards.nr_cells_per_metacells,
        metacell_key = config.get("metacell_key", "Metacell_Key"),

    conda:
        "../envs/scenicplus.yaml"

    output:
        config_dir = f"{scp_preprocessing_out_dir}/{{sample}}/{{nr_cells_per_metacells}}config.yaml"
    
    shell:
        r"""
        set -euo pipefail
        python scripts/create_scenic_config.py \
            --cisTopic_obj {input.cisTopic_obj} \
            --scrna_data {input.scrna_data} \
            --region_set_folder {input.region_set_folder} \
            --ctx_db {input.ctx_db} \
            --dem_db {input.dem_db} \
            --motif_annotations {params.motif_annotations} \
            --scp_out_dir {params.scp_out_dir} \
            --out_dir {output.config_dir} \
            --temp_dir {params.temp_dir} \
            --nr_cells_per_metacells {params.nr_cells_per_metacells} \
            --metacell_key {params.metacell_key}
        """

rule scenicplus:
    input:
        config_file = lambda wildcards: f"{scp_preprocessing_out_dir}/{wildcards.sample}/{wildcards.nr_cells_per_metacells}config.yaml",

    output:
        scplus_mudata = f"{scp_out_dir}/{{sample}}_metacell_{{nr_cells_per_metacells}}/scplusmdata.h5mu",
    
    params:
        scp_pipeline_dir = config["scp_pipeline_dir"],
        n_cpu = config.get("n_cpu", 1),
    
    conda:
        "../envs/scenicplus.yaml",

    shell:
        r"""
        set -euo pipefail
        cd {params.scp_pipeline_dir}
        snakemake \
        --latency-wait 30 \
        --cores {params.n_cpu} \
        --rerun-triggers input \
        --configfile {input.config_file}
        """