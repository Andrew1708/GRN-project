rule scdart_preprocessing:
    input:
        cistopic = lambda wildcards: f"{pycistopic_model_out_dir}/{wildcards.sample}_model_pycistopic.pkl",
        h5ad = lambda wildcards: f"{source_rna_dir}/{wildcards.sample}_filtered.h5ad",
    
    output:
        region2gene_mat = f"{config['scdart_temp_dir']}/{{sample}}/{{sample}}_region2gene_dense.csv",
        rna_counts = f"{config['scdart_temp_dir']}/{{sample}}/{{sample}}_RNA_counts_dense.csv",
        atac_counts = f"{config['scdart_temp_dir']}/{{sample}}/{{sample}}_ATAC_counts_dense.csv",

    params:
        out_dir = lambda wildcards: f"{config['scdart_temp_dir']}/{wildcards.sample}",
        chromsizes = config["chromsizes"],
        tss = config["tss"],
    
    conda:
        "../envs/scenicplus.yaml",

    shell:
        """
        python scripts/scdart_preprocessing.py \
            --cistopic_path {input.cistopic} \
            --h5ad_path {input.h5ad} \
            --chromsizes_path {params.chromsizes} \
            --tss_path {params.tss} \
            --out_dir {params.out_dir} 
        """

rule scdart:
    input:
        rna = lambda wildcards: f"{config['scdart_temp_dir']}/{wildcards.sample}/{wildcards.sample}_RNA_counts_dense.csv",
        atac = lambda wildcards: f"{config['scdart_temp_dir']}/{wildcards.sample}/{wildcards.sample}_ATAC_counts_dense.csv",
        reg = lambda wildcards: f"{config['scdart_temp_dir']}/{wildcards.sample}/{wildcards.sample}_region2gene_dense.csv",

    output:
        f"{config['scdart_out_dir']}/{{sample}}/matches.csv"
    
    params:
        out_dir = lambda wildcards: f"{config['scdart_out_dir']}/{wildcards.sample}",
        temp_dir = lambda wildcards: f"{config['scdart_temp_dir']}/{wildcards.sample}",
        scdart_dir= config["scdart_dir"],
        cuda = config.get("cuda", "cpu"),
    
    conda:
        "../envs/scdart.yaml",

    shell:
        """
        python scripts/scdart_integration.py \
            --project_name {wildcards.sample} \
            --rna {input.rna} \
            --atac {input.atac} \
            --reg {input.reg} \
            --temp_dir {params.temp_dir} \
            --out_dir {params.out_dir} \
            --cuda {params.cuda} \
            --scdart_dir {params.scdart_dir} \
            """

rule scbridge_preprocessing:
    input:
        cistopic = lambda wildcards: f"{pycistopic_model_out_dir}/{wildcards.sample}_model_pycistopic.pkl",
        h5ad = lambda wildcards: f"{source_rna_dir}/{wildcards.sample}_filtered.h5ad"

    output:
        rna = f"{config['scbridge_temp_dir']}/{{sample}}/rna_filtered.h5ad",
        act_matrix = f"{config['scbridge_temp_dir']}/{{sample}}/gene_act_matrix.h5ad",

    params:
        out_dir = lambda wildcards: f"{config['scbridge_temp_dir']}/{wildcards.sample}",
        celltype_col = config.get("celltype_col", "Classified_Celltype"),
        chromsizes = config["chromsizes"],
        tss = config["tss"],
    
    conda:
        "../envs/scenicplus.yaml",

    shell:
        """
        python scripts/scbridge_preprocessing.py \
            --cistopic_path {input.cistopic} \
            --rna_path {input.h5ad} \
            --chromsizes_path {params.chromsizes} \
            --tss_path {params.tss} \
            --celltype_col {params.celltype_col} \
            --out_dir {params.out_dir} \
        """

rule run_scbridge:
    input:
        rna_file = lambda wildcards: f"{config['scbridge_temp_dir']}/{wildcards.sample}/rna_filtered.h5ad",
        gene_act_matrix = lambda wildcards: f"{config['scbridge_temp_dir']}/{wildcards.sample}/gene_act_matrix.h5ad",
    output:
        matches= f"{config['scbridge_out_dir']}/{{sample}}/matches.csv",
    params:
        main_script = "/home/andrem/GRN-project/resources/scBridge/main.py",
        data_path = lambda wildcards: f"{config['scbridge_temp_dir']}/{wildcards.sample}/",
        gpu_device= config.get("cuda", "cpu"),
        rna_file = lambda wildcards: f"{config['scbridge_temp_dir']}/{wildcards.sample}/rna_filtered-integrated.h5ad",
        act_matrix = lambda wildcards: f"{config['scbridge_temp_dir']}/{wildcards.sample}/gene_act_matrix-integrated.h5ad",
        out_dir = lambda wildcards: f"{config['scbridge_out_dir']}/{wildcards.sample}",
    conda:
        "../envs/scbridge.yaml"
    shell:
        """
        CUDA_VISIBLE_DEVICES={params.gpu_device} \
        python {params.main_script} \
        --data_path {params.data_path} \
        --source_data rna_filtered.h5ad \
        --target_data gene_act_matrix.h5ad \

        python scripts/scbridge_integration.py \
            --rna_path {params.rna_file} \
            --atac_path {params.act_matrix} \
            --out_dir {params.out_dir} \
        """

rule multivi_integration:
    input:
        rna_path= lambda wildcards: f"{config['scrna_source_dir']}/{wildcards.sample}_filtered.h5ad",
        atac_fragments= lambda wildcards: f"{config['scatac_source_dir']}/{wildcards.sample}_filtered_peak_matrix.mtx",
        atac_peaks= lambda wildcards: f"{config['scatac_source_dir']}/{wildcards.sample}_peaks.bed",
        atac_barcodes= lambda wildcards: f"{config['scatac_source_dir']}/{wildcards.sample}_filtered_barcodes.tsv",
    output:
        matches= f"{config['multivi_out_dir']}/{{sample}}/matches.csv",
    params:
        out_dir= lambda wildcards: f"{config['multivi_out_dir']}/{wildcards.sample}",
    conda:
        "../envs/multivi.yaml"
    shell:
        """
        python scripts/multivi_integration.py \
            --rna_path {input.rna_path} \
            --atac_fragments {input.atac_fragments} \
            --atac_peaks {input.atac_peaks} \
            --atac_barcodes {input.atac_barcodes} \
            --out_dir {params.out_dir}
        """



def get_integration_tool_out_dir(tool):
    if tool == "scdart":
        return config["scdart_out_dir"]
    elif tool == "scbridge":
        return config["scbridge_out_dir"]
    else:
        raise ValueError(f"Unknown integration tool: {tool}")
        return

rule create_mo_data:
    input:
        match_file = lambda wildcards: f"{get_integration_tool_out_dir(wildcards.tool)}/{wildcards.sample}/matches.csv",
        rna_file = lambda wildcards: f"{config['scrna_source_dir']}/{wildcards.sample}_filtered.h5ad",
        cistopic_file = lambda wildcards: f"{config['pycistopic_model_out_dir']}/{wildcards.sample}_model_pycistopic.pkl",
    output:
        rna_mo_data= f"{config["scrna_source_dir"]}/{{sample}}_MO_{{tool}}_filtered.h5ad",
        atac_cistopic_mo=  f"{config["pycistopic_out_dir"]}/{{sample}}_MO_{{tool}}_pycistopic_obj.pkl",
        mo_fragment_matrix= f"{config["scatac_source_dir"]}/{{sample}}_MO_{{tool}}_filtered_peak_matrix.mtx",
        mo_bed_file= f"{config["scatac_source_dir"]}/{{sample}}_MO_{{tool}}_peaks.bed",
        mo_barcodes= f"{config["scatac_source_dir"]}/{{sample}}_MO_{{tool}}_filtered_barcodes.tsv",

    params:
         input_rna_dir = config["scrna_source_dir"],
         input_atac_dir = config["scatac_source_dir"],
         cistopic_no_model_dir = config["pycistopic_out_dir"],
         project = lambda wildcards: f"{wildcards.sample}_MO_{wildcards.tool}",
    
    conda:
        "../envs/scenicplus.yaml",
    
    shell:
        """
        python scripts/create_multiome.py \
            --match_file {input.match_file} \
            --rna_file {input.rna_file} \
            --cistopic_file {input.cistopic_file} \
            --input_rna_dir {params.input_rna_dir} \
            --input_atac_dir {params.input_atac_dir} \
            --cistopic_dir {params.cistopic_no_model_dir} \
            --project_name {params.project} \
        """