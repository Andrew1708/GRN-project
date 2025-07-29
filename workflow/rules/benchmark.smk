def get_input_dir(grn_tool, sample):
    if grn_tool == "scenicplus":
        return f"{config["input_scenicplus"]}/{sample}/scplusmdata.h5mu"
    elif grn_tool == "celloracle":
        return f"{config["co_out_dir"]}/{sample}/mdata.h5mu"
    elif grn_tool == "linger":
        return f"{config["linger_out_dir"]}/{sample}/mdata.h5mu"
    else:
        raise ValueError(f"Unknown GRN tool: {grn_tool}. Please specify a valid tool.")

rule benchmark:
    input:
        grn_path= lambda wildcards: get_input_dir(wildcards.grn_tool, wildcards.sample)
    output:
        f"{config['benchmark_out_dir']}/{{grn_tool}}/{{sample}}_benchmark_results.csv"
    
    params:
        tfb_golden= config["tfb_golden"],
        prt_golden= config["prt_golden"],
        frc_golden= config["frc_golden"],
        gst_golden= config["gst_golden"],
        tfm_golden= config["tfm_golden"],
        celltype_col = config.get("celltype_col", "Classified_Celltype"),
        output_dir= lambda wildcards: f"{config['benchmark_out_dir']}/{wildcards.grn_tool}",

    
    conda:
        "../envs/benchmark.yaml"
    
    shell:
        r"""
        set -euo pipefail
        python scripts/benchmark.py \
            --grn_path {input.grn_path} \
            --grn_tool {wildcards.grn_tool} \
            --tfb_golden {params.tfb_golden} \
            --prt_golden {params.prt_golden} \
            --frc_golden {params.frc_golden} \
            --gst_golden {params.gst_golden} \
            --tfm_golden {params.tfm_golden} \
            --project_name {wildcards.sample} \
            --celltype_col {params.celltype_col} \
            --output_dir {params.output_dir}
        """


rule fix_benchmark:
    input:
        grn_path = lambda wc: get_input_dir(wc.grn_tool, wc.sample),
        benchmark_table = lambda wc: f"{config['benchmark_out_dir']}/{wc.grn_tool}/{wc.sample}_benchmark_results.csv",
    output:
        corrected_table = f"{config['benchmark_out_dir']}/{{grn_tool}}/{{sample}}_benchmark_results_corrected.csv",
    params:
        tfb_golden = config["tfb_golden"],
        prt_golden = config["prt_golden"],
        frc_golden = config["frc_golden"],
        gst_golden = config["gst_golden"],
        tfm_golden = config["tfm_golden"],
        celltype_col = config.get("celltype_col", "Classified_Celltype"),
        test_list = ["gst"],
    conda:
        "../envs/benchmark.yaml"
    shell:
        r"""
        set -euo pipefail
        python scripts/fix_metric.py \
            --benchmark_table {input.benchmark_table} \
            --tests {params.test_list} \
            --grn_path {input.grn_path} \
            --grn_tool {wildcards.grn_tool} \
            --tfb_golden {params.tfb_golden} \
            --prt_golden {params.prt_golden} \
            --frc_golden {params.frc_golden} \
            --gst_golden {params.gst_golden} \
            --tfm_golden {params.tfm_golden} \
            --project_name {wildcards.sample} \
            --celltype_col {params.celltype_col} \
            --out_file {output.corrected_table}
        """

