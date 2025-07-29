rule aggregate_grns:
    input:
        # Point to the output GRNs from each tool
        scenicplus_dir=config["scp_out_dir"],
        celloracle_dir=config["co_out_dir"],
        linger_dir=config["linger_out_dir"],
    output:
        csv=f"../results/grn_aggregated_top{{n}}.csv"
    params:
        # Filenames inside each tool folder
        scenicplus_filename="scplusmdata.h5mu",
        celloracle_filename="mdata.h5mu",
        linger_filename="mdata.h5mu",
        scenicplus_filter="metacell_15",  # Optional pattern
        celloracle_filter="without_atac",  # Optional pattern
        linger_filter="scdart",      # Optional pattern
        top_n= lambda wildcards: wildcards.n,  # Number of top interactions to extract
    conda:
        "../envs/scenicplus.yaml"  # Optional: conda env with mudata, scanpy, etc.
    shell:
        """
        python scripts/results_get_top_interactions.py \
            --scenicplus_dir {input.scenicplus_dir} \
            --scenicplus_filename {params.scenicplus_filename} \
            --scenicplus_filter {params.scenicplus_filter} \
            --celloracle_dir {input.celloracle_dir} \
            --celloracle_filename {params.celloracle_filename} \
            --celloracle_filter {params.celloracle_filter} \
            --linger_dir {input.linger_dir} \
            --linger_filename {params.linger_filename} \
            --linger_filter {params.linger_filter} \
            --top_n {params.top_n} \
            --output_csv {output.csv}
        """

rule compute_overlap_matrix_by_dir:
    input:
        scenicplus_dir=config["scp_out_dir"],
        celloracle_dir=config["co_out_dir"],
        linger_dir=config["linger_out_dir"],
    output:
        matrix = f"../results/overlap_matrix.csv"
    params:
        scenicplus_filename = "scplusmdata.h5mu",
        celloracle_filename = "mdata.h5mu",
        linger_filename = "mdata.h5mu",
        use_region_flag = "" #--use_region
    conda:
        "../envs/benchmark.yaml"  # Optional: conda env with mudata, scanpy, etc.
    shell:
        """
        python scripts/results_overlap_matrix.py \
            --scenicplus_dir {input.scenicplus_dir} \
            --scenicplus_filename {params.scenicplus_filename} \
            --celloracle_dir {input.celloracle_dir} \
            --celloracle_filename {params.celloracle_filename} \
            --linger_dir {input.linger_dir} \
            --linger_filename {params.linger_filename} \
            {params.use_region_flag} \
            --output_csv {output.matrix}
        """
