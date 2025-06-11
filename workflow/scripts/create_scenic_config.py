import argparse
import os 
import yaml

N_CPU = 16
#metacell_key = "Metacell_Key"

def parse_args():
    parser = argparse.ArgumentParser(description="Create SCENIC configuration file.")
    
    parser.add_argument("--cisTopic_obj", type=str, required=True, help="Path to the Cistopic object.")
    parser.add_argument("--scrna_data", type=str, required=True, help="Path to the single-cell RNA data (AnnData).")
    parser.add_argument("--region_set_folder", type=str, required=True, help="Folder containing region sets.")
    parser.add_argument("--ctx_db", type=str, required=True, help="Path to the context database.")
    parser.add_argument("--dem_db", type=str, required=True, help="Path to the DEM database.")
    parser.add_argument("--motif_annotations", type=str, required=True, help="Path to motif annotations.")
    parser.add_argument("--scp_out_dir", type=str, required=True, help="Output directory for SCPlus results.")
    parser.add_argument("--out_dir", type=str, required=True, help="Output directory for results.")
    parser.add_argument("--temp_dir", type=str, default="/tmp/scenic_temp", help="Temporary directory for intermediate files.")
    parser.add_argument("--is_multiome", type=str, default="False", choices=["True", "False"], help="Whether the data is multiome or not.")
    parser.add_argument("--nr_cells_per_metacells", type=int, default=2, help="Number of cells per metacell.")
    parser.add_argument("--metacell_key", type=str, default="Metacell_Key", help="Key for metacells in the AnnData object.")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    nr_cells_per_metacells = int(args.nr_cells_per_metacells) 
    metacell_key = args.metacell_key

    is_multiome = True if nr_cells_per_metacells == 0 else False

    config = {
        "input_data": {
            "cisTopic_obj_fname": f"{args.cisTopic_obj}",
            "GEX_anndata_fname": f"{args.scrna_data}",
            "region_set_folder": f"{args.region_set_folder}",
            "ctx_db_fname": f"{args.ctx_db}",
            "dem_db_fname": f"{args.dem_db}",
            "path_to_motif_annotations": f"{args.motif_annotations}"
        },
        "output_data": {
            "combined_GEX_ACC_mudata": f"{args.scp_out_dir}/ACC_GEX.h5mu",
            "dem_result_fname": f"{args.scp_out_dir}/dem_results.hdf5",
            "ctx_result_fname": f"{args.scp_out_dir}/ctx_results.hdf5",
            "output_fname_dem_html": f"{args.scp_out_dir}/dem_results.html",
            "output_fname_ctx_html": f"{args.scp_out_dir}/ctx_results.html",
            "cistromes_direct": f"{args.scp_out_dir}/cistromes_direct.h5ad",
            "cistromes_extended": f"{args.scp_out_dir}/cistromes_extended.h5ad",
            "tf_names": f"{args.scp_out_dir}/tf_names.txt",
            "genome_annotation": f"{args.scp_out_dir}/genome_annotation.tsv",
            "chromsizes": f"{args.scp_out_dir}/chromsizes.tsv",
            "search_space": f"{args.scp_out_dir}/search_space.tsv",
            "tf_to_gene_adjacencies": f"{args.scp_out_dir}/tf_to_gene_adj.tsv",
            "region_to_gene_adjacencies": f"{args.scp_out_dir}/region_to_gene_adj.tsv",
            "eRegulons_direct": f"{args.scp_out_dir}/eRegulon_direct.tsv",
            "eRegulons_extended": f"{args.scp_out_dir}/eRegulon_extended.tsv",
            "AUCell_direct": f"{args.scp_out_dir}/AUCell_direct.h5mu",
            "AUCell_extended": f"{args.scp_out_dir}/AUCell_extended.h5mu",
            "scplus_mdata": f"{args.scp_out_dir}/scplusmdata.h5mu"
        },
        "params_general": {
            "temp_dir": f"{args.temp_dir}",
            "n_cpu": N_CPU,
            "seed": 666
        },
        "params_data_preparation": {
            "bc_transform_func": "\"lambda x: f'{x}'\"",
            "is_multiome": is_multiome,
            "key_to_group_by": metacell_key,
            "nr_cells_per_metacells": nr_cells_per_metacells,
            "direct_annotation": "Direct_annot",
            "extended_annotation": "Orthology_annot",
            "species": "hsapiens",
            "biomart_host": "http://ensembl.org/",
            "search_space_upstream": "1000 150000",
            "search_space_downstream": "1000 150000",
            "search_space_extend_tss": "10 10"
        },
        "params_motif_enrichment": {
            "species": "homo_sapiens",
            "annotation_version": "v10nr_clust",
            "motif_similarity_fdr": 0.001,
            "orthologous_identity_threshold": 0.0,
            "annotations_to_use": "Direct_annot Orthology_annot",
            "fraction_overlap_w_dem_database": 0.4,
            "dem_max_bg_regions": 500,
            "dem_balance_number_of_promoters": True,
            "dem_promoter_space": 1000,
            "dem_adj_pval_thr": 0.05,
            "dem_log2fc_thr": 1.0,
            "dem_mean_fg_thr": 0.0,
            "dem_motif_hit_thr": 3.0,
            "fraction_overlap_w_ctx_database": 0.4,
            "ctx_auc_threshold": 0.005,
            "ctx_nes_threshold": 3.0,
            "ctx_rank_threshold": 0.05
        },
        "params_inference": {
            "tf_to_gene_importance_method": "GBM",
            "region_to_gene_importance_method": "GBM",
            "region_to_gene_correlation_method": "SR",
            "order_regions_to_genes_by": "importance",
            "order_TFs_to_genes_by": "importance",
            "gsea_n_perm": 1000,
            "quantile_thresholds_region_to_gene": "0.85 0.90 0.95",
            "top_n_regionTogenes_per_gene": "5 10 15",
            "top_n_regionTogenes_per_region": "",
            "min_regions_per_gene": 0,
            "rho_threshold": 0.05,
            "min_target_genes": 10
        }
    }

    os.makedirs(os.path.dirname(args.out_dir), exist_ok=True)

    with open(args.out_dir, "w") as f:
        yaml.dump(config, f, sort_keys=False)

