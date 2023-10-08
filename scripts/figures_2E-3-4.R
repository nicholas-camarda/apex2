rm(list = ls())
library(ReactomePA)
library(STRINGdb)
library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)
library(readxl)
library(ggpubr)
library(GetoptLong)
library(patchwork)
library(progressr)
library(openxlsx)
library(ggprism)

source(file.path("scripts", "figure4_helper_functions.R"))
handlers(list(
    handler_progress(
        format   = ":spin :current/:total [:bar] :percent in :elapsed ETA: :eta", # | (:message)
        width    = 60,
        complete = "+"
    )
))

remove_fig_number <- function(char_vec) {
    gsub(x = char_vec, pattern = "^[0-9]+[A-Z]+-*[A-Z]*\\.\\s+", "")
}

p_friend <- function(char) {
    gsub(pattern = "\\.", replacement = "_", x = char)
}

make_vector_genes <- function(df) {
    entrez_ids_temp <- mapIds(org.Hs.eg.db,
        keys = df$gene_symbol,
        column = "ENTREZID",
        keytype = "SYMBOL",
        multiVals = "first"
    ) # 'first' will return the first match if multiple matches exist.

    entrez_ids_to_gene_symbol_mapping <- enframe(entrez_ids_temp, name = "gene_symbol", value = "entrez_id") %>%
        unnest(cols = c(entrez_id)) %>%
        filter(!is.na(entrez_id))

    # For use with GO and StringDB
    f_vector_genes <- entrez_ids_to_gene_symbol_mapping$gene_symbol %>% unique()
    # For use with KEGG
    f_entrez_ids <- entrez_ids_to_gene_symbol_mapping$entrez_id %>% unique()
    return(list(f_vector_genes, f_entrez_ids))
}


make_figure3 <- function(data, func_base_results_path, pathway_p_adjust_thresh = 0.1, go_pathway_fc_threshold = 2) {
    fig3_results_path <- file.path(func_base_results_path, "3A-D", "All data")
    fig3_sig_results_path <- file.path(dirname(fig3_results_path), "significant")

    cache_dir_path <- file.path("cache", "enrichment_rds_cache", basename(func_base_results_path), "All data")
    cache_dir_path_sig <- file.path(dirname(cache_dir_path), "significant")

    # all analytes
    all_gene_accession_mapping <- data %>%
        distinct(Accession, gene_symbol) %>%
        # critical step because need to assess all variants of gene name
        separate_longer_delim(gene_symbol, delim = "_")

    lst_res <- make_vector_genes(all_gene_accession_mapping)
    vector_genes <- lst_res[[1]]
    entrez_ids <- lst_res[[2]]

    # do the same but make it on significant analytes
    all_gene_accession_mapping_sig <- data %>%
        group_by(Accession) %>%
        mutate(median_log2_fc = median(log2_fc, na.rm = TRUE)) %>%
        ungroup() %>%
        filter(abs(median_log2_fc) > log2(go_pathway_fc_threshold) & adjusted_cmbd_fisher_p_val < pathway_p_adjust_thresh) %>%
        distinct(Accession, gene_symbol) %>%
        # critical step because need to assess all variants of gene name
        separate_longer_delim(gene_symbol, delim = "_")

    # even if accession and gene_symbol map isn't 1:1, the mapping to entrez_ids
    # inside make_vector_genes() function will ensure uniqueness
    lst_res_sig <- make_vector_genes(all_gene_accession_mapping_sig)
    vector_genes_sig <- lst_res_sig[[1]]
    entrez_ids_sig <- lst_res_sig[[2]]

    # Convert HUGO to Entrez
    message(qq("\nStarting pathway analysis on Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)..."))


    message(qq("Found @{length(vector_genes)} entrez ids mapped to gene symbols!"))
    message(qq("Found @{length(vector_genes_sig)} SIGNIFICANT entrez ids mapped to gene symbols!"))


    if (length(list.files(cache_dir_path_sig, recursive = TRUE, all.files = TRUE)) != 0) {
        message("Already did GO enrichment. Skipping..")
    } else {
        message("Didn't find any cached files... proceeding with enrichment analysis")
        dir.create(fig3_results_path, showWarnings = FALSE, recursive = TRUE)
        dir.create(fig3_sig_results_path, showWarnings = FALSE, recursive = TRUE)
        dir.create(cache_dir_path, showWarnings = FALSE, recursive = TRUE)
        dir.create(cache_dir_path_sig, showWarnings = FALSE, recursive = TRUE)

        my_dotplot_theme <- theme_prism(base_size = 14) +
            theme(
                strip.text.x = element_text(size = rel(1.25), face = "bold"),
                axis.text = element_text(size = rel(0.8)), # , angle = 45, hjust = 1, vjust = 1),
                panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = 0.5),
                panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = 0.25)
                # aspect.ratio = 5 / 4
            )

        enrich_name_lst <- c(
            "3A. GO enrichment Molecular Function", "3B. GO enrichment Cellular Compartment",
            "3C. GO enrichment Biological Process", "3ABC. GO enrichment MF-CC-BP", "3D. KEGG enrichment"
        )
        names(enrich_name_lst) <- c("MF", "CC", "BP", "ALL", "KEGG")
        with_progress({
            p <- progressr::progressor(steps = length(enrich_name_lst))
            enrichment_results <- lapply(seq_len(length(enrich_name_lst)), FUN = function(i) {
                full_name <- enrich_name_lst[i]
                message(full_name)
                if (names(full_name) != "KEGG") {
                    # Perform GO enrichment analysis
                    go_res <- enrichGO(
                        gene = vector_genes,
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        ont = names(full_name), # Change to "CC" for Cellular Component or "MF" for Molecular Function
                        pAdjustMethod = "fdr",
                        qvalueCutoff = pathway_p_adjust_thresh
                    )
                    # Visualize the results
                    num_genes <- length(unique(vector_genes))
                    go_res_plot <- dotplot(go_res, title = qq("@{full_name} for\nConserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)")) +
                        labs(caption = qq("Enrichment q-value cutoff = @{pathway_p_adjust_thresh}\nAdjust method = FDR\nUsing entire conserved gene set regardless of significance\nNum genes = @{num_genes}")) +
                        my_dotplot_theme

                    # Perform GO enrichment analysis on significant Accessions
                    go_res_sig <- enrichGO(
                        gene = vector_genes_sig,
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        ont = names(full_name), # Change to "CC" for Cellular Component or "MF" for Molecular Function
                        pAdjustMethod = "fdr",
                        qvalueCutoff = pathway_p_adjust_thresh
                    )
                    # Visualize the results
                    num_sig_genes <- length(unique(vector_genes_sig))
                    go_res_sig_plot <- dotplot(go_res_sig, title = qq("@{full_name} for\nSignificant Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)")) +
                        labs(caption = qq("Enrichment q-value cutoff = @{pathway_p_adjust_thresh}\nAdjust method = FDR\nUsing significant Q (FDR<@{pathway_p_adjust_thresh}) + FC (abs(FC)>@{go_pathway_fc_threshold})\nNum genes = @{num_sig_genes}")) +
                        my_dotplot_theme
                    p()
                    return(list(go_res, go_res_plot, go_res_sig, go_res_sig_plot) %>% set_names(c("go_res", "go_res_plot", "go_res_sig", "go_res_sig_plot")))
                } else {
                    kegg_res <- enrichKEGG(
                        gene = entrez_ids,
                        pAdjustMethod = "fdr",
                        qvalueCutoff = pathway_p_adjust_thresh,
                        organism = "hsa", # for Homo sapiens. Change according to your species.
                    )

                    # Visualize the results
                    num_entrez <- length(unique(entrez_ids))
                    kegg_res_plot <- dotplot(kegg_res, title = qq("@{full_name} for\nConserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)")) +
                        labs(caption = qq("Enrichment q-value cutoff\nq-value cutoff = @{pathway_p_adjust_thresh}\nAdjust method = FDR\nUsing entire conserved gene set regardless of significance\nNum genes = @{num_entrez}")) +
                        my_dotplot_theme

                    kegg_res_sig <- enrichKEGG(
                        gene = entrez_ids_sig,
                        pAdjustMethod = "fdr",
                        qvalueCutoff = pathway_p_adjust_thresh,
                        organism = "hsa", # for Homo sapiens. Change according to your species.
                    )

                    # Visualize the results
                    num_sig_entrez <- length(unique(entrez_ids_sig))

                    kegg_res_sig_plot <- dotplot(kegg_res, title = qq("@{full_name} for\nSignificant Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)")) +
                        labs(caption = qq("Enrichment q-value cutoff\nq-value cutoff = @{pathway_p_adjust_thresh}\nAdjust method = FDR\nUsing significant Q (FDR<@{pathway_p_adjust_thresh}) + FC (abs(FC)>@{go_pathway_fc_threshold})\nNum genes = @{num_sig_entrez}")) +
                        my_dotplot_theme
                    p()
                    return(list(kegg_res, kegg_res_plot, kegg_res_sig, kegg_res_sig_plot) %>% set_names(c("kegg_res", "kegg_res_plot", "kegg_res_sig", "kegg_res_sig_plot")))
                }
            }) %>%
                set_names(names(enrich_name_lst))
        })


        # str(enrichment_results, max.level = 2)
        go_res_mf_plot <- enrichment_results$MF$go_res_plot
        go_res_cc_plot <- enrichment_results$CC$go_res_plot
        go_res_bp_plot <- enrichment_results$BP$go_res_plot
        go_res_all_plot <- enrichment_results$ALL$go_res_plot
        kegg_res_plot <- enrichment_results$KEGG$kegg_res_plot

        go_res_mf_sig_plot <- enrichment_results$MF$go_res_sig_plot
        go_res_cc_sig_plot <- enrichment_results$CC$go_res_sig_plot
        go_res_bp_sig_plot <- enrichment_results$BP$go_res_sig_plot
        go_res_all_sig_plot <- enrichment_results$ALL$go_res_sig_plot
        kegg_res_sig_plot <- enrichment_results$KEGG$kegg_res_sig_plot

        figs3_a_d <- tribble(
            ~data, ~fn,
            list(go_res_mf_plot), file.path(fig3_results_path, qq("3A. Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GeneOntology_MF.png")),
            list(go_res_cc_plot), file.path(fig3_results_path, qq("3B. Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GeneOntology_CC.png")),
            list(go_res_bp_plot), file.path(fig3_results_path, qq("3C. Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GeneOntology_BP.png")),
            list(go_res_all_plot), file.path(fig3_results_path, qq("3ABC. Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GeneOntology_MF-CC-BP.png")),
            list(kegg_res_plot), file.path(fig3_results_path, qq("3D. Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-KEGG_pathway_enrichment.png")),
            ###
            list(go_res_mf_sig_plot), file.path(fig3_sig_results_path, qq("3A. Significant-Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GeneOntology_MF.png")),
            list(go_res_cc_sig_plot), file.path(fig3_sig_results_path, qq("3B. Significant-Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GeneOntology_CC.png")),
            list(go_res_bp_sig_plot), file.path(fig3_sig_results_path, qq("3C. Significant-Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GeneOntology_BP.png")),
            list(go_res_all_sig_plot), file.path(fig3_sig_results_path, qq("3ABC. Significant-Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GeneOntology_MF-CC-BP.png")),
            list(kegg_res_sig_plot), file.path(fig3_sig_results_path, qq("3D. Significant-Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-KEGG_pathway_enrichment.png"))
        )

        message("Saving figures 3A - D and table for 3ABC")
        walk2(.x = figs3_a_d$data, .y = figs3_a_d$fn, .f = function(plt, fn) {
            ggsave(
                plot = plt[[1]],
                filename = fn,
                width = 10, height = 8
            )
        })

        go_res_mf <- enrichment_results$MF$go_res
        go_res_cc <- enrichment_results$CC$go_res
        go_res_bp <- enrichment_results$BP$go_res
        go_res_all <- enrichment_results$ALL$go_res
        kegg_res <- enrichment_results$KEGG$kegg_res

        go_res_sig_mf <- enrichment_results$MF$go_res_sig
        go_res_sig_cc <- enrichment_results$CC$go_res_sig
        go_res_sig_bp <- enrichment_results$BP$go_res_sig
        go_res_sig_all <- enrichment_results$ALL$go_res_sig
        kegg_sig_res <- enrichment_results$KEGG$kegg_res_sig

        message("Caching...")
        figs3_a_d_cache <- tribble(
            ~data, ~fn_rds,
            list(go_res_mf), file.path(cache_dir_path, qq("3A. Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GeneOntology_MF.rds")),
            list(go_res_cc), file.path(cache_dir_path, qq("3B. Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GeneOntology_CC.rds")),
            list(go_res_bp), file.path(cache_dir_path, qq("3C. Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GeneOntology_BP.rds")),
            list(go_res_all), file.path(cache_dir_path, qq("3ABC. Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GeneOntology_MF-CC-BP.rds")),
            list(kegg_res), file.path(cache_dir_path, qq("3D. Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-KEGG_pathway_enrichment.rds")),
            list(vector_genes_sig), file.path(cache_dir_path, qq("3ABCD. Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-vector_genes.rds")),
            ##
            list(go_res_sig_mf), file.path(cache_dir_path_sig, qq("3A. Significant-Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GeneOntology_MF.rds")),
            list(go_res_sig_cc), file.path(cache_dir_path_sig, qq("3B. Significant-Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GeneOntology_CC.rds")),
            list(go_res_sig_bp), file.path(cache_dir_path_sig, qq("3C. Significant-Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GeneOntology_BP.rds")),
            list(go_res_sig_all), file.path(cache_dir_path_sig, qq("3ABC. Significant-Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GeneOntology_MF-CC-BP.rds")),
            list(kegg_sig_res), file.path(cache_dir_path_sig, qq("3D. Significant-Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-KEGG_pathway_enrichment.rds")),
            list(vector_genes_sig), file.path(cache_dir_path_sig, qq("3ABCD. Significant-Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-vector_genes.rds"))
        )
        walk2(.x = figs3_a_d_cache$data, .y = figs3_a_d_cache$fn_rds, .f = function(obj, fn) {
            write_rds(
                x = obj,
                file = fn
            )
        })

        # save the go_res_all and go_res_sig_all in 3ABC
        write.xlsx(go_res_all,
            file = file.path(fig3_results_path, qq("3ABC. Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GeneOntology_MF-CC-BP_table.xlsx"))
        )
        write.xlsx(go_res_sig_all,
            file = file.path(fig3_sig_results_path, qq("3ABC. Significant-Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GeneOntology_MF-CC-BP_table.xlsx"))
        )
    }

    message("Done!")
    return(list(vector_genes, vector_genes_sig))
}

do_stringdb_analysis <- function(data, func_base_results_path = "~/Downloads", string_db_q_thresh = 0.1, string_db_fc_threshold = 2, my_input_intervals = c(10, 20, 30, 40, 50)) {
    ### STRINGdB Section ###

    # for help with these functions
    # STRINGdb$help("plot_network")
    #  STRINGdb$help("get_png")

    message(qq("Starting StringDB analysis..."))

    dir.create(func_base_results_path, showWarnings = FALSE, recursive = TRUE)
    p <- progressr::progressor(steps = 12)


    message("\nConnecting to STRINGdb for species 9606, version 11.5, score threshold > 200..")
    # Create a STRINGdb object
    # score threshold default = 400 | a lot of people online use 200 threshold
    string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 200)

    # Create the mapped gene set to STRINGdb database
    my_mapped_proteins <- generate_my_mapped_proteins(
        my_dat = data,
        string_db = string_db
    )

    # Create the significant mapped proteins gene set by filtering down full one that was just mapped
    # then add color to it based on FC

    significant_mapped_proteins_df <- my_mapped_proteins %>%
        filter(qvalue < string_db_q_thresh & abs(median_logFC) > log2(string_db_fc_threshold)) %>%
        arrange(desc(median_logFC))
    significant_mapped_proteins_df

    ### STRINGdb$help("add_diff_exp_color")

    my_mapped_proteins_sig <- string_db$add_diff_exp_color(
        significant_mapped_proteins_df, ## we can afford to color these slightly differently than the strong thresholds
        logFcColStr = "median_logFC"
    )
    p()

    mapped_dfs <- list(my_mapped_proteins, my_mapped_proteins_sig)
    names(mapped_dfs) <- c("On All Data", "On Significant FC and Q")
    # walk through these lists
    walk2(mapped_dfs, as.list(names(mapped_dfs)), .f = function(x_mapped_df, subdir_name) {
        message(qq("\nStarting with @{subdir_name}"))
        # Get the interactions
        message("Gathering interactions")
        #     STRINGdb$help("get_interactions")
        x_interactions_df <- string_db$get_interactions(string_ids = x_mapped_df$STRING_id) %>%
            arrange(desc(as.numeric(combined_score))) %>%
            distinct()
        # linearize this df while maintaining the combined_score order
        # if a protein has a higher score in one column, it will always be positioned
        # relative to its highest score
        x_interactions_temp <- vector("character", nrow(x_interactions_df) * 2)
        x_interactions_temp[seq(1, by = 2, length.out = nrow(x_interactions_df))] <- x_interactions_df$from
        x_interactions_temp[seq(2, by = 2, length.out = nrow(x_interactions_df))] <- x_interactions_df$to
        x_interactions <- unique(x_interactions_temp)
        x_gene_names <- x_mapped_df %>%
            filter(STRING_id %in% x_interactions) %>%
            .$gene %>%
            unique()

        # discrepancy sometimes between mapped df and the proteins for which interactions have been found
        no_interactions <- x_mapped_df %>%
            filter(STRING_id %in% setdiff(x_mapped_df$STRING_id, x_interactions)) %>%
            distinct(gene) %>%
            .$gene %>%
            str_c(collapse = ", ")

        if (length(no_interactions) > 0) {
            message(qq("Num mapped proteins = @{nrow(x_mapped_df)}"))
            message(qq("No interactions found for these genes: @{no_interactions}"))
            message(qq("Total num proteins with interactions = @{length(x_interactions)}"))
        }
        p()
        # write_tsv(x_mapped_df, file.path("cache", "DEBUG_mapped_df.tsv"))
        # write_tsv(x_interactions, file.path("cache", "DEBUG_plot_interactions.tsv"))
        # x_interactions_df <- read_tsv(file.path("cache", "DEBUG_plot_interactions.tsv"))
        # x_mapped_df <- read_tsv(file.path("cache", "DEBUG_mapped_df.tsv"))


        # string_db$plot_network(x_interactions[1:25])

        results_3e_interaction_enrichments_subpath <- file.path(
            func_base_results_path,
            "3E",
            "3E Part 2. Enrichment on Significant STRINGdb PP Interactions (Tables)",
            subdir_name
        )
        dir.create(results_3e_interaction_enrichments_subpath, showWarnings = FALSE, recursive = TRUE)
        # make tables for enrichment analysis of top interactions
        # set_intervals <- generate_quant_intervals(tibble(x_interactions))
        set_intervals <- find_good_intervals(tibble(x_interactions), my_input_intervals)

        message("Generating tables of enrichments for top STRINGdb PP interactions")
        generate_interaction_enrichment_tables(
            interactions = x_interactions,
            output_directory_path = results_3e_interaction_enrichments_subpath,
            string_db = string_db,
            intervals = set_intervals
        )
        p()

        results_3e_supp_fig1_subpath <- file.path(
            func_base_results_path,
            "3E",
            "3E Part 1. Top STRINGdb PP Interactions",
            subdir_name
        )
        dir.create(results_3e_supp_fig1_subpath, showWarnings = FALSE, recursive = TRUE)
        # plot the interactions
        message("Plotting network interactions...")
        message("These are the strongest interactions as predicted by STRINGdb..")
        plot_stringdb_interactions(
            gene_names = x_gene_names,
            interactions = x_interactions,
            output_directory_path = results_3e_supp_fig1_subpath,
            string_db = string_db,
            intervals = set_intervals
        )
        p()

        results_3e_supp_fig2_subpath <- file.path(
            func_base_results_path,
            "3E",
            "3E (Supplement Fig 1). Top Scoring PP Interaction Clusters",
            subdir_name
        )
        dir.create(results_3e_supp_fig2_subpath, showWarnings = FALSE, recursive = TRUE)
        plot_top_scoring_pp_interaction_clusters(
            mapped_df = x_mapped_df,
            output_directory_path = results_3e_supp_fig2_subpath,
            string_db = string_db
        )
        p()
    })

    write_tsv(tibble(NOTE = "Plot network takes the TOP X proteins with the highest combined interaction scores and plots them."),
        file = file.path(func_base_results_path, "3E", "README.tsv")
    )
    message("Done with Figure 3E. ")

    # make figure 4a and init subdir
    results_4a_pp_interact_sig_hits_subpath <- file.path(
        func_base_results_path,
        "4A",
        "PP Interaction Network of (ONLY) Significant FC Hits"
    )
    dir.create(results_4a_pp_interact_sig_hits_subpath, showWarnings = FALSE, recursive = TRUE)

    message("Saving the significant proteins dataframe for later inspection if desired...")
    write.xlsx(my_mapped_proteins_sig,
        file = file.path(
            dirname(results_4a_pp_interact_sig_hits_subpath),
            "4A. Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-mapped_significant_proteins.xlsx"
        )
    )
    message("Done!")

    # use a logical number of intervals
    set_intervals <- find_good_intervals(my_mapped_proteins_sig, my_input_intervals)
    # generate networks for significant fc and q hits
    lst_of_sig_hits <- generate_significant_hits_interaction_network(
        mapped_df = my_mapped_proteins_sig,
        output_directory_path = results_4a_pp_interact_sig_hits_subpath,
        string_db = string_db,
        intervals = set_intervals
    )

    p()

    # GO on hits
    results_4a_enrich_sig_hits_subpath <- file.path(
        func_base_results_path,
        "4A",
        "Enrichment on Significant FC Hits"
    )
    dir.create(results_4a_enrich_sig_hits_subpath, showWarnings = FALSE, recursive = TRUE)
    generate_signif_hits_enrichment_tables(
        lst_of_sig_hits = lst_of_sig_hits,
        output_directory_path = results_4a_enrich_sig_hits_subpath,
        string_db = string_db
    )

    p()

    results_4a_perturbation_sig_hits_subpath <- file.path(
        func_base_results_path,
        "4A",
        "Regulated Pathways Perturbed in PP Interaction Networks of Significant FC Hits"
    )
    dir.create(results_4a_perturbation_sig_hits_subpath, showWarnings = FALSE, recursive = TRUE)
    generate_pert_pathways(
        mapped_df = my_mapped_proteins_sig,
        lst_of_sig_hits = lst_of_sig_hits,
        output_directory_path = results_4a_perturbation_sig_hits_subpath,
        string_db = string_db
    )
    p()

    message("Done!\n")
}

# My variables
all_data <- read_rds(file.path("cache", "FINAL-cache-including_kruse.rds"))
conserved_data <- all_data %>% filter(n_owners == 4)


# make figure 3
my_vector_genes_lst <- make_figure3(
    data = conserved_data,
    func_base_results_path = file.path("results", "plots", "Figures 3 and 4", qq("Strict Thresholds for Fig 3 and 4 (q0_1 and FC2)")),
    pathway_p_adjust_thresh = 0.1,
    go_pathway_fc_threshold = 2
)
my_vector_genes <- my_vector_genes_lst[[1]]

# make figures 3E and 4
with_progress(do_stringdb_analysis(
    data = conserved_data,
    func_base_results_path = file.path("results", "plots", "Figures 3 and 4", qq("Strict Thresholds for Fig 3 and 4 (q0_1 and FC2)")),
    string_db_q_thresh = 0.1,
    string_db_fc_threshold = 2,
    my_input_intervals = c(10, 20, 30, 40, 50)
))


# make figure 3
my_vector_genes_lst <- make_figure3(
    data = conserved_data,
    func_base_results_path = file.path("results", "plots", "Figures 3 and 4", qq("Relaxed Thresholds for Fig 3 and 4 (q0_2 and FC1_33)")),
    pathway_p_adjust_thresh = 0.2,
    go_pathway_fc_threshold = 1.33
)

# make figures 3E and 4
with_progress(do_stringdb_analysis(
    data = conserved_data,
    func_base_results_path = file.path("results", "plots", "Figures 3 and 4", qq("Relaxed Thresholds for Fig 3 and 4 (q0_2 and FC1_33)")),
    string_db_q_thresh = 0.2,
    string_db_fc_threshold = 1.33,
    my_input_intervals = c(10, 20, 30, 40, 50)
))




# DEBUG:
# set_and_print_variables <- function() {
#     verbose <<- my_verbosity
#     cat("verbose:", verbose, "\n")

#     # Assuming my_combined_dat is defined elsewhere
#     combined_dat <<- my_combined_dat
#     cat("combined_dat:", dim(combined_dat), "\n")

#     results_path <<- file.path("results", "plots", "Figure 2", "2B-E", "Including Kruse")
#     cat("results_path:", results_path, "\n")

#     # cache_name_helper <- "complete"
#     # cat("cache_name_helper:", cache_name_helper, "\n")
#     fc_cutoff <<- my_fc_cutoff
#     cat("fc_cutoff:", fc_cutoff, "\n")

#     stat_value_type <<- my_stat_value_type
#     cat("stat_value_type:", stat_value_type, "\n")

#     stat_cutoff <<- my_p_adjust_cutoff
#     cat("stat_cutoff:", stat_cutoff, "\n")
# }

# # Call the function to set and print the variables
# set_and_print_variables()
# set_and_print_variables()
