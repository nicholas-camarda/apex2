#' Generates intervals for plotting top hits, enrichments, etc., in data visualization.
#'
#' This function creates a sequence of intervals based on the row count of the provided dataframe.
#' It ensures that the intervals start from a minimum value, increment by a specified amount,
#' and do not exceed a maximum value.
#'
#' @param df A dataframe for which the intervals are to be generated.
#' @param min_int The minimum interval value (default is 25).
#' @param to_seq_by The increment value for the intervals (default is 50).
#' @param max_int The maximum interval value (default is 300).
#' @return A numeric vector of intervals.
#' @examples
#' df <- data.frame(matrix(runif(1000), ncol = 10))
#' intervals <- generate_intervals(df)
generate_intervals <- function(df, min_int = 25, to_seq_by = 50, max_int = 300) {
    # Determine the first interval value
    first_interval <- min(min_int, nrow(df))

    # Generate the intervals
    intervals <- seq(first_interval, pmin(nrow(df), max_int), by = to_seq_by)

    # Check if the max row number is included, if not, append it to intervals
    if (tail(intervals, n = 1) != pmin(nrow(df), max_int)) {
        intervals <- c(intervals, pmin(nrow(df), max_int))
    }
    return(intervals)
}

#' Identifies suitable intervals from a given set, based on the number of rows in a dataframe.
#'
#' This function examines a vector of intervals and selects those that are less than or equal to the number
#' of rows in the provided dataframe, plus the next higher interval in the vector.
#'
#' @param mapped_df A dataframe whose row count is used to determine the good intervals.
#' @param intervals A numeric vector of predefined intervals.
#' @return A numeric vector of selected intervals.
#' @examples
#' mapped_df <- data.frame(matrix(runif(80), ncol = 10))
#' intervals <- c(25, 50, 100)
#' good_intervals <- find_good_intervals(mapped_df, intervals)
find_good_intervals <- function(mapped_df, intervals) {
    # Find the position of nrow(mapped_df) within the intervals vector
    position <- findInterval(nrow(mapped_df), intervals)

    # Select all intervals up to and including the interval right above nrow(mapped_df)
    updated_intervals <- intervals[1:(position + 1)]

    # Create a logical vector indicating which intervals to do
    which_intervals_to_do <- intervals %in% updated_intervals

    # Finalize the interval numeric vector by subsetting the initial intervals vector with the new intervals to do
    final_intervals <- intervals[which_intervals_to_do]
}

#' Generates quantitative intervals for data segmentation.
#'
#' This function creates intervals for breaking a dataframe into chunks, based on the number of rows.
#' It dynamically adjusts the interval sizes and ensures they do not exceed specified maximum and minimum limits.
#'
#' @param mapped_df A dataframe to be segmented.
#' @param set_max The maximum limit for the intervals (default is 300).
#' @param set_min The minimum limit for the intervals (default is 5).
#' @return A numeric vector of quantitative intervals.
#' @examples
#' mapped_df <- data.frame(matrix(runif(300), ncol = 10))
#' quant_intervals <- generate_quant_intervals(mapped_df)
generate_quant_intervals <- function(mapped_df, set_max = 300, set_min = 5) {
    # # Example usage:
    # # Assuming your data frame is named mapped_df
    # mapped_df <- data.frame(matrix(runif(300 * 10), ncol = 10)) # Example data frame with 300 rows
    # intervals <- generate_quant_intervals(mapped_df)
    # print(intervals)

    n <- nrow(mapped_df)
    max_limit <- min(n, set_max) # Set the maximum limit to 300 or n, whichever is smaller

    # Start with a minimum group size of 10, or 5% of the number of rows, whichever is smaller
    intervals <- min(set_min, 0.05 * n)

    # Set the increment factor based on the size of mapped_df
    # increment_factor <- ifelse(n > 299, 3, 2) # Adjust these values as needed
    increment_factor <- 1.5

    # Generate intervals
    while (tail(intervals, n = 1) < max_limit) {
        next_interval <- tail(intervals, n = 1) * increment_factor
        if (next_interval > max_limit) {
            # If the next interval exceeds the maximum limit, use the maximum limit instead
            next_interval <- max_limit
        }
        intervals <- c(intervals, next_interval)
    }

    # Ensure the maximum size is included
    if (tail(intervals, n = 1) != max_limit) {
        intervals <- c(intervals, max_limit)
    }
    final_intervals <- map_int(unique(intervals), ceiling)
    final_intervals
    return(final_intervals)
}



#' Generates a dataframe of gene names mapped to STRINGdb IDs and summarizes data into median log2(FC).
#'
#' This function processes a given dataset to map gene symbols to STRINGdb identifiers and calculates
#' the median log fold change (log2(FC)) for each gene. It is specifically tailored for proteomic or genomic datasets.
#'
#' @param my_dat The input dataset containing gene information [owner, gene_symbol, adjusted_cmbd_fisher_p_val, log2_fc]
#' @param string_db The STRINGdb database used for mapping gene symbols.
#' @return A dataframe with mapped protein information.
#' @examples
#' my_dat <- read.csv("gene_data.csv")
#' string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 200, input_directory = ".")
#' mapped_proteins_df <- generate_my_mapped_proteins(my_dat, string_db)
generate_my_mapped_proteins <- function(my_dat, string_db) {
    gene_mapping <- my_dat %>%
        # filter(gene_symbol %in% vector_genes) %>% # no need to match to vector genes because Stringdb has it's own database that will do the mapping
        distinct(owner, gene_symbol, adjusted_cmbd_fisher_p_val, log2_fc) %>%
        mutate(og_gene_symbol = gene_symbol) %>%
        # need to do this so that we can identify a gene using all the variants of that gene name available to us
        separate_longer_delim(col = gene_symbol, delim = "_")

    string_db_plus_our_dat <- gene_mapping %>%
        rename(qvalue = adjusted_cmbd_fisher_p_val) %>%
        group_by(og_gene_symbol, qvalue) %>%
        summarize(
            median_logFC = median(log2_fc, na.rm = TRUE),
            gene_lst = list(unique(gene_symbol)), .groups = "drop"
        ) %>%
        unnest(c(gene_lst)) %>%
        distinct(gene_lst, median_logFC, qvalue) %>%
        rename(gene = gene_lst)

    gene_mapping_dist <- gene_mapping %>%
        distinct(gene_symbol, og_gene_symbol) %>%
        mutate(gene_symbol = toupper(gene_symbol)) %>%
        rename(gene = gene_symbol)

    message("Mapping to STRINGdB and sorting from most significant to least...")
    # Map your differentially expressed proteins to STRING and sort them by most significant to least (ascending order)
    mapped_df_temp <- string_db$map(as.data.frame(string_db_plus_our_dat), "gene", removeUnmappedRows = TRUE) %>%
        arrange(qvalue)

    mapped_df <- mapped_df_temp %>%
        left_join(gene_mapping_dist, by = "gene", relationship = "many-to-many") %>% # this is expected, we fix in the next step
        distinct(across(c(4, 5)), .keep_all = TRUE)

    return(mapped_df)
}

#' Saves data into an Excel file for further analysis or reporting.
#'
#' This convenience function takes a given dataset and writes it to an Excel file in a specified directory,
#' naming the file based on provided parameters.
#'
#' @param dat The dataset to be saved.
#' @param output_directory_path The directory path where the file will be saved.
#' @param my_str_name A string to be included in the filename (default is an empty string).
#' @param x_i An identifier to be included in the filename (default is NA).
#' @examples
#' dat <- data.frame(matrix(runif(100), ncol = 10))
#' make_tbl_data(dat, "path/to/output/directory", "dataset_example", 1)
make_tbl_data <- function(dat, output_directory_path, my_str_name = "", x_i = NA) {
    tbl_path <- file.path(output_directory_path, "tables")
    dir.create(tbl_path)
    data_fn <- file.path(tbl_path, qq("@{my_str_name}@{x_i}.xlsx"))
    write.xlsx(dat, data_fn)
}

#' @note generates plots of interactions
#' @param interactions
plot_stringdb_interactions <- function(gene_names, interactions, output_directory_path, string_db, intervals = c(25, 50, 100, 200, 300)) {
    walk(intervals, .f = function(i) {
        message(qq("Top @{i} interactions"))
        top_x_interactions <- interactions[1:i]
        make_tbl_data(dat = tibble(gene = gene_names) %>% slice(1:i), output_directory_path = output_directory_path, my_str_name = "stringdb_interactors_top", x_i = i)
        string_db_interaction_plot <- string_db$plot_network(top_x_interactions)
        # copy the plot from the viewport
        current_plot <- recordPlot()

        interaction_fn <- file.path(
            output_directory_path,
            qq("3E Part 1. Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-top_@{i}_stringdb_pp_interactions.png")
        )
        # Save the plot as a PNG file
        png(interaction_fn, width = 3840, height = 2160, res = 300)
        replayPlot(current_plot)
        dev.off()
    })
}

#' @note perform enrichmetn analysis on the top interactions
generate_interaction_enrichment_tables <- function(interactions, output_directory_path, string_db, intervals = c(25, 50, 100, 200, 300)) {
    walk(intervals, .f = function(i) {
        top_x_interactions <- interactions[1:i]
        message(qq("Enrichment analysis on the top @{i} interactions"))
        interaction_pp_enrichments <- string_db$get_enrichment(top_x_interactions)

        openxlsx::write.xlsx(
            x = interaction_pp_enrichments,
            file = file.path(output_directory_path, qq("3E Part 2. Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GO_and_KEGG-top_stringdb_pp-@{i}.xlsx"))
        )
    })
}


#' @Note from the significant FC and q-value hits plot the interaction network plot the top x
#' @note this function generates its own intervals, separate from the above functions!
#' @return list of top x hits dataframes
generate_significant_hits_interaction_network <- function(mapped_df, output_directory_path, string_db, intervals = c(25, 50, 100, 200, 300)) {
    message("Generating networks of interactions for significant q and FC hits!")

    my_lst_of_sig_hits <- lapply(intervals, FUN = function(x) {
        message(qq("Top @{x} hits"))
        hits <- mapped_df$STRING_id[1:x]
        x_gene <- mapped_df %>%
            distinct(gene) %>%
            slice(1:x)
        make_tbl_data(dat = x_gene, output_directory_path = output_directory_path, my_str_name = "significant_hits_interactors", x_i = x)

        string_db_plot1_idx <- string_db$plot_network(hits)
        # copy the plot from the viewport
        current_plot <- recordPlot()

        mapped_fn <- file.path(output_directory_path, qq("Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-fc_qvalue_top_hits_@{x}.png"))
        # Save the plot as a PNG file
        png(mapped_fn, width = 3840, height = 2160, res = 300)
        replayPlot(current_plot)
        dev.off()

        return(hits)
    }) %>%
        set_names(str_c("h", intervals))

    return(my_lst_of_sig_hits)
}

#' @note perform enrichmemt analysis on the top significnat hits
generate_signif_hits_enrichment_tables <- function(lst_of_sig_hits, output_directory_path, string_db) {
    message("Generating tables of enrichments for significant q and FC hits!")
    walk(seq_len(length(lst_of_sig_hits)), .f = function(i) {
        fn_name <- names(lst_of_sig_hits)[i]
        message(fn_name)
        all_enrichments <- string_db$get_enrichment(lst_of_sig_hits[[i]])

        openxlsx::write.xlsx(
            x = all_enrichments,
            file = file.path(output_directory_path, qq("Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-GO_and_KEGG-sig_fc_hits-@{fn_name}.xlsx"))
        )
    })
}



generate_pert_pathways <- function(mapped_df, lst_of_sig_hits, output_directory_path, string_db) {
    # post payload information to the STRING server
    payload_id <- string_db$post_payload(mapped_df$STRING_id,
        colors = mapped_df$color
    )
    message("Generating graphs with fold changes to visualize up-regulated / down-regulated pathway members for significant q and FC hits...")
    # display a STRING network png with the "halo"
    walk(seq_len(length(lst_of_sig_hits)), .f = function(i) {
        x <- lst_of_sig_hits[[i]]
        fn_name <- names(lst_of_sig_hits)[i]
        x_gene_pert <- mapped_df %>%
            filter(STRING_id %in% lst_of_sig_hits[[i]]) %>%
            distinct(gene)
        make_tbl_data(dat = x_gene_pert, output_directory_path = output_directory_path, my_str_name = "fc_regulated_pathways_top", x_i = i)
        colored_plot_hits <- string_db$plot_network(x, payload_id = payload_id)
        # copy the plot from the viewport
        current_plot <- recordPlot()

        pathways_fn <- file.path(output_directory_path, qq("Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-regulated-pathways_@{fn_name} (down-green up-red).png"))
        # Save the plot as a PNG file
        png(pathways_fn, width = 3840, height = 2160, res = 300)
        replayPlot(current_plot)
        dev.off()
    })
}


plot_top_scoring_pp_interaction_clusters <- function(mapped_df, output_directory_path, string_db) {
    # clusterings
    clustersList <- string_db$get_clusters(mapped_df$STRING_id, algorithm = "fastgreedy")

    # Fastgreedy Algorithm:
    # Approach: This algorithm is a hierarchical agglomerative clustering method. It starts by treating each node as a separate community and ]
    # then merges communities in a way that maximizes modularity (a measure of community structure) in a greedy manner.
    # Advantages:
    # It's efficient and faster compared to some other methods, making it suitable for large networks.
    # It often produces a good community structure with a high modularity score.
    # Disadvantages:
    # It's a deterministic algorithm, so it may get stuck in local optima.
    # It doesn’t allow for overlapping communities.

    # plot first 4 clusters
    # par(mfrow = c(1, 1))
    message("Clustering interactions and returning the top X most significant clusters..")
    for (i in seq_len(length(clustersList))) {
        clusters_fn <- file.path(output_directory_path, qq("3E Supplement. Conserved Data (Rockman-Kruse-Von Zastrow-Rajagopal)-top@{i}_clusters.png"))
        # get the gene names involved
        x_gene_clst <- mapped_df %>%
            filter(STRING_id %in% clustersList[[i]]) %>%
            distinct(gene)
        make_tbl_data(dat = x_gene_clst, output_directory_path = output_directory_path, my_str_name = "cluster_interactors", x_i = i)
        # plot the network using the string_ids
        clsuter_plot <- string_db$plot_network(clustersList[[i]])
        current_plot <- recordPlot()

        png(clusters_fn, width = 3840, height = 2160, res = 300)
        replayPlot(current_plot)
        dev.off()
    }
}

#' Removes figure numbers from a character vector.
#'
#' This function takes a character vector and removes figure number patterns
#' that match the format 'number + letter(s) + optional dash and letters' followed by a period.
#'
#' @param char_vec A character vector from which figure numbers will be removed.
#' @return A character vector with figure numbers removed.
#' @examples
#' char_vec <- c("1A. Figure", "2B-C. Another Figure")
#' remove_fig_number(char_vec)
remove_fig_number <- function(char_vec) {
    gsub(x = char_vec, pattern = "^[0-9]+[A-Z]+-*[A-Z]*\\.\\s+", "")
}

#' Replaces periods with underscores in a character string.
#'
#' This function takes a character string and replaces all periods ('.') with underscores ('_').
#'
#' @param char A character string to be processed.
#' @return A character string with periods replaced by underscores.
#' @examples
#' char <- "example.string"
#' p_friend(char)
p_friend <- function(char) {
    gsub(pattern = "\\.", replacement = "_", x = char)
}

#' Maps gene symbols to Entrez IDs and creates a vector of unique gene symbols.
#'
#' This function utilizes the org.Hs.eg.db package to map provided gene symbols to Entrez IDs.
#' It then creates a vector of unique gene symbols and Entrez IDs for further analysis.
#'
#' @param df A dataframe containing gene_symbol column.
#' @return A list containing two elements: a vector of unique gene symbols and a vector of unique Entrez IDs.
#' @examples
#' df <- data.frame(gene_symbol = c("BRCA1", "TP53"))
#' gene_lists <- make_vector_genes(df)
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

#' Generates enrichment analysis plots and saves results for a given dataset.
#'
#' This function performs GO and KEGG enrichment analyses on the given dataset.
#' It saves the plots, results, and caching data in specified directories.
#'
#' @param data The dataset for enrichment analysis.
#' @param func_base_results_path The base path for saving results.
#' @param pathway_p_adjust_thresh p-value adjustment threshold for pathway analysis.
#' @param go_pathway_fc_threshold Fold change threshold for GO pathway analysis.
#' @return A list of vectors containing gene symbols for all and significant data.
#' @examples
#' data <- read.csv("gene_data.csv")
#' make_figure3(data, "path/to/results")
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

#' Performs STRINGdb analysis on the provided dataset.
#'
#' This function processes the dataset for STRINGdb analysis, including interaction retrieval,
#' enrichment analysis, and plotting interactions and clusters. Results are saved in specified paths.
#'
#' @param data The dataset to be analyzed.
#' @param func_base_results_path The base directory for saving analysis results.
#' @param string_db_q_thresh Q-value threshold for STRINGdb analysis.
#' @param string_db_fc_threshold Fold change threshold for STRINGdb analysis.
#' @param my_input_intervals Intervals for STRINGdb analysis.
#' @return None, but generates and saves analysis results.
#' @examples
#' data <- read.csv("protein_data.csv")
#' do_stringdb_analysis(data)
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
        significant_mapped_proteins_df,
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

        # sort interactions by combined_score to later pluck 1:x top interactions
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
        print(head(x_interactions_df))

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
