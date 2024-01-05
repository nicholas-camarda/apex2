#' @Note this function generates the intervals for plotting top x hits, enrichments, etc.
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

#' @Note takes intervals and finds along the intervals vector which ones are less than nrow(mapped_df) and then the interval just above nrow(mapped_df)
#' @note for example, nrow(mapped_df) = 32. if intervals = c(25, 50, 100), then find_good_intervals returns c(25, 50)
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

#' @note this function takes a group size and your mapped_df, and breaks it into chunks
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



#' @note generates mapped_df df, which is a dataframe of gene names mapped to stringdb IDs
#' @note also summarizes my_dat into median log2(FC)
#' @return mapped proteins df
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

#' @note convenience function to save the data that goes into making the stringdb plots
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
