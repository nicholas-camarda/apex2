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

handlers(list(
    handler_progress(
        format   = ":spin :current/:total [:bar] :percent in :elapsed ETA: :eta", # | (:message)
        width    = 60,
        complete = "+"
    )
))

all_data <- read_rds(file.path("cache", "FINAL-all_data.rds"))
data_no_raja <- read_rds(file.path("cache", "FINAL-data_no_raja.rds"))

pathway_p_adjust_thresh <- 0.2
string_db_p_value_threshold <- 0.01
string_db_fc_threshold <- 0.5

base_results_path <- file.path("results", "plots", "Figures 3 and 4")
dir.create(base_results_path, showWarnings = FALSE, recursive = TRUE)


do_pathway_analysis <- function(data, data_type) {
    vector_genes <- unique(data$gene_symbol)
    # default


    results_path <- file.path(base_results_path, data_type)
    dir.create(results_path, showWarnings = FALSE, recursive = TRUE)

    p <- progressr::progressor(steps = 14)
    # Convert HUGO to Entrez
    message(qq("Starting pathway analysis on @{data_type}..."))
    entrez_ids <- mapIds(org.Hs.eg.db,
        keys = vector_genes,
        column = "ENTREZID",
        keytype = "SYMBOL",
        multiVals = "first"
    ) # 'first' will return the first match if multiple matches exist.

    p()
    message("GO enrichment MF")
    # Perform GO enrichment analysis
    go_res_mf <- enrichGO(
        gene = vector_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "MF", # Change to "CC" for Cellular Component or "MF" for Molecular Function
        pAdjustMethod = "fdr",
        qvalueCutoff = pathway_p_adjust_thresh
    )

    # # View the results
    # head(go_res_mf)

    # Visualize the results
    go_res_mf_plot <- dotplot(go_res_mf, title = qq("A. GO enrichment Molecular Function for @{data_type}")) +
        labs(caption = qq("Using Fisher's Combined p-value\nq-value cutoff = @{pathway_p_adjust_thresh}\nAdjust method = FDR"))
    go_res_mf_plot

    p()
    message("GO enrichment CC")

    # Perform GO enrichment analysis
    go_res_cc <- enrichGO(
        gene = vector_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "CC", # Change to "CC" for Cellular Component or "MF" for Molecular Function
        pAdjustMethod = "fdr",
        qvalueCutoff = pathway_p_adjust_thresh
    )

    # # View the results
    # head(go_res_cc)

    # Visualize the results
    go_res_cc_plot <- dotplot(go_res_cc, title = qq("B. GO enrichment Cellular Compartment for @{data_type}")) +
        labs(caption = qq("Using Fisher's Combined p-value\nq-value cutoff = @{pathway_p_adjust_thresh}\nAdjust method = FDR"))
    go_res_cc_plot

    p()

    message("GO enrichment BP")

    # Perform GO enrichment analysis
    go_res_bp <- enrichGO(
        gene = vector_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP", # Change to "CC" for Cellular Component or "MF" for Molecular Function
        pAdjustMethod = "fdr",
        qvalueCutoff = pathway_p_adjust_thresh
    )

    # # View the results
    # head(go_res_bp)

    # Visualize the results
    go_res_bp_plot <- dotplot(go_res_bp, title = qq("C. GO enrichment Biological Process for @{data_type}")) +
        labs(caption = qq("Using Fisher's Combined p-value\nq-value cutoff = @{pathway_p_adjust_thresh}\nAdjust method = FDR"))
    go_res_bp_plot

    p()

    message("GO enrichment ALL - MF CC BP")
    # Perform GO enrichment analysis
    go_res_all <- enrichGO(
        gene = vector_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "ALL", # Change to "CC" for Cellular Component or "MF" for Molecular Function
        pAdjustMethod = "fdr",
        qvalueCutoff = pathway_p_adjust_thresh
    )
    p()

    # Visualize the results
    go_res_all_plot <- dotplot(go_res_all, title = qq("C. GO enrichment MF-CC-BP for @{data_type}")) +
        labs(caption = qq("Using Fisher's Combined p-value\nq-value cutoff = @{pathway_p_adjust_thresh}\nAdjust method = FDR"))
    go_res_all_plot
    write.xlsx(go_res_all,
        file = file.path(results_path, qq("ABC--@{data_type}-GeneOntology_MF-CC-BP_table.xlsx"))
    )

    # # Perform Reactome pathway analysis
    # reactome_res <- enrichPathway(
    #     gene = entrez_ids,
    #     organism = "human",
    #     qvalueCutoff = pathway_p_adjust_thresh,
    #     pAdjustMethod = "fdr"
    # )

    # # View the results
    # head(reactome_res)

    # Perform KEGG pathway enrichment analysis
    # use default statistical cut offs
    message("KEGG enrichment")

    kegg_res <- enrichKEGG(
        gene = entrez_ids,
        pAdjustMethod = "fdr",
        qvalueCutoff = pathway_p_adjust_thresh,
        organism = "hsa", # for Homo sapiens. Change according to your species.
    )

    # Visualize the results
    kegg_res_plot <- dotplot(kegg_res, title = qq("D. KEGG enrichment for @{data_type}")) +
        labs(caption = qq("Using Fisher's Combined p-value\nq-value cutoff = @{pathway_p_adjust_thresh}\nAdjust method = FDR"))
    kegg_res_plot

    p()

    message("Saving figures 3A - D")
    figs3_a_d <- tribble(
        ~data, ~fn,
        list(go_res_mf_plot), file.path(results_path, qq("A--@{data_type}-GeneOntology_MF.png")),
        list(go_res_cc_plot), file.path(results_path, qq("B--@{data_type}-GeneOntology_CC.png")),
        list(go_res_bp_plot), file.path(results_path, qq("C--@{data_type}-GeneOntology_BP.png")),
        list(go_res_all_plot), file.path(results_path, qq("ABC--@{data_type}-GeneOntology_MF-CC-BP.png")),
        list(kegg_res_plot), file.path(results_path, qq("D--@{data_type}-KEGG_pathway_enrichment.png"))
    )
    walk2(.x = figs3_a_d$data, .y = figs3_a_d$fn, .f = function(plt, fn) {
        ggsave(
            plot = plt[[1]],
            filename = fn,
            width = 10, height = 8
        )
    })
    p()
    message("Done!")


    ### STRINGdB Section ###

    # for help with these functions
    # STRINGdb$help("plot_network")
    #  STRINGdb$help("get_png")
    message("Starting StringDB analysis...")

    stringdb_results_path <- file.path(results_path, "StringDB")
    dir.create(results_path, showWarnings = FALSE, recursive = TRUE)
    # Create a STRINGdb object
    # score threshold default = 400
    string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 200)

    string_db_plus_our_dat <- data %>%
        distinct(gene_symbol, cmbd_fisher_p_val, log2_fc) %>%
        rename(gene = gene_symbol, pvalue = cmbd_fisher_p_val) %>%
        group_by(gene, pvalue) %>%
        summarize(logFC = mean(log2_fc, na.rm = TRUE), .groups = "drop")

    message("Mapping to STRINGdB and sorting from most significant to least...")
    # Map your differentially expressed proteins to STRING and sort them by most significant to least
    my_mapped_proteins <- string_db$map(as.data.frame(string_db_plus_our_dat), "gene", removeUnmappedRows = TRUE) %>%
        arrange(pvalue)
    p()

    # Get the interactions
    message("Gathering interactions")

    interactions <- string_db$get_interactions(string_ids = my_mapped_proteins$STRING_id) %>%
        arrange(desc(combined_score))
    `Interaction Score` <- interactions$combined_score
    p()

    # hist(`Interaction Score`, main = "Histogram of interactions")
    # head(interactions %>% arrange(desc(combined_score)))
    # Visualize the network

    stringdb_interactions_results_path <- file.path(stringdb_results_path, "interactions")
    dir.create(stringdb_interactions_results_path, showWarnings = FALSE, recursive = TRUE)

    message("PLotting network interactions...")
    message("These are the strongest interactions as predicted by STRINGdb..")
    walk(c(25, 50, 100, 200, 300), .f = function(i) {
        message(qq("Top @{i} interactions"))
        top_25_interactions <- interactions[1:i, , drop = FALSE]
        string_db_interaction_plot <- string_db$plot_network(top_25_interactions)
        # copy the plot from the viewport
        current_plot <- recordPlot()


        interaction_fn <- file.path(
            stringdb_interactions_results_path,
            qq("E--@{data_type}-top_@{i}_stringdb_pp_interactions.png")
        )
        # Save the plot as a PNG file
        png(interaction_fn, width = 3840, height = 2160, res = 300)
        replayPlot(current_plot)
        dev.off()
    })

    p()

    openxlsx::write.xlsx(interactions,
        file = file.path(stringdb_interactions_results_path, "E--interactions_list (ENSMBL).xlsx")
    )

    stringdb_hits_results_path <- file.path(stringdb_results_path, "hits")
    dir.create(stringdb_hits_results_path, showWarnings = FALSE, recursive = TRUE)

    message("Generating graphs of top hits based on logFC pvalue")
    lst_of_hits <- lapply(c(25, 50, 100, 200, 300), FUN = function(x) {
        message(qq("Top @{x} hits"))
        hits <- my_mapped_proteins$STRING_id[1:x]
        string_db_plot1_idx <- string_db$plot_network(hits)
        # copy the plot from the viewport
        current_plot <- recordPlot()

        mapped_fn <- file.path(stringdb_hits_results_path, qq("@{data_type}-fc_pvalue_top_hits_@{x}.png"))
        # Save the plot as a PNG file
        png(mapped_fn, width = 3840, height = 2160, res = 300)
        replayPlot(current_plot)
        dev.off()

        return(hits)
    }) %>%
        set_names(str_c("h", c(25, 50, 100, 200, 300)))

    p()
    # head(subset(my_mapped_proteins, log10(pvalue) >= -log10(0.01) | abs(logFC) >= 0.5))

    # let's add color to up-regulated and down regulated genes
    my_mapped_proteins_sig <- string_db$add_diff_exp_color(
        subset(
            my_mapped_proteins,
            log10(pvalue) >= -log10(string_db_p_value_threshold) | abs(logFC) >= string_db_fc_threshold
        ),
        logFcColStr = "logFC"
    )

    stringdb_pathways_results_path <- file.path(stringdb_results_path, "regulated_pathways")
    dir.create(stringdb_pathways_results_path, showWarnings = FALSE, recursive = TRUE)

    # post payload information to the STRING server
    payload_id <- string_db$post_payload(my_mapped_proteins_sig$STRING_id,
        colors = my_mapped_proteins_sig$color
    )
    message("Generating graphs with fold changes to visualize up-regulated / down-regulated pathway members...")
    # display a STRING network png with the "halo"
    walk(seq_len(length(lst_of_hits)), .f = function(i) {
        x <- lst_of_hits[[i]]
        fn_name <- names(lst_of_hits)[i]
        colored_plot_hits <- string_db$plot_network(x, payload_id = payload_id)
        # copy the plot from the viewport
        current_plot <- recordPlot()

        pathways_fn <- file.path(stringdb_pathways_results_path, qq("@{data_type}-regulated-pathways_@{fn_name}.png"))
        # Save the plot as a PNG file
        png(pathways_fn, width = 3840, height = 2160, res = 300)
        replayPlot(current_plot)
        dev.off()
    })

    p()

    stringdb_clusters_results_path <- file.path(stringdb_results_path, "interaction_clusters")
    dir.create(stringdb_clusters_results_path, showWarnings = FALSE, recursive = TRUE)
    # clusterings
    clustersList <- string_db$get_clusters(my_mapped_proteins$STRING_id[1:300])
    # plot first 4 clusters
    # par(mfrow = c(1, 1))
    message("Clustering interactions and returning the top 4 most significant clusters..")
    for (i in seq_len(4)) {
        clusters_fn <- file.path(stringdb_clusters_results_path, qq("@{data_type}-top@{i}_clusters.png"))
        clsuter_plot <- string_db$plot_network(clustersList[[i]])
        current_plot <- recordPlot()

        png(clusters_fn, width = 3840, height = 2160, res = 300)
        replayPlot(current_plot)
        dev.off()
    }

    p()

    # GO on hits
    stringdb_enrichment_results_path <- file.path(stringdb_results_path, "enrichment_on_hits")
    dir.create(stringdb_enrichment_results_path, showWarnings = FALSE, recursive = TRUE)

    message("Generating tables of enrichments for top hits")
    walk(seq_len(length(lst_of_hits)), .f = function(i) {
        fn_name <- names(lst_of_hits)[i]
        message(fn_name)
        all_enrichments <- string_db$get_enrichment(lst_of_hits[[i]])

        openxlsx::write.xlsx(
            x = all_enrichments,
            file = file.path(stringdb_enrichment_results_path, qq("@{data_type}-GO_and_KEGG-hits@{fn_name}.xlsx"))
        )
    })
    p()

    message("Done!")

    # # because i'm lazy
    # # Assign the comma-separated values to a string
    # text <- "entrez_ids, kegg_res, kegg_res_plot, go_res_bp, go_res_bp_plot, go_res_cc, go_res_cc_plot, go_res_mf, go_res_mf_plot, reactome_res, reactome_plot_res, interactions, mapped_proteins, string_db_plot"
    # # Split the string into a list of strings
    # list_of_strings <- strsplit(text, ", ")[[1]]

    # return(list(
    #     entrez_ids, kegg_res, kegg_res_plot, go_res_bp, go_res_bp_plot,
    #     go_res_cc, go_res_cc_plot, go_res_mf, go_res_mf_plot,
    #     reactome_res, reactome_plot_res, interactions, mapped_proteins,
    #     string_db_plot
    # ) %>% set_names(c(list_of_strings)))
}


with_progress(do_pathway_analysis(data = all_data, data_type = "All_data"))
with_progress(do_pathway_analysis(data = all_data %>% filter(num_accession == 3), data_type = "conserved_data"))
with_progress(do_pathway_analysis(data = data_no_raja, data_type = "No_Raja"))
