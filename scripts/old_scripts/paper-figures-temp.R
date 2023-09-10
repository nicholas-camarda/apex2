library(ReactomePA)
library(STRINGdb)
library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)
library(readxl)
library(ggpubr)
library(GetoptLong)
library(patchwork)

cmbd_p_vals <- read_rds("results/cmbd_p_vals-everything.rds")
pathway_p_adjust_thresh <- 0.25
# fisher_to_plot_abbrev
diff_expressed_proteins_all_four_df <- cmbd_p_vals %>%
  filter(n == 4) %>%
  filter(adj_cmbd_fisher_p_val < pathway_p_adjust_thresh)

# FIX there's somethign wrong here
diff_expressed_proteins_all_four_df %>% filter(!is.na(reg_fc))

diff_expressed_proteins_all_four <- diff_expressed_proteins_all_four_df %>%
  .$gene_symbol %>%
  unique()
length(diff_expressed_proteins_all_four)

diff_expressed_proteins_no_raj_df <- cmbd_p_vals %>%
  mutate(not_raj = map2_lgl(
    .x = "Rajagopal", .y = labs_lst,
    .f = function(x, y) !(x %in% y[[1]])
  )) %>%
  filter(n >= 2, not_raj, !is.na(gene_symbol)) %>%
  # filter(adj_cmbd_fisher_p_val < pathway_p_adjust_thresh) %>%
  unnest(data)

diff_expressed_proteins_no_raj <- diff_expressed_proteins_no_raj_df %>%
  .$gene_symbol %>%
  unique()

nrow(diff_expressed_proteins_no_raj)

to_analyze <- list(diff_expressed_proteins_all_four, diff_expressed_proteins_no_raj)


lst_result <- lapply(to_analyze, FUN = function(vector_genes) {
  # Convert HUGO to Entrez
  entrez_ids <- mapIds(org.Hs.eg.db,
    keys = vector_genes,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  ) # 'first' will return the first match if multiple matches exist.


  # Perform KEGG pathway enrichment analysis
  kegg_res <- enrichKEGG(
    gene = entrez_ids,
    pAdjustMethod = "fdr",
    organism = "hsa", # for Homo sapiens. Change according to your species.
    qvalueCutoff = pathway_p_adjust_thresh
  )

  # View the results
  head(kegg_res)

  # Visualize the results
  kegg_res_plot <- dotplot(kegg_res, title = qq("KEGG enrichment (Fisher's) for all datasets\nq = @{pathway_p_adjust_thresh}"))


  # Perform GO enrichment analysis
  go_res_bp <- enrichGO(
    gene = vector_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP", # Change to "CC" for Cellular Component or "MF" for Molecular Function
    pAdjustMethod = "fdr",
    qvalueCutoff = pathway_p_adjust_thresh
  )

  # View the results
  head(go_res_bp)

  # Visualize the results
  go_res_bp_plot <- dotplot(go_res_bp, title = qq("GO enrichment BP (Fisher's) for all datasets\nq = @{pathway_p_adjust_thresh}"))

  # Perform GO enrichment analysis
  go_res_cc <- enrichGO(
    gene = vector_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "CC", # Change to "CC" for Cellular Component or "MF" for Molecular Function
    pAdjustMethod = "fdr",
    qvalueCutoff = pathway_p_adjust_thresh
  )

  # View the results
  head(go_res_cc)

  # Visualize the results
  go_res_cc_plot <- dotplot(go_res_cc, title = qq("GO enrichment CC (Fisher's) for all datasets\nq = @{pathway_p_adjust_thresh}"))

  # Perform GO enrichment analysis
  go_res_mf <- enrichGO(
    gene = vector_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "MF", # Change to "CC" for Cellular Component or "MF" for Molecular Function
    pAdjustMethod = "fdr",
    qvalueCutoff = pathway_p_adjust_thresh
  )

  # View the results
  head(go_res_mf)

  # Visualize the results
  go_res_mf_plot <- dotplot(go_res_mf, title = qq("GO enrichment MF (Fisher's) for all datasets\nq = @{pathway_p_adjust_thresh}"))


  # Perform Reactome pathway analysis
  reactome_res <- enrichPathway(
    gene = entrez_ids,
    organism = "human",
    qvalueCutoff = pathway_p_adjust_thresh,
    pAdjustMethod = "fdr"
  )

  # View the results
  head(reactome_res)

  # Visualize the results
  reactome_plot_res <- dotplot(reactome_res, title = qq("Reactome (Fisher's) for all datasets\nq = @{pathway_p_adjust_thresh}"))



  # Create a STRINGdb object
  string_db <- STRINGdb$new(version = "11.0b", species = 9606, score_threshold = 400)

  # Map your differentially expressed proteins to STRING
  mapped_proteins <- string_db$map(data.frame(gene = vector_genes), "gene", removeUnmappedRows = TRUE)

  # Get the interactions
  interactions <- string_db$get_interactions(string_ids = mapped_proteins$STRING_id)

  # Visualize the network
  string_db_plot <- string_db$plot_network(interactions)

  # because i'm lazy
  # Assign the comma-separated values to a string
  text <- "entrez_ids, kegg_res, kegg_res_plot, go_res_bp, go_res_bp_plot, go_res_cc, go_res_cc_plot, go_res_mf, go_res_mf_plot, reactome_res, reactome_plot_res, interactions, mapped_proteins, string_db_plot"
  # Split the string into a list of strings
  list_of_strings <- strsplit(text, ", ")[[1]]

  return(list(
    entrez_ids, kegg_res, kegg_res_plot, go_res_bp, go_res_bp_plot,
    go_res_cc, go_res_cc_plot, go_res_mf, go_res_mf_plot,
    reactome_res, reactome_plot_res, interactions, mapped_proteins,
    string_db_plot
  ) %>% set_names(c(list_of_strings)))
})

lst_result
