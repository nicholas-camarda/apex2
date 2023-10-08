# rm(list = ls())
library(tidyverse)
library(readxl)
library(ggpubr)
library(GetoptLong)
library(metap)
library(rstatix)
library(ggprism)
library(latex2exp)
library(ggrepel)
library(openxlsx)
library(RVenn)
# library(ggVennDiagram)
library(ggvenn)

source(file.path("scripts", "helper_functions.R"))
my_sig_q_val <- 0.1
my_fc_cutoff <- 2

# TODO #1
## table
# hits that are significant in the other 3
# hits that are significant in Rajagopal
# what is shared
# what is not shared

all_data <- read_rds(file.path("cache", "FINAL-cache-including_kruse.rds")) %>%
    group_by(Accession) %>%
    mutate(median_log2_fc = median(log2_fc, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(significant_fc = abs(median_log2_fc) > log2(my_fc_cutoff)) %>%
    dplyr::select(-starts_with(c("use", "stat", "significant_stat_val", "len", "row_number"))) %>%
    mutate(
        list_gene_symbol = map_chr(list_gene_symbol, .f = function(x) str_c(x, collapse = ", ")),
        gene_symbol = map_chr(gene_symbol, .f = function(x) str_split(x, pattern = "_", simplify = TRUE)[, 1])
    ) %>%
    rename(
        original_p_val = p_val
    ) %>%
    dplyr::select(
        owner, Accession, gene_symbol, log2_fc,
        cmbd_fisher_p_val, adjusted_cmbd_fisher_p_val, tx,
        everything()
    ) %>%
    arrange(desc(log2_fc), adjusted_cmbd_fisher_p_val)

# all_data %>% arrange(desc(Accession))
significant_data <- all_data %>%
    filter(
        significant_fc,
        adjusted_cmbd_fisher_p_val < my_sig_q_val # this takes care of pvalue significance across all three
    )
significant_data

make_venn <- function(df, variable = "list_gene_symbol") {
    venns_df <- lapply(df %>%
        distinct(owner, !!sym(variable)) %>%
        split(.$owner), FUN = function(x) x %>% pluck(variable))

    pair_differences <- RVenn::discern_pairs(RVenn::Venn(venns_df))
    pair_intersections <- RVenn::overlap_pairs(RVenn::Venn(venns_df))
    name_title <- ifelse(variable == "Accession", variable, "Gene Symbol")
    venn_diagram <- ggvenn::ggvenn(
        data = venns_df,
        fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
        stroke_size = 0.5, set_name_size = 5
    ) +
        ggtitle(qq("Venn Diagram of Significant @{name_title} Between Datasets")) +
        # theme(plot.background = element_rect(fill = "white")) +
        theme_prism(base_size = 18) +
        theme(
            axis.line = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_blank(), axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(), legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        )
    return(list(venn_diagram, pair_differences, pair_intersections) %>% set_names(c("ggvenn", "pair_diff", "pair_intersect")))
}

save_results <- function(lst, var_name = "accession", subdir = "all_pairs") {
    dir.create(qq("results/tables/@{subdir}"), showWarnings = FALSE, recursive = TRUE)
    ggsave(
        plot = lst$ggvenn,
        filename = file.path("results", "tables", subdir, qq("@{var_name}_venn_diagram.png")),
        width = 10, height = 10
    )

    openxlsx::write.xlsx(lst$pair_diff, file = qq("results/tables/@{subdir}/@{var_name}-pairwise_set_differences.xlsx"))
    openxlsx::write.xlsx(lst$pair_intersect, file = qq("results/tables/@{subdir}/@{var_name}-pairwise_set_intersections.xlsx"))
}

## All pairs
venn_gene_lst <- make_venn(
    df = significant_data,
    variable = "list_gene_symbol"
)
save_results(venn_gene_lst, var_name = "gene_symbol")

venn_gene_lst_all_others <- make_venn(
    df = significant_data %>%
        mutate(owner = ifelse(owner != "Rajagopal", "Other", owner)),
    variable = "list_gene_symbol"
)
save_results(
    lst = venn_gene_lst_all_others,
    var_name = "raja_vs_all_others-gene_symbol",
    subdir = "raja_vs_all_others"
)
