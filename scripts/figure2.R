## ----read--------------------------------------------------------------------------------------------------------------------
# Load the library
library(tidyverse)
library(readxl)
library(ggpubr)
library(GetoptLong)
library(metap)
library(rstatix)
library(ggvenn)
library(ggprism)
library(latex2exp)
library(ggrepel)


# read in the data
my_fns <- tribble(
    ~fn_name, ~sheet,
    "1-s2.0-S0092867417303471-mmc1.xlsx", 2, # This is the sheet we want to use (use the 2 minute time points)
    "5909_SupplementalData_020822.xlsx", 5,
    "Summary_90s_log2_1v1comp.xlsm", 2,
    "Summary_10m_log2_1v1comp.xlsm", 2,
    "1-s2.0-S0092867417303021-mmc3.xlsx", 2
) # use 1, 3, and 10 min time points

fns_temp_tbl <- tibble(fn_path = dir("data", recursive = T, pattern = "xls", full.names = TRUE)) %>%
    mutate(
        owner = str_trim(string = gsub(x = basename(dirname(fn_path)), pattern = "Data", replacement = ""), side = "both"),
        fn_name = basename(fn_path), .before = 1
    )


fns_tbl <- fns_temp_tbl %>%
    inner_join(my_fns, by = "fn_name") %>%
    mutate(data = map2(.x = fn_path, .y = sheet, .f = function(x, y) {
        read_excel(path = x, sheet = y, col_types = "guess")
    }))
fns_tbl

# rajagopal

# is this fold change or difference??
rajagopal_cxcl2_full <- fns_tbl %>%
    filter(owner == "Rajagopal") %>%
    dplyr::select(owner, data) %>%
    unnest(data) %>%
    dplyr::select(owner, Accession, CXCL12_v_Veh_FC, CXCL12_v_Veh_pval, Stats) %>%
    mutate(
        tx = "veh vs cxcl12",
        timing = "3min"
    ) %>%
    group_by(tx) %>%
    rename(log2_fc = CXCL12_v_Veh_FC, p_val = CXCL12_v_Veh_pval) %>%
    mutate(
        adjusted_p_val = p.adjust(p_val, method = "fdr"),
        reg_fc = 2^log2_fc
    )

rajagopal_cxcl2 <- rajagopal_cxcl2_full %>%
    dplyr::select(-Stats)
head(rajagopal_cxcl2)


fn <- file.path("results", "plots", "data_stats", "Rajagopal_3min_data_hist.png")
png(fn, width = 1920, height = 1920, res = 300)
hist(rajagopal_cxcl2$log2_fc, main = "Rajagopal (3min)", xlab = "Log2 Fold Change")
dev.off()

find_na <- function(data, column = "log2_fc") {
    idx <- data %>%
        ungroup() %>%
        pull(column) %>%
        is.na() %>%
        which()

    if (length(idx) == 0) {
        return(tibble())
    }
    return(data %>% mutate(row_number = row_number(), .before = 1) %>% slice(idx))
}

find_na(rajagopal_cxcl2, column = "log2_fc")

## ----kruse-------------------------------------------------------------------------------------------------------------------
# kruse
kruse_ang2_2min <- fns_tbl %>%
    filter(owner == "Kruse") %>%
    dplyr::select(owner, data) %>%
    unnest(data) %>%
    dplyr::select(owner, Accession, `Gene Symbol`, contains("2 min") & contains("Ang II")) %>%
    rename(
        rep1_2min = "2 min Ang II biolgical replicate 1",
        rep2_2min = "2 min Ang II biolgical replicate 2",
        fc_2min = "2 min Ang II enrichment relative to ligand-free control"
    ) %>%
    mutate(timing = "2min") %>%
    rename(rep1 = rep1_2min, rep2 = rep2_2min, reg_fc = fc_2min)

kruse_ang2_21min <- fns_tbl %>%
    filter(owner == "Kruse") %>%
    dplyr::select(owner, data) %>%
    unnest(data) %>%
    dplyr::select(owner, Accession, `Gene Symbol`, contains("21 min") & contains("Ang II")) %>%
    rename(
        rep1_21min = "21 min Ang II biolgical replicate 1",
        rep2_21min = "21 min Ang II biological replicate 2",
        fc_21min = "21 min Ang II enrichment relative to ligand-free control"
    ) %>%
    mutate(timing = "21min") %>%
    rename(rep1 = rep1_21min, rep2 = rep2_21min, reg_fc = fc_21min)

# we're only going to use the 2min time point
kruse_ang2_2min_21min_full <- bind_rows(kruse_ang2_2min, kruse_ang2_21min) %>%
    mutate(tx = "veh vs angII") %>%
    rename(gene_symbol = `Gene Symbol`) %>%
    mutate(log2_fc = log2(reg_fc))

kruse_ang2_2min_21min <- kruse_ang2_2min_21min_full %>%
    dplyr::select(-rep1, -rep2)

head(kruse_ang2_2min_21min)
fn <- file.path("results", "plots", "data_stats", "Kruse_2min_data_hist.png")
png(fn, width = 1920, height = 1920, res = 300)
hist(kruse_ang2_2min_21min %>% filter(timing == "2min") %>% .$log2_fc, main = "Kruse (2min)", xlab = "Log2 Fold Change")
dev.off()

find_na(kruse_ang2_2min_21min, column = "log2_fc")


## ----rockman-----------------------------------------------------------------------------------------------------------------
# Rockman

# Rockman indices == 3 and 4
rockman_10m_90s_full <- bind_rows(
    fns_tbl %>%
        slice(3) %>%
        dplyr::select(owner, data) %>%
        unnest(data) %>%
        mutate(timing = "10min"),
    fns_tbl %>%
        slice(4) %>%
        dplyr::select(owner, data) %>%
        unnest(data) %>%
        mutate(timing = "90sec")
) %>%
    dplyr::select(owner, Accession, `Gene Name`, Ang, NS, `Ang-NS`, `Ang-NS pval`, `Ang-NSpval (adj)`, timing) %>%
    rename(gene_symbol = `Gene Name`) %>%
    rename(reg_fc = `Ang-NS`, p_val = `Ang-NS pval`, adjusted_p_val = `Ang-NSpval (adj)`) %>%
    mutate(log2_fc = reg_fc) %>%
    # I think Rockman is already log transforme, because we see negative values
    # mutate(log2_fc = log2(reg_fc)) %>%
    mutate(tx = "ns vs angII")

rockman_10m <- rockman_10m_90s_full %>% filter(timing == "10min")
rockman_90s <- rockman_10m_90s_full %>% filter(timing == "90sec")

rockman_10m_90s <- rockman_10m_90s_full %>%
    dplyr::select(-Ang, -NS)
head(rockman_10m_90s)

fn <- file.path("results", "plots", "data_stats", "Rockman_90sec_data_hist.png")
png(fn, width = 1920, height = 1920, res = 300)
hist(rockman_10m_90s %>% filter(timing == "90sec") %>% .$log2_fc, main = "Rockman (90sec)", xlab = "Log2 Fold Change")
dev.off()

fn <- file.path("results", "plots", "data_stats", "Rockman_10min_data_hist.png")
png(fn, width = 1920, height = 1920, res = 300)
hist(rockman_10m_90s %>% filter(timing == "10min") %>% .$log2_fc, main = "Rockman (10min)", xlab = "Log2 Fold Change")
dev.off()

find_na(rockman_10m_90s, column = "log2_fc")
find_na(rockman_10m_90s, column = "p_val")

## ----zastrow-----------------------------------------------------------------------------------------------------------------
# Von Zastrow
Zastro_all <- fns_tbl %>%
    filter(owner == "Von Zastrow") %>%
    dplyr::select(owner, data) %>%
    unnest(data) %>%
    dplyr::select(owner, `Uniprot accession`, `Uniprot name`, `Gene symbol`, FC, FDR, Timing) %>%
    rename(
        Accession = `Uniprot accession`, gene_symbol = `Gene symbol`, timing = `Timing`, reg_fc = FC,
        adjusted_p_val = FDR
    ) %>%
    # hacking this, making the unadjusted p-value = to FDR. this is probably wrong...
    mutate(
        p_val = adjusted_p_val,
        # probably the rest of the data is log2 transformed
        log2_fc = log2(reg_fc)
    ) %>%
    mutate(
        timing = str_split(timing, pattern = "_", simplify = TRUE)[, 3],
        tx = "veh vs DADLE"
    ) %>%
    filter(timing %in% c("1min", "3min", "10min"))

zastro_1min_3min_10min <- Zastro_all %>%
    mutate(gene_symbol = str_split(`Uniprot name`, pattern = "_", simplify = TRUE)[, 1]) %>%
    dplyr::select(-`Uniprot name`) %>%
    filter(!is.na(reg_fc) & !is.na(adjusted_p_val))


# i think zastro is in fold change, and the rest are in log fold change
# zastro_1min_3min_10min %>% filter(log2_fc < 0)
head(zastro_1min_3min_10min)
# hist(zastro_1min_3min_10min$log2_fc, main = "Von Zastrow", xlab = "Log2 Fold Change")

fn <- file.path("results", "plots", "data_stats", "VonZastrow_3min_data_hist.png")
png(fn, width = 1920, height = 1920, res = 300)
hist(zastro_1min_3min_10min %>% filter(timing == "3min") %>% .$log2_fc, main = "Von Zastrow (3min)", xlab = "Log2 Fold Change")
dev.off()
find_na(zastro_1min_3min_10min, column = "adjusted_p_val")

# get mapping of gene symbols, https://www.uniprot.org/id-mapping
# accessions_to_convert <- tibble(accession = combined_dat$Accession %>% unique())
# write_tsv(file = "extra/uniprot_accessions.tsv", x = accessions_to_convert)

gene_mapping_tbl <- read_tsv(
    file = "extra/uniprot-download_true_fields_accession_2Creviewed_2Cid_2Cprotein_nam-2022.08.23-15.54.46.84.tsv",
    show_col_types = FALSE
)

my_gene_map_temp <- gene_mapping_tbl %>%
    dplyr::select(Entry, `Gene Names`) %>%
    mutate(
        gene_symbol_all = map(.x = `Gene Names`, .f = function(x) as.list(str_split(x, pattern = " ", simplify = TRUE))),
        gene_symbol = map_chr(.x = gene_symbol_all, .f = function(x) x[[1]])
    ) %>%
    rename(Accession = Entry) %>%
    distinct()

helper_accessions <- kruse_ang2_2min_21min %>%
    distinct(Accession, gene_symbol) %>%
    arrange(gene_symbol)

helper_accessions2 <- zastro_1min_3min_10min %>%
    distinct(Accession, gene_symbol) %>%
    arrange(gene_symbol)

my_gene_map <- bind_rows(my_gene_map_temp, helper_accessions, helper_accessions2) %>%
    distinct(Accession, gene_symbol) %>%
    arrange(gene_symbol) %>%
    filter(gene_symbol != "-") %>%
    group_by(Accession) %>%
    summarize(list_gene_symbol = list(gene_symbol), chr_gs = str_c(gene_symbol, collapse = ",")) %>%
    mutate(len = map_int(list_gene_symbol, .f = length))

my_gene_map %>% filter(len > 1)


# Full data
combined_dat_full <- bind_rows(
    rajagopal_cxcl2,
    kruse_ang2_2min_21min,
    rockman_10m_90s,
    zastro_1min_3min_10min
) %>%
    ungroup() %>%
    dplyr::select(-gene_symbol) %>%
    # seems like i need to do this, acecssions with dashses aren't processed properly
    mutate(Accession = str_split(Accession, "-", simplify = TRUE)[, 1]) %>%
    left_join(my_gene_map, by = "Accession") %>%
    mutate(gene_symbol = map2_chr(list_gene_symbol, len, .f = function(x, y) {
        gs <- ifelse(!is.na(y), x[[1]], NA)
        return(gs)
    }))

# combined_dat_full %>% slice(14)
# my_gene_map %>% filter(Accession == "Q9P273")

combined_dat_temp <- combined_dat_full %>%
    filter((owner == "Rajagopal" & timing == "3min") |
        (owner == "Kruse" & timing == "2min") |
        (owner == "Rockman" & timing == "90sec") |
        owner == "Von Zastrow" & timing == "3min") %>%
    filter(!is.na(log2_fc)) %>% # c("1min", "3min", "10min") %>%
    filter(!is.na(Accession)) %>%
    mutate(row_number = row_number(), .before = 1)

remove_dups <- combined_dat_temp %>%
    # there are some duplicates so we keep the one with the lowest p-value
    group_by(Accession, gene_symbol, tx, owner) %>%
    mutate(dups = max(n())) %>%
    filter(dups >= 2) %>%
    group_by(Accession, owner) %>%
    arrange(Accession, owner) %>%
    mutate(which_best = ifelse(!is.na(p_val),
        min(p_val) == p_val,
        max(abs(log2_fc)) == log2_fc
    )) %>%
    filter(!which_best)
remove_dups

combined_dat <- combined_dat_temp %>%
    filter(!(row_number %in% remove_dups$row_number))

message(qq("Original dataset size: @{nrow(combined_dat_full)}"))
message(qq("Removed @{nrow(combined_dat_full) - nrow(combined_dat)} rows with NA accession or NA log2_fc"))
message(qq("Dataset size with outliers: @{nrow(combined_dat)}"))

dim(combined_dat)

fn <- file.path("results", "plots", "data_stats", "combined_data_hist.png")
png(fn, width = 1920, height = 1920, res = 300)
hist(combined_dat$log2_fc, main = "Combined Data", xlab = "Log2 Fold Change")
dev.off()

# detect outliers by group
extreme_outliers_by_owner <- combined_dat %>%
    group_by(owner) %>%
    rstatix::identify_outliers(variable = "log2_fc") %>%
    ungroup() %>%
    filter(is.extreme) %>%
    arrange(desc(log2_fc))
# Outliers
print(extreme_outliers_by_owner, n = Inf)

extreme_outliers_all <- combined_dat %>%
    # group_by(owner) %>%
    rstatix::identify_outliers(variable = "log2_fc") %>%
    ungroup() %>%
    filter(is.extreme) %>%
    arrange(desc(log2_fc))

print(extreme_outliers_all, n = Inf)

# remove all the outliers
combined_dat_full_no_outliers <- combined_dat %>%
    filter(!(row_number %in% extreme_outliers_all$row_number))
message(qq("Final Dataset size without outliers: @{nrow(combined_dat_full_no_outliers)}"))

fn <- file.path("results", "plots", "data_stats", "combined_no_outliers_data_hist.png")
png(fn, width = 1920, height = 1920, res = 300)
hist(combined_dat_full_no_outliers$log2_fc, main = "Combined Data (no outliers)", xlab = "Log2 Fold Change")
dev.off()

# Let's go with the the full data
combined_dat_no_outliers <- combined_dat_full_no_outliers

## ----metamsd, echo = FALSE---------------------------------------------------------------------------------------------------
# Code to generate datasets usable in MetaMSD
#### START #######
to_write <- tibble(datasets = map(.x = combined_dat_no_outliers %>% split(.$owner), .f = function(d) {
    d %>%
        dplyr::select(Protein = Accession, Sign = log2_fc, Pvalue = p_val) %>%
        filter(!is.na(Protein))
})) %>%
    mutate(
        name_data = names(datasets),
        output_fns = file.path(getwd(), "../MetaMSD/MetaMSD/input", str_c(name_data, ".txt"))
    )

find_nas <- to_write %>%
    unnest(datasets)

find_na(find_nas, column = "Sign")
check_na_proteins <- find_na(find_nas, column = "Pvalue")

check_df1 <- find_nas %>%
    dplyr::select(name_data, Protein, Sign, Pvalue) %>%
    distinct() %>%
    pivot_wider(id_cols = c(Protein), names_from = name_data, values_from = Pvalue)

to_write_no_na <- tibble(datasets = map(.x = combined_dat %>% split(.$owner), .f = function(d) {
    d %>%
        dplyr::select(Protein = Accession, Sign = log2_fc, Pvalue = p_val) %>%
        filter(!(Protein %in% check_na_proteins$Protein))
})) %>%
    mutate(
        name_data = names(datasets),
        output_fns = file.path(getwd(), "../MetaMSD/MetaMSD/input", str_c(name_data, "_no_na.txt"))
    )

walk2(.x = as.list(to_write$datasets), .y = as.list(to_write$output_fns), .f = function(x, y) write_tsv(x = x, file = y))
walk2(.x = as.list(to_write_no_na$datasets), .y = as.list(to_write_no_na$output_fns), .f = function(x, y) write_tsv(x = x, file = y))
#### END #######

fn <- file.path("results", "plots", "data_stats", "p_value_distr", "raj_pval_distr.png")
png(fn, width = 1920, height = 1920, res = 300)
raja_dstr_plot <- plotp(combined_dat_no_outliers %>% filter(owner == "Rajagopal") %>% .$p_val,
    main = "Distribution of p-values in Rajagopal dataset (no outliers)"
)
dev.off()


save_pval_distr <- function(data) {
    message("Saving p-value distributions...")
    fn <- file.path("results", "plots", "data_stats", "p_value_distr", "rockman_pval_distr.png")
    png(fn, width = 1920, height = 1920, res = 300)
    rockman_dstr_plot <- plotp(data %>% filter(owner == "Rockman") %>% .$p_val,
        main = "Distribution of p-values in Rockman dataset (no outliers)"
    )
    dev.off()
    #' @note no pvalues from Kruse
    # plotp(combined_dat_no_outliers %>% filter(owner == "Kruse") %>% .$p_val,
    #       main = "Distribution of p-values in Kruse dataset (no outliers)")


    fn <- file.path("results", "plots", "data_stats", "p_value_distr", "vonzastrow_pval_distr.png")
    png(fn, width = 1920, height = 1920, res = 300)
    vonzastrow_dstr_plot <- plotp(data %>% filter(owner == "Von Zastrow") %>% .$p_val,
        main = "Distribution of p-values in Von Zastrow dataset (no outliers)"
    )
    dev.off()


    fn <- file.path("results", "plots", "data_stats", "p_value_distr", "all-together_pval_distr.png")
    png(fn, width = 1920, height = 1920, res = 300)
    # All together
    all_distr_plot <- plotp(data$p_val,
        main = "Distribution of p-values in Combined dataset (no outliers)"
    )
    dev.off()
    message("Done!")
}

save_pval_distr(data = combined_dat_no_outliers)
# large sample sizes here makes test for normality meaningless here due to CLT

generate_cmbd_fisher_pval_data <- function(data) {
    input_to_stat <- data %>%
        dplyr::select(-chr_gs, -adjusted_p_val, -reg_fc) %>%
        mutate(
            two_sidedness = TRUE,
            # toinvert will depend on the result of the statistical test, based on which 'side' you were testing
            # i'm assuming that all tests are two sided and none need to be inverted
            toinvert = ifelse(two_sidedness, FALSE, TRUE)
        ) %>%
        # remove entries with na pvalues - they won't contribute
        filter(!is.na(p_val))
    # nest(data = c(owner, log2_fc, p_val, tx, timing, two_sidedness, toinvert))
    message("Removing NA fisher p-vals\n(i.e. n < 2 or any NA p-values used to construct the combined p-value)")
    onesided_fisher <- two2one(p = input_to_stat$p_val, two = input_to_stat$two_sidedness, invert = input_to_stat$toinvert)
    final_fisher_df <- input_to_stat %>%
        bind_cols(onesided_fisher = tibble(onesided_fisher)) %>%
        group_by(Accession) %>%
        mutate(cmbd_fisher_p_val = metap::sumlog(p = onesided_fisher)$p) %>%
        filter(!is.na(cmbd_fisher_p_val)) %>%
        dplyr::select(-c(two_sidedness, toinvert, onesided_fisher)) %>%
        ungroup() %>%
        mutate(
            adjusted_cmbd_fisher_p_val = p.adjust(cmbd_fisher_p_val, method = "fdr")
        ) %>%
        ungroup() %>%
        bind_rows(data %>%
            filter(owner == "Kruse") %>%
            dplyr::select(-chr_gs, -adjusted_p_val, -reg_fc)) %>%
        group_by(Accession) %>%
        mutate(num_accession = max(n())) %>%
        filter(num_accession > 1) %>%
        ungroup() %>%
        arrange(Accession) %>%
        mutate(rown_number = row_number())

    message("Note: Generated new row numbers for the final dataset.\nThese will NOT match previous row numbers!!\nThese are just for idx purposes, no real functional use")
    # get all the log2_fc data for those analytes, since we removed some rows with NA p-values above
    return(final_fisher_df)
}

message("Generating combined p-values!")
cmbd_p_vals_all_completed_data <- generate_cmbd_fisher_pval_data(data = combined_dat_no_outliers)
cmbd_p_vals_all_completed_data_no_raja <- generate_cmbd_fisher_pval_data(data = combined_dat_no_outliers %>% filter(owner != "Rajagopal"))


## P-value cutoffs
# p_cutoff <- 0.001
p_adjust_cutoff <- 0.2
results_path <- file.path("results", "plots", "Figure 2")
dir.create(results_path, showWarnings = FALSE, recursive = TRUE)

# sig_fisher_df <- cmbd_p_vals_all_completed_data %>%
#     filter(adjusted_cmbd_fisher_p_val < p_adjust_cutoff | owner == "Kruse")
sig_fisher_df <- cmbd_p_vals_all_completed_data %>%
    filter(cmbd_fisher_p_val < p_cutoff | owner == "Kruse")

naming_tribble <- tribble(
    ~data, ~gvar, ~title, ~path,
    list(sig_fisher_df), "Accession", "Waterfall plot of significant analytes by Fisher p-value", "All_data-accession_fc_waterfall.png",
    list(sig_fisher_df %>% filter(num_accession == 3)), "Accession", "Waterfall plot of significant analytes by Fisher p-value\n(Conserved hits)", "All_data-accession_fc_waterfall-conserved.png",
    list(sig_fisher_df), "gene_symbol", "Waterfall plot of significant analytes by Fisher p-value", "All_data-gene_symbol_fc_waterfall.png",
    list(sig_fisher_df %>% filter(num_accession == 3)), "gene_symbol", "Waterfall plot of significant analytes by Fisher p-value\n(Conserved hits)", "All_data-gene_symbol_fc_waterfall-conserved.png"
)
message("Generating water fall plots!")

#' @note plotting watefall
plot_fc <- function(df, grouping_variable = "Accession", p_cutoff_val = NA, title_ = NA) {
    grouping_var_sym <- sym(grouping_variable)

    p1_plot_dat <- df %>%
        group_by(Accession, gene_symbol) %>%
        summarize(
            mean_ = mean(log2_fc, na.rm = TRUE),
            log2_fc = list(log2_fc),
            owner = list(owner), .groups = "drop"
        ) %>%
        arrange(desc(mean_)) %>%
        mutate(order = row_number()) %>%
        mutate(
            gene_symbol = factor(gene_symbol, levels = gene_symbol[order]),
            Accession = factor(Accession, levels = Accession[order]),
            up_down = factor(ifelse(mean_ > 0, "++", "--"),
                levels = c("++", "--")
            )
        )

    if (grouping_variable == "gene_symbol") grouping_variable <- "Gene Symbol"

    ggplot(p1_plot_dat %>% unnest(log2_fc),
        mapping = aes(
            x = !!grouping_var_sym, y = log2_fc,
            color = up_down, group = !!grouping_var_sym
        )
    ) +
        geom_boxplot() +
        geom_jitter(width = 0.1) +
        theme_prism(base_size = 16) +
        theme(
            strip.text.x = element_text(size = rel(1.25), face = "bold"),
            axis.text.x = element_text(size = rel(0.75)), # , angle = 45, hjust = 1, vjust = 1),
            panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = 0.5),
            panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = 0.25)
        ) +
        scale_color_manual(values = c("black", "red")) +
        ggtitle(qq("@{title_}")) +
        labs(caption = qq("p-val cuttoff: @{p_cutoff_val}")) +
        ylab(TeX("$\\bf{log_2(FC)}$")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        geom_hline(yintercept = 0, color = "orange") +
        labs(x = grouping_variable)
}

for (i in seq_len(nrow(naming_tribble))) {
    waterfall_plot <- plot_fc(
        df = naming_tribble$data[[i]][[1]],
        grouping_variable = naming_tribble$gvar[[i]],
        p_cutoff_val = p_cutoff,
        title_ = naming_tribble$title[[i]]
    )
    ggsave(
        plot = waterfall_plot,
        filename = file.path(results_path, naming_tribble$path[[i]]),
        width = 10, height = 8
    )
}

plot_volcano <- function(df, p_cutoff_val = 0.001, label_variable = "Accession", fc_thresh = 2, title_ = "Volcano plot") {
    label_var_sym <- sym(label_variable)
    plot_df <- df %>%
        group_by(!!label_var_sym) %>%
        mutate(
            negative_log_10_cmbd_fisher_p_val = -log10(cmbd_fisher_p_val),
            mean_log2_fc = mean(log2_fc)
        ) %>%
        mutate(significant = ifelse(abs(mean_log2_fc) > fc_thresh & negative_log_10_cmbd_fisher_p_val > -log10(p_cutoff_val), TRUE, FALSE)) %>%
        mutate(significant_alpha = ifelse(abs(mean_log2_fc) > fc_thresh & negative_log_10_cmbd_fisher_p_val > -log10(p_cutoff_val), 1, 0.01))

    label_df <- plot_df %>%
        filter(significant == TRUE) %>%
        distinct(
            .data[[label_variable]], mean_log2_fc,
            negative_log_10_cmbd_fisher_p_val,
            significant, significant_alpha
        )
    if (label_variable == "gene_symbol") label_variable <- "Gene Symbol"
    # Assuming df is a data frame with columns logFC (log fold change) and pval (p-value)
    ggplot(plot_df, aes(
        x = mean_log2_fc,
        y = negative_log_10_cmbd_fisher_p_val
    ), color = "black", ) +
        geom_hline(yintercept = p_cutoff_val, color = "orange") +
        geom_vline(xintercept = c(-fc_thresh, fc_thresh), color = "blue") +
        geom_point(
            mapping = aes(
                fill = significant,
                alpha = significant_alpha
            ),
            pch = 21,
            show.legend = FALSE, size = 2
        ) +
        xlab(TeX("$\\bf{log_2(FC)}$")) +
        ylab(TeX("$\\bf{-log_10(p)}$")) +
        scale_fill_manual(values = c("black", "red")) +
        geom_label_repel(label_df,
            mapping = aes(label = !!label_var_sym), show.legend = FALSE
        ) +
        theme_prism(base_size = 16) +
        theme(
            strip.text.x = element_text(size = rel(1.25), face = "bold"),
            axis.text.x = element_text(size = rel(0.75)), # , angle = 45, hjust = 1, vjust = 1),
            panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = 0.5),
            panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = 0.25)
        ) +
        ggtitle(qq("@{title_}")) +
        labs(caption = qq("p-val cuttoff: @{p_cutoff_val}\nFold Change cutoff: @{fc_thresh}")) +
        ylab(TeX("$\\bf{log_2(FC)}$")) +
        labs(x = label_variable)
}


naming_tribble_volcano <- tribble(
    ~data, ~gvar, ~title, ~path,
    list(cmbd_p_vals_all_completed_data), "Accession", "Volcano plot of significant analytes by Fisher p-value", "All_data-accession_volcano.png",
    list(cmbd_p_vals_all_completed_data %>% filter(num_accession == 3)), "Accession", "Volcano plot of significant analytes by Fisher p-value\n(Conserved hits)", "All_data-accession_volcano-conserved.png",
    list(cmbd_p_vals_all_completed_data_no_raja), "Accession", "Volcano plot of significant analytes by Fisher p-value\n(No Rajagopal data)", "No_Raj-accession_volcano.png",
    #####
    list(cmbd_p_vals_all_completed_data), "gene_symbol", "Volcano plot of significant analytes by Fisher p-value", "All_data-gene_symbol_volcano.png",
    list(cmbd_p_vals_all_completed_data %>% filter(num_accession == 3)), "gene_symbol", "Volcano plot of significant analytes by Fisher p-value\n(Conserved hits)", "All_data-gene_symbol_volcano-conserved.png",
    list(cmbd_p_vals_all_completed_data_no_raja), "gene_symbol", "Waterfall plot of significant analytes by Fisher p-value\n(No Rajagopal data)", "No_Raj-gene_symbol_volcano.png"
)

for (i in seq_len(nrow(naming_tribble_volcano))) {
    volcano_plot <- plot_volcano(
        df = naming_tribble_volcano$data[[i]][[1]],
        p_cutoff_val = p_cutoff,
        fc_thresh = 2,
        label_variable = naming_tribble_volcano$gvar[[i]],
        title_ = naming_tribble_volcano$title[[i]]
    )
    ggsave(
        plot = volcano_plot,
        filename = file.path(results_path, naming_tribble_volcano$path[[i]]),
        width = 8, height = 8
    )
}

message("Generating volcano plots!")

# no rajagopal
accession_volcano <- plot_volcano(
    df = cmbd_p_vals_all_completed_data_no_raja,
    p_cutoff_val = p_cutoff,
    fc_thresh = 2,
    label_variable = "Accession",
    title_ = "Volcano plot of significant analytes by Fisher p-value\n(No Rajagopal)"
)

gene_symbol_volcano <- plot_volcano(
    df = cmbd_p_vals_all_completed_data_no_raja,
    p_cutoff_val = p_cutoff,
    fc_thresh = 2,
    label_variable = "gene_symbol",
    title_ = "Volcano plot of significant analytes by Fisher p-value\n(No Rajagopal)"
)

acc_volcano_fn <- file.path(results_path, "No_Raj-accession_volcano.png")
ggsave(acc_volcano_fn, accession_volcano, width = 8, height = 8)
gs_volcano_fn <- file.path(results_path, "No_Raj-gene_symbol_volcano.png")
ggsave(gs_volcano_fn, gene_symbol_volcano, width = 8, height = 8)


message("Caching data for later...")

write_rds(cmbd_p_vals_all_completed_data, file = file.path("cache", "FINAL-all_data.rds"))
write_rds(cmbd_p_vals_all_completed_data_no_raja, file = file.path("cache", "FINAL-data_no_raja.rds"))

message("All done!")





# # Install the ReactomePA package
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ReactomePA")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("STRINGdb")

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# BiocManager::install("multtest")
# Install the clusterProfiler package
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("clusterProfiler")

# Install the package if you haven't already
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db")


# OLD

# n_tally_all <- combined_dat_no_outliers %>%
#     dplyr::select(owner, Accession) %>%
#     distinct() %>% # we are counting unique occurrences of a protein per lab (owner)
#     group_by(Accession) %>%
#     summarize(
#         n = n(),
#         labs = str_c(sort(unique(owner)), collapse = ","),
#         labs_lst = list(owner), .groups = "keep"
#     ) %>%
#     arrange(desc(n))

# four_n_tally_all <- n_tally_all %>%
#     filter(n == 4)

# three_n_all <- n_tally_all %>%
#     filter(n >= 3) %>%
#     nrow()
# three_n_only_raj <- n_tally_all %>%
#     filter(Accession %in% raj_proteins$Accession) %>%
#     filter(n >= 3) %>%
#     nrow()

# message(qq("Number of proteins conserved in any 3 or more datasets: @{three_n_all}\nNumber of proteins conserved in at least 3 datasets, using only proteins from Raj: @{three_n_only_raj}"))


# message(qq("Number of shared proteins across all datasets, whether or not they are in Rajagopal's dataset: @{nrow(four_n_tally_all)}"))


# raj_proteins <- combined_dat_no_outliers %>%
#     filter(owner == "Rajagopal")

# four_n_tally_only_raj <- n_tally_all %>%
#     filter(Accession %in% raj_proteins$Accession)


# message(qq("Number of shared proteins across all datasets, using only proteins in Rajagopal's dataset: @{nrow(four_n_tally_only_raj)}"))

## conserved hits without Raj, at least in one other dataset
# n_tally_no_raj_temp <- n_tally_all %>%
#     # check if raj is in these
#     mutate(not_raj = map2_lgl(
#         # find accessions that weren't present in Raj dataset
#         .x = "Rajagopal", .y = labs_lst,
#         .f = function(x, y) !(x %in% y[[1]])
#     ))

# # filter based on the number of times this accession occurs in other datasets
# n_tally_no_raj <- n_tally_no_raj_temp %>%
#     filter(n >= 2, not_raj)
# n_tally_no_raj
# message(qq("Number of proteins conserved in at least 2 datasets and not in Rajagopal's dataset : @{nrow(n_tally_no_raj)}"))

## no raja strict set difference
# accession_volcano <- plot_volcano(
#   df = cmbd_p_vals_all_completed_data_no_raja_set_difference,
#   p_cutoff_val = p_cutoff,
#   label_variable = "Accession",
#   fc_thresh = 2,
#   title_ = "Volcano plot of significant analytes by Fisher p-value\n(No Rajagopal Set Difference)"
# )

# gene_symbol_volcano <- plot_volcano(
#   df = cmbd_p_vals_all_completed_data_no_raja_set_difference,
#   p_cutoff_val = p_cutoff,
#   label_variable = "gene_symbol",
#   fc_thresh = 2,
#   title_ = "Volcano plot of significant analytes by Fisher p-value\n(No Rajagopal Set Difference)"
# )

# no_raja_set_difference <- combined_dat_no_outliers %>%
#     filter(Accession %in% n_tally_no_raj$Accession)
# cmbd_p_vals_all_completed_data_no_raja_set_difference <- generate_cmbd_fisher_pval_data(data = no_raja_set_difference)
# write_rds(cmbd_p_vals_all_completed_data, "results/cmbd_p_vals-everything.rds")


# acc_volcano_fn <- file.path(results_path, "No_Raj_set_diff-accession_volcano.png")
# ggsave(acc_volcano_fn, accession_volcano, width = 8, height = 8)
# gs_volcano_fn <- file.path(results_path, "No_Raj_set_diff-gene_symbol_volcano.png")
# ggsave(gs_volcano_fn, gene_symbol_volcano, width = 8, height = 8)

# to use later
