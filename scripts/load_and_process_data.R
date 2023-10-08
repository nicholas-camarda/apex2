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
# selected these sheets based off of discussions with Dylan
my_fns <- tribble(
    ~fn_name, ~sheet,
    "1-s2.0-S0092867417303471-mmc1.xlsx", 2, # This is the sheet we want to use (use the 2 minute time points)
    "5909_SupplementalData_020822.xlsx", 5,
    "Summary_90s_log2_1v1comp.xlsm", 2,
    "Summary_10m_log2_1v1comp.xlsm", 2,
    "1-s2.0-S0092867417303021-mmc3.xlsx", 2
)

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

# ## ----rajagopal-------------------------------------------------------------------------------------------------------------------


extract_gn_from_description_raja <- Vectorize(function(string) {
    # First, use str_detect to identify if the pattern exists
    if (str_detect(string, "GN=[A-Za-z0-9]+")) {
        # Then, use str_extract to extract the actual matched part
        result <- str_extract(string, "GN=[A-Za-z0-9]+")
        # Further extract the part after 'GN='
        gene_name <- str_extract(result, "[A-Za-z0-9]+$")
        return(gene_name)
    } else {
        return("Pattern not found")
    }
})

rajagopal_cxcl2 <- fns_tbl %>%
    filter(owner == "Rajagopal") %>%
    dplyr::select(owner, data) %>%
    unnest(data) %>%
    dplyr::select(owner, Accession, Description, CXCL12_v_Veh_FC, CXCL12_v_Veh_pval) %>%
    mutate(
        gene_symbol = extract_gn_from_description_raja(Description),
        gene_symbol = ifelse(gene_symbol == "Pattern not found", Accession, gene_symbol)
    ) %>%
    mutate(
        tx = "veh vs cxcl12",
        timing = "3min"
    ) %>%
    group_by(tx) %>%
    rename(log2_fc = CXCL12_v_Veh_FC, p_val = CXCL12_v_Veh_pval) %>%
    mutate(
        adjusted_p_val = p.adjust(p_val, method = "fdr"),
        reg_fc = 2^log2_fc
    ) %>%
    dplyr::select(-Description) %>%
    dplyr::select(owner, Accession, gene_symbol, everything())

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
    separate_longer_delim(`Gene Names`, " ")
# mutate(
#     gene_symbol_all = map(.x = `Gene Names`, .f = function(x) as.list(str_split(x, pattern = " ", simplify = TRUE))),
#     gene_symbol = map_chr(.x = gene_symbol_all, .f = function(x) x[[1]])
# ) %>%
# rename(Accession = Entry) %>%
#     distinct() %>%
#     arrange(Accession) %>%
#     mutate(owner = "UniProt")

helper_accessions <- kruse_ang2_2min_21min %>%
    distinct(Accession, gene_symbol, owner) %>%
    arrange(gene_symbol)

helper_accessions2 <- zastro_1min_3min_10min %>%
    distinct(Accession, gene_symbol, owner) %>%
    arrange(gene_symbol)

helper_accessions3 <- rajagopal_cxcl2 %>%
    distinct(Accession, gene_symbol, owner) %>%
    arrange(gene_symbol)

helper_accessions4 <- zastro_1min_3min_10min %>%
    distinct(Accession, gene_symbol, owner) %>%
    arrange(gene_symbol)


my_gene_map <- bind_rows(my_gene_map_temp, helper_accessions, helper_accessions2, helper_accessions3, helper_accessions4) %>%
    distinct(Accession, gene_symbol) %>%
    arrange(gene_symbol) %>%
    filter(gene_symbol != "-" & !is.na(gene_symbol)) %>%
    group_by(Accession) %>%
    summarize(list_gene_symbol = list(gene_symbol), chr_gs = str_c(gene_symbol, collapse = "_")) %>%
    mutate(len = map_int(list_gene_symbol, .f = length)) %>%
    mutate(gene_symbol = chr_gs) %>%
    mutate(gene_symbol = factor(gene_symbol), Accession = factor(Accession))

# my_gene_map %>% filter(gene_symbol == "AP3D1")
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
    distinct()

combined_dat_temp <- combined_dat_full %>%
    filter((owner == "Rajagopal" & timing == "3min") |
        (owner == "Kruse" & timing == "2min") |
        (owner == "Rockman" & timing == "90sec") |
        owner == "Von Zastrow" & timing == "3min") %>%
    filter(!is.na(log2_fc)) %>% # c("1min", "3min", "10min") %>%
    filter(!is.na(Accession)) %>%
    mutate(row_number = row_number(), .before = 1)

remove_dups_accession <- combined_dat_temp %>%
    # there are some duplicates so we keep the one with the lowest p-value
    group_by(Accession, tx, owner) %>%
    mutate(dups = max(n())) %>%
    filter(dups >= 2) %>%
    group_by(Accession, owner) %>%
    arrange(Accession) %>%
    mutate(which_best = ifelse(!is.na(p_val),
        min(p_val) == p_val,
        max(abs(log2_fc)) == log2_fc
    )) %>%
    filter(!which_best)
remove_dups_accession

combined_dat <- combined_dat_temp %>%
    filter(!(row_number %in% remove_dups_accession$row_number))

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
combined_dat_full_no_extreme_outliers <- combined_dat %>%
    filter(!(row_number %in% extreme_outliers_all$row_number))
message(qq("Final Dataset size without EXTREME outliers: @{nrow(combined_dat_full_no_extreme_outliers)}"))

fn <- file.path("results", "plots", "data_stats", "combined_dat_full_no_extreme_outliers_hist.png")
png(fn, width = 1920, height = 1920, res = 300)
hist(combined_dat_full_no_extreme_outliers$log2_fc, main = "Combined Data (no EXTREME outliers)", xlab = "Log2 Fold Change")
dev.off()

# Let's go with the the full data
message("Using combined_dat_full_no_extreme_outliers")
combined_dat <- combined_dat_full_no_extreme_outliers

fn <- file.path("results", "plots", "data_stats", "p_value_distr", "raj_pval_distr.png")
png(fn, width = 1920, height = 1920, res = 300)
raja_dstr_plot <- plotp(combined_dat %>% filter(owner == "Rajagopal") %>% .$p_val,
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
    # plotp(combined_dat %>% filter(owner == "Kruse") %>% .$p_val,
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

save_pval_distr(data = combined_dat)

message("Writing initial dataset to cache dir...")
write_rds(combined_dat, file = file.path("cache", "initial-dataset_combined-dat.rds"))
