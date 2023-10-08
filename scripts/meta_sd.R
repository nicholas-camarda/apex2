# ----metamsd, echo = FALSE---------------------------------------------------------------------------------------------------
# Code to generate datasets usable in MetaMSD
### START #######
combined_dat <- read_rds(file.path("cache", "initial-dataset_combined-dat.rds"))
{
    to_write <- tibble(datasets = map(.x = combined_dat %>% split(.$owner), .f = function(d) {
        d %>%
            dplyr::select(Protein = Accession, Sign = log2_fc, Pvalue = adjusted_p_val) %>%
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
            dplyr::select(Protein = Accession, Sign = log2_fc, Pvalue = adjusted_p_val) %>%
            filter(!(Protein %in% check_na_proteins$Protein))
    })) %>%
        mutate(
            name_data = names(datasets),
            output_fns = file.path(getwd(), "../MetaMSD/MetaMSD/input", str_c(name_data, "_no_na.txt"))
        )

    walk2(.x = as.list(to_write$datasets), .y = as.list(to_write$output_fns), .f = function(x, y) write_tsv(x = x, file = y))
    walk2(.x = as.list(to_write_no_na$datasets), .y = as.list(to_write_no_na$output_fns), .f = function(x, y) write_tsv(x = x, file = y))
}
### END #######
