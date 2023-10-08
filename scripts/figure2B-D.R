# rm(list = ls())
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

source(file.path("scripts", "figure2_helper_functions.R"))

make_figure2_plots <- function(combined_dat, results_path = "", stat_value_type = "q_value", stat_cutoff = 0.1, fc_cutoff = 2, verbose = FALSE) {
    dir.create(results_path, showWarnings = FALSE, recursive = TRUE)
    helper_fn_title <- basename(results_path)
    helper_fn_name <- convert_to_filename_friendly(helper_fn_title)


    max_datasets <- length(unique(combined_dat$owner))
    message(qq("Working with @{max_datasets} datasets."))
    cat(unique(combined_dat$owner), sep = " | ")
    cat("\n")

    message("Generating combined p-values and correcting for multiple hypotheses...")
    cmbd_p_vals_all_completed_data <- generate_cmbd_fisher_pval_data(
        data = combined_dat
    )
    message("Generating separate pvalues without Rajagopal contribution...")
    cmbd_p_vals_all_completed_data_no_raja <- generate_cmbd_fisher_pval_data(
        data = combined_dat %>% filter(owner != "Rajagopal")
    )
    message("Done!")

    # for (stat_value_type in c("p_value", "q_value")) {
    # results_path <- file.path(results_path, stat_value_type)
    # results_path <- file.path(results_path)
    # dir.create(results_path, showWarnings = FALSE, recursive = TRUE)

    lst_res <- make_sig_dfs(
        data1 = cmbd_p_vals_all_completed_data,
        data2 = cmbd_p_vals_all_completed_data_no_raja,
        my_stat_value_type = stat_value_type,
        my_stat_cutoff = stat_cutoff
    )
    sig_fisher_df <- lst_res[[1]]
    sig_fisher_no_raja_df <- lst_res[[2]]
    stat_cutoff <- lst_res[[3]]

    if (verbose) {
        message("Sig DFs")
        message("Full")
        print(sig_fisher_df)
        message("No Raj")
        print(sig_fisher_no_raja_df)
        message("Stat cutoff:")
        print(stat_cutoff)
    }


    message(qq("Init waterfall `full` df size: @{dim(sig_fisher_df)[1]}, @{dim(sig_fisher_df)[2]}"))
    message(qq("Init waterfall df NO RAJA size: @{dim(sig_fisher_no_raja_df)[1]}, @{dim(sig_fisher_no_raja_df)[2]}"))

    waterfall_plot_name_lst <- make_naming_tribble_waterfall(
        data = sig_fisher_df,
        data2 = sig_fisher_no_raja_df,
        my_stat_value_type = stat_value_type,
        num_full_set = 2,
        num_conserved = max_datasets,
        helper_title = helper_fn_title
    )


    message("Generating water fall plots!")
    wf_path <- file.path(results_path, "2CE. waterfall")
    dir.create(wf_path, showWarnings = FALSE, recursive = TRUE)

    for (i in seq_len(nrow(waterfall_plot_name_lst))) {
        message(i)
        if (verbose) print(waterfall_plot_name_lst[i, ])
        waterfall_plot_lst_info <- plot_fc(
            df = waterfall_plot_name_lst$data[[i]][[1]],
            grouping_variable = waterfall_plot_name_lst$gvar[[i]],
            fc_cutoff = fc_cutoff,
            stat_cutoff_val = stat_cutoff,
            my_stat_value_type = stat_value_type,
            title_ = waterfall_plot_name_lst$title[[i]],
            found_in_title = waterfall_plot_name_lst$found_in[[i]]
        )
        num_analytes_plt <- waterfall_plot_lst_info[[2]]

        num_increments <- num_analytes_plt %/% 50
        my_width <- 11 * (1 + 0.5 * num_increments)
        my_height <- 8.5 * (1 + 0.15 * num_increments)

        fn_name <- add_string_to_filename_basic(waterfall_plot_name_lst$path[[i]], string_to_add = helper_fn_name)
        ggsave(
            plot = waterfall_plot_lst_info[[1]],
            filename = file.path(wf_path, fn_name),
            width = my_width, height = my_height
        )
    }

    volcano_plot_name_lst <- make_naming_tribble_volcano(
        data1 = sig_fisher_df,
        data2 = sig_fisher_no_raja_df,
        my_stat_value_type = stat_value_type,
        num_full_set = 2,
        num_conserved = max_datasets,
        helper_title = helper_fn_title
    )

    message("Generating volcano plots!")
    volcano_path <- file.path(results_path, "2BD. volcano")
    dir.create(volcano_path, showWarnings = FALSE, recursive = TRUE)
    for (i in seq_len(nrow(volcano_plot_name_lst))) {
        message(i)
        if (verbose) print(volcano_plot_name_lst[i, ])
        volcano_plot_lst_info <- plot_volcano(
            df = volcano_plot_name_lst$data[[i]][[1]],
            stat_cutoff_val = stat_cutoff,
            my_stat_value_type = stat_value_type,
            fc_thresh = fc_cutoff,
            label_variable = volcano_plot_name_lst$gvar[[i]],
            title_ = volcano_plot_name_lst$title[[i]],
            found_in_title = volcano_plot_name_lst$found_in[[i]]
        )

        num_analytes_vplt <- volcano_plot_lst_info[[2]]

        num_v_increments <- num_analytes_vplt %/% 50
        my_v_width <- 10 * (1 + 0.15 * num_v_increments)
        my_v_height <- 10 * (1 + 0.15 * num_v_increments)
        fn_name <- add_string_to_filename_basic(volcano_plot_name_lst$path[[i]], string_to_add = helper_fn_name)

        ggsave(
            plot = volcano_plot_lst_info[[1]],
            filename = file.path(volcano_path, fn_name),
            width = my_v_width, height = my_v_height
        )
    }
    message("All done!")
    # }

    message("Caching data for later...")
    complete_cache_fn_name <- file.path("cache", add_string_to_filename_basic("FINAL-cache.rds", string_to_add = helper_fn_name, file_ext = "rds"))
    no_raja_cache_fn_name <- file.path("cache", add_string_to_filename_basic("FINAL-cache_no_raja.rds", string_to_add = helper_fn_name, file_ext = "rds"))

    write_rds(sig_fisher_df %>% mutate(stat_value_type, stat_value_cutoff = stat_cutoff), file = complete_cache_fn_name)
    write_rds(sig_fisher_no_raja_df %>% mutate(stat_value_type, stat_value_cutoff = stat_cutoff), file = no_raja_cache_fn_name)
    cat("\n\n")
}

my_combined_dat <- read_rds(file.path("cache", "initial-dataset_combined-dat.rds"))

make_figure2_plots(
    combined_dat = my_combined_dat,
    results_path = file.path("results", "plots", "Figure 2", "Strict Thresholds for Fig 2 (q0_1 and FC2)", "2B-E", "Including Kruse"),
    stat_value_type = "q_value",
    stat_cutoff = 0.1,
    fc_cutoff = 2,
    verbose = FALSE
)

make_figure2_plots(
    combined_dat = my_combined_dat %>% filter(owner != "Kruse"),
    results_path = file.path("results", "plots", "Figure 2", "Strict Thresholds for Fig 2 (q0_1 and FC2)", "2B-E", "Excluding Kruse"),
    stat_value_type = "q_value",
    stat_cutoff = 0.1,
    fc_cutoff = 2,
    verbose = FALSE
)

make_figure2_plots(
    combined_dat = my_combined_dat,
    results_path = file.path("results", "plots", "Figure 2", "Relaxed Thresholds for Fig 2 (q0_2 and FC1_33)", "2B-E", "Including Kruse"),
    stat_value_type = "q_value",
    stat_cutoff = 0.2,
    fc_cutoff = 1.33,
    verbose = FALSE
)

make_figure2_plots(
    combined_dat = my_combined_dat %>% filter(owner != "Kruse"),
    results_path = file.path("results", "plots", "Figure 2", "Relaxed Thresholds for Fig 2 (q0_2 and FC1_33)", "2B-E", "Excluding Kruse"),
    stat_value_type = "q_value",
    stat_cutoff = 0.2,
    fc_cutoff = 1.33,
    verbose = FALSE
)


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
