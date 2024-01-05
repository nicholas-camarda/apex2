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


# Raja only
raja_result_path <- file.path("results", "plots", "Figure 2", "Strict Thresholds for Fig 2 (q0_1 and FC2)", "2B-E", "Rajagopal Only")
waterfall_plot_easy_lst <- easy_plot_fc(
    dat1 = my_combined_dat %>% filter(owner == "Rajagopal"),
    fc_cutoff = 2,
    q_val_cutoff = 0.1,
    title_ = "Rajagopal only"
)
ggsave(
    filename = file.path(raja_result_path, "Rajagopal_waterfall.png"),
    plot = waterfall_plot_easy_lst[[1]],
    width = 20,
    height = 12
)
openxlsx::write.xlsx(waterfall_plot_easy_lst[[2]] %>% dplyr::select(-owners, -n_owners), file = file.path(raja_result_path, "raja_only_waterfall.xlsx"))


volcano_plot_easy_lst <- easy_plot_volcano(
    df = my_combined_dat %>% filter(owner == "Rajagopal"),
    fc_thresh = 2,
    q_val_cutoff = 0.1,
    title_ = "Rajagopal only"
)
ggsave(
    filename = file.path(raja_result_path, "Rajagopal_volcano.png"),
    plot = volcano_plot_easy_lst[[1]],
    width = 20,
    height = 12
)
openxlsx::write.xlsx(volcano_plot_easy_lst[[2]] %>% dplyr::select(-owners, -n_owners), file = file.path(raja_result_path, "raja_only_volcano-all.xlsx"))



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
