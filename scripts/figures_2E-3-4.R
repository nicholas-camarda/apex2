rm(list = ls())
# library(ReactomePA)
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
library(ggprism)

source(file.path("scripts", "figure4_helper_functions.R"))
handlers(list(
    handler_progress(
        format   = ":spin :current/:total [:bar] :percent in :elapsed ETA: :eta", # | (:message)
        width    = 60,
        complete = "+"
    )
))


# My variables
all_data <- read_rds(file.path("cache", "FINAL-cache-including_kruse.rds"))
conserved_data <- all_data %>% filter(n_owners == 4)


# STRICT
# make figure 3
my_vector_genes_lst <- make_figure3(
    data = conserved_data,
    func_base_results_path = file.path("results", "plots", "Figures 3 and 4", qq("Strict Thresholds for Fig 3 and 4 (q0_1 and FC2)")),
    pathway_p_adjust_thresh = 0.1,
    go_pathway_fc_threshold = 2
)
my_vector_genes <- my_vector_genes_lst[[1]]

# make figures 3E and 4
with_progress(do_stringdb_analysis(
    data = conserved_data,
    func_base_results_path = file.path("results", "plots", "Figures 3 and 4", qq("Strict Thresholds for Fig 3 and 4 (q0_1 and FC2)")),
    string_db_q_thresh = 0.1,
    string_db_fc_threshold = 2,
    my_input_intervals = c(10, 20, 30, 40, 50)
))


# make figure 3
# relaxed
my_vector_genes_lst <- make_figure3(
    data = conserved_data,
    func_base_results_path = file.path("results", "plots", "Figures 3 and 4", qq("Relaxed Thresholds for Fig 3 and 4 (q0_2 and FC1_33)")),
    pathway_p_adjust_thresh = 0.2,
    go_pathway_fc_threshold = 1.33
)

# make figures 3E and 4
with_progress(do_stringdb_analysis(
    data = conserved_data,
    func_base_results_path = file.path("results", "plots", "Figures 3 and 4", qq("Relaxed Thresholds for Fig 3 and 4 (q0_2 and FC1_33)")),
    string_db_q_thresh = 0.2,
    string_db_fc_threshold = 1.33,
    my_input_intervals = c(10, 20, 30, 40, 50)
))




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
# set_and_print_variables()
