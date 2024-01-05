#' Generates combined Fisher p-values for a dataset.
#'
#' This function computes combined Fisher p-values for each unique 'Accession' in the dataset.
#' It handles two-sidedness and includes data adjustments using the 'metap' package.
#'
#' @param data The dataset for which combined Fisher p-values are to be calculated.
#' @param verbose A logical indicating whether to print additional messages (default is FALSE).
#' @return A dataframe with combined Fisher p-values and adjusted values.
#' @examples
#' data <- read.csv("path/to/data.csv")
#' generate_cmbd_fisher_pval_data(data)
generate_cmbd_fisher_pval_data <- function(data, verbose = FALSE) {
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

    # this operation is always performed on unadjusted p-values
    onesided_fisher <- two2one(p = input_to_stat$p_val, two = input_to_stat$two_sidedness, invert = input_to_stat$toinvert)
    if (verbose) {
        message("Removing NA fisher p-vals\n(i.e. n < 2 or any NA p-values used to construct the combined p-value)")
        message("Adding in Kruse data with no p-values after Fisher combined-pvalue step")
    }

    fisher_df <- input_to_stat %>%
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
        mutate(
            cmbd_fisher_p_val = ifelse(is.na(cmbd_fisher_p_val), dplyr::first(cmbd_fisher_p_val), cmbd_fisher_p_val),
            adjusted_cmbd_fisher_p_val = ifelse(is.na(adjusted_cmbd_fisher_p_val), dplyr::first(adjusted_cmbd_fisher_p_val), adjusted_cmbd_fisher_p_val)
        ) %>%
        ungroup() %>%
        suppressWarnings()

    # get all the log2_fc data for those analytes, since we removed some rows with NA p-values above
    return(fisher_df)
}

#' Counts the number of unique 'owners' for each Accession in the dataset.
#'
#' This function counts the occurrences of each unique Accession across different 'owners' in the dataset.
#'
#' @param sig_data The dataset containing Accession and owner columns.
#' @return A dataframe with each Accession, the number of owners, and a concatenated string of owners.
#' @examples
#' sig_data <- read.csv("significant_data.csv")
#' accession_counts <- count_owners(sig_data)
count_owners <- function(sig_data) {
    count_accessions <- sig_data %>%
        distinct(Accession, owner) %>%
        group_by(Accession) %>%
        summarize(
            n_owners = max(n()),
            owners = str_c(sort(owner), collapse = "|"), .groups = "drop"
        ) %>%
        arrange(desc(n_owners))

    return(count_accessions)
}


#' Creates dataframes for significant findings based on statistical values.
#'
#' This function processes two datasets to flag significant findings based on a specified statistical value
#' (e.g., q-value or p-value) and its cutoff. It returns dataframes indicating significant findings.
#'
#' @param data1 First dataset for analysis.
#' @param data2 Second dataset for analysis.
#' @param my_stat_value_type Type of statistical value used for significance ('q_value' or 'p_value').
#' @param my_stat_cutoff Cutoff value for determining significance.
#' @return A list of two dataframes with significant findings and the statistical cutoff used.
#' @examples
#' data1 <- read.csv("dataset1.csv")
#' data2 <- read.csv("dataset2.csv")
#' sig_dfs <- make_sig_dfs(data1, data2, "q_value", 0.05)
make_sig_dfs <- function(data1, data2, my_stat_value_type, my_stat_cutoff = NA) {
    if (my_stat_value_type == "q_value") {
        which_p <- sym("adjusted_cmbd_fisher_p_val")
    } else {
        which_p <- sym("cmbd_fisher_p_val")
    }

    sig_fisher_df_temp <- data1 %>%
        mutate(
            use_stat_val = !!which_p,
            use_neg_log10_stat_val = -log10(use_stat_val),
            significant_stat_val = ifelse(!is.na(use_stat_val), use_stat_val < my_stat_cutoff, FALSE)
        )

    sig_fisher_df_no_raja_temp <- data2 %>%
        mutate(
            use_stat_val = !!which_p,
            use_neg_log10_stat_val = -log10(use_stat_val),
            significant_stat_val = ifelse(!is.na(use_stat_val), use_stat_val < my_stat_cutoff, FALSE)
        )
    stat_cutoff <- my_stat_cutoff

    # print(sig_fisher_df_temp)
    # print(sig_fisher_df_no_raja_temp)
    # print(stat_cutoff)

    counts_sig_df <- count_owners(sig_fisher_df_temp)
    counts_no_raja_sig_df <- count_owners(sig_fisher_df_no_raja_temp)

    sig_fisher_df <- sig_fisher_df_temp %>% left_join(counts_sig_df, by = join_by(Accession))
    sig_fisher_no_raja_df <- sig_fisher_df_no_raja_temp %>% left_join(counts_no_raja_sig_df, by = join_by(Accession))

    # print(sig_fisher_df)
    # print(count_accessions)
    # print(sig_fisher_no_raja_df)
    # print(count_accessions_no_raja)
    return(list(sig_fisher_df, sig_fisher_no_raja_df, stat_cutoff))
}

#' Creates a tribble for naming waterfall plots.
#'
#' This function prepares a tibble (tribble) with configurations for generating waterfall plots.
#' It includes specific settings for different subsets of the data, such as full dataset or conserved hits,
#' and differentiates between Accession and gene_symbol as grouping variables.
#'
#' @param data The primary dataset used for waterfall plotting.
#' @param data2 An alternative dataset used for specific waterfall plots.
#' @param my_stat_value_type The statistical value type ('p_value' or 'q_value') used in the plots.
#' @param num_full_set The number of full datasets to consider for plotting (default is 2).
#' @param num_conserved The number for conserved datasets to consider (default is 4).
#' @param helper_title An additional title to add to the plot titles (default is an empty string).
#' @return A tribble with configurations for different waterfall plots.
#' @examples
#' data <- read.csv("path/to/data.csv")
#' data2 <- read.csv("path/to/data2.csv")
#' tribble_waterfall <- make_naming_tribble_waterfall(data, data2, "q_value")
make_naming_tribble_waterfall <- function(data, data2, my_stat_value_type, num_full_set = 2, num_conserved = 4, helper_title = "") {
    naming_tribble <- tribble(
        ~data, ~found_in, ~gvar, ~title, ~path,
        list(data %>% filter(n_owners >= num_full_set)), qq("num_ds >= @{num_full_set}"), "Accession", qq("Waterfall plot of significant analytes by Fisher @{my_stat_value_type}\n@{helper_title}"), qq("All_data-accession_fc_waterfall-@{my_stat_value_type}.png"),
        list(data %>% filter(n_owners >= num_full_set)), qq("num_ds >= @{num_full_set}"), "gene_symbol", qq("Waterfall plot of significant analytes by Fisher @{my_stat_value_type}\n@{helper_title}"), qq("All_data-gene_symbol_fc_waterfall-@{my_stat_value_type}.png"),
        list(data %>% filter(n_owners == num_conserved)), qq("num_ds == @{num_conserved}"), "Accession", qq("2C. Waterfall plot of significant analytes by Fisher @{my_stat_value_type}\n@{helper_title} (Conserved hits)"), qq("2C. Conserved_data-accession_fc_waterfall-@{my_stat_value_type}.png"),
        list(data %>% filter(n_owners == num_conserved)), qq("num_ds == @{num_conserved}"), "gene_symbol", qq("2C. Waterfall plot of significant analytes by Fisher @{my_stat_value_type}\n@{helper_title} (Conserved hits)"), qq("2C. Conserved_data-gene_symbol_fc_waterfall-@{my_stat_value_type}.png"),
        list(data2), qq("num_ds == @{num_conserved-1}"), "Accession", qq("2E. Waterfall plot of significant analytes by Fisher @{my_stat_value_type}\n@{helper_title} and Excluding Rajagopal)"), qq("2E. No_Raj-accession_waterfall-@{my_stat_value_type}.png"),
        list(data2), qq("num_ds == @{num_conserved-1}"), "gene_symbol", qq("2E. Waterfall plot of significant analytes by Fisher @{my_stat_value_type}\n@{helper_title} and Excluding Rajagopal)"), qq("2E. No_Raj-gene_symbol_waterfall-@{my_stat_value_type}.png")
    )
    return(naming_tribble)
}

#' Creates a tribble for naming volcano plots.
#'
#' This function generates a tibble (tribble) with configurations for creating volcano plots.
#' It sets up parameters for different data subsets and grouping variables, allowing for flexible
#' and targeted volcano plot generation.
#'
#' @param data1 The first dataset for volcano plot generation.
#' @param data2 The second dataset for specific volcano plot settings.
#' @param my_stat_value_type The statistical value type ('p_value' or 'q_value') used in the plots.
#' @param num_full_set The number of full datasets to include (default is 2).
#' @param num_conserved The number of conserved datasets to include (default is 4).
#' @param helper_title Additional title text for the plots (default is an empty string).
#' @return A tribble with configurations for different volcano plots.
#' @examples
#' data1 <- read.csv("path/to/data1.csv")
#' data2 <- read.csv("path/to/data2.csv")
#' tribble_volcano <- make_naming_tribble_volcano(data1, data2, "q_value")
make_naming_tribble_volcano <- function(data1, data2, my_stat_value_type, num_full_set = 2, num_conserved = 4, helper_title = "") {
    naming_tribble <- tribble(
        ~data, ~found_in, ~gvar, ~title, ~path,
        list(data1 %>% filter(n_owners >= num_full_set)), qq("num_ds >= @{num_full_set}"), "Accession", qq("Volcano plot of significant analytes by Fisher @{my_stat_value_type}\n@{helper_title}"), qq("All_data-accession_volcano-@{my_stat_value_type}.png"),
        list(data1 %>% filter(n_owners >= num_full_set)), qq("num_ds >= @{num_full_set}"), "gene_symbol", qq("Volcano plot of significant analytes by Fisher @{my_stat_value_type}\n@{helper_title}"), qq("All_data-gene_symbol_volcano-@{my_stat_value_type}.png"),
        list(data1 %>% filter(n_owners == num_conserved)), qq("num_ds == @{num_conserved}"), "Accession", qq("2B. Volcano plot of significant analytes by Fisher @{my_stat_value_type}\n@{helper_title} (Conserved hits)"), qq("2B. Conserved_data-accession_volcano-@{my_stat_value_type}.png"),
        list(data1 %>% filter(n_owners == num_conserved)), qq("num_ds == @{num_conserved}"), "gene_symbol", qq("2B. Volcano plot of significant analytes by Fisher @{my_stat_value_type}\n@{helper_title} (Conserved hits)"), qq("2B. Conserved_data-gene_symbol_volcano-@{my_stat_value_type}.png"),
        list(data2), qq("num_ds == @{num_conserved-1}"), "Accession", qq("2D. Volcano plot of significant analytes by Fisher @{my_stat_value_type}\n(@{helper_title} and Excluding Rajagopal)"), qq("2D. No_Raj-accession_volcano-@{my_stat_value_type}.png"),
        list(data2), qq("num_ds == @{num_conserved-1}"), "gene_symbol", qq("2D. Volcano plot of significant analytes by Fisher @{my_stat_value_type}\n(@{helper_title} and Excluding Rajagopal)"), qq("2D. No_Raj-gene_symbol_volcano-@{my_stat_value_type}.png")
    )
    return(naming_tribble)
}


#' Creates and returns a waterfall plot for fold change analysis.
#'
#' This function generates a waterfall plot based on the fold change and statistical significance of
#' each gene or Accession in the provided dataset.
#'
#' @param df The dataset for plotting.
#' @param grouping_variable The variable to group data by (default is "Accession").
#' @param fc_cutoff Fold change cutoff for significance.
#' @param stat_cutoff_val Statistical value cutoff for significance.
#' @param my_stat_value_type Type of statistical value ('q_value' or 'p_value').
#' @param title_ Title of the plot.
#' @param found_in_title Part of the caption indicating dataset inclusion criteria.
#' @return A ggplot object representing the waterfall plot.
#' @examples
#' df <- read.csv("gene_data.csv")
#' waterfall_plot <- plot_fc(df)
plot_fc <- function(df, grouping_variable = "Accession", fc_cutoff = 2, stat_cutoff_val = 0.001, my_stat_value_type = "q_value", title_ = "Waterfall plot", found_in_title = ">= 2") {
    # DEBUG
    # i = 1; df = waterfall_plot_name_lst$data[[i]][[1]];grouping_variable = waterfall_plot_name_lst$gvar[[i]]; stat_cutoff_val = stat_cutoff; title_ = waterfall_plot_name_lst$title[[i]]; fc_cutoff = my_fc_cutoff
    grouping_var_sym <- sym(grouping_variable)

    p1_plot_dat <- df %>%
        group_by(!!grouping_var_sym) %>%
        reframe(
            # mean_log2_fc = mean(log2_fc, na.rm = TRUE),
            median_log2_fc = median(log2_fc, na.rm = TRUE),
            log2_fc = list(log2_fc),
            use_stat_val = unique(use_stat_val),
            significant = unique(significant_stat_val),
            n_owners = unique(n_owners),
            owners = unique(owners)
        ) %>%
        ## threshold on MEAN FC also
        filter(abs(median_log2_fc) > log2(fc_cutoff) & significant) %>%
        arrange(!!grouping_var_sym) %>%
        # these two arrange calls HAVE to be separate
        arrange(desc(median_log2_fc)) %>%
        mutate(order = row_number())


    if (grouping_variable == "gene_symbol") {
        p1_plot_dat <- p1_plot_dat %>%
            mutate(
                gene_symbol = factor(gene_symbol, levels = gene_symbol[order]),
                up_down = factor(ifelse(median_log2_fc > 0, "++", "--"),
                    levels = c("++", "--")
                )
            )
        grouping_variable <- "Gene Symbol"
    } else {
        p1_plot_dat <- p1_plot_dat %>%
            mutate(
                Accession = factor(Accession, levels = Accession[order]),
                up_down = factor(ifelse(median_log2_fc > 0, "++", "--"),
                    levels = c("++", "--")
                )
            )
    }

    owners_labeller_df <- df %>%
        left_join(
            p1_plot_dat %>%
                dplyr::select(!!grouping_var_sym, median_log2_fc, order, up_down, significant),
            by = join_by(!!grouping_var_sym)
        ) %>%
        ## threshold on FC also
        filter(abs(median_log2_fc) > log2(fc_cutoff) & significant) %>%
        dplyr::select(!!grouping_var_sym, owner, log2_fc, median_log2_fc, order, up_down) %>%
        arrange(order)

    num_analytes <- nrow(p1_plot_dat)

    plot_dat <- p1_plot_dat %>% unnest(log2_fc)

    helper_found_in_title_caption <- str_split(found_in_title, "num_ds ", simplify = TRUE)[, 2]
    wf_plot_with_kruse <- ggplot(plot_dat,
        mapping = aes(
            x = !!grouping_var_sym, y = log2_fc,
            color = up_down, group = !!grouping_var_sym
        )
    ) +
        geom_hline(yintercept = 0, color = "orange") +
        geom_hline(yintercept = c(-log2(fc_cutoff), log2(fc_cutoff)), color = "#9700e8", alpha = 0.6, linewidth = 1.1, linetype = 4) +
        # boxplot first
        geom_boxplot() +
        # then jitters
        geom_jitter(owners_labeller_df,
            mapping = aes(x = !!grouping_var_sym, y = log2_fc, shape = owner, fill = up_down),
            color = "black",
            width = 0.1, size = rel(2),
            inherit.aes = FALSE
        ) +
        theme_prism(base_size = 20) +
        theme(
            strip.text.x = element_text(size = rel(1.25), face = "bold"),
            axis.text.x = element_text(size = rel(0.7), angle = 90, hjust = 1), # , angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = rel(0.7)),
            panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = 0.5),
            panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = 0.25)
        ) +
        scale_fill_manual(values = c("#474747", "red")) +
        scale_color_manual(values = c("#474747", "red")) +
        scale_shape_manual(values = c("Rajagopal" = 21, "Kruse" = 22, "Von Zastrow" = 23, "Rockman" = 24)) +
        ggtitle(qq("@{title_}")) +
        labs(caption = qq("@{my_stat_value_type} cuttoff: @{stat_cutoff_val}\n abs(raw median(FC)) > @{fc_cutoff}\nFound in @{helper_found_in_title_caption} dataset(s)\nN = @{num_analytes}")) +
        ylab(TeX("$\\bf{log_2(FC)}$")) +
        labs(x = grouping_variable)

    wf_plot_with_kruse

    return(list(wf_plot_with_kruse, num_analytes))
}

#' Generates a volcano plot for the given dataset.
#'
#' This function creates a volcano plot based on log fold change and negative log10 statistical values.
#' Significant points are highlighted and labeled.
#'
#' @param df The dataset for the volcano plot.
#' @param stat_cutoff_val Statistical cutoff value for significance.
#' @param label_variable The variable used for labeling significant points.
#' @param my_stat_value_type Type of statistical value ('q_value' or 'p_value').
#' @param fc_thresh Fold change threshold for significance.
#' @param found_in_title Part of the caption indicating dataset inclusion criteria.
#' @param title_ Title of the volcano plot.
#' @return A ggplot object representing the volcano plot.
#' @examples
#' df <- read.csv("gene_data.csv")
#' volcano_plot <- plot_volcano(df, 0.001, "Accession", "q_value", 2)
plot_volcano <- function(df, stat_cutoff_val = 0.001, label_variable = "Accession", my_stat_value_type = "q_value", fc_thresh = 2, found_in_title = "num_ds >= 2", title_ = "Volcano plot") {
    # DEBUG
    # i = 6; df = volcano_plot_name_lst$data[[i]][[1]];label_variable = volcano_plot_name_lst$gvar[[i]];fc_thresh = my_fc_cutoff;stat_cutoff_val = stat_cutoff;my_stat_value_type = stat_value_type;title_ = volcano_plot_name_lst$title[[i]];found_in_title = volcano_plot_name_lst$found_in[[i]]

    label_var_sym <- sym(label_variable)
    counted_owners_df <- count_owners(df)

    if (my_stat_value_type == "q_value") {
        my_ylab <- TeX("$\\bf{-log_{10}(Q)}$")
    } else {
        my_ylab <- TeX("$\\bf{-log_{10}(P)}$")
    }


    plot_df_temp <- df %>%
        # left_join(counted_owners_df, by = join_by(Accession)) %>%
        group_by(!!label_var_sym) %>%
        mutate(
            # mean_log2_fc = mean(log2_fc, na.rm = TRUE),
            median_log2_fc = median(log2_fc, na.rm = TRUE),
            significant_stat = unique(significant_stat_val),
        ) %>%
        mutate(significant = ifelse(abs(median_log2_fc) > log2(fc_thresh) & significant_stat, TRUE, FALSE))
    max_for_alpha <- max(sqrt(plot_df_temp$median_log2_fc^2 + plot_df_temp$use_neg_log10_stat_val^2))

    plot_df <- plot_df_temp %>%
        mutate(
            res = sqrt(median_log2_fc^2 + use_neg_log10_stat_val^2),
            significant_alpha = ifelse(significant,
                (res / max_for_alpha)^1.5,
                0
            ), .before = 1 # 1, 0.01)
        ) %>%
        distinct(
            .data[[label_variable]], median_log2_fc,
            use_neg_log10_stat_val,
            significant, significant_alpha, n_owners, owners
        )

    label_df <- plot_df %>%
        filter(significant)

    num_analytes <- nrow(label_df %>% distinct(.data[[label_variable]]))
    if (label_variable == "gene_symbol") label_variable <- "Gene Symbol"
    helper_found_in_title_caption <- str_split(found_in_title, "num_ds ", simplify = TRUE)[, 2]


    volcano_plot <- ggplot(plot_df, aes(
        x = median_log2_fc,
        y = use_neg_log10_stat_val,
    ), color = "black") +
        geom_hline(yintercept = -log10(stat_cutoff_val), color = "orange") +
        geom_vline(xintercept = c(-log2(fc_thresh), log2(fc_thresh)), color = "#370afc", alpha = 0.6, linewidth = 1.1, linetype = 4) +
        geom_point(
            mapping = aes(
                fill = significant,
                alpha = significant_alpha
            ),
            pch = 21,
            show.legend = FALSE, size = 2.5
        ) +
        xlab(TeX("$\\bf{median(log_2(FC))}$")) +
        ylab(my_ylab) +
        scale_fill_manual(values = c("black", "red")) +
        geom_label_repel(label_df,
            mapping = aes(label = !!label_var_sym), force = 10, max.overlaps = 20, min.segment.length = 0, show.legend = FALSE
        ) +
        theme_prism(base_size = 20) +
        theme(
            strip.text.x = element_text(size = rel(1.25), face = "bold"),
            axis.text = element_text(size = rel(0.8)), # , angle = 45, hjust = 1, vjust = 1),
            panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = 0.5),
            panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = 0.25)
        ) +
        ggtitle(qq("@{title_}")) +
        labs(caption = qq("@{my_stat_value_type} cuttoff: @{stat_cutoff_val}\n abs(raw median(FC)) > @{fc_thresh}\nFound in @{helper_found_in_title_caption} dataset(s)\nN = @{num_analytes}"))
    volcano_plot
    return(list(volcano_plot, num_analytes))
}

#' Adds a specified string to a filename, maintaining the file extension.
#'
#' This function modifies a given filename by inserting an additional string before the file extension.
#'
#' @param filename The original filename.
#' @param string_to_add The string to be added to the filename.
#' @return A new filename with the added string.
#' @examples
#' filename <- "data_analysis.csv"
#' new_filename <- add_string_to_filename(filename, "updated")
add_string_to_filename <- function(filename, string_to_add) {
    # Usage example:
    # add_string_to_filename("myfile.txt", "_new")

    # Split the filename into the base name and the extension
    file_parts <- str_split(filename, "\\.", simplify = TRUE)[1, ]

    # Check if the filename has an extension
    if (length(file_parts) > 1) {
        # Combine the base name and the string to add, and then add the extension back
        new_filename <- paste0(file_parts[1], "-", string_to_add, ".", file_parts[2])
    } else {
        # If there's no extension, just add the string to the end of the filename
        new_filename <- paste0(filename, "-", string_to_add)
    }

    return(new_filename)
}

#' Adds a specified string to a filename, maintaining the file extension, for less complicated paths.
#' Less functionality, but simpler.
#'
#' This function modifies a given filename by inserting an additional string before the file extension.
#'
#' @param filename The original filename.
#' @param string_to_add The string to be added to the filename.
#' @return A new filename with the added string.
#' @examples
#' filename <- "data_analysis.csv"
#' new_filename <- add_string_to_filename(filename, "updated")
add_string_to_filename_basic <- function(filename, string_to_add, file_ext = "png") {
    new_str <- str_split(filename, qq("\\.@{file_ext}"), simplify = TRUE)[, 1]
    new_filename <- paste0(new_str, "-", string_to_add, ".", file_ext)
}


#' Converts a string to a filename-friendly format.
#'
#' This function transforms an input string to a lowercase, underscore-separated format,
#' suitable for filenames.
#'
#' @param input_string The string to be converted.
#' @return A filename-friendly version of the input string.
#' @examples
#' friendly_name <- convert_to_filename_friendly("Example String")
convert_to_filename_friendly <- function(input_string) {
    # Usage example:
    # convert_to_filename_friendly("Include Kruse data")
    # Convert the string to lowercase
    lowercase_string <- tolower(input_string)

    # Replace spaces with underscores
    filename_friendly_string <- gsub(" ", "_", lowercase_string)

    # Remove any other non-alphanumeric characters (except underscores)
    filename_friendly_string <- gsub("[^a-z0-9_]", "", filename_friendly_string)

    return(filename_friendly_string)
}

#' Creates a directory name based on a statistical value.
#'
#' This function generates a directory name using a provided statistical value and type.
#'
#' @param sval The statistical value.
#' @param type The type of statistical value ('q' for q-value, etc.).
#' @return A directory name string.
#' @examples
#' dir_name <- make_stat_directory(0.05, "p")
make_stat_directory <- function(sval, type = "q") {
    sval_str <- as.character(sval)
    sval_str <- gsub("\\.", "-", sval_str)
    dir_name <- paste0(type, sval_str)

    return(dir_name)
}



#' Creates a waterfall plot with less restrictive filtering.
#'
#' This function facilitates the generation of a waterfall plot for protein expression data.
#' It applies less stringent criteria for filtering, focusing on fold change and q-value thresholds.
#' The function also focusing only on data grouping and labeling based on gene symbols, rather than also on Accessions.
#'
#' @param dat1 The dataset containing gene expression data.
#' @param q_val_cutoff The q-value cutoff for significance (default is 0.1).
#' @param fc_cutoff The fold change cutoff for significance (default is 2).
#' @param title_ The title for the plot (default is an empty string).
#' @return A list containing the ggplot object for the waterfall plot and the final dataset used for plotting.
#' @examples
#' data <- read.csv("path/to/expression_data.csv")
#' waterfall_plot <- easy_plot_fc(data, 0.05, 2, "Gene Expression Waterfall Plot")
easy_plot_fc <- function(dat1, q_val_cutoff = 0.1, fc_cutoff = 2, title_ = "") {
    ds <- dat1 %>%
        group_by(gene_symbol) %>%
        left_join(dat1 %>% count_owners()) %>%
        reframe(
            # mean_log2_fc = mean(log2_fc, na.rm = TRUE),
            median_log2_fc = median(log2_fc, na.rm = TRUE),
            log2_fc = list(log2_fc),
            use_neg_log10_stat_val = -log10(adjusted_p_val),
            significant_stat_val = ifelse(!is.na(adjusted_p_val), adjusted_p_val < q_val_cutoff, FALSE),
            significant = unique(significant_stat_val),
            n_owners = unique(n_owners),
            owners = unique(owners)
        ) %>%
        ## threshold on MEAN FC also
        filter(abs(median_log2_fc) > log2(fc_cutoff) & significant) %>%
        arrange(gene_symbol) %>%
        # these two arrange calls HAVE to be separate
        arrange(desc(median_log2_fc)) %>%
        mutate(order = row_number()) %>%
        mutate(
            gene_symbol = factor(gene_symbol, levels = gene_symbol[order]),
            up_down = factor(ifelse(median_log2_fc > 0, "++", "--"),
                levels = c("++", "--")
            )
        )
    owners_labeller_df <- dat1 %>%
        left_join(
            ds %>%
                dplyr::select(gene_symbol, median_log2_fc, order, up_down, significant),
            by = join_by(gene_symbol)
        ) %>%
        ## threshold on FC also
        filter(abs(median_log2_fc) > log2(fc_cutoff) & significant) %>%
        dplyr::select(gene_symbol, owner, log2_fc, median_log2_fc, order, up_down) %>%
        arrange(order)

    num_analytes <- nrow(ds)

    ds_final <- ds %>% unnest(log2_fc)

    wf_plot <- ggplot(ds_final,
        mapping = aes(
            x = gene_symbol, y = log2_fc,
            color = up_down, group = gene_symbol
        )
    ) +
        geom_hline(yintercept = 0, color = "orange") +
        geom_hline(yintercept = c(-log2(fc_cutoff), log2(fc_cutoff)), color = "#9700e8", alpha = 0.6, linewidth = 1.1, linetype = 4) +
        # boxplot first
        geom_boxplot() +
        # then jitters
        geom_point(owners_labeller_df,
            mapping = aes(x = gene_symbol, y = log2_fc, shape = owner, fill = up_down),
            color = "black",
            size = rel(4),
            inherit.aes = FALSE
        ) +
        theme_prism(base_size = 20) +
        theme(
            strip.text.x = element_text(size = rel(1.25), face = "bold"),
            axis.text.x = element_text(size = rel(0.7), angle = 90, hjust = 1), # , angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = rel(0.9)),
            panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = 0.5),
            panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = 0.25)
        ) +
        scale_fill_manual(values = c("#474747", "red")) +
        scale_color_manual(values = c("#474747", "red")) +
        scale_shape_manual(values = c("Rajagopal" = 21, "Kruse" = 22, "Von Zastrow" = 23, "Rockman" = 24)) +
        ggtitle(qq("@{title_}")) +
        labs(caption = qq("q value cuttoff: @{q_val_cutoff}\n abs(raw median(FC)) > @{fc_cutoff}\nN = @{num_analytes}")) +
        ylab(TeX("$\\bf{log_2(FC)}$")) +
        labs(x = "Gene Symbol")

    return(list(wf_plot, ds_final %>% dplyr::select(-median_log2_fc, -significant_stat_val) %>% rename(neg_log10_q_val = use_neg_log10_stat_val)))
}

#' Generates a volcano plot with simplified filtering.
#'
#' This function creates a volcano plot for gene expression data, applying basic criteria
#' for significance based on fold change and q-value thresholds. It includes labeling and aesthetic adjustments
#' for better visualization.
#'
#' @param df The dataset to be used for the volcano plot.
#' @param fc_thresh The fold change threshold for significance.
#' @param q_val_cutoff The q-value cutoff for significance.
#' @param title_ The title for the plot.
#' @return A list containing the ggplot object for the volcano plot and the processed dataset.
#' @examples
#' data <- read.csv("path/to/expression_data.csv")
#' volcano_plot <- easy_plot_volcano(data, 2, 0.05, "Gene Expression Volcano Plot")
easy_plot_volcano <- function(df, fc_thresh, q_val_cutoff, title_) {
    counted_owners_df <- count_owners(df)

    my_ylab <- TeX("$\\bf{-log_{10}(Q)}$")

    plot_df_temp <- df %>%
        group_by(gene_symbol) %>%
        mutate(
            # mean_log2_fc = mean(log2_fc, na.rm = TRUE),
            use_neg_log10_stat_val = -log10(adjusted_p_val),
            median_log2_fc = median(log2_fc, na.rm = TRUE),
        ) %>%
        mutate(significant = ifelse((abs(median_log2_fc) > log2(fc_thresh)) & (adjusted_p_val < q_val_cutoff), TRUE, FALSE)) %>%
        left_join(counted_owners_df)
    max_for_alpha <- max(sqrt(plot_df_temp$median_log2_fc^2 + plot_df_temp$use_neg_log10_stat_val^2))

    plot_df <- plot_df_temp %>%
        mutate(
            res = sqrt(median_log2_fc^2 + use_neg_log10_stat_val^2),
            significant_alpha = ifelse(significant,
                (res / max_for_alpha)^1.5,
                0
            ), .before = 1 # 1, 0.01)
        ) %>%
        distinct(
            gene_symbol, log2_fc,
            use_neg_log10_stat_val,
            significant, significant_alpha, n_owners, owners
        )

    label_df <- plot_df %>%
        filter(significant)

    num_analytes <- nrow(label_df %>% distinct(gene_symbol))
    label_variable <- "Gene Symbol"

    volcano_plot <- ggplot(plot_df, aes(
        x = log2_fc,
        y = use_neg_log10_stat_val,
    ), color = "black") +
        geom_hline(yintercept = -log10(q_val_cutoff), color = "orange") +
        geom_vline(xintercept = c(-log2(fc_thresh), log2(fc_thresh)), color = "#370afc", alpha = 0.6, linewidth = 1.1, linetype = 4) +
        geom_point(
            mapping = aes(
                fill = significant,
                alpha = significant_alpha
            ),
            pch = 21,
            show.legend = FALSE, size = 2.5
        ) +
        ylim(c(0, max(plot_df$use_neg_log10_stat_val, na.rm = TRUE) + 0.25)) +
        xlab(TeX("$\\bf{log_2(FC)}$")) +
        ylab(my_ylab) +
        scale_fill_manual(values = c("black", "red")) +
        geom_label_repel(label_df,
            mapping = aes(label = gene_symbol),
            force = 10, max.overlaps = 10,
            min.segment.length = 0.2, show.legend = FALSE
        ) +
        theme_prism(base_size = 20) +
        theme(
            strip.text.x = element_text(size = rel(1.25), face = "bold"),
            axis.text = element_text(size = rel(0.8)), # , angle = 45, hjust = 1, vjust = 1),
            panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = 0.5),
            panel.grid.minor = element_line(colour = "gray", linetype = 2, linewidth = 0.25)
        ) +
        ggtitle(qq("@{title_}")) +
        labs(caption = qq("q value cuttoff: @{q_val_cutoff}\n abs(raw median(FC)) > @{fc_thresh}\nN = @{num_analytes}"))
    return(list(volcano_plot, plot_df %>%
        dplyr::select(-significant_alpha) %>%
        arrange(desc(log2_fc), use_neg_log10_stat_val) %>%
        rename(neg_log10_q_val = use_neg_log10_stat_val)))
}


#' Generates waterfall and volcano plots for combined data (Figure 2 plots).
#'
#' This function creates waterfall and volcano plots based on combined data from multiple datasets.
#' It involves generating combined p-values, creating significant data frames, and plotting.
#'
#' @param combined_dat The dataset containing combined information from multiple studies.
#' @param results_path Path to save the generated plots.
#' @param stat_value_type Type of statistical value ('p_value' or 'q_value') for significance testing.
#' @param stat_cutoff Cutoff threshold for statistical significance.
#' @param fc_cutoff Fold change cutoff for significance.
#' @param verbose Boolean to control the display of detailed messages during processing.
#' @return Saves waterfall and volcano plots to the specified results path.
#' @examples
#' combined_dat <- read.csv("path/to/combined_data.csv")
#' make_figure2_plots(combined_dat, "path/to/results")
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
