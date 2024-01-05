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


#' @note plotting watefall
plot_fc <- function(df, grouping_variable = "Accession", fc_cutoff = 2, stat_cutoff_val = 0.001, my_stat_value_type = "q_value", title_ = "Waterfall plot", found_in_title = ">= 2") {
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

#' @note make volcano plots
plot_volcano <- function(df, stat_cutoff_val = 0.001, label_variable = "Accession", my_stat_value_type = "q_value", fc_thresh = 2, found_in_title = "num_ds >= 2", title_ = "Volcano plot") {
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

add_string_to_filename_basic <- function(filename, string_to_add, file_ext = "png") {
    new_str <- str_split(filename, qq("\\.@{file_ext}"), simplify = TRUE)[, 1]
    new_filename <- paste0(new_str, "-", string_to_add, ".", file_ext)
}

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

make_stat_directory <- function(sval, type = "q") {
    sval_str <- as.character(sval)
    sval_str <- gsub("\\.", "-", sval_str)
    dir_name <- paste0(type, sval_str)

    return(dir_name)
}



#' @note less restrictions for easy plotting of waterfall plot
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

#' @note easy volcano
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
