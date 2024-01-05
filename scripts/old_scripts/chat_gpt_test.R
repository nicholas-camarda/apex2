# chat_gpt_test.R

library(ggplot2)

# Load the data
data_set1 <- read.csv("~/Downloads/final_data_set1_new.csv") %>%
    as_tibble() %>%
    arrange(desc(neg_log10_adjust_p_val))
data_set2 <- read.csv("~/Downloads/final_data_set2_new.csv") %>%
    as_tibble() %>%
    filter(abs(log2_fc_mean) > log2(2) & neg_log10_adjust_p_val > -log10(0.1)) %>%
    arrange(desc(log2_fc_mean))

## Create a function to plot a volcano plot
plot_volcano <- function(data, title) {
    ggplot(data, aes(x = log2_fc_mean, y = neg_log10_adjust_p_val)) +
        geom_point(alpha = 0.6) +
        geom_hline(yintercept = -log10(0.1), color = "red", linetype = "dashed") +
        geom_vline(xintercept = log2(2), color = "blue", linetype = "dashed") +
        geom_vline(xintercept = -log2(2), color = "blue", linetype = "dashed") +
        labs(title = title, x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
        theme_minimal() +
        theme(
            legend.position = "none"
        )
}

# Plot the volcano plots
plot1 <- plot_volcano(data_set1, "Volcano Plot (Rockman and Von Zastrow)")
plot2 <- plot_volcano(data_set2, "Volcano Plot (Rockman, Von Zastrow, and Kruse)")

# Print the plots
print(plot1)
print(plot2)
