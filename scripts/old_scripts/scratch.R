library(tidyverse)

init_data <- read_rds("cache/initial-dataset_combined-dat.rds")

zastrow_rockman <- init_data %>%
    filter(owner %in% c("Von Zastrow", "Rockman"))

# power calculation for this -- what's the deal with so little hits when we look at von zastrow and rockman alone, and including kruse?
zastrow_rockman %>%
    group_by(Accession) %>%
    summarize(mean_log2_fc = mean(log2_fc, na.rm = TRUE)) %>%
    filter(abs(mean_log2_fc) > 1)



## rockman kruse and von zastrow --> what is significant?
## compare these hits to rajagopal
## is there significant overlap here?

# TODO #1 - done
## table
# hits that are significant in the other 3
# hits that are significant in Rajagopal
# what is shared
# what is not shared

# TODO #2
## all possible combinations (three +) for waterfalls and volcanos (including and excluding kruse and raja)
# for the pathway stuff, strictly conserved across all 4

# TODO #3  - DONE
## Clarify file and folder names -- label explicitly

# TODO #4
## send to everyone with explanation

# TODO #5
# write methods!
