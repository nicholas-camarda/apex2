library(tidyverse)
library(readxl)

# read in the data
my_fns <- tribble(~fn_name, ~sheet, 
                  "1-s2.0-S0092867417303471-mmc1.xlsx", 2,  # This is the sheet we want to use (use the 2 minute time points)
                  "5909_SupplementalData_020822.xlsx", 5,
                  "Summary_90s_log2_1v1comp.xlsm", 2, 
                  "Summary_10m_log2_1v1comp.xlsm", 2, 
                  "1-s2.0-S0092867417303021-mmc3.xlsx", 2) # use 1, 3, and 10 min time points

fns_tbl <- tibble(fn_path = dir("data", recursive = T, pattern = "xls", full.names = TRUE)) %>%
  mutate(owner = str_trim(string = gsub(x = basename(dirname(fn_path)), pattern = "Data", replacement = ""), side = "both"), 
         fn_name = basename(fn_path), .before = 1) %>%
  inner_join(my_fns) %>%
  mutate(data = map2(.x = fn_path, .y = sheet, .f = function(x,y){
    read_excel(path = x, sheet = y)
  })); fns_tbl

# rajagopal

# is this fold change or difference??
rajagopal_cxcl2 <- fns_tbl$data[[2]] %>%
  dplyr::select(Accession, Peptides, Unique, CXCL12_v_Veh_FC, CXCL12_v_Veh_pval, Stats) %>%
  mutate(tx = "veh vs cxcl12") %>%
  rename(fc = CXCL12_v_Veh_FC, p_val = CXCL12_v_Veh_pval); rajagopal_cxcl2

colnames(rajagopal_cxcl2)


# kruse
kruse_ang2_2min <- fns_tbl$data[[1]] %>%
  dplyr::select(Accession, `Gene Symbol`, contains("2 min") & contains("Ang II")) %>%
  rename(rep1_2min = "2 min Ang II biolgical replicate 1",
         rep2_2min = "2 min Ang II biolgical replicate 2",
         fc_2min = "2 min Ang II enrichment relative to ligand-free control") %>%
  mutate(timing = "2min") %>%
  rename(rep1 = rep1_2min, rep2 = rep2_2min, fc = fc_2min); 

kruse_ang2_21min <- fns_tbl$data[[1]] %>%
  dplyr::select(Accession, `Gene Symbol`, contains("21 min") & contains("Ang II")) %>%
  rename(rep1_21min = "21 min Ang II biolgical replicate 1",
         rep2_21min = "21 min Ang II biological replicate 2",
         fc_21min = "21 min Ang II enrichment relative to ligand-free control") %>%
  mutate(timing = "21min") %>%
  rename(rep1 = rep1_21min, rep2 = rep2_21min, fc = fc_21min); 

kruse_ang2_2min_21min <- bind_rows(kruse_ang2_2min, kruse_ang2_21min) %>%
  mutate(tx = "veh vs angII"); kruse_ang2_2min_21min

colnames(kruse_ang2_2min_21min)

# Rockman

# is fc really a fold change, or just a difference
rockman_10m_90s <- bind_rows(fns_tbl$data[[3]] %>% 
                               mutate(timing = "10min"), 
                             fns_tbl$data[[4]] %>%
                               mutate(timing = "90sec")) %>%
  dplyr::select(Accession, `Gene Name`, Ang , NS, `Ang-NS`, `Ang-NS pval`, `Ang-NSpval (adj)`, timing) %>%
  rename(gene_symbol = `Gene Name`) %>%
  # this is just difference... not fold change..?
  rename(fc = `Ang-NS`, p_val = `Ang-NS pval`, adjusted_p_val = `Ang-NSpval (adj)`) %>%
  mutate(tx = "ns vs angII"); rockman_10m_90s

colnames(rockman_10m_90s)

rockman_10m <- rockman_10m_90s %>% filter(timing == "10min")
rockman_90s <- rockman_10m_90s %>% filter(timing == "90sec")

# Von Zastrow
Zastro_all <- fns_tbl$data[[5]] %>%
  dplyr::select(`Uniprot accession`, `Uniprot name`, `Gene symbol`, FC, FDR, Timing) %>%
  rename(Accession = `Uniprot accession`, gene_symbol = `Gene symbol`, timing = `Timing`) %>%
  mutate(timing = str_split(timing, pattern = "_", simplify = TRUE)[,3]) %>%
  filter(timing %in% c("1min", "3min", "10min")); vZastro_all
