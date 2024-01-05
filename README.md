# Identification of G Protein and β-Arrestin Independent GPCR Signaling Partners

This codebase accompanies the manuscript titled as above. It includes a comprehensive set of R functions for data preprocessing, analysis, and visualization in the context of GPCR signaling research.

## Data Preprocessing and Identifier Mapping

Data from Rajagopal, Kruse, Rockman, and Von Zastrow datasets were preprocessed and mapped to human-readable gene names for further analysis. This involved filtering, log2 transformation, and outlier removal. Protein accessions were mapped using UniProtKB and Entrez Gene identifiers for subsequent analyses in Gene Ontology/KEGG and STRINGdb.

## Combined Fisher's P-Value Calculation

We employed Fisher's method for combined p-value calculation using the "metap" R package. This involved converting two-sided p-values to one-sided, grouping by accession, and correcting the resultant combined Fisher p-values for multiple hypothesis testing.

## Gene Ontology/KEGG Enrichment and STRINGdb Analysis

The "clusterProfiler" and "STRINGdb" R packages facilitated GO/KEGG enrichment and STRINGdb protein-protein interaction, network visualization, and clustering analyses. We set specific thresholds for enrichment and interaction scores and used the FastGreedy algorithm for clustering.

## Custom R Functions

The repository includes several custom R functions to aid in data analysis and visualization:

### Visualization Functions

- `plot_fc`: For generating waterfall plots, showcasing protein expression data.
- `plot_volcano`: For creating volcano plots, highlighting significant protein expressions.

### Data Analysis Functions

- `generate_intervals`, `find_good_intervals`, `generate_quant_intervals`: Functions for interval generation in data analysis.
- `generate_my_mapped_proteins`: Maps gene names to STRINGdb IDs and summarizes data.
- `make_tbl_data`: Saves data into Excel files for further analysis.
- `plot_stringdb_interactions`, `generate_interaction_enrichment_tables`: Functions for STRINGdb interaction analysis and visualization.

### Utility Functions

- `remove_fig_number`, `p_friend`: Assist in text and data formatting.
- `make_vector_genes`: Creates a list of unique gene symbols and corresponding Entrez IDs.
- `make_naming_tribble_waterfall`, `make_naming_tribble_volcano`: Utility functions for naming and organizing plot data.
- `generate_significant_hits_interaction_network`, `generate_signif_hits_enrichment_tables`, `generate_pert_pathways`: For analyzing and visualizing significant interaction networks.

### Easy Plotting Functions

- `easy_plot_fc`, `easy_plot_volcano`: Simplified functions for quick and efficient generation of waterfall and volcano plots.

## Installation and Usage

To use these functions, clone this repository and source the functions in your R environment. Example usage for each function is provided within the function documentation. The usage order of the main executing scripts was:

- `load_and_process_data.R`
- `figure2B-D.R`
- `figures_2E-3-4.R`
- `table_1.R`

## Contributions and Feedback

Contributions to this codebase are welcome. Please feel free to fork the repository, make changes, and submit pull requests. For bugs, questions, or suggestions, please open an issue in the repository.

---

This project represents a collaborative effort to advance the understanding of GPCR signaling pathways. We hope this codebase will be a valuable resource for researchers in the field.
