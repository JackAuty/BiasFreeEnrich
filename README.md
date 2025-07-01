
# BiasFreeEnrich

**BiasFreeEnrich** is an R package designed for permutation-based protein signature enrichment analysis. Traditional enrichment tools often overestimate significance by comparing to whole-genome backgrounds. This leads to biologically misleading results, especially in tissue-specific or small-scale proteomics datasets. **BiasFreeEnrich** corrects for this bias by using permutation-based null distributions derived from user-defined control protein lists, enabling more accurate and testable hypotheses.

## 🔍 Overview

Enrichment analysis is widely used to interpret proteomics and transcriptomics data by identifying over-represented pathways or ontologies. However, many methods do not account for sample-specific background protein expression. **BiasFreeEnrich** addresses this by:

- Performing pathway enrichment against background gene lists provided by the user.
- Using permutation tests to simulate the null distribution of overlapping genes.
- Returning adjusted p-values using both Benjamini-Hochberg (BH) and discrete q-value methods (Chen et al.).
- Supporting multiple databases: GO (2025), KEGG, and WikiPathways.
- You can perform a standard enrichment analysis against the entire genome by simply omitting a background dataset. In this case, a   Benjamini-Hochberg (BH) false discovery rate (FDR) correction will be applied by default.


## 🚀 Installation

```r
# Install dependencies
install.packages(c("stringr", "DiscreteQvalue", "remotes", "enrichR"))
remotes::install_github("wjawaid/enrichR")

# Install BiasFreeEnrich
remotes::install_github("JackAuty/BiasFreeEnrich")
```

📦 Package Structure
Main function: BiasFreeEnrich() or alias BFE()

Input: Vector of significant genes and control genes

Output: A list of enriched pathways across selected databases, adjusted for sample bias

🧪 Example
```r

library(BiasFreeEnrich)

sig_gene <- c("IL1B", "IL1A", "IL6")
ctrl_gene <- c("IL1B", "IL1A", "IL6", "IL10", "TGFB", "IL4")

results <- BiasFreeEnrich(
  sig_gene,
  ctrl_gene,
  databases = c("GO_Biological_Process_2025", "KEGG_2021_Human")
)

# View results
print(results)

```
📚 Available Databases
Gene Ontology (GO):

GO_Biological_Process_2025

GO_Cellular_Component_2025

GO_Molecular_Function_2025

KEGG:

KEGG_2021_Human

KEGG_2019_Mouse

WikiPathways:

WikiPathways_2024_Human

WikiPathways_2024_Mouse

Alternatively, you may specify species = "Human" or species = "Mouse" to include all relevant databases automatically.

⚠️ Do not specify both databases and species.

🧠 Why Use BiasFreeEnrich?
Traditional gene enrichment methods compare against the whole genome, which leads to:

False positives due to tissue-specific protein expression.

Excessively conservative corrections assuming continuous p-value distributions.

BiasFreeEnrich uses:

Permutation tests to generate null models tailored to your dataset.

Discrete FDR correction (Chen et al.), more appropriate for gene count data.

Custom background gene sets, essential for unbiased enrichment.

🗂️ Output Format
It returns a list containing each database searched.
Each database returns a dataframe with:

Term: Pathway or ontology name

Significant_genes: Enriched gene symbols

P_value_background_adjusted: Permutation-derived p-value

P_value_background_adjusted_BH_FDR: BH-adjusted FDR

P_value_background_adjusted_CHEN_FDR: Discrete q-value FDR

P_value_whole_genome: Original enrichR p-value

Sig_gen_number: Number of significant genes in pathway

Ctrl_gene_number: Number of control genes in pathway

