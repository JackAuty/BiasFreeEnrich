% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BFE.R
\name{BiasFreeEnrich}
\alias{BiasFreeEnrich}
\alias{BFE}
\title{BiasFreeEnrich: Permutation-based Protein Signature Analysis}
\usage{
BiasFreeEnrich(
  sig_gene,
  ctrl_gene,
  min_protein_count = 2,
  repeats = 10000,
  databases = c("GO_Biological_Process_2025", "GO_Cellular_Component_2025",
    "GO_Molecular_Function_2025", "KEGG_2021_Human", "KEGG_2019_Mouse",
    "WikiPathways_2024_Mouse", "WikiPathways_2024_Human"),
  species = c("Human", "Mouse")
)

BFE(
  sig_gene,
  ctrl_gene,
  min_protein_count = 2,
  repeats = 10000,
  databases = c("GO_Biological_Process_2025", "GO_Cellular_Component_2025",
    "GO_Molecular_Function_2025", "KEGG_2021_Human", "KEGG_2019_Mouse",
    "WikiPathways_2024_Mouse", "WikiPathways_2024_Human"),
  species = c("Human", "Mouse")
)
}
\arguments{
\item{sig_gene}{A vector of significant gene symbols corresponding to proteins of interest.}

\item{ctrl_gene}{A vector of gene symbols from the control group corresponding to proteins.}

\item{min_protein_count}{Minimum number of proteins in a pathway to be considered meaningful (default = 2).}

\item{repeats}{Number of permutations to run for each ontology. Higher values increase precision but slow execution (default = 10000).}

\item{databases}{Vector of databases to compare against. Options include:
"GO_Biological_Process_2025", "GO_Cellular_Component_2025", "GO_Molecular_Function_2025",
"KEGG_2021_Human", "KEGG_2019_Mouse", "WikiPathways_2024_Mouse", "WikiPathways_2024_Human" (default = all of them)}

\item{species}{INSTEAD of selecting the databases you can choose a species. It will then enrich on all databases from that species: "mouse" or "human". Do not attempt to select databases AND species.}
}
\value{
A list of significantly enriched pathways for each selected database.
}
\description{
BiasFreeEnrich: Permutation-based Protein Signature Analysis
}
\examples{
sig_gene <- c("IL1B", "IL1A", "IL6")
ctrl_gene <- c("IL1B", "IL1A", "IL6", "IL10", "TGFB", "IL4")
results <- BiasFreeEnrich(sig_gene, ctrl_gene, databases = c("GO_Biological_Process_2025", "KEGG_2021_Human"))
}
