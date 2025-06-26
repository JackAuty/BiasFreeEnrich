
#' BiasFreeEnrich: Permutation-based Protein Signature Analysis
#'
#' @param sig_gene A vector of significant gene symbols corresponding to proteins of interest.
#' @param ctrl_gene A vector of gene symbols from the control group corresponding to proteins.
#' @param min_protein_count Minimum number of proteins in a pathway to be considered meaningful (default = 2).
#' @param repeats Number of permutations to run for each ontology. Higher values increase precision but slow execution (default = 10000).
#' @param databases Vector of databases to compare against. Options include:
#'   "GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023",
#'   "KEGG_2021_Human", "KEGG_2019_Mouse", "WikiPathways_2024_Mouse", "WikiPathways_2024_Human" (default = all of them)
#' @param species INSTEAD of selecting the databases you can choose a species. It will then enrich on all databases from that species: "mouse" or "human". Do not attempt to select databases AND species.
#'
#' @return A list of significantly enriched pathways for each selected database.
#' @export
#'
#' @examples
#' sig_gene <- c("IL1B", "IL1A", "IL6")
#' ctrl_gene <- c("IL1B", "IL1A", "IL6", "IL10", "TGFB", "IL4")
#' results <- BiasFreeEnrich(sig_gene, ctrl_gene, databases = c("GO_Biological_Process_2025", "KEGG_2025_Human"))
#' @name BiasFreeEnrich
BiasFreeEnrich<-function(sig_gene,
                         ctrl_gene,
                         min_protein_count=2,
                         repeats = 10000,
                         databases = c("GO_Biological_Process_2025",
                                       "GO_Cellular_Component_2025",
                                       "GO_Molecular_Function_2025",
                                       "KEGG_2021_Human",
                                       "KEGG_2019_Mouse",
                                       "WikiPathways_2024_Mouse",
                                       "WikiPathways_2024_Human"),
                         species = c("Human", "Mouse")){

  pb <- txtProgressBar(min = 0, max = (5+2*length(databases)), style = 3)
  prog<-1
  setTxtProgressBar(pb, prog)
  cat("\n")
  if(!missing(ctrl_gene)){
    if(length(sig_gene) > length(ctrl_gene)){
      stop("Error: Significant gene list is
           larger than the Control gene list")}
  }

  if(is.numeric(databases)){
    db_text<- c("GO_Biological_Process_2025",
                "GO_Cellular_Component_2025",
                "GO_Molecular_Function_2025",
                "KEGG_2021_Human",
                "KEGG_2019_Mouse",
                "WikiPathways_2024_Mouse",
                "WikiPathways_2024_Human")[databases]
  } else {
    db_text<-databases
  }
  species <- c(species, "GO")
  spec<-paste(species, collapse = "|")
  db<-db_text[grepl(spec, db_text)]

  if(!missing(ctrl_gene)){

    ctrl_gene <- toupper(unique(as.character(ctrl_gene)))
  }


  sig_gene <- toupper(unique(as.character(sig_gene)))

  results_list<-list()
  prog<-2
  setTxtProgressBar(pb, prog)
  cat("\n")

  if(missing(ctrl_gene)){
    cat("Warning no background data provided")
    #Loop through database
    for(k in 1:length(db)){
      prog <- prog +1
      setTxtProgressBar(pb, prog)
      Significant_only <- enrichr(sig_gene, databases = db[k])

      Significant_only_df <- data.frame(Significant_only[1])

      colselecnames <- c("Term", "Overlap", "Adjusted", "Genes")
      colselect<-paste(colselecnames, collapse = "|")

      Significant_only_df <- Significant_only_df[,grepl(
        colselect, colnames(Significant_only_df))]
      Significant_only_df<-Significant_only_df[,-which(grepl(
        "Old", colnames(Significant_only_df)))]
      colnames(Significant_only_df)<-colselecnames
      Significant_only_df$Database <- names(Significant_only)[1]


      prog <- prog +1
      setTxtProgressBar(pb, prog)
      cat("\n")
      if(nrow(Significant_only_df) == 0){
        cat(paste("No significant Enrichment found in", db[k]))
        next
      }
      colnames(Significant_only_df)[colnames(Significant_only_df)=="Genes"] <- "Significant_genes"
      colnames(Significant_only_df)[colnames(Significant_only_df)=="Adjusted"] <- "P_value_BH_FDR"
      Significant_only_df$Significant_genes <- str_replace_all(
        Significant_only_df$Significant_genes,
        ";", ", ")
      Significant_only_df$Warning<-"No background dataset provided"
      results_list[[db[k]]] <- Significant_only_df

    }
    setTxtProgressBar(pb, prog+3)
    cat("\n")
    return(results_list)
  }



  #Loop through database
  for(k in 1:length(db)){
    prog <- prog +1
    setTxtProgressBar(pb, prog)
    cat("\n")


    ctrl_pathways <- enrichr(ctrl_gene, databases = db[k])

    ctrl_pathways_result <- data.frame(ctrl_pathways[1])

    colselecnames <- c("Term", "Overlap", "Adjusted", "Genes")
    colselect<-paste(colselecnames, collapse = "|")

    ctrl_pathways_result <- ctrl_pathways_result[,grepl(
      colselect, colnames(ctrl_pathways_result))]
    ctrl_pathways_result<-ctrl_pathways_result[,-which(grepl(
      "Old", colnames(ctrl_pathways_result)))]
    colnames(ctrl_pathways_result)<-colselecnames
    ctrl_pathways_result$Database <- names(ctrl_pathways)[1]


    prog <- prog +1
    setTxtProgressBar(pb, prog)
    cat("\n")
    sig_pathways <- enrichr(sig_gene, databases = db[k])

    sig_pathways_result <- data.frame(sig_pathways[1])

    sig_pathways_result <- sig_pathways_result[,grepl(
      colselect, colnames(sig_pathways_result))]

    sig_pathways_result<-sig_pathways_result[,-which(grepl(
      "Old", colnames(sig_pathways_result)))]

    colnames(sig_pathways_result)<-colselecnames
    sig_pathways_result$Database <- names(sig_pathways)[1]




    ctrl_pathways_result$gene_number<-  str_count(
      ctrl_pathways_result$Genes, ";")+1
    ctrl_pathways_result<-ctrl_pathways_result[
      ctrl_pathways_result$gene_number>min_protein_count,]

    sig_pathways_result$gene_number<-  str_count(
      sig_pathways_result$Genes, ";")+1

    sig_pathways_result<-sig_pathways_result[
      sig_pathways_result$gene_number>min_protein_count,]
    if(nrow(sig_pathways_result) == 0){
      cat(paste("No significant Enrichment found in", db[k]))
      next
    }

    sig_pathways_result$P_value_background_adjusted<-NA
    sig_pathways_result$Ctrl_gene_number <- NA



    sig_pathways_result$Term_id<-paste(
      sig_pathways_result$Term,
      sig_pathways_result$Database, sep =" ")
    sig_pathways_result$Term_id<-str_remove_all(
      sig_pathways_result$Term_id, "Term")

    ctrl_pathways_result$Term_id<-paste(
      ctrl_pathways_result$Term,
      ctrl_pathways_result$Database, sep =" ")
    ctrl_pathways_result$Term_id<-str_remove_all(
      ctrl_pathways_result$Term_id, "Term")

    for(i in 1:nrow(sig_pathways_result)){

      if(sig_pathways_result$Term_id[i]
         %in% ctrl_pathways_result$Term_id ){

        sample_ctrl <- c(
          rep(1, ctrl_pathways_result$gene_number[which(
            ctrl_pathways_result$Term_id ==
              sig_pathways_result$Term_id[i])]),
          rep(0, length(ctrl_gene)-
                ctrl_pathways_result$gene_number[which(
                  ctrl_pathways_result$Term_id ==
                    sig_pathways_result$Term_id[i])]))

        n = rep(9999, repeats)
        for (j in 1:repeats) {
          n[j] <- sum(sample(sample_ctrl, length(sig_gene)))
        }
        sig_interection_total <-
          sig_pathways_result$gene_number[i]
        sig_pathways_result$P_value_background_adjusted [i] <-
          sum(n >= sig_interection_total) / repeats
        sig_pathways_result$Ctrl_gene_number[i]<-ctrl_pathways_result$
          gene_number[which(
            ctrl_pathways_result$Term_id ==
              sig_pathways_result$Term_id[i])]

      } else{
        sig_pathways_result$P_value_background_adjusted [i] <- 0
      }

    }

    sig_pathways_result$P_value_whole_genome<-
      sig_pathways_result$Adjusted

    sig_pathways_result$P_value_whole_genome<-round(
      sig_pathways_result$P_value_whole_genome, 4)



    sig_pathways_result<-sig_pathways_result[order(
      sig_pathways_result$P_value_background_adjusted),]



    sig_pathways_result <- sig_pathways_result[,-which(
      colnames(sig_pathways_result)=="Adjusted")]

    sig_pathways_result_return<-sig_pathways_result[,c(1,8,4,3,5,7,9,6)]




    sig_pathways_result_return$P_value_background_adjusted_BH_FDR <- p.adjust(
      sig_pathways_result$P_value_background_adjusted,
      method="BH"
    )
    sig_pathways_result_return$P_value_background_adjusted_CHEN_FDR <-
      DQ(sig_pathways_result$P_value_background_adjusted,
         ss=sort(unique(sig_pathways_result$P_value_background_adjusted)),
         method = "Chen")$q.value
    length(db)
    cat(k)
    cat(db[k])


    colnames(sig_pathways_result_return)[
      which(colnames(sig_pathways_result_return) == "Genes"|
              colnames(sig_pathways_result_return) == "gene_number")] <-
      c("Significant_genes", "Sig_gen_number")

    sig_pathways_result_return<-sig_pathways_result_return[,
                                                           -which(colnames(
                                                             sig_pathways_result_return
                                                           ) == "Term_id")]
    sig_pathways_result_return$Significant_genes <- str_replace_all(
      sig_pathways_result_return$Significant_genes,
      ";", ", ")

    results_list[[db[k]]] <- sig_pathways_result_return

  }
  setTxtProgressBar(pb, prog+3)
  cat("\n")
  # Load the dataset (it will create an object called gene_list in your workspace)
  data("gene_list", package = "BiasFreeEnrich")
  gene_list$x <- toupper(gene_list$x)

  genes_provided <- toupper(c(sig_gene, ctrl_gene))
  genes_in_list <- genes_provided[genes_provided %in% gene_list$x]
  genes_not_recognised <- setdiff(genes_provided, genes_in_list)

  # cat number of recognised genes
  cat(paste("Number of total genes recognised:", length(genes_in_list)))
  cat("\n")
  # cat number of unrecognised genes
  cat(paste("Number of genes not recognised:", length(genes_not_recognised)))
  cat("\n")
  # Show up to 10 unrecognised genes
  if (length(genes_not_recognised) > 0) {
    to_show <- head(genes_not_recognised, 10)
    cat("Genes not recognised (up to 10 shown):")
    cat(to_show)

    if (length(genes_not_recognised) > 10) {
      cat("Note: There are more than 10 unrecognised genes.\n")
    }
  }
  return(results_list)
}



# Alias for BiasFreeEnrich
#' @rdname BiasFreeEnrich
#' @export
BFE <- BiasFreeEnrich

