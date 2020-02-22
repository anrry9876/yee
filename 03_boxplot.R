setwd("/Users/simon/GoogleDrive/NTUH/Projects/Survival")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(beeswarm)
library(ggpubr)

KIRC <- readRDS("TCGA-KIRC_RNAseq_count.rds")
KIRP <- readRDS("TCGA-KIRP_RNAseq_count.rds")
BLCA <- readRDS("TCGA-BLCA_RNAseq_count.rds")

## Interested gene
query.gene <- "HMGB1"
query.id <- rowData(LGG)$ensembl_gene_id[rowData(LGG)$external_gene_name == query.gene]

# Get table
get_table <- function(data, gene){

  # Exp data
  d.g <- data.frame(assay(data)[gene,])
  colnames(d.g) <- 'expression'
  d.g <- tibble::rownames_to_column(d.g, "barcode")

  # Clinical data
  aa <- colData(data)
  bb <- cbind(aa$barcode, aa$definition)
  bb <- data.frame(bb)
  colnames(bb) <- c("barcode", "definition")
  bb$tumor <- ifelse(bb$definition == "Solid Tissue Normal", "Normal", "Tumor")

  # Merge exp & clinical
  output <- merge(bb, d.g, by = 1)

  return(output)
}
# Plot
get_boxplot <- function(data, data_name, gene_name){

  amax <- max(data$expression)
  amin <- min(data$expression)
  m <- (amax-amin)/10

  p <- compare_means(expression ~ tumor, data, method = "wilcox.test")
  if ( p$p < 0.0001){
    p$p <- 'p < 0.0001'
  }else{
    p$p <- paste('p =', formatC(p$p,format = "e", digits = 2))
  }

  beeswarm(expression ~ tumor, data = data,
           pch = 16, col = rainbow(8),
           corral = 'wrap', main = paste(gene_name, 'in', data_name),
           add = FALSE,
           ylim = c(amin-m, amax+m),
           ylab = 'expression', xlab = ""
  )
  boxplot(expression ~ tumor, data = data,
          outline=FALSE, col = rgb(0.1,0.1,0.1,0),
          add = TRUE
  )

  text(x = 1.5, y = amax+0.9*m, labels = formatC(p$p,format = "e", digits = 2))
  #lines(x = c(1,1,2,2), y = c(amax+2,amax+4,amax+4,amax+2))

}
get_triple_boxplot <- function(d1, d2, d3, gene_name){
  png(paste(gene_name, 'boxplot', 'png', sep = '.'), width = 30, height = 10, res = 600, units = "cm")
  par(mar = c(4, 4, 2, 2), cex = 0.8, mfrow=c(1,3))
  get_boxplot(d1, 'KIRC', gene_name)
  get_boxplot(d2, 'KIRP', gene_name)
  get_boxplot(d3, 'BLCA', gene_name)
  dev.off()
}

# OPA1
query.gene <- "OPA1"
query.id <- rowData(LGG)$ensembl_gene_id[rowData(LGG)$external_gene_name == query.gene]
op.c <- get_table(KIRC, query.id)
op.p <- get_table(KIRP, query.id)
op.b <- get_table(BLCA, query.id)
get_triple_boxplot(op.c, op.p, op.b, query.gene)
# MFN1
query.gene <- "MFN1"
query.id <- rowData(LGG)$ensembl_gene_id[rowData(LGG)$external_gene_name == query.gene]
m1.c <- get_table(KIRC, query.id)
m1.p <- get_table(KIRP, query.id)
m1.b <- get_table(BLCA, query.id)
get_triple_boxplot(m1.c, m1.p, m1.b, query.gene)
# MFN2
query.gene <- "MFN2"
query.id <- rowData(LGG)$ensembl_gene_id[rowData(LGG)$external_gene_name == query.gene]
m2.c <- get_table(KIRC, query.id)
m2.p <- get_table(KIRP, query.id)
m2.b <- get_table(BLCA, query.id)
get_triple_boxplot(m2.c, m2.p, m2.b, query.gene)
# LONP1
query.gene <- "LONP1"
query.id <- rowData(LGG)$ensembl_gene_id[rowData(LGG)$external_gene_name == query.gene]
lo.c <- get_table(KIRC, query.id)
lo.p <- get_table(KIRP, query.id)
lo.b <- get_table(BLCA, query.id)
get_triple_boxplot(lo.c, lo.p, lo.b, query.gene)
# Prohibitin
pr.c <- get_table(KIRC, "ENSG00000167085")
pr.p <- get_table(KIRP, "ENSG00000167085")
pr.b <- get_table(BLCA, "ENSG00000167085")
# PGC-1alpha
pg.c <- get_table(KIRC, "ENSG00000109819")
pg.p <- get_table(KIRP, "ENSG00000109819")
pg.b <- get_table(BLCA, "ENSG00000109819")
# Drp1
dr.c <- get_table(KIRC, "ENSG00000087470")
dr.p <- get_table(KIRP, "ENSG00000087470")
dr.b <- get_table(BLCA, "ENSG00000087470")
# Binp3
bi.c <- get_table(KIRC, "ENSG00000176171")
bi.p <- get_table(KIRP, "ENSG00000176171")
bi.b <- get_table(BLCA, "ENSG00000176171")
# PKM2
pk.c <- get_table(KIRC, "ENSG00000067225")
pk.p <- get_table(KIRP, "ENSG00000067225")
pk.b <- get_table(BLCA, "ENSG00000067225")
# FoxO1
fo.c <- get_table(KIRC, "ENSG00000150907")
fo.p <- get_table(KIRP, "ENSG00000150907")
fo.b <- get_table(BLCA, "ENSG00000150907")
# Glutaminase
gl.c <- get_table(KIRC, "ENSG00000115419")
gl.p <- get_table(KIRP, "ENSG00000115419")
gl.b <- get_table(BLCA, "ENSG00000115419")
# HMGB1
hm.c <- get_table(KIRC, "ENSG00000189403")
hm.p <- get_table(KIRP, "ENSG00000189403")
hm.b <- get_table(BLCA, "ENSG00000189403")
# CDK5
cd.c <- get_table(KIRC, "ENSG00000164885")
cd.p <- get_table(KIRP, "ENSG00000164885")
cd.b <- get_table(BLCA, "ENSG00000164885")




get_triple_boxplot(pr.c, pr.p, pr.b, 'Prohibitin')
get_triple_boxplot(pg.c, pg.p, pg.b, 'PGC1alpha')
get_triple_boxplot(dr.c, dr.p, dr.b, 'Drp1')
get_triple_boxplot(bi.c, bi.p, bi.b, 'Binp3')
get_triple_boxplot(pk.c, pk.p, pk.b, 'PKM2')
get_triple_boxplot(fo.c, fo.p, fo.b, 'FoxO1')
get_triple_boxplot(gl.c, gl.p, gl.b, 'Glutaminase')
get_triple_boxplot(hm.c, hm.p, hm.b, 'HMGB1')
get_triple_boxplot(cd.c, cd.p, cd.b, 'CDK5')










#png("example.png", width = 30, height = 10, res = 600, units = "cm")
#par(mar = c(4, 4, 2, 2), cex = 0.8, mfrow=c(1,3))
#get_boxplot(pr.c, 'KIRC', 'Prohibitin')
#get_boxplot(pr.p, 'KIRP', 'Prohibitin')
#get_boxplot(pr.b, 'BLCA', 'Prohibitin')
#dev.off()
