##
## Linear regression on interaction effects in single-cell data
##
library(pbapply)
library(ggplot2)
require(broom)

##
## Paths to files (change if needed).
##
genes <- read.table("~/Documents/scRNA-seq/data/lane_1/genes.tsv")
genotypes <- read.table("~/Documents/scRNA-seq/R/data/genotypes/maf_10.genotypes.txt", check.names = F)
dir <- list.files(path="~/Documents/scRNA-seq/data/th_cells_per_sample/magic_output/", pattern="*.tsv", full.names=T, recursive=FALSE)

##
## Read in data.
##
samples <- list()
cell.counts <- c()
sample.names <- vector()

i <- 1
for (file in dir) {
  sample.names <- c(sample.names, tools::file_path_sans_ext(basename(file)))
  sample <- read.csv(file = file, sep = "\t", row.names = 1)
  rownames(sample) <- substr(rownames(sample), start = 7, stop = 21)
  sample <- t(sample)
  cell.counts <- c(cell.counts, nrow(sample))
  
  samples[[i]] <- sample
  i <- i + 1
}

genotypes <- genotypes[,match(sample.names, colnames(genotypes))]

##
## Return a vector with factor for the given snp.id per sample.
##
get.snp <- function(snp.id) {
  snp <- unlist(genotypes[snp.id,])
  snp[snp == "C/T"] <- "T/C"
  snp[snp == "C/G"] <- "G/C"
  snp[snp == "C/A"] <- "A/C"
  snp[snp == "T/G"] <- "G/T"
  snp[snp == "T/A"] <- "A/T"
  snp[snp == "T/C"] <- "C/T"
  snp[snp == "G/A"] <- "A/G"
  snp[snp == "G/C"] <- "C/G"
  snp[snp == "G/T"] <- "T/G"
  return(droplevels(snp))
}

##
## Returns the correlation matrix for the given eqtl gene against all other genes for every sample
##
## Correlation to eQTL gene:
##          Sample 1  Sample 2  Sample 3
## Gene 1     0.5       -0.7      0.5
## Gene 2     0.9       1         1
## Gene 3     0         0         0
##
create.cor.matrix <- function(exp.matrices, gene.name, cor.method = "spearman") {
  samples <- exp.matrices
  
  cor.vectors <- list()
  for (i in 1:length(samples)) {
    if (gene.name %in% colnames(samples[[i]])) {
      samp.cor <- cor(samples[[i]][,gene.name], samples[[i]], method = cor.method)
      cor.vectors[[i]] <- t(samp.cor)
    } else {
      cor.vectors[[i]] <- matrix(NA, nrow = ncol(samples[[i]]), ncol = 1, dimnames = list(colnames(samples[[i]]), NA))
    }
  }
  
  genes <- unique(unlist(lapply(cor.vectors, rownames)))
  cor.matrix <- matrix(NA, nrow = length(genes), ncol = length(cor.vectors), 
                       dimnames = list(genes, sample.names))
  for (i in seq_along(cor.vectors)) {
    cor.matrix[rownames(cor.vectors[[i]]), i] <- cor.vectors[[i]]
  }
  
  cor.matrix[is.na(cor.matrix)] <- 0
  cor.matrix <- cor.matrix[-which(rownames(cor.matrix) == gene.name),]
  
  return(cor.matrix)
}

##
## Writes the correlation matrix to a csv file
##
cor.matrix.to.csv <- function(cor.matrix, snp, path) {
  
  mean.matrix <- matrix(NA, nrow = nrow(cor.matrix), ncol=3, dimnames = list(rownames(cor.matrix), levels(snp)))
  for (genotype in levels(snp)) {
    mean.matrix[,genotype] <- rowMeans(cor.matrix[,snp == genotype], na.rm = F)
  }

  filename <- paste0(gene.name, "_", snp.name, ".csv")
  filepath <- file.path(output.dir, filename)

  write.table( as.matrix( t(as.character(snp)) ), file = filepath, sep = "\t", quote = F, col.names = F)
  write.table(cbind(cor.matrix, mean.matrix), file = filepath, sep = "\t", quote = F, col.names=NA, append = T)
  
}

##
## Applies the linear model correlation~snp overthe correlation of all genes.
##
## Returns a matrix with the summary statistics of the model for every gene.
##
interaction.regression <- function(cor.matrix, eqtl.gene, snp.name) {
  
  snp <- as.numeric(get.snp(snp.name))
  
  cor.statistics <- do.call("rbind", apply(cor.matrix, 1, function(x) {
    model <- lm(formula = x~snp, weights = sqrt(cell.counts))
    #model <- lm(formula = x~snp)
    return(glance(model))
  }))
  
  return(cor.statistics)
}

##
## Plots the interaction effect of the genotype on the co-expression of the given genes
##
## TODO: Add zero values to outputData if absent
plot.correlations.qtl <- function(gene.1, gene.2, snp.id, exp.matrices, xlim = NULL, p.value, r) {
  snp <- unlist(genotypes[snp.id,])
  snp[snp == "C/T"] <- "T/C"
  snp[snp == "C/G"] <- "G/C"
  snp[snp == "C/A"] <- "A/C"
  snp[snp == "T/G"] <- "G/T"
  snp[snp == "T/A"] <- "A/T"
  snp[snp == "T/C"] <- "C/T"
  snp[snp == "G/A"] <- "A/G"
  snp[snp == "G/C"] <- "C/G"
  snp[snp == "G/T"] <- "T/G"
  snp <- droplevels(snp)
  
  i <- 1
  outputData = data.frame(x=numeric(0), y=numeric(0), genotype=character(0))
  for (sample in samples) {
    if (gene.2 %in% colnames(sample) & gene.1 %in% colnames(sample)) {
      outputData = rbind(outputData, data.frame(x=sample[,gene.2], y=sample[,gene.1], genotype=snp[sample.names[[i]]]))
    }
    
    i <- i + 1
  }
  
  plot <- ggplot(outputData, aes(x=x, y=y, colour=genotype)) + 
    geom_point(alpha = 0.4) + 
    geom_smooth(method = "lm", fullrange = T, se=F) + 
    ggtitle(paste(genes[genes$V1 == gene.1,]$V2, snp.id, "R-squared:", r, "P-value:", p.value)) + 
    xlab(gene.2) + ylab(gene.1)
  
  if (!is.null(xlim)) {
    plot <- plot + xlim(xlim)
  }
  plot
}

##
## Returns a matrix with the linear model r-squared for every gene (rows) per eQTL gene (columns)
##
create.inter.matrix <- function(eqtl.file, exp.matrices , sig.thresh = 0.05, output.dir) {
  
  r.squared.list <- pbapply(eqtl.file, 1, function(eqtl) {
    
    cor.matrix <- create.cor.matrix(exp.matrices = exp.matrices, gene.name = eqtl["ProbeName"])
    interaction.statistics <- interaction.regression(cor.matrix = cor.matrix, eqtl.gene = eqtl["ProbeName"], snp.name = eqtl["SNPName"])
    
    return(data.frame(r.squard = interaction.statistics$r.squared, row.names = rownames(cor.matrix)))
  })
  
  #colnames(r.squared.list) <- eqtl.file$ProbeName
  return(r.squared.list)
}

##
## Plots all the significant interations effects of the eQTLs in an eQTL file.
##
plot.interactions.eqtl <- function(eqtl.file, exp.matrices , sig.thresh = 0.05, output.dir) {
  
  pbapply(eqtl.file, 1, function(eqtl) {
    
    cor.matrix <- create.cor.matrix(exp.matrices = exp.matrices, gene.name = eqtl["ProbeName"])
    interaction.statistics <- interaction.regression(cor.matrix = cor.matrix, eqtl.gene = eqtl["ProbeName"], snp.name = eqtl["SNPName"])
    
    interaction.statistics$adjusted.p.value <-  p.adjust(interaction.statistics$p.value, method = "bonf")
    significant.interactions <- interaction.statistics[interaction.statistics$adjusted.p.value <= sig.thresh,]
    if (nrow(significant.interactions) < 1) return(NULL)
    significant.interactions <- significant.interactions[order(significant.interactions$adjusted.p.value),]


    pdf(paste0(output.dir,"/" , eqtl["ProbeName"], ".pdf"), onefile = TRUE)

    print(significant.interactions)

    lapply(rownames(significant.interactions), function(interaction.gene) {

      print(plot.correlations.qtl(eqtl["ProbeName"],
                                  interaction.gene,
                                  snp.id = eqtl["SNPName"],
                                  exp.matrices = exp.matrices,
                                  r = significant.interactions[interaction.gene, "r.squared"],
                                  p.value = significant.interactions[interaction.gene, "adjusted.p.value"]))
    })
    dev.off()
  })
}

r.squared.list <- plot.interactions.eqtl(th.eqtls, samples, output.dir = "")
r.squared.list
r.squared.table <- do.call("cbind")



