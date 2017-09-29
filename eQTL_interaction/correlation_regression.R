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
th.eqtls <- read.table("~/Documents/scRNA-seq/R/data/eqtl/cis_100K_MAF_10/th-cells.txt", stringsAsFactors = F, header = T)

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
interaction.regression <- function(cor.matrix, eqtl.gene, snp) {
  
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
create.inter.matrix <- function(eqtl.file, exp.matrices, output.dir, permutations = F, n.perm = 10, fdr.thresh = 0.05) {
  
  if (permutations) {
    perm.sample.orders <- list()
    for (i in 1:n.perm) {
      perm.sample.orders[[i]] <- sample(1:length(exp.matrices), length(exp.matrices), replace = F)
    }
  }
  
  r.squared.matrix <- NULL
  p.value.matrix <- NULL
  
  r.permuted <- list()
  p.value.permuted <- list()
  p.value.tresholds <- NULL
  
  pb <- txtProgressBar(min = 0, max = nrow(eqtl.file), style = 1)
  setTxtProgressBar(pb, 0)
  
  for (i in 1:nrow(eqtl.file)) {
    
    eqtl <- eqtl.file[i,]
    eqtl <- unlist(eqtl)

    cor.matrix <- create.cor.matrix(exp.matrices = exp.matrices, gene.name = eqtl["ProbeName"])
    snp <- as.numeric(get.snp(eqtl["SNPName"]))
    
    interaction.statistics <- interaction.regression(cor.matrix = cor.matrix, eqtl.gene = eqtl["ProbeName"], snp = snp)
    r.squared.matrix <- cbind(r.squared.matrix, interaction.statistics$statistic / sqrt(length(snp) -2 + interaction.statistics$statistic ** 2))
    p.value.matrix <- cbind(p.value.matrix, interaction.statistics$p.value)
    
    print("Starting permutations")
    if (permutations) {
      for (current.perm in 1:n.perm) {
        
        permutated.snp <- snp[perm.sample.orders[[current.perm]]]
        perm.interaction.statistics <- interaction.regression(cor.matrix = cor.matrix, eqtl.gene = eqtl["ProbeName"], snp = permutated.snp)
        r.perm <- perm.interaction.statistics$statistic / sqrt(length(snp) -2 + perm.interaction.statistics$statistic ** 2)
        if (length(r.permuted) < current.perm){
          r.permuted[[current.perm]] = r.perm
          p.value.permuted[[current.perm]] = perm.interaction.statistics$p.value
        } else {
          r.permuted[[current.perm]] <- cbind(r.permuted[[current.perm]], r.perm)
          p.value.permuted[[current.perm]] <- cbind(p.value.permuted[[current.perm]], perm.interaction.statistics$p.value)
        }
      }
      
      p.value.tresh <- NA
      for (current.p.value.thresh in unique(sort(interaction.statistics$p.value, decreasing=T))){
        signif.interactions <- length(which(interaction.statistics$p.value <= current.p.value.thresh))
        permuted.signif.interactions <- c()
        for (current.perm in 1:n.perm){
          permuted.signif.interactions <- c(permuted.signif.interactions, length(which(perm.interaction.statistics$p.value <= current.p.value.thresh)))
        }
        if (mean(permuted.signif.interactions)/signif.interactions <= fdr.thresh){
          p.value.thresh <- current.p.value.thresh
          break
        }
      }
      p.value.tresholds <- c(p.value.tresholds, p.value.thresh)
    
    }
    setTxtProgressBar(pb, i)
    
  }
  
  p.value.matrix <- rbind(p.value.tresholds, p.value.matrix)
  if (permutations){
    rownames(p.value.matrix)[1] = "significance_threshold"
  }
  
  interaction.list <- list(r.squared.matrix, p.value.matrix)
  
  close(pb)
  
  return(interaction.list)
  # r.squared.list <- pbapply(eqtl.file, 1, function(eqtl) {
  #   snp <- as.numeric(get.snp(eqtl["SNPName"]))
  
  #   cor.matrix <- create.cor.matrix(exp.matrices = exp.matrices, gene.name = eqtl["ProbeName"])
  #   interaction.statistics <- interaction.regression(cor.matrix = cor.matrix, eqtl.gene = eqtl["ProbeName"], snp = snp )
  
  
  
  #   return(data.frame(r.squard = interaction.statistics$r.squared, row.names = rownames(cor.matrix)))
  # })
  
  #colnames(r.squared.list) <- eqtl.file$ProbeName
  # return(r.squared.list)
}

inter.matrix <- create.inter.matrix(eqtl.file = th.eqtls[1:3,], exp.matrices = samples,output.dir = "~/Desktop/interactions/", permutations = T)

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
r.squared.table <- do.call("cbind", r.squared.list)
colnames(r.squared.table) <- th.eqtls$ProbeName

r.squared.table[is.na(r.squared.table)] <- 0
length(which(as.numeric(as.matrix(r.squared.table)) == 0))

write.table(r.squared.table, file = "~/Desktop/interactions/interaction_table_t_cd4_eqtls_with_NA.tsv", sep = "\t",  quote = F, col.names = NA)
hist(apply(r.squared.table, 2, max), 100)

