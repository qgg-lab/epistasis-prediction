# simulate phenotypes
# ============================================================

args <- commandArgs(TRUE)

# function to read tped genotypes
# assume that the genotypes have been coded by --recode12
# ============================================================

read.tped <- function(file.prefix) {
  # check file existence
  if (!file.exists(paste(file.prefix, ".tped", sep = "")) || !file.exists(paste(file.prefix, ".tfam", sep = ""))) {
    stop("cannot find tped/tfam files!")
  }
  
  # read genotypes
  tfam <- read.table(paste(file.prefix, ".tfam", sep = ""), as.is = TRUE, header = FALSE)
  n.ind <- nrow(tfam)
  tped <- matrix(scan(paste(file.prefix, ".tped", sep = ""), what = ""), ncol = 4 + 2*n.ind, byrow = TRUE)

  snp.name <- tped[, 2]
  snp.loc <- tped[, c(1, 4)]
  geno.code <- matrix(4 - (as.numeric(tped[, seq(from = 5, length = n.ind, by = 2)]) + as.numeric(tped[, seq(from = 6, length = n.ind, by = 2)])), ncol = n.ind)
  geno.code[geno.code == 4] <- NA
  
  rownames(geno.code) <- snp.name
  colnames(geno.code) <- tfam[, 1]
  
  return(geno.code)
  
}

# simulate phenotype
# ============================================================

geno <- t(read.tped(args[1]))
H2 <- as.numeric(args[2])
out.file1 <- args[3]
out.file2 <- args[4]

# additive phenotype
p <- ncol(geno); n <- nrow(geno)
geno <- geno[, sample(1:p)]
add.beta <- rnorm(p, 0, 1)
add.g <- geno %*% add.beta
add.pheno <- add.g + rnorm(n, 0, sqrt(var(add.g)*(1 - H2)/H2))

# dominance phenotype
dom.geno <- geno
dom.geno[dom.geno == 2] <- 0
dom.g <- add.g + dom.geno %*% (add.beta*sample(c(-1, 1), size = p, replace = T))
dom.pheno <- dom.g + rnorm(n, 0, sqrt(var(dom.g)*(1 - H2)/H2))

# overdominance phenotype
od.beta <- rnorm(p, 0, 1)
od.g <- dom.geno %*% od.beta
od.pheno <- od.g + rnorm(n, 0, sqrt(var(od.g)*(1 - H2)/H2))

# axa phenotype
aa.geno <- (geno[, 1:floor(p/2)] - 1)*(geno[, (floor(p/2) + 1):(2*floor(p/2))] - 1)
aa.beta <- rnorm(floor(p/2), 0, 1)
aa.g <- aa.geno %*% aa.beta
aa.pheno <- aa.g + rnorm(n, 0, sqrt(var(aa.g)*(1 - H2)/H2))

# axd phenotype
ad.geno <- (geno[, 1:floor(p/2)] - 1)*(dom.geno[, (floor(p/2) + 1):(2*floor(p/2))])
ad.beta <- rnorm(floor(p/2), 0, 1)
ad.g <- ad.geno %*% ad.beta
ad.pheno <- ad.g + rnorm(n, 0, sqrt(var(ad.g)*(1 - H2)/H2))

# dxd phenotype
dd.geno <- (dom.geno[, 1:floor(p/2)])*(dom.geno[, (floor(p/2) + 1):(2*floor(p/2))])
dd.beta <- rnorm(floor(p/2), 0, 1)
dd.g <- dd.geno %*% dd.beta
dd.pheno <- dd.g + rnorm(n, 0, sqrt(var(dd.g)*(1 - H2)/H2))

# random phenotype
rand.pheno <- rnorm(n, 0, 1)

# dd causal tped file
# genotype for dxd
dd.tped <- matrix(ncol = nrow(dd.geno)*2 + 4, nrow = floor(p/2))
dd.tped[, 1] <- rep(1, nrow(dd.tped))
dd.tped[, 2] <- paste("1_", 1:nrow(dd.tped), sep = "")
dd.tped[, 3] <- rep(0, nrow(dd.tped))
dd.tped[, 4] <- 1:nrow(dd.tped)

dd.tped[, seq(5, ncol(dd.tped), 2)] <- t(dd.geno) + 1
dd.tped[, seq(6, ncol(dd.tped), 2)] <- t(dd.geno) + 1

write.table(dd.tped, file = out.file2, sep = " ", col.names = F, row.names = F, quote = F)

write.table(cbind(rownames(geno), add.pheno, dom.pheno, od.pheno, aa.pheno, ad.pheno, dd.pheno, rand.pheno), file = out.file1, sep = " ", col.names = F, row.names = F, quote = F)
