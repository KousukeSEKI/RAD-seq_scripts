sink("Path to /CIM_trait.txt")

library(qtl)
cross <- read.cross(format="csv", file="Path to/csv.csv", na.strings=c("-","NA","N"), genotypes=c("A","H","B"), alleles=c("A","B"), estimate.map=F)
summary(cross)


interval <- 1
cross <- calc.genoprob(cross, step=1, error.prob=0.001, map.function=c("kosambi"))
cross <- sim.geno(cross, step=1, n.draw=1000, error.prob=0.001, map.function=c("kosambi"))
trait <- #trait_number
n.covar <- # the number of covariate markers
window.size <- 10 
outcim.hk <- cim(cross, pheno.col = trait, method = "hk", n.marcovar = n.covar, window = window.size)
plot(outcim.hk)
opermcim.hk <- cim(cross, pheno.col = trait, method = "hk", n.marcovar = n.covar, window = window.size, n.perm = 10000)

summary(opermcim.hk, alpha = 0.01)
summary(opermcim.hk, alpha = 0.05)

temp1 <- summary(outcim.hk, perms = opermcim.hk, alpha = 0.01) 
print(temp1)
temp5 <- summary(outcim.hk, perms = opermcim.hk, alpha = 0.05) 
print(temp5)

qtl1 <- makeqtl(cross, chr = temp1$chr, pos = temp1$pos, what = "prob")
res1 <-fitqtl(cross, pheno.col=trait,qtl=qtl1,method=c("hk"),model=c("normal"), dropone=TRUE, run.checks=TRUE, maxit=1000,get.ests=TRUE)
summary(res1)

qtl5 <- makeqtl(cross, chr = temp5$chr, pos = temp5$pos, what = "prob")
res5 <-fitqtl(cross, pheno.col=trait,qtl=qtl5,method=c("hk"),model=c("normal"), dropone=TRUE, run.checks=TRUE, maxit=1000,get.ests=TRUE)
summary(res5)

sink()
