#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(survival))

CTL_genes = c('CD8A', 'CD8B', 'GZMA', 'GZMB', 'PRF1')

readmat = function(mat) as.matrix(read.table(mat, sep='\t', header=T, check.names=F, quote=NULL))

writemat = function(result, output)
{
  result = result[!is.na(result[,"p"]),,drop=F]
  FDR = p.adjust(result[,"p"], method="fdr")
  result = cbind(result, FDR)
  write.table(result, file=output, sep='\t', quote=F)
}

commands = commandArgs(trailingOnly=T)

# three inputs here
mat = t(readmat(commands[1])) # expression, mutation, or RPPA matrix
pivot = t(readmat(commands[2]))  # CD8A, CD8B, GZMA, GZMB, PRF1 expression
survival = readmat(commands[3]) # clinical information

# remove negative survival values
survival = survival[survival[,1] > 0,]

output = commands[4]

common = Reduce(intersect, list(rownames(mat),rownames(survival),rownames(pivot)))
print(paste(length(common), "samples"))

# stop at low sample numbers
if(length(common) < 20) stop("two few samples")

not_included = setdiff(CTL_genes, colnames(pivot))

if(length(not_included) > 0) stop(paste(c("pivot genes", not_included, "are not included"), ' '))

pivot = rowMeans(pivot[common,CTL_genes])
mat = mat[common,,drop=F]
survival = survival[common,,drop=F]

# stop at low death rate
death_rate = sum(survival[,2])/dim(survival)[1]
if(death_rate < 0.1) stop(paste("death rate", death_rate, "is too low"))

# split up survival and background
surv = Surv(survival[,1], survival[,2])

if(dim(survival)[2] > 2){
  B = survival[,3:dim(survival)[2], drop=F]
}else{
  B = survival[,c(), drop=F]
}

# build up regression data space
B = cbind(B, pivot, pivot, pivot)
B = as.data.frame(B)
N_B = dim(B)[2]

# test pivot effect first
coxph.pivot = summary(coxph(surv~., data=B[,c(-N_B, -(N_B-1)), drop=F]))$coef
write.table(coxph.pivot, file=paste0(output, ".pivot"), sep='\t', quote=F)

colnames(B)[N_B-1] = "partner"
colnames(B)[N_B] = "Interaction"

# iterate over features
features = colnames(mat)
N = length(features)

result_interaction = matrix(nrow=N, ncol=2, dimnames=list(features, c('z','p')))
result_main = matrix(nrow=N, ncol=2, dimnames=list(features, c('z','p')))
result_partial = matrix(nrow=N, ncol=2, dimnames=list(features, c('z','p')))
result_base = matrix(nrow=N, ncol=2, dimnames=list(features, c('z','p')))

for (i in 1:N)
{
  title = features[i]
  partner = mat[,i]
  
  B[,N_B-1] = partner
  B[,N_B] = partner * pivot
  
  # region 1: model with interaction
  errflag = F
  coxph.fit = tryCatch(coxph(surv~., data=B),
                       error = function(e) errflag <<- T,
                       warning = function(w) errflag <<- T)
  
  if(!errflag)
  {
    reg.summary = summary(coxph.fit)$coef
    result_interaction[i,] = reg.summary["Interaction", c("z", "Pr(>|z|)")]
    result_main[i,] = reg.summary["partner", c("z", "Pr(>|z|)")]
  }
  
  # region 2: model without interaction
  errflag = F
  coxph.fit = tryCatch(coxph(surv~., data=B[,-N_B]),
                       error = function(e) errflag <<- T,
                       warning = function(w) errflag <<- T)
  
  if(!errflag)
  {
    reg.summary = summary(coxph.fit)$coef
    result_partial[i,] = reg.summary["partner", c("z", "Pr(>|z|)")]
  }
  
  # region 3: model with only base effect
  errflag = F
  coxph.fit = tryCatch(coxph(surv~., data=B[,c(-N_B, -(N_B-2)), drop=F]),
                       error = function(e) errflag <<- T,
                       warning = function(w) errflag <<- T)
  
  if(!errflag)
  {
    reg.summary = summary(coxph.fit)$coef
    result_base[i,] = reg.summary["partner", c("z", "Pr(>|z|)")]
  }
}

writemat(result_interaction, paste0(output, ".interaction"))
writemat(result_main, paste0(output, ".main"))
writemat(result_partial, paste0(output, ".partial"))
writemat(result_base, paste0(output, ".base"))
