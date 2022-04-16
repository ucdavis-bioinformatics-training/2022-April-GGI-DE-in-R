## R code for selected quiz answers

Note: There are multiple ways of doing just about anything in R, the code shown below is one of many ways.

### DE

#### Quiz 1
* How many genes are in the counts table?: `nrow(counts)`
* How many samples are in the counts table?: `ncol(counts)`
* What is the total count across all genes for sample SRR8245064?: `colSums(counts)["SRR8245064"]`
* How many genes have a count of 0 in every sample?: `sum(rowSums(counts) == 0)`

#### Quiz 2
* Which sample has the largest normalization factor?: `rownames(d0$samples)[which.max(d0$samples$norm.factors)]`
* Is the sample with the largest normalization factor the sample with the smallest total counts?: `rownames(d0$samples)[which.min(d0$samples$lib.size)] == rownames(d0$samples)[which.max(d0$samples$norm.factors)]`
* Make an MDS plot of the unfiltered data.: `plotMDS(d0, col = as.numeric(factor(metadata$simplified_cell_type)))`

### Quiz 3
* Based on the above model, how many genes are significantly differentially expressed between naive-like and memory-like?: `length(which(top.table$adj.P.Val < 0.05))`
* Based on the above model, and without taking significance into account, how many genes have higher expression in naive-like than in memory-like?: `length(which(top.table$logFC > 0))`
* 

### Enrichment

#### Quiz 1
