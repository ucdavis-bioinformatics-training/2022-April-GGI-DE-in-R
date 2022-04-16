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

#### Quiz 3
* Based on the above model, how many genes are significantly differentially expressed between naive-like and memory-like?: `length(which(top.table$adj.P.Val < 0.05))`
* Based on the above model, and without taking significance into account, how many genes have higher expression in naive-like than in memory-like?: `length(which(top.table$logFC > 0))`
* How many genes have an _unadjusted_ p-value less than 0.05 for the comparison of naive to memory-like in the above model?: `length(which(top.table$P.Value < 0.05))`
* What is the adjusted p-value for the last gene with unadjusted P < 0.05?: `top.table$adj.P.Val[max(which(top.table$P.Value < 0.05))]`

#### Quiz 4
* For the model ~0 + simplified_cell_type, how many genes are differentially expressed between naive-like and effector-like?: 
```
mm.1 <- model.matrix(~0 + simplified_cell_type, data = metadata)
y.1 <- voom(d, mm.1); fit.1 <- lmFit(y.1, mm.1)
tmp.1 <- contrasts.fit(fit.1, contrasts = makeContrasts(simplified_cell_typenaive_like - simplified_cell_typeeffector_like, levels = colnames(coef(fit.1))))
tmp.1 <- eBayes(tmp.1)
length(which(topTable(tmp.1, n = Inf, sort.by = "P")$adj.P.Val < 0.05))
```

### Enrichment
#### Quiz 1
* Rerun the KS test analysis using the molecular function (MF) ontology.  What is the top GO term listed?:
```
GOdata.1 <- new("topGOdata", ontology = "MF", allGenes = geneList, geneSelectionFun = function(x)x, annot = annFUN.org , mapping = "org.Mm.eg.db")
resultKS.1 <- runTest(GOdata.1, algorithm = "weight01", statistic = "ks")
tab.1 <- GenTable(GOdata.1, raw.p.value = resultKS.1, topNodes = length(resultKS.1@score), numChar = 120)
tab.1[1, "Term"]
```

* How many genes from the top table are annotated with the term 'protein serine/threonine kinase activity'?: `tab.1[which(tab.1[,'Term'] == 'protein serine/threonine kinase activity'),'Annotated']`

#### Quiz 2
* How many pathways have a p-value less than 0.05?: `length(which(outdat$p.value < 0.05))`
* Which pathway has the most genes annotated to it (excluding genes not in the top table)?: `outdat$pathway.name[which.max(outdat$Annotated)]`
* Make a pathview diagram for mmu04740: 
```
foldChangeList <- tmp$logFC
xx <- as.list(org.Mm.egENSEMBL2EG)
names(foldChangeList) <- xx[tmp$Gene.stable.ID]
head(foldChangeList)

mmu04740 <- pathview(gene.data  = foldChangeList,
                     pathway.id = "mmu04740",
                     species    = "mmu",
                     limit      = list(gene=max(abs(foldChangeList)), cpd=1))
```
