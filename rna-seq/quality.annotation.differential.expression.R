library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(org.Hs.eg.db)

# read count matrix
expr = load(...)

myCPM <- cpm(expr)
thresh <- myCPM > 0.5
table(rowSums(thresh))
keep <- rowSums(thresh) >= 2
counts.keep <- expr[keep,]

plot(myCPM[,1],expr[,1])
abline(v=0.5)
y <- DGEList(counts.keep)
plotMDS(y)
y <- calcNormFactors(y)
barplot(y$samples$lib.size, names=colnames(y), las=2)
logcounts <- cpm(y, log=TRUE)

# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")


# define here your groups of interest
gpe1 = c(...)
gpe2 = c(...)

design = model.matrix(~ 0 + gpe1 + gpe2)
v <- voom(y, design, plot = TRUE)
fit <- lmFit(v)

cont.matrix <- makeContrasts(Gpe1vsGpe2=gpe1-gpe2, levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)

# annotation
ann <- select(org.Hs.eg.db, keys=rownames(fit.cont), columns=c("SYMBOL"), keytype="ENTREZID")
fit.cont$genes <- ann

# result
topTable(fit.cont,coef="CD4NaiveVsMemory",sort.by="logFC", n=50)
volcanoplot(fit.cont, coef=1, highlight=100, names=fit.cont$genes$SYMBOL)