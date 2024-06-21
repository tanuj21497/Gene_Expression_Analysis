
#1----------------------------------------

library(GEOquery)
gset <- getGEO("GSE25097")

library(limma)
library(ggplot2)
library(enrichR)
library(clusterProfiler)
gset <- gset[[1]]




fvarLabels(gset) <- make.names(fvarLabels(gset))


gsms <- paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "00000000000000000000000000000000000000000000000000",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "11111111111111111111111111111111111111111111111111",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

ex <- exprs(gset)
ex[which(ex <= 0)] <- NaN




#2  and 3 -----------------------------------


exprs(gset) <- log2(ex) # log2 transform

exprs(gset) <- normalizeBetweenArrays(exprs(gset)) # normalize data

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("A","B"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

nall <- nrow(gset)
gset <- gset[complete.cases(exprs(gset)), ]

boxplot(exprs(gset))

v <- vooma(gset, design, plot=T)
plotMA(exprs(gset)[,c(1,3)])
plot(x=exprs(gset)[,1],y=exprs(gset)[,3])

v$genes <- fData(gset)
fit  <- lmFit(v)


fd <- fData(gset)
metaData = pData(gset)

#4--------------------------------------------
m <- c()
p <- c()
subset1_metaData <- metaData[metaData$source_name_ch1 == "non_tumor", ]

subset2_metaData <- metaData[metaData$source_name_ch1 == "tumor", ]
n <- c()


exp1 <- exprs(gset)[, colnames(gset) %in% rownames(subset1_metaData)]
exp2 <- exprs(gset)[, colnames(gset) %in% rownames(subset2_metaData)]

colo <- c()


for(i in (1: nrow(exp1))){
  
  t<- t.test(exp1[i, ], exp2[i, ])
  n <- c(n, t$p.value)
  
  
}
p <- p.adjust(n, method = "holm")
genelist <- c()
for(i in (1:length(p))){
  if(p[i]<0.05){colo <- c(colo, "black")
  genelist <- c(genelist, rownames(exp1)[i])
  }
  else{
    colo <- c(colo, "red")
  }
}

for(i in (1: nrow(exp1))){
  
  u <- mean(exp1[i, ]) - mean(exp2[i, ])
  m <- c(m, u)
  
}



results <- data.frame(
  log_fc = m,
  p_value = p,
  co = colo
)

r <- results[results$co=="black", ]


# Draw a volcano plot
ggplot(results, aes(x = log_fc, y = -log10(p_value), color=co)) +
  geom_point() +
  ggtitle("Volcano Plot") +
  xlab("Log Fold Change") +
  ylab("-log10(P-value)") +
  theme_minimal()


genl <- c()

for(i in (1:length(p))){
  if(p[i]<0.05 & abs(m[i]) >= 1){
    genl <- c(genl, rownames(exp1)[i])
  }
}

#5----------------------------------------------------------------------------
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

dT <- decideTests(fit2, adjust.method="holm", p.value=0.05)

colnames(fit2) 
volcanoplot(fit2, coef=1, main=colnames(fit2)[1], pch=5,
            highlight=length(which(dT[,1]!=0)), names=rep('+', nrow(fit2)))

#7, 8, 9-------------------------------------------
genelist_new <- c()

genelist_new <- c()  # create empty vector for storing matched values
i <- 1  # initialize counter for iterating through fd$ID and fd$EntrezGeneID

while (i <= length(fd[ , 1])) {
  if (fd$ID[i] %in% genl) {
    genelist_new <- c(genelist_new, fd$EntrezGeneID[i])
  }
  
  i <- i + 1
  
}
genelist_new1 <- c()
for (i in 1:length(genl)) {
  if (!is.na(genelist_new[i])) {
    genelist_new1 <- c(genelist_new1, genelist_new[i])
  }
}


# res <- enrichr(genelist_new1, "hsa")
res2 <- enrichKEGG(gene=genelist_new1, organism='hsa', keyType='kegg', pvalueCutoff=0.05, qvalueCutoff=0.2)
dotplot(res2)
# Using the 'enrich' object generated in the previous example
cnetplot(res2,  layoutType='circle', nodeSize='Count')

heatplot(res2, showCategory = 10, pvalue=0.05)
