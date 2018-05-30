

library(rcompanion)
library(vcd)
library(DescTools)

rm(list=ls())
gc()

# Begin our data set
# Load data
setwd("C:/Northwell/Research/Projects/CaGMFM")
source("gentab.R")
dat <- as.data.frame(readxl::read_xlsx(path="FM-Peds mutaton project.xlsx"))
nogs <- 29 # 26 genes + 3 statistics
coln <- colnames(dat) # extract column names

# Define the variables for data processing
gname <- NULL # vector of gene name 
pos1 <- 2 #position of element in the vector
pos2 <- 2 #position of column
mdat <- vector(mode="list",length=nogs)
rowlen <- dim(dat)[1] # the number of row
for (i in 1:nogs) {
  mdat[[i]] <- dat[2:rowlen, pos2:(pos2+2)] # subset the data
  for (j in 1:3) # Convert "chr" into "num"
    mdat[[i]][,j] <- as.numeric(mdat[[i]][,j])
  rownames(mdat[[i]]) <- dat[2:rowlen,1] # assign the row name with cancer types
  colnames(mdat[[i]]) <- c("PEDI","AYA","ADULT")
  gname <- c(gname, coln[pos1])
  pos1 <- pos1+3; pos2 <- pos2+3
}
names(mdat) <- gname #Name the row of each dataframes in the list
TT <- mdat$`TOTAL SAMPLES`

# Statistical testing the significances of mutations among different patients groups
cal.pval <- function(dt) {
  # number of row in the data size from the new "list" object
  rowdim <- dim(dt)[1] 
  pv <- vector(mode="numeric",length=(rowdim-1))
  for (k in 1:(rowdim-1)) {
    job <- dt[c(k,rowdim),]
    job <- rbind(job, (job[2,]-job[1,])) # fine the number of non-cancer samples
    job <- job[c(1,3),] #switch rows between 2nd and 3rd
    #Contingency Tables - Fisher's Exact Test
    #need 2 or more non-zero column marginals
    if (sum(apply(job,2,sum)!=0)>=2 && sum(job[1,]>=1)) { 
      ft <- fisher.test(t(job), simulate.p.value = TRUE, B = 1e5)
      pv[k] <- ft$p.value
    } else {
      pv[k] <- 1
    }
  }
  return(pv)
}

# calculate Cochran-Mantel-Haenszel test for repeated tests of independence
cal.CMH <- function(dt, total, CMH=TRUE) {
  # number of row in the data size from the new "list" object
  rowdim <- dim(dt)[1]; coldim <- dim(dt)[2]
  pv <- vector(mode="numeric",length=(rowdim-1))
  CMH.dat <- vector(mode="list",length=coldim)
  # derive the number of non-Cancer
  nonD <- dt; for (i in 1:coldim) nonD[,i] <- nonD[rowdim,i]-nonD[,i]
  # derive the number of non-Gene
  nonG <- dt; nonG <- total - nonG
  # derive the number of both non-Cancer and non-Gene
  nonT <- dt; for (i in 1:coldim) nonT[,i] <- total[rowdim,i] - nonG[,i] - nonD[,i] - dt[,i]
  job <- vector(mode="list",4)
  job[[1]] <- rep(c("PEDS","AYA","ADULT"), each=4)
  job[[2]] <- rep(rep(c("caner","non-cancer"),each=2),3)
  job[[3]] <- rep(c("Gene","non-Gene"),6)
  for (k in 1:(rowdim-1)) {
    #Contingency Tables
    cct <- NULL
    for (h in 1:coldim) {
      cct <- c(cct, dt[k,h], nonG[k,h], nonD[k,h], nonT[k,h])
    }
    job[[4]] <- cct
    job <- as.data.frame(job)
    colnames(job) <- c("Group", "Disease", "Gene", "Count")
    Tab <- xtabs(Count ~  Disease + Gene + Group, data=job)
    
    #need 2 or more non-zero column marginals
    if (CMH) { 
      PT <- mantelhaen.test(Tab, exact = TRUE)
      pv[k] <- PT$p.value
    } else {
      #BDT <- BreslowDayTest(Tab, correct=TRUE)
      #if (BDT$p.value > 0.01) pv[k] <- PT$p.value
      WT <- woolf_test(Tab)
      pv[k] <- WT$p.value
    }
  }
  return(pv)
}

# Main code for data analyses
#dout <- lapply(mdat[1:26],cal.pval)
CMH <- lapply(mdat[1:26], cal.CMH, TT)
WT <- lapply(mdat[1:26], cal.CMH, TT, CMH=FALSE)
CMH.df<- as.data.frame(CMH) # convert list into data frame
WT.df <- as.data.frame(WT)
rownames(CMH.df) <- dat[2:(rowlen-1),1] #Name the row of each dataframes
rownames(WT.df) <- dat[2:(rowlen-1),1]

if (sum(CMH.df==0) >=1) CMH.df[CMH.df==0] <- 0.001
if (sum(WT.df==0) >=1) WT.df[WT.df==0] <- 0.001
# Adjust P-value by FDR
#adj.p.df <- apply(CMH.df,2,p.adjust,method="fdr")
WT.log.df <- -log10(WT.df)
CMH.log.df <- -log10(CMH.df)

## PCA
library(ggfortify)
autoplot(prcomp(dout.df))
autoplot(prcomp(dout.df), data=dout.df, label=TRUE, label.size=3,
         loadings = TRUE, loadings.label = TRUE, loadings.label.size  = 3)

## Heatmep
library(ComplexHeatmap)
library(circlize)
# collect total numer of patient for each gene
SS <- mdat[[1]][72,]
for (i in 2:26) SS <- rbind(SS, mdat[[i]][72,])
rownames(SS) <- names(mdat)[1:26]
# Plot
ha1 = HeatmapAnnotation(dist1 = anno_barplot(TT[1:71,], bar_width = 1, gp = gpar(col = NA, fill=c("red","yellow","green")), 
                                             border = FALSE, axis = TRUE))
ha2 = rowAnnotation(dist2 = anno_barplot(SS, bar_width = 1, gp = gpar(col = NA, fill = c("red","yellow","green")), 
                                         border = FALSE, which = "row", axis = TRUE), width = unit(1, "cm"))
ha_column = HeatmapAnnotation(cn = function(index) {
  cancer = as.numeric(colnames(mat))
  which_decade = which(year %% 10 == 0)
  grid.text(year[which_decade], which_decade/length(year), 1, just = c("center", "top"))
})

Heatmap(t(CMH.log.df), name = "-log10(P-value)", col = colorRamp2(c(0, 5, 10, 15), c("white", "cornflowerblue", "yellow", "red")),
        show_row_dend = FALSE, rect_gp = gpar(col= "white"), show_column_names = TRUE,
        row_names_side = "left", row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 10),
        column_title = 'The ACMG genes associate with pan-cancer cross age groups',
        top_annotation = ha1, top_annotation_height = unit(1, "cm")) + ha2

Heatmap(t(WT.log.df), name = "-log10(P-value)", col = colorRamp2(c(0, 1, 2, 3), c("white", "cornflowerblue", "yellow", "red")),
         show_row_dend = FALSE, rect_gp = gpar(col= "white"), show_column_names = TRUE,
        row_names_side = "left", row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 10),
        column_title = 'The ACMG genes associate with pan-cancer by varing three age groups',
        top_annotation = ha1, top_annotation_height = unit(1, "cm")) + ha2

decorate_heatmap_body("cases", {
  i = which(colnames(mat) == "1961")
  x = i/ncol(mat)
  grid.lines(c(x, x), c(0, 1), gp = gpar(lwd = 2))
  grid.text("Vaccine introduced", x, unit(1, "npc") + unit(5, "mm"))
})





## Examples and testing
# Example for WT1
Job <- matrix(c(5,9, 4,19, 8, 176), 2, 3,
              dimnames = list(cancer = c("AML", "non-AML"),
                              WT1 = c("PEDI", "AYA", "ADULT")))
FT <- fisher.test(t(Job), simulate.p.value = TRUE, 
                  B = 1e5, hybrid = TRUE) # larger than 2x2 table
PT = pairwiseNominalIndependence(FT, fisher = TRUE, gtest  = FALSE, 
                                 chisq  = FALSE, digits = 3)

library(DescTools)
Input =("
        Frequency  WT1       non-WT1
        cp.adult   3         3316
        gp.all     48        138584
        ")

ddt = as.matrix(read.table(textConnection(Input), 
                           header=TRUE, row.names=1))

GTest(ddt, correct="william") 
chisq.test(ddt, correct = TRUE)
# The input for C-M-H test
Input =(
  "                  Group PEDS AYA ADULT
  GMCancer             Gene
  AML     WT1            5      4      8
          non-WT1        194    152    933
  non-AML WT1            9      19     176
          non-WT1        3111   11038  124017
  ")

Tab <- as.table(read.ftable(textConnection(Input)))
mantelhaen.test(Tab, exact = TRUE)
# Mutiple comparisons
n <- dim(Tab)[3]
for(i in 1:n){
  Name = dimnames(Tab)[3]$Group[i]
  P.value = GTest(Tab[,,i])$p.value
  cat(Name, "\n")
  cat("Fisher test p-value: ", P.value, "\n")
  cat("\n")
}

Input =(
  "                  Group PEDS AYA ADULT
GMCancer      Gene
adrenal           WT1        1      0      0
              non-WT1        154    87    243
non-adrenal       WT1        13     23    184
              non-WT1        3151   11103 124716
")




# plotting data for C-M-H test
library(ggplot2)
library(grid)
library(dplyr)
# convering flat table into data frame
ftab <- ftable(Tab)
data <- vector(mode="list",5)
data[[1]] <- rep(c("PEDS","AYA","ADULT"), each=2)
data[[2]] <- rep(c("AML","non-AML"),3)
data[[3]] <- as.integer(as.vector(ftab[c(1,3),])) #as.vector(ftab[c(1,3),])
data[[4]] <- as.integer(as.vector(ftab[c(2,4),])) #as.vector(ftab[c(2,4),])
data[[5]] <- data[[3]]/data[[4]]
Data <- as.data.frame(data)
colnames(Data) <- c("Group", "Cancer", "Count", "Total", "Proportion")


### Add confidence intervals
Fun.low = function (x){
  binom.test(x["Count"], x["Total"],
             0.5)$ conf.int[1]
}

Fun.up = function (x){
  binom.test(x["Count"], x["Total"],
             0.5)$ conf.int[2]
}

Data =
  mutate(Data,
         low.ci = apply(Data[c("Count", "Total")], 1, Fun.low),
         upper.ci = apply(Data[c("Count", "Total")], 1, Fun.up)
  )

### Plotting
ggplot(Data, 
       aes(x = Group, y = Proportion, fill = Cancer, 
           ymax=upper.ci, ymin=low.ci))  +
  geom_bar(stat="identity", position = "dodge", width = 0.7) +
  geom_bar(stat="identity", position = "dodge", 
           colour = "black", width = 0.7, 
           show_guide = FALSE)  +
  scale_y_continuous(breaks = seq(0, 0.80, 0.01), 
                     limits = c(0, 0.10), 
                     expand = c(0, 0))  +
  scale_fill_manual(name = "Count type" , 
                    values = c('grey80', 'grey30'), 
                    labels = c("AML", 
                               "non-AML"))  +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  labs(x = "Age group", 
       y = "Proportion")  +
  ## ggtitle("Main title") + 
  theme_bw()  +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey50"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.key = element_rect(fill = "black"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5)
  )


# population simulation
library(coala)
mr <- 3.11666666667e-05 #*4*10000
model <- coal_model(sample_size = rep(10000,100), loci_number = 1, loci_length = 47856, ploidy = 1)
model <- model + feat_mutation(rate = mr, model = "IFS") +
  feat_recombination(5) +
  feat_migration(0.5, symmetric = TRUE) +
  sumstat_sfs(population = "all")
model
stats <- simulate(model, seed = 20)
