#!/usr/bin/env Rscript
####
# Tingfen Yan and Narisu Narisu
# 03/10/2024
####

args=commandArgs(trailingOnly = T)
infile <- args[1]
outfile <- args[2]

.libPaths(c("/cluster/ifs/users/narisu/anaconda/envs/rna/lib/R/library"))
library(ggplot2)
library(tidyr)
library(ggrepel)
library(gridExtra)
library(cowplot)

options(stringsAsFactors = F)

f=read.delim(infile)
# Calculate -log10(p value) for T2D and ASI
f$credible_set_ASI[!is.na(f$credible_set_ASI)]
f$sum_stat= -log10(f$T2D_pvalue)
f$credible_set = f$ppa
f$credible_set_ASI = -log10(f$asi)

f1=f[c(3:4,8:10)]
head(f1)
# Rename the column name
colnames(f1) <- c("end","SNP","T2D association","Credible set","Credible set ASI")
data_long <- gather(f1, condition, measurement, "T2D association":"Credible set ASI", factor_key=TRUE)

p <- ggplot(data = data_long, aes(x = end, y = measurement,col=condition)) + geom_point(size=1) 
# Subset T2D association
t2d.mtx <-data_long[data_long$condition=="T2D association",]
# Subset credible set
cs.mtx <- data_long[data_long$condition=="Credible set",]
# Subset credible set ASI
asi.mtx=data_long[data_long$condition=="Credible set ASI",]
# Subeset data with -log10(pvalue) > -log10(0.05)=1.3
asi.mtx1=asi.mtx[ !is.na(asi.mtx$measurement) & asi.mtx$measurement >4,]

# Subset t2d.mtx and cs.mtx by the rsId got from asi.mtx1
t2d.mtx1 <- t2d.mtx[t2d.mtx$SNP %in% asi.mtx1$SNP,]
cs.mtx1 <- cs.mtx[cs.mtx$SNP %in% asi.mtx1$SNP,]

t2d= ggplot(data_long[data_long$condition=="T2D association",],aes(x = end, y = measurement,col=condition)) + 
    geom_point(size=0.5,col="orange") +
    geom_point(data=t2d.mtx1,aes(x = end, y = measurement),col="red",size=0.5) +
    facet_grid(condition ~ ., ) + 
    theme(panel.background = element_rect(fill = 'white', colour = 'black'), legend.position = "none") +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'), legend.position = "none") +
    labs(x="", y=expression('-log'[10]*'(P value)')) +
    geom_vline(xintercept=asi.mtx1$end, linetype="dotted") +
    theme(
        plot.margin = unit(c(0.2, 0.2, -0.45, 0.2), "cm"),
        axis.title.y = element_text(hjust = 0.25, size=4),
        axis.text=element_text(size=6),
        strip.text.y = element_text( size = 4))

max_ppa=max(data_long[data_long$condition=="Credible set",]$measurement) + 0.05
cs= ggplot(data_long[data_long$condition=="Credible set",],aes(x = end, y = measurement,col=condition)) + 
    geom_point(size=0.5,col="blue") +
    geom_point(data=cs.mtx1,aes(x = end, y = measurement),col="red",size=0.5) +
    facet_grid(condition ~ ., ) + 
    coord_cartesian(ylim=c(-0.01,max_ppa)) +
    theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'), legend.position = "none") +
    labs(x="", y="    PPA") +
    geom_vline(xintercept=asi.mtx1$end, linetype="dotted") +
    theme(
        plot.margin = unit(c(0.1, 0.2, -0.35, 0.2), "cm"),
        axis.title.y = element_text(hjust = 0.25, size=4),
        axis.text=element_text(size=6),
        strip.text.y = element_text(
        size = 4))
asi= ggplot(data_long[data_long$condition=="Credible set ASI",],aes(x = end, y = measurement)) + 
    geom_point(size=0.5) +
    geom_point(data=asi.mtx1,aes(x = end, y = measurement),col="red",size=0.5) +
    geom_text_repel(data=asi.mtx1,size=1.5, aes(label=SNP),hjust = 1) +
    facet_grid(condition ~ ., ) + 
    theme(panel.background = element_rect(fill = 'white', colour = 'black'), legend.position = "none") +
    labs(x="", y=expression('-log'[10]*'(P value)'))  +
    geom_vline(xintercept=asi.mtx1$end, linetype="dotted") +
    theme(
        plot.margin = unit(c(0, 0.2, -0.35, 0.2), "cm"),
        axis.title.y = element_text(hjust = 0.25, size=4),
        axis.text=element_text(size=6),
        strip.text.y = element_text(
        size = 4))
    
  pdf(outfile,width=3, height=2)
  plot_grid(t2d, cs, asi, align = "v", nrow = 3)
  dev.off()

  
