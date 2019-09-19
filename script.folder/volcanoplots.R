#Step 1: installation of of relevant software 
 install.packages("ggplot2")
 install.packages("ggrepel")
 #Install ggrepel package if needed
 install.packages("devtools")
 devtools::install_github("slowkow/ggrepel")
 install.packages("dplyr")
 install.packages("extrafont")

# Load packages
library(dplyr)
library(ggplot2)
library(ggrepel)
library(devtools)
library(extrafont)

#Install the desired system fonts (at the moment only TrueType fonts):
font_import(pattern="[C/c]omic")
font_import(pattern="[A/a]rial")
font_import(pattern="[T/t]imes New Roman")

#Reset R's Brain.
rm(list=ls())
#Set working directory.
setwd("/home/icipe/Documents/careenwork/objectivethree/knownmiRNA/volcanoplots /")
#loading data 
data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/knownmiRNA/edgeranalysis/Teneralmale_vs_Adultmale_Analysis.csv", sep = "," , row.names = 1, header = T)
#reshape data 
Teneralfemale_vs_Adultfemale_volcanoplots=(data [ ,c(1, 2, 5, 6)])
head(Teneralfemale_vs_Adultfemale_volcanoplots)
write.csv(Teneralfemale_vs_Adultfemale_volcanoplots, file = "/home/icipe/Documents/careenwork/objectivethree/knownmiRNA/volcanoplots /Teneralmale_vs_Adultmale_volcanoplots")
X<-Teneralfemale_vs_Adultfemale_volcanoplots
rm("data")

# Load Tab delimited data file  from excel into R environment (read.delim = data is tab delimited, 
#header=TRUE = data has header information, ""= path to the  tab delimited data file 
#rna_seq_data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novel_miRNA/volcanoplotfiles/Teneralfemale_vs_teneralmales_Analysis.csv_volcanoplot", header=TRUE) 

# view/confirm the table/data arrangment in the imported datasheet
head(X)

# Create/add a column "sig" on the datasheet for identification of differentially expressed transcripts
# Genes are significantly expressed  if 1) FDR corrected p Value  > 0.05 AND there a 2 fold chmage in expression (log2 2 =1)
X = mutate(X, sig=ifelse(X$FDR <0.05 & abs(logFC)>1, "FDR p value < 0.05", "Not Sig"))

# Create a Volcano plot of log2 fold change vs -log10 pvalue, using the 'sig' column to determine the differentially expressed spots as 'red' differential and 'black' not changing. 
# keep the points at - .3 size big. 

p = ggplot(X, aes(logFC, -log10(PValue))) + geom_point(aes(col=sig),size = 1) + theme_bw(base_size = 9) + coord_cartesian(ylim=c(0,300))+ coord_cartesian(xlim=c(-10,15))  + scale_color_manual(values=c("red","black")) + ggtitle("Differentilaly expressed miRNA in Adultmale compared to Teneralmale ")

# View the plot 
p
# add labels to the most differentially expressed FDR <0.05  and abs(Log2.Fold_Change)>1
par(lwd=0.5)
p+geom_text_repel(data=filter(X, FDR <0.05 & abs(logFC)>1),
                  family = "Times New Roman", aes(label=Row.names),size = 2, 
                  arrow = arrow(length = unit(0.01, 'npc')), 
                  force = 7,box.padding = unit(0.4, "lines"), 
                  point.padding = unit(0.3, "lines"))
                 

