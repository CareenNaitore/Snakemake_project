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
setwd("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/")
data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/edgeranalysis/teneralmale_vs_teneralfemale_Analysis.csv", sep = "," , row.names = 1, header = T)

#reshape data 
Teneralstage=(data [ ,c(1, 2, 5, 6)])
head(Teneralstage)
write.csv(Teneralstage, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/teneralmale_vs_teneralfemale_Analysis.csv")
Teneralstage <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/teneralmale_vs_teneralfemale_Analysis.csv", sep = "," , row.names = 1, header = T)
# view/confirm the table/data arrangment in the imported datasheet
rm("data")
head(Teneralstage)

# Create/add a column "sig" on the datasheet for identification of differentially expressed transcripts
# Genes are significantly expressed  if 1) FDR corrected p Value  > 0.05 AND there a 2 fold change in expression (log2 2 =1)
Teneralstage = mutate(Teneralstage, sig=ifelse(Teneralstage$FDR <0.05 & abs(logFC)>1, "FDR p value < 0.05", "Not Sig"))

setwd("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/")
data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/edgeranalysis/Nonteneralfemale_vs_Nonteneralmale_Analysis.csv", sep = "," , row.names = 1, header = T)

#reshape data 
Nonteneralstage=(data [ ,c(1, 2, 5, 6)])
head(Nonteneralstage)
write.csv(Nonteneralstage, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/Nonteneralfemale_vs_Nonteneralmale_Analysis.csv")
Nonteneralstage<- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/Nonteneralfemale_vs_Nonteneralmale_Analysis.csv", sep = "," , row.names = 1, header = T)
# view/confirm the table/data arrangment in the imported datasheet
rm("data")
head(Nonteneralstage)

Nonteneralstage = mutate(Nonteneralstage, sig=ifelse(Nonteneralstage$FDR <0.05 & abs(logFC)>1, "FDR p value < 0.05", "Not Sig"))
head(Nonteneralstage)

# Create/add a column "sig" on the datasheet for identification of differentially expressed transcripts
# Genes are significantly expressed  if 1) FDR corrected p Value  > 0.05 AND there a 2 fold change in expression (log2 2 =1) 

attach(Teneralstage)
p1 = ggplot(Teneralstage, aes(logFC, -log10(PValue))) + geom_point(aes(col=sig),size = 4) + theme_bw(base_size = 9) + coord_cartesian(ylim=c(0,300))+ coord_cartesian(xlim=c(-15,15))  + scale_color_manual(values=c("red","black")) 
p1 = ggplot(Teneralstage, aes(logFC, -log10(PValue))) + geom_point(aes(col=sig),size = 4) + theme_bw(base_size = 9) + coord_cartesian(ylim=c(0,300))+ coord_cartesian(xlim=c(-15,15))  + scale_color_manual(values=c("red","black"))  
#View plot
p1
# add labels to the most differentially expressed FDR <0.05  and abs(Log2.Fold_Change)>1
par(lwd=0.5)
p1 <-p1+geom_text_repel(data=filter(Teneralstage, FDR <0.05 & abs(logFC)>1),
                        family = "Times New Roman", aes(label=Row.names),size = 9, 
                        arrow = arrow(length = unit(0.01, 'npc')), 
                        force = 7,box.padding = unit(0.4, "lines"), 
                        point.padding = unit(0.3, "lines"))+labs(x="log2 fold change for male compared to female",title="Teneralfemale vs teneral male",
                                                                 y="-log10(PValue)")+theme(plot.title=element_text(size=12, hjust=0.5, face="bold", colour="black", vjust=-1))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme( plot.title=element_text(size=50),
                                                                                                                                                                                                                                                                                                                                                     axis.text=element_text(size=50),
                                                                                                                                                                                                                                                                                                                                                     axis.title=element_text(size=50))+theme(legend.title = element_text( size = 20),legend.text = element_text(size = 20))

p1
                                                                                                                                  
detach(Teneralstage)
attach(Nonteneralstage)
p2 = ggplot(Nonteneralstage, aes(logFC, -log10(PValue))) + geom_point(aes(col=sig),size = 4) + theme_bw(base_size = 9) + coord_cartesian(ylim=c(0,300))+ coord_cartesian(xlim=c(-15,15))  + scale_color_manual(values=c("red","black")) 
p2 = ggplot(Nonteneralstage, aes(logFC, -log10(PValue))) + geom_point(aes(col=sig),size = 4) + theme_bw(base_size = 9) + coord_cartesian(ylim=c(0,300))+ coord_cartesian(xlim=c(-15,15))  + scale_color_manual(values=c("red","black"))  
#View plot
p2
# add labels to the most differentially expressed FDR <0.05  and abs(Log2.Fold_Change)>1
par(lwd=0.5)
p2 <-p2+geom_text_repel(data=filter(Nonteneralstage, FDR <0.05 & abs(logFC)>1),
                        family = "Times New Roman", aes(label=Row.names),size = 9, 
                        arrow = arrow(length = unit(0.01, 'npc')), 
                        force = 7,box.padding = unit(0.4, "lines"), 
                        point.padding = unit(0.3, "lines"))+labs(x="log2 fold change for male compared to female",title="Non-teneral female vs non-teneral male",
                                                                 y="-log10(PValue)")+theme(plot.title=element_text(size=12, hjust=0.5, face="bold", colour="black", vjust=-1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme( plot.title=element_text(size=50),
                                                                                                                                                                                                                                                                                                                                                   axis.text=element_text(size=50),
                                                                                                                                                                                                                                                                                                                                                   axis.title=element_text(size=50))+theme(legend.title = element_text( size = 20),legend.text = element_text(size = 20))
p2
#combining the graphs
library(gridExtra)
grid.arrange(p1, p2,nrow=1,top=textGrob("Differentially expressed sex-specific miRNAs", gp=gpar(fontsize=50, font=25)))


