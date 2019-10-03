library(dplyr)
library(ggplot2)
library(ggrepel)
library(devtools)
library(extrafont)

#Install the desired system fonts (at the moment only TrueType fonts):
#font_import(pattern="[C/c]omic")
#font_import(pattern="[A/a]rial")
#font_import(pattern="[T/t]imes New Roman")
#Reset R's Brain.
rm(list=ls())
setwd("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/")
data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/edgeranalysis/Larvae_vs_pupae_Analysis.csv", sep = "," , row.names = 1, header = T)

#reshape data 
Larvae_vs_Pupae_volcanoplots=(data [ ,c(1, 2, 5, 6)])
head(Larvae_vs_Pupae_volcanoplots)
write.csv(Larvae_vs_Pupae_volcanoplots, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/Larvae_vs_pupae_Analysis.csv")
larvae_vs_pupae <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/Larvae_vs_pupae_Analysis.csv", sep = "," , row.names = 1, header = T)
# view/confirm the table/data arrangment in the imported datasheet
rm("data")
head(larvae_vs_pupae)

# Create/add a column "sig" on the datasheet for identification of differentially expressed transcripts
# Genes are significantly expressed  if 1) FDR corrected p Value  > 0.05 AND there a 2 fold change in expression (log2 2 =1
larvae_vs_pupae = mutate(larvae_vs_pupae, sig=ifelse(larvae_vs_pupae$FDR <0.05 & abs(logFC)>1, "FDR p value < 0.05", "Not Sig"))

setwd("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/")
data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/edgeranalysis/Pupae_vs_teneralfemale_Analysis.csv", sep = "," , row.names = 1, header = T)

#reshape data 
Pupae_vs_teneralfemale_volcanoplots=(data [ ,c(1, 2, 5, 6)])
head(Pupae_vs_teneralfemale_volcanoplots)
write.csv(Pupae_vs_teneralfemale_volcanoplots, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/pupae_vs_teneralfemale_Analysis.csv")
Pupae_vs_teneralfemale <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/pupae_vs_teneralfemale_Analysis.csv", sep = "," , row.names = 1, header = T)
# view/confirm the table/data arrangment in the imported datasheet
rm("data")
head(Pupae_vs_teneralfemale)

# Create/add a column "sig" on the datasheet for identification of differentially expressed transcripts
# Genes are significantly expressed  if 1) FDR corrected p Value  > 0.05 AND there a 2 fold change in expression (log2 2 =1)
Pupae_vs_teneralfemale = mutate(Pupae_vs_teneralfemale, sig=ifelse(Pupae_vs_teneralfemale$FDR <0.05 & abs(logFC)>1, "FDR p value < 0.05", "Not Sig"))

setwd("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/")
data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/edgeranalysis/Pupae_vs_teneralmale_Analysis.csv", sep = "," , row.names = 1, header = T)

#reshape data 
Pupae_vs_teneralmale_volcanoplots=(data [ ,c(1, 2, 5, 6)])
head(Pupae_vs_teneralmale_volcanoplots)
write.csv(Pupae_vs_teneralmale_volcanoplots, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/Pupae_vs_teneralmale_Analysis.csv")
Pupae_vs_teneralmale <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/Pupae_vs_teneralmale_Analysis.csv", sep = "," , row.names = 1, header = T)
# view/confirm the table/data arrangment in the imported datasheet
rm("data")
head(Pupae_vs_teneralmale)

# Create/add a column "sig" on the datasheet for identification of differentially expressed transcripts
# Genes are significantly expressed  if 1) FDR corrected p Value  > 0.05 AND there a 2 fold change in expression (log2 2 =1)
Pupae_vs_teneralmale = mutate(Pupae_vs_teneralmale, sig=ifelse(Pupae_vs_teneralmale$FDR <0.05 & abs(logFC)>1, "FDR p value < 0.05", "Not Sig"))

setwd("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/")
data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/edgeranalysis/teneralfemale_vs_Adultfemale_Analysis.csv", sep = "," , row.names = 1, header = T)

#reshape data 
teneralfemale_vs_adultfemalevolcanoplots=(data [ ,c(1, 2, 5, 6)])
head(teneralfemale_vs_adultfemalevolcanoplots)
write.csv(teneralfemale_vs_adultfemalevolcanoplots, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/teneralfemale_vs_adulfemale_Analysis.csv")
teneralfemale_vs_adultfemale <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/teneralfemale_vs_adulfemale_Analysis.csv", sep = "," , row.names = 1, header = T)
# view/confirm the table/data arrangment in the imported datasheet
rm("data")

head(teneralfemale_vs_adultfemale)

# Create/add a column "sig" on the datasheet for identification of differentially expressed transcripts
# Genes are significantly expressed  if 1) FDR corrected p Value  > 0.05 AND there a 2 fold change in expression (log2 2 =1)
teneralfemale_vs_adultfemale = mutate(teneralfemale_vs_adultfemale, sig=ifelse(teneralfemale_vs_adultfemale$FDR <0.05 & abs(logFC)>1, "FDR p value < 0.05", "Not Sig"))



setwd("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/")
data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/edgeranalysis/teneralmale_vs_Adultmale_Analysis.csv", sep = "," , row.names = 1, header = T)

#reshape data 
teneralmale_vs_Adultmalevolcanoplots=(data [ ,c(1, 2, 5, 6)])
head(teneralmale_vs_Adultmalevolcanoplots)
write.csv(teneralmale_vs_Adultmalevolcanoplots, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/teneralmale_vs_Adultmale_Analysis.csv")
teneralmale_vs_adultmale <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/teneralmale_vs_Adultmale_Analysis.csv", sep = "," , row.names = 1, header = T)
# view/confirm the table/data arrangment in the imported datasheet
rm("data")
head(teneralmale_vs_adultmale)

# Create/add a column "sig" on the datasheet for identification of differentially expressed transcripts
# Genes are significantly expressed  if 1) FDR corrected p Value  > 0.05 AND there a 2 fold change in expression (log2 2 =1)
teneralmale_vs_adultmale  = mutate(teneralmale_vs_adultmale , sig=ifelse(teneralmale_vs_adultmale$FDR <0.05 & abs(logFC)>1, "FDR p value < 0.05", "Not Sig"))

setwd("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/")
data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/edgeranalysis/nongravid_vs_gravid.csv", sep = "," , row.names = 1, header = T)

#reshape data 
nongravid_vs_gravid_volcanoplots=(data [ ,c(1, 2, 5, 6)])
head(nongravid_vs_gravid_volcanoplots)
write.csv(nongravid_vs_gravid_volcanoplots, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/nongravid_vs_gravid_Analysis.csv")
# view/confirm the table/data arrangment in the imported datasheet
rm("data")
nongravid_vs_gravid_Analysis.csv
nongravid_vs_gravid<- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Volcanoplots/nongravid_vs_gravid_Analysis.csv", sep = "," , row.names = 1, header = T)
# view/confirm the table/data arrangment in the imported datasheet
head(nongravid_vs_gravid)

# Create/add a column "sig" on the datasheet for identification of differentially expressed transcripts
# Genes are significantly expressed  if 1) FDR corrected p Value  > 0.05 AND there a 2 fold change in expression (log2 2 =1)
nongravid_vs_gravid = mutate(nongravid_vs_gravid , sig=ifelse(nongravid_vs_gravid$FDR <0.05 & abs(logFC)>1, "FDR p value < 0.05", "Not Sig"))


attach(larvae_vs_pupae)
p1 = ggplot(larvae_vs_pupae, aes(logFC, -log10(PValue))) + geom_point(aes(col=sig),size = 4) + theme_bw(base_size = 9) + coord_cartesian(ylim=c(0,300))+ coord_cartesian(xlim=c(-15,15))  + scale_color_manual(values=c("red","black")) 
p1 = ggplot(larvae_vs_pupae, aes(logFC, -log10(PValue))) + geom_point(aes(col=sig),size = 4) + theme_bw(base_size = 9) + coord_cartesian(ylim=c(0,300))+ coord_cartesian(xlim=c(-15,15))  + scale_color_manual(values=c("red","black"))  
  #View plot
p1
# add labels to the most differentially expressed FDR <0.05  and abs(Log2.Fold_Change)>1
par(lwd=0.5)
p1 <-p1+geom_text_repel(data=filter(larvae_vs_pupae, FDR <0.05 & abs(logFC)>1),
                      family = "", aes(label=Row.names),size = 9, 
                      arrow = arrow(length = unit(0.00, 'npc')), 
                      force = 7,box.padding = unit(0.2, "lines"), 
                      point.padding = unit(0.01, "lines"))+labs(x="log2 fold change",title="C. Larvae vs. pupae",
                                                               y="-log10(PValue)")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(plot.title=element_text(size=50, face = "bold"),
                                                                                                                                                                                                                                                        axis.text=element_text(size=50),
                                                                                                                                                                                                                                                        axis.title=element_text(size=50))+theme(legend.title = element_text( size = 20),legend.text = element_text(size = 20))+ theme(legend.position = "none") 

p1
detach(larvae_vs_pupae)
attach(Pupae_vs_teneralfemale)
p2 = ggplot(Pupae_vs_teneralfemale, aes(logFC, -log10(PValue))) + geom_point(aes(col=sig),size = 4) + theme_bw(base_size = 9) + coord_cartesian(ylim=c(0,300))+ coord_cartesian(xlim=c(-15,15))  + scale_color_manual(values=c("red","black")) 
p2 = ggplot(Pupae_vs_teneralfemale, aes(logFC, -log10(PValue))) + geom_point(aes(col=sig),size = 4) + theme_bw(base_size = 9) + coord_cartesian(ylim=c(0,300))+ coord_cartesian(xlim=c(-15,15))  + scale_color_manual(values=c("red","black"))
#View plot

p2
# add labels to the most differentially expressed FDR <0.05  and abs(Log2.Fold_Change)>1
par(lwd=0.5)
p2 <-p2+geom_text_repel(data=filter(Pupae_vs_teneralfemale, FDR <0.05 & abs(logFC)>1),
                       family = "", aes(label=Row.names),size = 9, 
                       arrow = arrow(length = unit(0.00, 'npc')), 
                       force = 7,box.padding = unit(0.2, "lines"), 
                       point.padding = unit(0.01, "lines"))+labs(x="log2 fold change",
                                                                 y = NULL,title="D. Pupae vs. teneral females")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(plot.title=element_text(size=50, face = "bold"),
                                                                                                                                                                                                                                                                                                        axis.text=element_text(size=50),
                                                                                                                                                                                                                                                                                                        axis.title=element_text(size=50)) +theme(legend.title = element_text( size = 20),legend.text = element_text(size = 20))+ theme(legend.position = "none")



detach(Pupae_vs_teneralfemale)
attach(Pupae_vs_teneralmale)
p3 = ggplot(Pupae_vs_teneralmale, aes(logFC, -log10(PValue))) + geom_point(aes(col=sig),size = 4) + theme_bw(base_size = 9) + coord_cartesian(ylim=c(0,300))+ coord_cartesian(xlim=c(-15,15))  + scale_color_manual(values=c("red","black")) 
p3 = ggplot(Pupae_vs_teneralmale, aes(logFC, -log10(PValue))) + geom_point(aes(col=sig),size = 4) + theme_bw(base_size = 9) + coord_cartesian(ylim=c(0,300))+ coord_cartesian(xlim=c(-15,15))  + scale_color_manual(values=c("red","black"))  
#View plot
p3
# add labels to the most differentially expressed FDR <0.05  and abs(Log2.Fold_Change)>1
par(lwd=0.5)
p3 <-p3+geom_text_repel(data=filter(Pupae_vs_teneralmale, FDR <0.05 & abs(logFC)>1),
                        family = "", aes(label=Row.names),size = 9, 
                        arrow = arrow(length = unit(0.00, 'npc')), 
                        force = 7,box.padding = unit(0.2, "lines"), 
                        point.padding = unit(0.01, "lines"))+labs(x="log2 fold change",title="E. Pupae vs. teneral males",
                                                                 y="-log10(PValue)") +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme( plot.title=element_text(size=50, face="bold"),
                                                                                                                                                                                                                                                          axis.text=element_text(size=50),
                                                                                                                                                                                                                                                          axis.title=element_text(size=50))+theme(legend.title = element_text( size = 20),legend.text = element_text(size = 20))+ theme(legend.position = "none")

detach(Pupae_vs_teneralmale)
attach(teneralfemale_vs_adultfemale)
p4 = ggplot(teneralfemale_vs_adultfemale, aes(logFC, -log10(PValue))) + geom_point(aes(col=sig),size = 4) + theme_bw(base_size = 9) + coord_cartesian(ylim=c(0,300))+ coord_cartesian(xlim=c(-15,15))  + scale_color_manual(values=c("red","black")) 
p4 = ggplot(teneralfemale_vs_adultfemale, aes(logFC, -log10(PValue))) + geom_point(aes(col=sig),size = 4) + theme_bw(base_size = 9) + coord_cartesian(ylim=c(0,300))+ coord_cartesian(xlim=c(-15,15))  + scale_color_manual(values=c("red","black"))  
#View plot
p4
# add labels to the most differentially expressed FDR <0.05  and abs(Log2.Fold_Change)>1
par(lwd=0.5)
p4 <-p4+geom_text_repel(data=filter(teneralfemale_vs_adultfemale, FDR <0.05 & abs(logFC)>1),
                        family = "", aes(label=Row.names),size = 9, 
                        arrow = arrow(length = unit(0.00, 'npc')), 
                        force = 7,box.padding = unit(0.2, "lines"), 
                        point.padding = unit(0.01, "lines"))+labs(x="log2 fold change",title= "F. Teneral vs. non-teneral females",
                                                                  y= NULL) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme( plot.title=element_text(size=50, face = "bold"),
                                                                                                                                                                                                                                                         axis.text=element_text(size=50),
                                                                                                                                                                                                                                                         axis.title=element_text(size=50))+theme(legend.title = element_text( size = 20),legend.text = element_text(size = 20))+ theme(legend.position = "none")
detach(teneralfemale_vs_adultfemale)
attach(teneralmale_vs_adultmale)
p5 = ggplot(teneralmale_vs_adultmale, aes(logFC, -log10(PValue))) + geom_point(aes(col=sig),size = 4) + theme_bw(base_size = 9) + coord_cartesian(ylim=c(0,300))+ coord_cartesian(xlim=c(-15,15))  + scale_color_manual(values=c("red","black")) 
p5 = ggplot(teneralmale_vs_adultmale, aes(logFC, -log10(PValue))) + geom_point(aes(col=sig),size = 4) + theme_bw(base_size = 9) + coord_cartesian(ylim=c(0,300))+ coord_cartesian(xlim=c(-15,15))  + scale_color_manual(values=c("red","black"))  
#View plot
p5
# add labels to the most differentially expressed FDR <0.05  and abs(Log2.Fold_Change)>1
par(lwd=0.5)
p5 <-p5+geom_text_repel(data=filter(teneralmale_vs_adultmale, FDR <0.05 & abs(logFC)>1),
                        family = "", aes(label=Row.names),size = 9, 
                        arrow = arrow(length = unit(0.01, 'npc')), 
                        force = 7,box.padding = unit(0.2, "lines"), 
                        point.padding = unit(0.01, "lines"))+labs(x="log2 fold change ",title="G. Teneral vs non-teneral males",
                                                                 y="-log10(PValue)" )+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme( plot.title=element_text(size=50 , face = "bold"),
                                                                                                                                                                                                                                                         axis.text=element_text(size=50),
                                                                                                                                                                                                                                                         axis.title=element_text(size=50))+theme(legend.title = element_text( size = 20),legend.text = element_text(size = 20))+ theme(legend.position = "none") 
p5
detach(teneralmale_vs_adultmale)
attach(nongravid_vs_gravid)
p6 = ggplot(nongravid_vs_gravid, aes(logFC, -log10(PValue))) + geom_point(aes(col=sig),size = 4)  + theme_bw(base_size = 9) + coord_cartesian(ylim=c(0,300))+ coord_cartesian(xlim=c(-15,15))  + scale_color_manual(values=c("red","black")) 
#View plot
p6
# add labels to the most differentially expressed FDR <0.05  and abs(Log2.Fold_Change)>1
par(lwd=0.5)
p6 <-p6+geom_text_repel(data=filter(nongravid_vs_gravid, FDR <0.05 & abs(logFC)>1),
                        family = "", aes(label=Row.names),size = 9, 
                        arrow = arrow(length = unit(0.01, 'npc')), 
                        force = 7,box.padding = unit(0.2, "lines"), 
                        point.padding = unit(0.01, "lines"))+labs(x="log2 fold change",title="H. Gravid vs non-gravid females",
                                                                 y= NULL) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+ theme(legend.title = element_text( size = 20),legend.text = element_text(size = 20))+theme(plot.title=element_text(size=50, face = "bold"),
                                                                                                                                                                                                                                                                                                                                                axis.text=element_text(size=50),
                                                                                                                                                                                                                                                                                                                                              axis.title=element_text(size=50))+ theme(legend.position = "none")

p6
#combining the graphs
library(gridExtra)
library(grid)
grid.arrange(p1, p2, p3, p4, p5 ,p6)
g <- arrangeGrob(p1, p2, p3, p4, p5, p6, nrow=2, top = "Differentilaly expressed stage-specific miRNAs " ) #generates g
library(devEMF)
emf('imPic.emf')
print(g)
dev.off()

#save

ggsave(file="whatever.pdf", g) #saves g