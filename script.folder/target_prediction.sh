#Downloaded from vectorbase  
#extraction of 3 utr region for Glossina morsitans 
#grep -B1 "[^ATGCatgc]" -v Glossina_morsitans.txt2| sed -e 's/--//g' >targetgenes2.txt2

#target prediction by miRanda

#for i in *.txt; do echo "$i" |uniq ; done # identifies the file that the command is work on 
#for i in *.txt;
#do
#a=$(basename "$i"  .txt) #identifies the index for for each library
#miranda "$i" targetgenes2.txt2  -sc 80 -en -14 -keyval  -out miranda/"$a"_miRandaresults
#done

#target prediction by RNAhybrid
#for i in *.txt; do echo "$i" |uniq ; done # identifies the file that the command is work on 
#for i in *.txt;
#do
#a=$(basename "$i"  .txt) #identifies the index for for each library
#RNAhybrid -e -20 -p 0.05 -s 3utr_fly -m 200000 -t targetgenes2.txt2  -q "$i"  -c > rnahybrid/"$a"_rnahybridresults
#done

#Extracting information from miRanda 
#for i in miranda/*_miRandaresults; do echo "$i" |uniq ; done # identifies the file that the command is work on 
#for i in miranda/*_miRandaresults;
#do
#a=$(basename "$i"  _miRandaresults) #identifies the index for for each library

 #grep -A 2 "Score for this Scan:" "$i" | sort | grep '>>' > miranda/data/"$a".miRandatxt

#done 

#Extracting information from miRanda
 
#for i in miranda/data/*.miRandatxt; do echo "$i" |uniq ; done # identifies the file that the command is work on 
#for i in miranda/data/*.miRandatxt;
#do
#a=$(basename "$i"  .miRandatxt) #identifies the index for for each library

 #awk '{print $1,$2}'  "$i" |  awk -F '>>' '{print $1,$2}' | awk -F "|" '{print $1,$2;}' | awk '{print $1,$2}' > miranda/data/venn/"$a".miRandatxt_venndiagram2

#done
 
#Extracting information from RNAHybrid 

for i in rnahybrid/*_rnahybridresults; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in rnahybrid/*_rnahybridresults;
do
a=$(basename "$i"  _rnahybridresults) #identifies the index for for each library

 cut -d ':' -f1,3 "$i" | awk -F ':' '{print $2,$1;}' | awk -F "|" '{print $1,$2;}'| awk '{print $1,$2}' > rnahybrid/venn/"$a".rnahybrid_venndiagram2

done   

# editting cytoscape data 

#awk -F "|" '{print $1,$2;}' UntitledDocument1| awk '{print $1,$2}' >teneralfemale_vs_pupaedownregulatedmiRNAtxt.cytoscape
# here we edit the data generated from interveen diagram to remove an extra name (transcript id) from the target gene any name after the "|" is removed and to reshape data, only use when you need to reshape data 
#for example
#converting this data :
#let-7-3p GMOY000080|GMOY000080-RA to let-7-3p GMOY000080
#for easier analysis 

