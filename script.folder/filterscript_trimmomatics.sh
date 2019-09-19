cd  .././results.folders
cd  ./filtered_sequences_trimmomatics
cd  ./miRNA_filtered
#for i in *.trim ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
#for i in *.trim ;
#do
#fastq2fasta.pl "$i" > "$i".fa2 #converts from fastq to fasta file
#done
#for i in *.fa2 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
#for i in *.fa2 ;
#do
#collapse_reads_md.pl  "$i" seq > "$i".fasta #collapses the reads to unique
#done
#for i in *.fasta ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
#for i in *.fasta ;
#do
 #blastn -query "$i" -db miRNA -out "$i".miRna.blast.outfmt6 -num_threads 16  -task blastn-short -max_target_seqs 1 -qcov_hsp_perc 80 -evalue 0.1  -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs qcovhsp" # blast nucleotide 
#done
#for i in *.outfmt6 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
#for i in *.outfmt6 ;
#do
 #cat "$i" | wc -l > $i.miRNAcount #count number for rRNA 
#done
#for i in *.outfmt6 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
#for i in *.outfmt6 ;
#do
#cut -f1 "$i" > "$i".filtered
#done
#for i in *.fasta; do echo "$i" |uniq ; done # identifies the file that the command is work on 
#for i in *.fasta;
#do
#a=$(basename "$i" .fasta)
# grep  -f "$a".fasta.miRna.blast.outfmt6.filtered -A1 "$i" |awk '/>/{n=1}; n {n--; next}; 1'|sed "s/--//g" |sed "/^\s*$/d"> "$i".miRNA
#done
cd  .././rRNA_filtered

#for i in *.trim ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
#for i in *.trim ;
#do
#fastq2fasta.pl "$i" > "$i".fa2 #converts from fastq to fasta file
#done
#for i in *.fa2 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
#for i in *.fa2 ;
#do
#collapse_reads_md.pl  "$i" seq > "$i".fasta #collapses the reads to unique
#done
#for i in *.fasta ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
#for i in *.fasta ;
#do
 #blastn -query "$i" -db rRNA -out "$i".rRNA.blast.outfmt6 -num_threads 16  -task blastn-short -max_target_seqs 1 -qcov_hsp_perc 80 -evalue 0.1  -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs qcovhsp" #blast command for rRNA
#done
#for i in *.outfmt6 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
#for i in *.outfmt6 ;
#do 
 #cat "$i" | wc -l > $i.rRNAcount #count number for rRNA 
#done  
for i in *.outfmt6 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.outfmt6 ;
do
cut -f1 "$i" > "$i".rRNA.filtered
done
for i in *.fasta; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fasta;
do
a=$(basename "$i" .fasta)
 grep  -weF  "$a".fasta.rRNA.blast.outfmt6.rRNA.filtered -A1 "$i" |awk '/>/{n=1}; n {n--; next}; 1'|sed "s/--//g" |sed "/^\s*$/d"> "$i".rRNA #filter sequences that are others 
done
cd  .././tRNA_filtered
#for i in *.trim ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
#for i in *.trim ;
#do
#fastq2fasta.pl "$i" > "$i".fa2 #converts from fastq to fasta file
#done
#for i in *.fa2 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
#for i in *.fa2 ;
#do
#collapse_reads_md.pl  "$i" seq > "$i".fasta #collapses the reads to unique
#done
for i in *.fasta ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fasta ;
do
 blastn -query "$i" -db tRNA -out "$i".tRNA.blast.outfmt6 -num_threads 16  -task blastn-short -max_target_seqs 1 -qcov_hsp_perc 80 -evalue 0.1  -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs qcovhsp" # count number for tRNA
 done
for i in *.outfmt6 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.outfmt6 ;
do
cut -f1 "$i" > "$i".tRNA.filtered
done
for i in *.outfmt6 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.outfmt6 ;
do
  cat "$i" | wc-l > $i.tRNAcount #count number for rRNA 
done 
for i in *.fasta; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fasta;
do
a=$(basename "$i" .fasta)
 grep  -wEf  "$a".fasta.tRNA.blast.outfmt6.tRNA.filtered -A1 "$i" |awk '/>/{n=1}; n {n--; next}; 1'|sed "s/--//g" |sed "/^\s*$/d"> "$i".tRNA
done
cd  .././others_filtered
for i in *.trim ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.trim ;
do
fastq2fasta.pl "$i" > "$i".fa2 #converts from fastq to fasta file
done
for i in *.fa2 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fa2 ;
do
collapse_reads_md.pl  "$i" seq > "$i".fasta #collapses the reads to unique
done

for i in *.fasta ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fasta ;
do
  blastn -query "$i" -db othersfiltered -out "$i".blast.outfmt6 -num_threads 16  -task blastn-short -max_target_seqs 1 -qcov_hsp_perc 80 -evalue 0.1  -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs qcovhsp" #blast command to identify 
done
for i in *.outfmt6; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.outfmt6 ;
do
cat "$i" | wc -l > "$i".otherscount #count number for others
done 
for i in *.outfmt6 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.outfmt6;
do
cut -f1 "$i" > "$i".others.filtered # identifies the file that the command 
done
for i in *.fasta; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fasta;
do
a=$(basename "$i" .fasta)
 grep -wEf  "$a".fasta.blast.outfmt6.others.filtered -A1 "$i" |awk '/>/{n=1}; n {n--; next}; 1'|sed "s/--//g" |sed "/^\s*$/d"> "$i".others  
done
#for i in *.others ;
#do
#a=$(basename "$i" .others)
 #echo "$i" 
 #; done # identifies the file that the command is work on 
#for i in `cat *.others|sort -u`;
#do awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (seq != "$i") {print header, seq, qheader, qseq}}' "$a".fastq   > "$a".others.filt;  done
# awk '{ print "\""$0"\""}' read62.others |awk 'BEGIN{FS=OFS="\""} {gsub(/[[:space:]]/,"",$2)} 1' # adding double quotation and removing space                          

