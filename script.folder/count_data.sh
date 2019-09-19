#awk '!($1 in sum) {f[n++] = $1}
 #    {sum[$1] += $2}
  #   END {for (i = 0; i < n; i++) print f[i], counts[$1], sum[f[i]]}' < Adfem6.arf> Adfem
cd .././results.folders
cd  ./filtered_sequences
cd  ./miRNA_filtered

for i in Tenmal*.fasta Wtmal*.fasta ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fasta ;
do
 blastn -query "$i" -db miRNA -out "$i".miRna.blast.outfmt6 -num_threads 16  -task blastn-short -max_target_seqs 1 -qcov_hsp_perc 80 -evalue 0.1  -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs qcovhsp" # blast nucleotide 
done
for in Tenmal*.outfmt6 Tenmal*.ou; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.outfmt6 ;
do
  grep -c "$i" > $i.miRNAcount #count number for rRNA 
done
for i in *.outfmt6 ; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.outfmt6 ;
do
cut -f1 "$i" > "$i".filtered
done

for i in *.fasta; do echo "$i" |uniq ; done # identifies the file that the command is work on 
for i in *.fasta;
do
a=$(basename "$i" .fasta)
 grep -A1 -f "$a".fa2.fasta.miRNA.blast.outfmt6.filtered "$i" |awk '/>/{n=1}; n {n--; next}; 1'|sed "s/--//g" |sed "/^\s*$/d"> "$i".miRNA
done

