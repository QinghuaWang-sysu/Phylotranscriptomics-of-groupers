## 1. Data collection
Raw data newly sequenced in this study can be download from NCBI Sequence Read Archive (SRA) with accession number PRJNA1344424.

Raw data sequenced in our previous study can be download from NCBI SRA with accession number SRP483481 and PRJCA002372.

Raw data collected from the previous study can be download from NCBI SRA with accession number PRJNA609274.

## 2. Transcriptome assembly, ORFs and CDS predict, and completeness assessment
### 2.1. Quality control
The raw sequencing RNA-seq reads were visualized and quality-checked using FastQC v0.11.9.
```
# conda install -c bioconda fastqc
for i in Species; do fastqc ${i}_1.fq.gz ; done
for i in Species; do fastqc ${i}_2.fq.gz ; done
```

### 2.2 Filtered using Trimmomatic
Raw reads were first filtered using Trimmomatic v.0.39 with a minimum reads length (MINLEN) of 100 to obtain clean reads.
```
for i in Species ; do java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 ${i}_1.fq.gz ${i}_2.fq.gz ${i}_1_paired.fq.gz ${i}_1_unpaired.fq.gz ${i}_2_paired.fq.gz ${i}_2_unpaired.fq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100 ; done
```

### 2.3 Transcriptome assembly
Transcriptome assembly was accomplished using Trinity v.2.8.5 (Grabherr et al., 2011) with --CPU set to 20, --min_kmer_cov set to 2, and all other parameters set to default.
```
for i in Species ; do Trinity --seqType fq --max_memory 50G --left ${i}_1.paired.fq.gz --right ${i}_2.paired.fq.gz --CPU 10 --min_kmer_cov 2 --min_contig_length 200 ; done
```

### 2.4 Preliminary redundancy reduction
Preliminary redundancy reduction was implemented with CD-HIT-EST v.4.6.
```
for i in Species; do /opt/biosoft/anaconda3_package/bin/cd-hit-est -i ${i}_total_trinity.fasta -o ${i}_cdhitest.fasta -c 0.95 -n 10 ; done
```

### 
The longest transcripts were selected using the “get_longest_isoform_seq_per_trinity_gene.pl” implemented in Trinity.
```
for i in Species; do /opt/biosoft/anaconda3_package/bin/misc/get_longest_isoform_seq_per_trinity_gene.pl ${i}_cdhitest.fasta > ${i}_longest.fasta ; done
```

Longest isoforms were used to predict open reading frames (ORFs) and coding sequences (CDS) using TransDecoder v.5.5.0 against the NCBI Non-Redundant Protein Sequence (Nr), UniProt Knowledge Base (UniProtKB/Swiss-Prot), and Protein Family (Pfam) databases with an E-value threshold of 1e-5.
```
### Longest isoforms were used to predict open reading frames (ORFs) and coding sequences (CDS) using TransDecoder v.5.5.0
for i in Species; do /TransDecoder-v5.5.0/TransDecoder.LongOrfs -t ${i}_longest.fasta ; done
## mkdir transdecoder_tri_Nr, transdecoder_tri_Swiss, and transdecoder_tri_Pfam files
mkdir transdecoder_tri_Nr
mkdir transdecoder_tri_Swiss
mkdir transdecoder_tri_Pfam
# cp pep and cds files to transdecoder_tri_Nr, transdecoder_tri_Swiss, and transdecoder_tri_Pfam files
for i in Species; do cp -r /${i}_longest.fasta.transdecoder_dir /transdecoder_tri_Nr ; done
for i in Species; do cp -r /${i}_longest.fasta.transdecoder_dir /transdecoder_tri_Swiss ; done
for i in Species; do cp -r /${i}_longest.fasta.transdecoder_dir /transdecoder_tri_Pfam ; done
for i in Species; do cp -r /${i}_longest.fasta.transdecoder_dir.__checkpoints_longorfs /transdecoder_tri_Nr ; done
for i in Species; do cp -r /${i}_longest.fasta.transdecoder_dir.__checkpoints_longorfs /transdecoder_tri_Swiss ; done
for i in Species; do cp -r /${i}_longest.fasta.transdecoder_dir.__checkpoints_longorfs /transdecoder_tri_Pfam ; done

## (1) Against the NCBI Non-Redundant Protein Sequence (Nr)
# Download RefSeq from NR database, and cat all RefSeq to 'nr'
# Diamond makedb the 'nr' for the 'nr_db'
diamond makedb --in nr -d nr_db
# Diamond blastp the nr_db for the Nr_blastp.outfmt6
for i in Species; do diamond blastp -d nr_db -q /${i}_longest.fasta.transdecoder_dir/longest_orfs.pep --evalue 1e-5 --max-target-seqs 1 > /transdecoder_tri_Nr/${i}_Nr_blastp.outfmt6 ; done
# Predict ORFs and CDS from NR database
for i in Species; do /TransDecoder-v5.5.0/TransDecoder.Predict -t ${i}_longest.fasta --retain_blastp_hits /transdecoder_tri_Nr/${i}_Nr_blastp.outfmt6 ; done

## (2) Against the UniProt Knowledge Base (UniProtKB/Swiss-Prot)
# Download uniprot_sprot from Swiss-Prot database
# Diamond makedb the 'uniprot_sprot' for the 'uniprot_sprot_db'
diamond makedb --in uniprot_sprot -d uniprot_sprot_db
# Diamond blastp the uniprot_sprot.fasta for the Swiss_blastp.outfmt6
for i in Species; do diamond blastp -d uniprot_sprot_db -q /transdecoder_tri_Swiss/${i}_longest.fasta.transdecoder_dir/longest_orfs.pep --evalue 1e-5 --max-target-seqs 1 > /transdecoder_tri_Swiss/${i}_Swiss_blastp.outfmt6 ; done
# Predict ORFs and CDS from NR database
for i in Species; do /TransDecoder-v5.5.0/TransDecoder.Predict -t ${i}_longest.fasta --retain_blastp_hits /transdecoder_tri_Swiss/${i}_Swiss_blastp.outfmt6 ; done

## (3) Against the Protein Family (Pfam)
# Download Trinotate-3.2.1 # https://github.com/Trinotate/Trinotate/releases
# Hmmscan for the new.out.tbl and the new.pfam.domtblout
for i in Species; do hmmscan --cpu 2 --tblout /transdecoder_tri_Pfam/${i}.new.out.tbl --domtblout /transdecoder_tri_Pfam/${i}.new.pfam.domtblout -E 1e-5 Trinotate/Pfam-A.hmm ${i}.fasta.transdecoder_dir/longest_orfs.pep ; done
# Predict ORFs and CDS from Pfam database
for i in Species; do TransDecoder-v5.5.0/TransDecoder.Predict -t ${i}_longest.fasta --retain_pfam_hits /transdecoder_tri_Pfam/${i}.new.pfam.domtblout ; done
```






