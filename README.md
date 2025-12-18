# Phylotranscriptomics of Groupers

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

### 2.2. Filtered using Trimmomatic
Raw reads were first filtered using Trimmomatic v.0.39 with a minimum reads length (MINLEN) of 100 to obtain clean reads.
```
for i in Species ; do java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 ${i}_1.fq.gz ${i}_2.fq.gz ${i}_1_paired.fq.gz ${i}_1_unpaired.fq.gz ${i}_2_paired.fq.gz ${i}_2_unpaired.fq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100 ; done
```

### 2.3. Transcriptome assembly
Transcriptome assembly was accomplished using Trinity v.2.8.5 (Grabherr et al., 2011) with --CPU set to 20, --min_kmer_cov set to 2, and all other parameters set to default.
```
for i in Species ; do Trinity --seqType fq --max_memory 50G --left ${i}_1.paired.fq.gz --right ${i}_2.paired.fq.gz --CPU 10 --min_kmer_cov 2 --min_contig_length 200 --output trinity_out_dir_${i} ; done
```

### 2.4. Merge transcriptome
All assemblies of each species were merged into one transcriptome.
```
for i in Species ; do cat trinity_out_dir_${i}1/Trinity.fasta trinity_out_dir_${i}2/Trinity.fasta trinity_out_dir_${i}3/Trinity.fasta > ${i}_total_trinity.fasta ; done
```

### 2.5. Preliminary redundancy reduction
Preliminary redundancy reduction was implemented with CD-HIT-EST v.4.6, implementing a sequence identity threshold of 0.95 and applying a word length of 10.
```
for i in Species; do /opt/biosoft/anaconda3_package/bin/cd-hit-est -i ${i}_total_trinity.fasta -o ${i}_cdhitest.fasta -c 0.95 -n 10 ; done
```

### 2.6. Obtaining the longest transcripts
The longest transcripts were selected using the “get_longest_isoform_seq_per_trinity_gene.pl” implemented in Trinity.
```
for i in Species; do /opt/biosoft/anaconda3_package/bin/misc/get_longest_isoform_seq_per_trinity_gene.pl ${i}_cdhitest.fasta > ${i}_longest.fasta ; done
```

### 2.7. ORFs and CDS predict
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

### 2.8. Local BLASTX searches
Local BLASTX searches were conducted for the predicted CDS to further reduce redundancy using DIAMOND against the NCBI RefSeq of all Actinopterygii datasets (189 species).
```
## mkdir transdecoder_tri_cds and transdecoder_tri_pep files
mkdir transdecoder_tri_cds
mkdir transdecoder_tri_pep

## Cat all cds files of the same species from three databses for Total.transdecoder.cds
for i in Species ; do cat /transdecoder_tri_Nr/${i}_longest.fasta.transdecoder.cds /transdecoder_tri_Pfam/${i}_longest.fasta.transdecoder.cds /transdecoder_tri_Swiss/${i}_longest.fasta.transdecoder.cds > /transdecoder_tri_cds/${i}_Total.transdecoder.cds ; done &

## Cat all pep files of the same species from three databses for Total.transdecoder.pep
for i in Species ; do cat /transdecoder_tri_Nr/${i}_longest.fasta.transdecoder.pep /transdecoder_tri_Pfam/${i}_longest.fasta.transdecoder.pep /transdecoder_tri_Swiss/${i}_longest.fasta.transdecoder.pep > /transdecoder_tri_pep/${i}_Total.transdecoder.pep ; done &

## Download all Actinopterygii datasets (189 species)
# Diamond makedb the 'fish_refseq' for the 'fish_refseq_db'
diamond makedb --in fish_refseq --db fish_refseq_db

# Local BLASTX searches for the refseqblastx.outfmt6
for i in Species ; do diamond blastx -q /transdecoder_tri_cds/${i}_Total.transdecoder.cds -d fish_refseq_db --evalue 1e-5 --max-target-seqs 1 >> ${i}.refseqblastx.outfmt6 ; done &

### Python process_data_${i}.py for the unique_protein
for i in Species ; do python process_data_${i}.py ; done

### Get the id.txt, seqkit grep for the refseq.cds.fa, and sed the sequence names
for i in Species ; do cut -f 1 unique_protein_${i}.refseqblastx.outfmt6 > ${i}_id.txt ; done
for i in Species ; do cat /transdecoder_tri_cds/${i}_Total.transdecoder.cds | seqkit grep -f ${i}_id.txt > ${i}.refseq.cds.fa ; done
for i in Species ; do cat /transdecoder_tri_pep/${i}_Total.transdecoder.pep | seqkit grep -f ${i}_id.txt > ${i}.refseq.pep ; done
for i in Species ; do sed -i "s/TRINITY/${i}/g" ${i}.refseq.cds.fa ; done
for i in Species ; do sed -i "s/TRINITY/${i}/g" ${i}.refseq.pep ; done
for i in Species ; do sed -i "s/TRINITY/${i}/g" ${i}_id.txt ; done

### The final ${i}.refseq.cds.fa and ${i}.refseq.pep were used to Ortholog identification
```

### 2.9. BUSCO completeness assessment
Transcriptome completeness was assessed using BUSCO v.5.5.0 according to the actinopterygii_odb10 database (highly conserved fish orthologs, https://busco.ezlab.org).
```
### The final ${i}.refseq.cds.fa and ${i}.refseq.pep were used to BUSCO assessment
## Download actinopterygii_odb10 from https://busco.ezlab.org
# Obtaing Docker for BUSCO assessment
for i in Species ; do docker run -v $...path:/work -w /work ezlabgva/busco:v5.5.0_cv1 busco -f -i ${i}.refseq.cds.fa -o ${i}_refseq_busco -l actinopterygii_odb10 -m transcriptome -c 8 ; done
```

## 3. Ortholog identification and sequence alignment
### 3.1. Ortholog identification
Orthogroups were identified using OrthoFinder v.2.5.5, applying an all-against-all BLAST algorithm to implement reciprocal best-hit BLASTs, using DIAMOND as the sequence similarity search engine, with the normalization of transcript length.

The phylogenetic framework of Epinephelidae was resolved using transcriptomic data from 32 species of the family and one outgroup from Serranidae (_Centropristis striata_; NCBI accession number: GCF_030273125).
```
## mkdir Orthofinder_33_Speceis
mkdir Orthofinder_33_Speceis
# Download RefSeq of outgroup Centropristis striata (NCBI accession number: GCF_030273125)
# cp Centropristis_striata_GCF_030273125.1_protein.faa to Orthofinder_33_Speceis file
# Ortholog identification, total 33 species, including 32 groupers and 1 outgroup (C. striata)
orthofinder -f /Orthofinder_33_Speceis/ -S diamond 

### Finally, '/Orthofinder_33_Speceis/OrthoFinder/Results_date/Single_Copy_Orthologue_Sequences' file obtain 108 OGs
```

### 3.2. Sequence alignment
Through ortholog identification, a total of 108 one-to-one single-copy orthologous genes (OGs) were identified and used for downstream analysis. All OGs were aligned using MAFFT v.7.453 with the L-INS-I algorithm. Subsequently, aligned OGs were trimmed using trimAl v.1.4 with the --automated1 parameter.
```
## mkdir Sequence_alignment, 01_108OGs_AA, 02mafft, 03trimal, and 04_AA_final
mkdir Sequence_alignment
mkdir Sequence_alignment/01_108OGs_AA
mkdir Sequence_alignment/02mafft
mkdir Sequence_alignment/03trimal
mkdir Sequence_alignment/04_AA_final

# cp 108 OGs to Sequence_alignment/01_108OGs_AA file
cp /Orthofinder_33_Speceis/OrthoFinder/Results_date/Single_Copy_Orthologue_Sequences/*.fa /Sequence_alignment/01_108OGs_AA

## MAFFT alignment
cd Sequence_alignment/01_108OGs_AA
export TMPDIR=/Sequence_alignment/01_108OGs_AA
for i in *.fa; do mafft --localpair --maxiterate 1000 $i > /Sequence_alignment/02mafft/$i ; done &
# rename the *.fa
rename 's/$/\.mafft/'  /Sequence_alignment/02mafft/*
# only retain the sequences of *.fa using AA_retain_sequence.py
python AA_retain_sequence.py

## trimAl trimmed
cd Sequence_alignment/02mafft
for i in *.mafft ; do trimal -in $i -out /Sequence_alignment/03trimal/${i}.trimal -htmlout output.html -automated1 ; done
cd Sequence_alignment/03trimal
for i in *.trimal ; do cat $i | seqkit replace -p "\s.+" >> $i.cat ; done
for i in *.cat ; do sed 's/^\(>.*\)/\1\t/' $i | tr '\n' ' ' | sed -e 's/ $/\n/' -e 's/ >/\n>/g' -e 's/ //g' -e 's/\t/\n/g' > $i.sed ; done

#### only retain the species name for each OGs in the 04fasta_AA file
```

## 4. Phylogenetic analyses
Phylogenetic trees of family Epinephelidae were reconstructed using two different dataset types (nucleotide and amino acid), five distinct phylogenetic inference methods (ML, BI, MSC, NJ, and ME), and two gene integration strategies (concatenation and coalescence). 

### 4.1. Phylogenetic analyses using two different dataset types
#### 4.1.1. Obtaining the corresponding nucleotide CDS
The corresponding nucleotide CDS were obtained and aligned using PAL2NAL v.14.0 against the amino acids of orthologs.
```
## For 32 Epinephelidae species
# mkdir Single_Copy_Orthologue_Sequences
mkdir Single_Copy_Orthologue_Sequences
# copy the 108 Single_Copy_Orthologue to Single_Copy_Orthologue_Sequences file
cp date/Single_Copy_Orthologue_Sequences/*fa Single_Copy_Orthologue_Sequences
# Generate a mapping table linking ortho sequence file names to their corresponding sequence identifiers
for i in `ls /Single_Copy_Orthologue_Sequences/*fa`; do
grep ">" $i | while read id; do
echo $i $id; done
done
# cat all *.refseq.cds.fa into a merge file
cat *.refseq.cds.fa > merge.ortho.original.cds
# Remove the suffixes from sequence identifiers and convert multi-line FASTA sequences into a single-line format
cat merge.ortho.original.cds | seqkit replace -p "\s.+" >> merge.ortho1.cds
sed 's/^\(>.*\)/\1\t/' merge.ortho1.cds | tr '\n' ' ' | sed -e 's/ $/\n/' -e 's/ >/\n>/g' -e 's/ //g' -e 's/\t/\n/g' > merge.ortho.cds
# obtaining the corresponding nucleotide CDS
cat output_cda.tab|while read id;do
	arr=($id);
	file=`basename ${arr[0]}`;
	name=${file%.*};
	grep -A1 ${arr[1]} merge.ortho.cds >> ./Orthogroup_cds_Sequences/${name}.fa;
done
echo "split cds done"
```
```
## For outgroup
# getting the cds of outgroup (C. striata: GCF_030273125)
# Remove the suffixes from sequence identifiers and convert multi-line FASTA sequences into a single-line format (outgroup)
cat cds_GCF_030273125.1.cds | seqkit replace -p "\s.+" >> cds_GCF_merge.ortho1.cds
sed 's/^\(>.*\)/\1\t/' cds_GCF_merge.ortho1.cds | tr '\n' ' ' | sed -e 's/ $/\n/' -e 's/ >/\n>/g' -e 's/ //g' -e 's/\t/\n/g' > cds_GCF_merge.ortho.cds
# obtaining the corresponding nucleotide CDS (outgroup)
cat output_cda.tab|while read id;do
	arr=($id);
	file=`basename ${arr[0]}`;
	name=${file%.*};
	grep -A1 ${arr[1]} cds_GCF_030273125.1.cds >> ./Orthogroup_cds_Sequences/${name}.fa;
done
echo "split cds done"
```
The CDS were obtained, and Sequence alignment was performed following the description in Section 3.2.
```
## mkdir Sequence_alignment, 11_108OGs_nuc, 12mafft, 13trimal, and 14_nuc_final
mkdir Sequence_alignment
mkdir Sequence_alignment/11_108OGs_nuc
mkdir Sequence_alignment/12mafft
mkdir Sequence_alignment/13trimal
mkdir Sequence_alignment/14_nuc_final

# cp the nuc sequence of 108 OGs to Sequence_alignment/11_108OGs_nuc file

## MAFFT alignment
cd /Sequence_alignment/11_108OGs_nuc
export TMPDIR=/Sequence_alignment/11_108OGs_nuc
for i in *.fa; do mafft --localpair --maxiterate 1000 $i > /Sequence_alignment/12mafft/$i ; done &
# rename the *.fa
rename 's/$/\.mafft/'  /Sequence_alignment/12mafft/*
# only retain the sequences of *.fa using nuc_retain_sequence.py
python nuc_retain_sequence.py

## trimAl trimmed
cd /Sequence_alignment/12mafft
for i in *.mafft ; do trimal -in $i -out /Sequence_alignment/13trimal/${i}.trimal -htmlout output.html -automated1 ; done
cd /Sequence_alignment/13trimal
for i in *.trimal ; do cat $i | seqkit replace -p "\s.+" >> $i.cat ; done
for i in *.cat ; do sed 's/^\(>.*\)/\1\t/' $i | tr '\n' ' ' | sed -e 's/ $/\n/' -e 's/ >/\n>/g' -e 's/ //g' -e 's/\t/\n/g' > $i.sed ; done

#### only retain the species name for each OGs in the 14fasta_nuc file
```

#### 4.1.2. The average pairwise distance between species
The average pairwise distance between species was calculated to further quantify the relationship within major clades for both nucleotide and amino acid datasets.
```
#### Rscript p-distance.R

#!/usr/bin/env Rscript

library(ape)

### 0. Read species order
species_order <- scan("id.txt", what = "character")

### 1. Read nucleotide sequences
nuc <- read.dna("Nuc108.fa", format = "fasta", as.character = FALSE)

### Reorder nucleotide sequences
nuc <- nuc[species_order, ]

### 2. Nucleotide p-distance
dist_nuc <- dist.dna(nuc, model = "raw", pairwise.deletion = TRUE)

### Convert dist to full matrix and reorder according to species order
dist_nuc_mat <- as.matrix(dist_nuc)
dist_nuc_mat <- dist_nuc_mat[species_order, species_order]

### 3. Read protein sequences
prot <- read.FASTA("Prot108.fa")

### Reorder protein sequences
prot <- prot[species_order]

### Convert protein sequences to equal-length matrix
prot_list <- lapply(prot, function(x) unlist(strsplit(as.character(x), "")))
max_len <- max(sapply(prot_list, length))

prot_mat <- do.call(rbind, lapply(prot_list, function(x) {
  length(x) <- max_len
  x[is.na(x)] <- "-"
  return(x)
}))

rownames(prot_mat) <- species_order

### 4. Compute protein p-distance manually
n <- nrow(prot_mat)
dist_prot <- matrix(0, n, n)
rownames(dist_prot) <- colnames(dist_prot) <- species_order

for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    s1 <- prot_mat[i, ]
    s2 <- prot_mat[j, ]
    valid <- (s1 != "-" & s2 != "-")
    if (sum(valid) > 0) {
      d <- sum(s1[valid] != s2[valid]) / sum(valid)
    } else {
      d <- NA
    }
    dist_prot[i, j] <- dist_prot[j, i] <- d
  }
}

### Get upper triangle values for summary
prot_values <- dist_prot[upper.tri(dist_prot)]

### 5. Output results
cat("Mean nucleotide p-distance:", mean(dist_nuc_mat), "\n")
cat("Mean protein p-distance:", mean(prot_values, na.rm = TRUE), "\n\n")

cat("Nucleotide summary:\n")
print(summary(as.vector(dist_nuc_mat)))

cat("\nProtein summary:\n")
print(summary(prot_values))

### 6. Save matrices in species_order
write.csv(dist_nuc_mat, "Nucleotide_pdistance_matrix.csv")
write.csv(dist_prot, "Protein_pdistance_matrix.csv")
```
```
#### Rscript plot_heatmap.R

#!/usr/bin/env Rscript

library(ggplot2)
library(tidyr)
library(dplyr)

plot_heatmap <- function(input_file, output_file_base, title_text, id_file) {

    # Read species order
    species_order <- read.table(id_file, stringsAsFactors = FALSE)[,1]

    # Read matrix
    df <- read.csv(input_file, row.names = 1, check.names = FALSE)

    # Keep 3 decimal places as character (retain trailing zeros)
    df <- as.data.frame(lapply(df, function(x) sprintf("%.3f", x)), row.names = rownames(df))

    # Sort rows and columns based on id.txt
    df <- df[species_order, species_order]

    # Convert to long format
    df_long <- df %>%
        mutate(Var1 = rownames(df)) %>%
        pivot_longer(
            cols = -Var1,
            names_to = "Var2",
            values_to = "Distance"
        )

    # Set Var1 and Var2 as factors to fix plotting order
    df_long$Var1 <- factor(df_long$Var1, levels = species_order)
    df_long$Var2 <- factor(df_long$Var2, levels = species_order)

    # Plot (single blue gradient + display values with trailing zeros)
    p <- ggplot(df_long, aes(x = Var1, y = Var2, fill = as.numeric(Distance))) +
        geom_tile(color = "white") +
        geom_text(aes(label = Distance), size = 2) +  # show value as string
        scale_fill_gradient(
            low = "#c6dbef",   # light blue
            high = "#08306b",  # dark blue
            name = "p-distance"
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_text(size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(size = 14, hjust = 0.5)
        ) +
        ggtitle(title_text)

    # Save PNG
    ggsave(paste0(output_file_base, ".png"), p, width = 10, height = 8, dpi = 300)

    # Save SVG
    ggsave(paste0(output_file_base, ".svg"), p, width = 10, height = 8)
}

# Nucleotide distance heatmap
plot_heatmap("Nucleotide_pdistance_matrix.csv",
             "Nucleotide_pdistance_heatmap_values_round3_fixed",
             "Nucleotide p-distance Heatmap",
             "id.txt")

# Protein distance heatmap
plot_heatmap("Protein_pdistance_matrix.csv",
             "Protein_pdistance_heatmap_values_round3_fixed",
             "Protein p-distance Heatmap",
             "id.txt")

```

### 4.2. Phylogenetic analyses using five distinct inference methods
#### 4.2.1. Maximum likelihood (ML) tree
ML phylogenetic inference was performed using IQ-TREE v.2.1.2 with 10,000 bootstrap replicates under the best-fit model determined by ModelFinder, which applied the ModelFinder Plus (MFP) option for -m parameter. Node support was evaluated by the SH-like approximate likelihood ratio test (SH-aLRT) with 10,000 replicates.
```
iqtree2 -s Grouper108OGs.fas -m MFP -nt 40 -bb 10000 -bnni -alrt 10000 -redo
```
#### 4.2.2. Bayesian inference (BI) tree
BI tree was inferred by MrBayes v.3.2.7 with the parameter set as follows: Four Markov Chain Monte Carlo (MCMC) chains (including three heated and one cold chains) were run independently for 10 million generations and sampled every 1,000 generations.
```
mpirun -np 4 --cpus-per-proc 8 mb Grouper108OGs.nex -mcmcappend yes
```
#### 4.2.3. Multispecies coalescent (MSC) tree and individual gene tree
Species tree inference was conducted with the MSC model implemented in ASTRAL-Pro.
```
### Extracting the names of all files in the Orthogroup_cds_Sequences directory, and saving them to PMLid.txt
for file in /Orthogroup_cds_Sequences/*; do
  echo $(basename "$file") >> PMLid.txt
done

### Copy aligned OGs from /Sequence_alignment/14fasta_nuc file to the corresponding Single_tree file
mkdir Single_tree
for i in `cat PMLid.txt`; 
do mkdir /Single_tree/${i}; cp /Sequence_alignment/14fasta_nuc/${i}.fasta /Single_tree/${i}/ ; done

### Reconstrcting the single tree of each OGs
for i in `cat /PMLid.txt`; 
do cd /Single_tree/${i}; iqtree -s ${i}.fasta -m GTR+F+R5 -nt 16 -bb 10000 -bnni -redo > ${i}.log; 
echo "ML tree for ${i} finished"; done

### Copy treefiles
mkdir singletreescopy
for i in `cat /PMLid.txt`; 
do cp /Single_tree/${i}/${i}.fasta.treefile /singletreescopy; 
done

### Rename .fasta.treefile
cd /singletreescopy
rename "s/.fasta.treefile//" *

### Use the pxrr program from Phyx to re-root all tree files (root = C_striata)
ls /singletreescopy/ > tree_id.txt
for i in `cat /tree_id.txt` ; do /phyx/src/pxrr -t /singletreescopy/${i} -g C_striata > /reroot/${i} ; done

### Use the pxcolt program from Phyx to collapse branches with low support
for i in `cat /tree_id.txt` ; do /phyx/src/pxcolt -t /reroot/${i} -l 0.1 > /collapsed/${i} ; done

### Merge all tree files
for i in `cat /tree_id.txt`; do cat /ASTRAL-master/Astral.5.7.8/collapsed/${i} >> /ASTRAL-master/Astral.5.7.8/collapse_genetrees.tre; done

### Individual gene tree
java -jar astral.5.7.8.jar -i collapse_genetrees.tre -o output_species_tree.tre 2> running.log

```

#### 4.2.4. Neighbour-joining (NJ) and Minimum-evolution (ME) tree
NJ and ME analyses were implemented using MEGA v.11.0 with 10,000 pseudo replicates.

### 4.3. Phylogenetic analyses using two gene integration strategies
The files of the concatenated dataset were obtained with the software PhyloSuite v.1.2.3.

The concordance and conflict among gene trees across datasets were evaluated by mapping individual gene trees onto the coalescent-based tree using PhyParts.
```
java -jar phyparts/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar \
-a 1 -v \
-d ASTRAL-master/Astral.5.7.8/collapsed \
-m species.tre \
-o phyparts/ASTRALtest/out
```

The pie chart of each internal branch (node) was summarized and visualized using the phypartspiecharts.py script (https://github.com/mossmatters/phyloscripts/tree/master/phypartspiecharts). 
```
conda activate py37
cd phyparts/ASTRALtest
python phypartspiecharts.py out.concon.tre out 108
```

## 5. Divergence time estimation 
Divergence time estimation was implemented using MCMCtree program in PAML v.4.10.7.

build Sequence108OGs.tree
```
### Sequence108OGs.tree
 33  1

(C_striata,((P_maculatus,P_leopardus),(V_louti,((((C_argus,A_rogaa),C_boenak),((C_sexmaculata,C_miniata),C_urodeta)),((((((E_undulosus,(E_flavocaeruleus,E_cyanopodus)),(E_maculatus,E_bleekeri)),(E_merra,E_fasciatus)),(E_quoyanus,E_trimaculatus)),((E_akaara,(E_fasciatomaculosus,E_awoara)),E_sexfasciatus)),(E_marginatus,(E_bruneus,(((A_leucogrammicus,(E_polyphekadion,E_fuscoguttatus)),((E_corallicola,E_coeruleopunctatus),C_altivelis)),(E_lanceolatus,E_coioides))))'>0.13'))'>0.267<0.547')'>0.223<0.623'));
```

build mcmctree1.ctl file
```
### mcmctree1.ctl
          seed = -2
       seqfile = Sequence108OGs.phy
      treefile = Sequence108OGs.tree
      mcmcfile = Sequence108OGs.txt
       outfile = Sequence108OGs.txt

         ndata = 1
       seqtype = 0    * 0: nucleotides; 1:codons; 2:AAs
       usedata = 3    * 0: no data; 1:seq like; 2:normal approximation; 3:out.BV (in.BV)
         clock = 3    * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge = '<0.65'  * safe constraint on root age, used if no fossil for root.

         model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0    * alpha for gamma rates at sites
         ncatG = 5    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0  * birth, death, sampling
   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

   rgene_gamma = 1 2   * gammaDir prior for rate for genes
  sigma2_gamma = 1 10   * gammaDir prior for sigma^2     (for clock=2 or 3)

      finetune = 1: 0.1  0.1  0.1  0.01 .5 * auto (0 or 1): times, musigma2, rates, mixing, paras, FossilErr

         print = 1   * 0: no mcmc sample; 1: everything except branch rates 2: everything
        burnin = 2000000
      sampfreq = 100
       nsample = 10000000

*** Note: Make your window wider (100 columns) before running the program.
```

copy the tmp0001.ctl with codeml1.ctl, and revise the parameter
```
### codeml1.ctl
seqfile = tmp0001.txt
treefile = codeml1.input.treefile
outfile = tmp0002.out
noisy = 3
model = 0
Small_Diff = 0.1e-6
getSE = 0
method = 1
clock = 1
```

copy the /PAML/paml-4.10.7/dat/wag.dat ./ and copy the tmp0001.trees with codeml1.input.treefile
```
### codeml1.input.treefile
  1

((P_maculatus, P_leopardus), (V_louti, ((((C_argus, A_rogaa), C_boenak), ((C_sexmaculata, C_miniata), C_urodeta)), ((((((E_undulosus, (E_flavocaeruleus, E_cyanopodus)), (E_maculatus, E_bleekeri)), (E_merra, E_fasciatus)), (E_quoyanus, E_trimaculatus)), ((E_akaara, (E_fasciatomaculosus, E_awoara)), E_sexfasciatus)), (E_marginatus, (E_bruneus, (((A_leucogrammicus, (E_polyphekadion, E_fuscoguttatus)), ((E_corallicola, E_coeruleopunctatus), C_altivelis)), (E_lanceolatus, E_coioides))))))'@.36')'@.42', C_striata);
```
codeml codeml1.ctl

copy mcmctree1.ctl with mcmctree2.ctl
```
### mcmctree2.ctl
          seed = -2
       seqfile = Sequence108OGs.phy
      treefile = Sequence108OGs.tree
      mcmcfile = Sequence108OGs.txt
       outfile = Sequence108OGs.txt

         ndata = 1
       seqtype = 0    * 0: nucleotides; 1:codons; 2:AAs
       usedata = 2    * 0: no data; 1:seq like; 2:normal approximation; 3:out.BV (in.BV)
         clock = 3    * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge = '<0.65'  * safe constraint on root age, used if no fossil for root.

         model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0    * alpha for gamma rates at sites
         ncatG = 5    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0  * birth, death, sampling
   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

   rgene_gamma = 1 14   * gammaDir prior for rate for genes
  sigma2_gamma = 1 10   * gammaDir prior for sigma^2     (for clock=2 or 3)

      finetune = 1: 0.1  0.1  0.1  0.01 .5 * auto (0 or 1): times, musigma2, rates, mixing, paras, FossilErr

         print = 1   * 0: no mcmc sample; 1: everything except branch rates 2: everything
        burnin = 2000000
      sampfreq = 100
       nsample = 10000000

*** Note: Make your window wider (100 columns) before running the program.
```
The final treefile was visualized using FigTree v.1.4.4.

## 6. Ancestral area reconstruction
Ancestral area reconstruction and biogeographic inference were performed on the time tree excluding the outgroup using the R package BioGeoBEARS.
Given that E. marginatus is the only Atlantic species represented in the samples of this study, its inclusion could potentially bias biogeographic inference.
Finally, ancestral area reconstruction was performed with BioGeoBEARS using a nucleotide dataset from 136 OGs of 31 Epinephelidae species excluding the outgroup and E. marginatus. 

The biogeographic distribution file was prepared based on Supplementary Table S4 (PAMLtreefile136OGs.phy)
The final PAML treefile was transfered to PAMLtreefile136OGs.new

All codes used for ancestral area reconstruction were sourced from the official BioGeoBEARS GitHub repository (https://github.com/nmatzke/BioGeoBEARS).
Here, showing the BAYAREALIKE model in R
```
# install.packages("magic")
# install.packages("geometry")
# install.packages("vegan")
# install.packages("permute")
# install.packages("optimx")
# install.packages("GenSA")
# install.packages("FD")
# install.packages("snow")
# install.packages("parallel")
# install.packages("rexpokit")
# install.packages("cladoRcpp")
# install.packages("devtools")

library(optimx)
library(GenSA)
library(FD)
library(snow)
library(parallel)
library(rexpokit)
library(cladoRcpp)

library(devtools)
# devtools::install_github(repo="nmatzke/BioGeoBEARS")
library(BioGeoBEARS)

calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
getwd()

# Step 1. Prepare the tree file
# Tree file: the divergence tree obtained from PAML is used
trfn = "PAMLtreefile136OGs.new"  # Input tree file (must be strictly bifurcating; branch lengths must be non-negative; extremely short branches may cause issues; species names must not contain spaces)
moref(trfn)  # Display detailed information of the tree
pdffn = "treeDD.pdf"
pdf(file = pdffn, width = 15, height = 20)  # Export the tree as a PDF file
tr = read.tree(trfn)
tr
plot(tr, cex = 0.5)
title("treeDD.pdf")
axisPhylo()  # Plot the time scale
mtext("Millions of years ago (Ma)", side = 1, line = 2)
dev.off()
cmdstr = paste0("open", pdffn)
system(cmdstr)


# Step 2. Prepare the geographic distribution file
# Geographic information
geogfn = "PAMLtreefile136OGs.phy"  # Input geographic distribution file
moref(geogfn)  # Show geographic data
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn = geogfn)  # Note: the delimiter between numeric and letter-coded areas must be a TAB, not a space
tipranges
max(rowSums(dfnums_to_numeric(tipranges@df)))  # Maximum number of areas occupied simultaneously by any taxon
max_range_size = 5  # Set the maximum ancestral range size (must not exceed the total number of defined areas)
numstates_from_numareas(numareas = 4, maxareas = 4, include_null_range = TRUE)  # See official documentation for alternative range-combination settings


#BAYAREALIKE

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE

BioGeoBEARS_run_object$on_NaN_error = -1e50
BioGeoBEARS_run_object$speedup = TRUE
BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE

BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Psychotria_BAYAREALIKE_M0_unconstrained_v1.Rdata"
res = bears_optim_run(BioGeoBEARS_run_object)
save(res, file=resfn)
resBAYAREALIKE = res

#load(resfn)
#resBAYAREALIKE = res

# Plot the result figures

# Only the visualization of the BAYAREALIKE results is shown here
pdffn = "Psychotria_M0_unconstrained_v1.pdf"
pdf(pdffn, width = 6, height = 6)

analysis_titletxt = "BioGeoBEARS BAYAREALIKE on Psychotria M0_unconstrained"
results_object = resBAYAREALIKE
scriptdir = np(system.file("extdata/a_scripts", package = "BioGeoBEARS"))  # Adjust plot margins according to the tree

res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.25, tipcex=0.3, statecex=0.3, splitcex=0.2, titlecex=0.4, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.25, tipcex=0.3, statecex=0.3, splitcex=0.2, titlecex=0.4, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)
```
