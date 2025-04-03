# Tracing-non-fungal-eukaryotic-diversity-via-shotgun-metagenomes-in-the-complex-mudflat-intertides

R code and some supplementary information for manuscript:
*"Tracing non-fungal eukaryotic diversity via shotgun metagenomes in the complex mudflat intertides"*

## Pipeline

An example workflow for the integrative approach,which combined different approaches for recovering eukaryotic sequences from shotgun metagenomes.

### Requirements:
* EukRep
* Tiara
* KMA
* CCMetagen
* SeqKit
* Kraken2
* SAMtools
* CoverM
* TaxonKit
* Csvtk

The tools can be installed via conda.

### Classify with EukRep, Tiara and CCMetagen

* Run EukRep on the pre-assembled shotgun metagenomic sample
```bash
EukRep -i input_contig.fa -o output-eukrep.fa
```
* Run Tiara on the pre-assembled shotgun metagenomic sample
```bash
tiara -i input_contig.fa -o output-tiara.txt --tf all -t 30 -p 0.65 0.60 --probabilities
```
* Run KMA and CCMetagen on the clean shotgun metagenomic reads
```bash
kma -ipe input_1.clean.fq.gz input_2.clean.fq.gz -o output -t_db cc_ncbi_nt/ncbi_nt_no_env_11jun2019 -t 4 -1t1 -mem_mode -and -apm f -ef
CCMetagen.py -i output.res -o ccmetagen -map output.mapstat -ef y
CCMetagen_merge.py -i ccmetagen -t Species -kr r -l Superkingdom -tlist Bacteria,Archaea,Virus -o output_cc.csv
```

### Annotating contigs with Kraken2

* Before that, the contigs obtained from EukRep and Tiara need to be de-duplicated and merged
```bash
awk -F " " '{print $1}' output-eukrep.fa > output-eukrep-clean.fa
awk -F " " '{print $1}' eukarya_tiara.fa > output-tiara-clean.fa
seqkit common output-eukrep-clean.fa output-tiara-clean.fa -o et-common.fa
grep ">" et-common.fa |sort|uniq|sed 's/>//g' > et-common-id.txt
seqkit grep -v -f et-common-id.txt output-tiara-clean.fa > output-tiara-new.fa
cat output-tiara-new.fa LYG.output-eukrep-clean.fa > output-et.fa
```
* Annotation with Kraken2
```bash
kraken2 --db kraken2/nt output-et.fa --threads 10 --output output-kranken2.txt  --memory-mapping --report output-kreport.txt --use-mpa-style --use-names sort
```

### Reads mapping

* Add the sample name to the contigs id to avoid confusion and combine all sample contigs as reference
```bash
grep ">" output-et.fa  |sort|uniq|sed 's/>//g' > id.txt
seqkit replace --ignore-case --kv-file rename.txt --pattern "^(\w+)" --replacement "{kv}" output-et.fa -o output-et-new.fa
cat *.fa > total.fa
```
* Getting reads mapping information with CoverM
```bash
coverm contig -1 input_1.clean.fq.gz -2 input_2.clean.fq.gz -r total.fa -p bwa-mem --min-read-percent-identity 95 --min-read-aligned-percent 90 -m count --min-covered-fraction 0 -t 20 --bam-file-cache-directory /coverm --discard-unmapped -o coverm.tsv
```
* De-duplication and merging with eukaryotic reads obtained from CCMetagen
```bash
samtools view -h input.bam > output.sam
awk '$1 !~ /^@/ {print $1, $3}' output.sam > output-id.txt
for i in *-id.txt; do cat $i | sort | uniq > $i-id-uniq.txt;done
for i in *.bam;do samtools fasta $i > $i.fa;done
seqkit common output-et.fa output-cc.fa -s -i -o output-common.fa
grep ">" output-common.fa |sort|uniq|sed 's/>//g' > common-id.txt
seqkit grep -n -v -f common-id.txt output-cc.fa > uncommon.fa
awk -F " " '{print $1}' uncommon.fa  > uncommon-clean.fa
grep ">" uncommon-clean.fa |sort|uniq|sed 's/>//g' > uncommon-id.txt
grep -wf uncommon-id.txt output.frag > uncommon-taxa.txt
```
* Processing the data to obtain a table of taxonomic information
```R
data <- read.csv("id.csv")
result <- aggregate(.~samples, data=data, sum)
write.table (result, file ="id.csv",sep =",", quote =FALSE)//R
```
```bash
cat TAXID.txt | taxonkit lineage -i 1 -r | tee lineage.txt
```
```bash
cat name.txt | taxonkit name2taxid | taxonkit lineage -i 2 -r -L > taxid.txt
```
* Get classification level information based on taxid (final canonical name)
```bash
cat taxid.txt \
> | taxonkit lineage \
>     | taxonkit reformat -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" -F -P \
>     | csvtk cut -t -f -2 \
>     | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species \
>     | csvtk pretty -t > reformat.txt
```
> **Note:** For inquiries, please [open a GitHub issue](https://github.com/HeHan-hub/Tracing-non-fungal-eukaryotic-diversity-via-shotgun-metagenomes-in-the-complex-mudflat-intertides/issues) or contact [me via email](mailto:hanhe0606@foxmail.com).  <dr>
> R scripts were written and run for R v4.4.0.
