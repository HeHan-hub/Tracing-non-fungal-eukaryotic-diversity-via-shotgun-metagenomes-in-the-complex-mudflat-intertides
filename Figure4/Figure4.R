library(devtools)
suppressWarnings(suppressMessages(library(amplicon)))
library(vegan)
library(dada2)
library(ggplot2)
library(picante)
setwd("E:/18S")

#dada2
list.files(path)
R1files_paths <- sort(list.files(path, pattern=".raw_1.fq",full.names = TRUE))
sample.names <- sapply(strsplit(basename(R1files_paths), "\\."), '[', 1)
R1files_noprimer <- file.path(path, "noprimer", 
                              paste0(sample.names, ".raw_1.noprimer.fastq"))

for (i in 1:85) {
  removePrimers(R1files_paths[i], R1files_noprimer[i], 
                primer.fwd="GCGGTAATTCCAGCTCCAA", 
                max.mismatch=2,orient=TRUE, trim.fwd=TRUE, 
                compress=FALSE, verbose=TRUE)
}

R2files_paths <- sort(list.files(path, pattern=".raw_2.fq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(R2files_paths),"\\."),'[',1)
R2files_noprimer <- file.path(path, "noprimer", 
                              paste0(sample.names, ".raw_2.noprimer.fastq"))

for (i in 1:85) {
  removePrimers(R2files_paths[i], R2files_noprimer[i], 
                primer.fwd="AATCCRAGAATTTCACCTCT", 
                max.mismatch=2,orient=TRUE, trim.fwd=TRUE, 
                compress=FALSE, verbose=TRUE)
}

#cd /mnt/e/18S
#perl matchseq.pl

fnFs <- sort(list.files(path, pattern=".raw_1.noprimer.fastq",full.names = TRUE))
fnRs <- sort(list.files (path,pattern=".raw_2.noprimer.fastq",full.names=TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "\\."), '[', 1)
sample.names

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(0,0),
                     maxN=0,maxEE=c(1.5,1.5),truncQ=2,rm.phix=TRUE,
                     compress=TRUE,multithread=FALSE,matchIDs = TRUE)
head(out)

errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
dadaFs <- dada(derepFs, err=errF, multithread=FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=FALSE)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

taxa <- assignTaxonomy(seqtab.nochim, paste0(path, "/pr2_version_5.0.0_SSU_dada2.fasta.gz"))
taxa.print <- taxa 
rownames(taxa.print) <- NULL
head(taxa.print)

setwd(path)
seqtable.taxa.plus <- cbind('#seq'=rownames(taxa), t(seqtab.nochim), taxa)
seqtab.nochim <- t(seqtab.nochim)

write.csv(seqtab.nochim, "dada2_counts.csv", quote=F, row.names = T)
write.csv(taxa,"18S-taxa.csv")
write.table(track , "dada2_track.txt", sep="\t", quote=F, row.names = F)


#Figure 4c-d
otu = read.csv('dada2_counts.csv', header=T, row.names=1) 
colSums(otu)
otu_Flattening = as.data.frame(t(rrarefy(t(otu), min(colSums(otu)))))
colSums(otu_Flattening)
write.table (otu_Flattening, file ="dada2_counts-flattening.csv",sep =",", quote =FALSE) #结果导出

otu <- read.csv("eukgenus-time-normalised.csv", row.names = 1)
#otu <- read.csv("dada2_counts_genus.csv", row.names = 1)

otu <- as.data.frame(t(otu))
richness <- rowSums(otu > 0)
richness
Shannon <- diversity(otu, index = 'shannon', base = exp(1))
Shannon
pielou <- Shannon/log(richness , base = exp(1))
pielou
diversity <- data.frame(richness,Shannon,pielou, Gini_simpson,goods_coverage)
write.csv(diversity, 'alpha_diversity-genus.csv', quote = FALSE)

a = read.csv('dada2_counts-flattening.csv', header=T, row.names=1)
a = t(a)
write.table (a, file ="dada2_counts-flattening-t.csv",sep =",", quote =FALSE)

metadata = read.table("group2.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)
head(metadata, n = 3)

alpha_div = read.csv('Richness.csv', header=T, row.names=1)
head(alpha_div, n = 3)
metadata$group = factor(metadata$group, levels = c("Metagenome","Amplicon(18S)"))
(p = alpha_boxplot(alpha_div, metadata, "Richness", "group"))

alpha_div = read.table("pielou.txt", header=T, row.names=1, sep="\t", comment.char="")
head(alpha_div, n = 3)
metadata$group = factor(metadata$group, levels = c("Metagenome","Amplicon(18S)"))
(p = alpha_boxplot(alpha_div, metadata, "pielou", "group"))


#Figure 4e
otu <- read.csv("eukspecies-time-normalised.csv", row.names = 1) 
otu <- t(otu)
dis <- vegan::vegdist(otu, method = 'bray')
dis <- as.matrix(dis)
write.table(dis, 'Bray-curtis.txt', sep = '\t', col.names = NA, quote = FALSE)
dis <- read.delim('Bray-curtis-1.txt', row.names = 1)

dis1 <- read.delim('Bray-curtis-meta.txt', row.names = 1)
dis2 <- read.delim('Bray-curtis-18S.txt', row.names = 1)
dis <- dplyr::bind_rows(dis1, dis2)
group <- read.delim('group2.txt', stringsAsFactors = FALSE)

env1 <- subset(group, group == 'Metagenome')$samples
dis_env1 <- dis[env1,env1]
env2 <- subset(group, group == 'Amplicon(18S)')$samples
dis_env2 <- dis[env2,env2]
dis_env1 <- as.vector(as.dist(dis_env1))
dis_env2 <- as.vector(as.dist(dis_env2))

dat <- data.frame(
  dis = c(dis_env1, dis_env2),
  group = factor(c(
    rep('Metagenome', length(dis_env1)), 
    rep('Amplicon(18S)', length(dis_env2))
  ), levels = c('Metagenome', 'Amplicon(18S)'))
)

p <- ggplot(dat, aes(group, dis)) +
  geom_boxplot(aes(fill = group), width = 0.4) +
  scale_fill_manual(values = c('#b03525', '#1b6baa')) +
  geom_jitter(size=0.4, alpha=0.9) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'), legend.position = 'none') +
  labs(x = NULL, y = 'Bray-Curtis dissimilarity\n')
p