### Call variants in a genome

#cd ~/Downloads/
#~/Downloads/sratoolkit.2.10.9-mac64/bin/fasterq-dump-orig.2.10.9 SRR5943504
#gzip SRR5943504_1.fastq
#gzip SRR5943504_2.fastq
#export PATH=$PATH:~/Downloads/bowtie2-2.3.5.1-macos-x86_64/
#~/Downloads/breseq-0.35.5-MacOSX-10.9+/bin/breseq -r ~/Downloads/genome_assemblies_genome_gb-4/ncbi-genomes-2021-03-06/GCF_002247485.1_ASM224748v1_genomic.gbff ~/Downloads/SRR5943504_1.fastq.gz ~/Downloads/SRR5943504_2.fastq.gz -o ~/Downloads/sonnei -j 12

### Figure 4ab

setwd("~/Desktop/Sequencing/Assembly/scripts/") # path to the main directory

library(DECIPHER)

# RA	NZ_NQBD01000011	253	C→T	pseudogene (31/273 nt)	CI730_RS03385 ←	IS1 family transposase
file <- "./breseq_sonnei/Pseudo_253.txt"
revC <- TRUE
readingFrame <- 2

# RA	1,352,127	(T)7→6	pseudogene (38/1657 nt)	ubiB →	ubiquinone biosynthesis regulatory protein kinase UbiB
#file <- "./breseq_noatunensis/Pseudo_1352127.txt"
#revC <- FALSE
#readingFrame <- 1

colors <- c(`-`="white",
	`A`="#3A6919",
	`C`="#EBD62F",
	`G`="#211819",
	`T`="#834184")

x <- readLines(file)
x <- gsub(" *(<|>).*$", "", x)
x <- gsub(" |‑", "", x)
x <- DNAStringSet(x[-2])
if (revC) {
	x <- reverseComplement(x)
	x[-1] <- rev(x[-1])
}
X <- AlignSeqs(x)
t <- TerminalChar(X)
X <- subseq(X,
	min(t[-1, 1]) + 1,
	-(min(t[-1, 2]) + 1))
s <- strsplit(as.character(X), "", fixed=TRUE)

plot(NA,
	xlim=c(0, width(X)[1]),
	ylim=c(0, length(X)),
	bty="n",
	xlab="",
	ylab="",
	xaxt="n",
	yaxt="n")
for (i in seq_along(s)) {
	if (i == 1L) {
		rect(seq_len(length(s[[i]])) - 1,
			length(X) - i + 1,
			seq_len(length(s[[i]])),
			length(X) - i,
			col=colors[s[[i]]],
			border=NA)
	} else {
		rect(seq_len(length(s[[i]])) - 1,
			length(X) - i + 1,
			seq_len(length(s[[i]])),
			length(X) - i,
			col=ifelse(s[[i]] != s[[1L]],
				colors[s[[i]]],
				"gray"),
			border=NA)
	}
}

t <- translate(subseq(RemoveGaps(X[1]),
	readingFrame))
t <- strsplit(as.character(t), "", fixed=TRUE)[[1]]
W <- which(t == "*")
pos <- which(s[[1]] != "-")
if (length(W) > 0) {
	text(x=pos[3*(W - 1) + 1 + readingFrame] - 0.5,
		y=length(X) + strheight("*") - 1,
		"*",
		xpd=TRUE,
		font=2)
}
