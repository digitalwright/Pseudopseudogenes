setwd("~/Downloads/SupplementalData/") # path to the main directory

data_path <- "./all-prokaryotes_v2.RData" # where data table is located
stats_path <- "./stats_v1.RData" # where to store the stats data
anno_path <- "./anno_v1.RData" # where to store the annotation data
colors <- c(total="#00000033",
	lines="green4",
	frame="#B0130033",
	stop="#00169933",
	points="gray",
	transposase="#FF0000",
	transporter="#40E0D0",
	`hypothetical protein`="#9932CC",
	other="#FFDEAD")

#r <- read.csv(data_path, quote='"', comment="") # for data downloaded from NCBI
load(data_path) # load the variable `r` from compressed file

# subset to those with PGAP annotations
w <- which(as.numeric(r$tRNA) > 0 & as.numeric(r$rRNA) > 0 & as.numeric(r$CDS) > 0 & as.numeric(r$Size.Mb.) > 0 & as.numeric(r$GC.) > 0 & as.numeric(r$Scaffolds) > 0)

### Download assembly statistics

if (file.exists(stats_path)) {
	load(stats_path)
} else {
	ftps <- r$GenBank.FTP
	https <- paste0("https", substring(ftps, 4))
	https <- paste0(https,
		"/",
		sapply(strsplit(https, "/", fixed=TRUE), tail, n=1),
		"_assembly_stats.txt")
	pBar <- txtProgressBar(style=3)
	stats <- vector("list", length(https))
	options(timeout=2, warn=1)
	tf <- tempfile()
	for (i in w[lengths(stats[w])==0]) {
		t <- try(download.file(https[i], tf, quiet=TRUE), silent=TRUE)
		if (!is(t, "try-error"))
			t <- try(stats[[i]] <- readLines(tf),
				silent=TRUE)
		setTxtProgressBar(pBar, i/length(https))
	}
	
	save(stats, file=stats_path, compress="xz")
}

### Download annotation data

if (file.exists(anno_path)) {
	load(anno_path)
} else {
	pBar <- txtProgressBar(style=3)
	anno <- vector("list", nrow(r))
	options(timeout=30, warn=1)
	library(parallel)
	missing <- w[lengths(anno[w])==0]
	batches <- floor(seq(1, length(missing), length.out=1000))
	for (i in seq_along(batches)[-1]) {
		anno[missing[batches[i - 1]:batches[i]]] <- mclapply(paste0("https://www.ncbi.nlm.nih.gov/assembly/", r$Assembly[missing[batches[i - 1]:batches[i]]]),
			function(x) {
				tf <- tempfile()
				t <- try(download.file(x, tf, quiet=TRUE),
					silent=TRUE)
				if (!is(t, "try-error"))
					t <- try(readLines(tf, warn=FALSE),
						silent=TRUE)
				if (!is(t, "try-error")) {
					grep("Pseudo Genes ", t, value=TRUE)
				} else {
					character()
				}
			},
			mc.cores=detectCores())
		setTxtProgressBar(pBar, i/length(batches))
	}
	
	save(anno, file=anno_path, compress="xz")
}

cds_tot <- as.numeric(sapply(sapply(sapply(anno, gsub, pattern=".*(CDSs \\(with protein\\)|Genes \\(coding\\))</span><span>::</span><span>([0-9,]+)</span>.*", replacement="\\2"), gsub, pattern=",", replacement=""),`[`, 1L))
pseudo_tot <- as.numeric(sapply(sapply(sapply(anno, gsub, pattern=".*Pseudo Genes \\(total\\)</span><span>::</span><span>([0-9,]+)</span>.*", replacement="\\1"), gsub, pattern=",", replacement=""),`[`, 1L))
pseudo_ambig <- as.numeric(sapply(sapply(sapply(anno, gsub, pattern=".*Pseudo Genes \\(ambiguous residues\\)</span><span>::</span><span>([0-9,]+) of .*", replacement="\\1"), gsub, pattern=",", replacement=""),`[`, 1L))
pseudo_frame <- as.numeric(sapply(sapply(sapply(anno, gsub, pattern=".*Pseudo Genes \\(frameshifted\\)</span><span>::</span><span>([0-9,]+) of .*", replacement="\\1"), gsub, pattern=",", replacement=""),`[`, 1L))
pseudo_incom <- as.numeric(sapply(sapply(sapply(anno, gsub, pattern=".*Pseudo Genes \\(incomplete\\)</span><span>::</span><span>([0-9,]+) of .*", replacement="\\1"), gsub, pattern=",", replacement=""),`[`, 1L))
pseudo_stop <- as.numeric(sapply(sapply(sapply(anno, gsub, pattern=".*Pseudo Genes \\(internal stop\\)</span><span>::</span><span>([0-9,]+) of .*", replacement="\\1"), gsub, pattern=",", replacement=""),`[`, 1L))
pseudo_multi <- as.numeric(sapply(sapply(sapply(anno, gsub, pattern=".*Pseudo Genes \\(multiple problems\\)</span><span>::</span><span>([0-9,]+) of .*", replacement="\\1"), gsub, pattern=",", replacement=""),`[`, 1L))

### Split violin plot function

split_violin <- function(X, Y1, Y2, colors, YMAX=max(Y1, Y2, na.rm=TRUE), adjust=1, las=2, ...) {
	u_x <- sort(unique(X))
	m <- match(X, u_x)
	
	plot(NA,
		xlim=c(0.5, length(u_x) + 0.5),
		ylim=c(0, YMAX),
		xaxt="n",
		...)
	axis(1, seq_along(u_x), u_x, las=las)
	#abline(v=seq_along(u_x))
	
	for (i in seq_along(u_x)) {
		w <- which(m==i & !is.na(Y1))
		if (length(w) > 1) {
			d1 <- density(Y1[w], adjust=adjust, from=min(Y1[w]), to=max(Y1[w]))
			x1 <- d1$x
			y1 <- d1$y
			y1 <- y1/sum(y1)*length(y1)
			y1 <- y1/max(y1)/2 + i
			polygon(c(y1, i, i),
				c(x1, max(x1), min(x1)),
				col=colors[1],
				border=substring(colors[1], 1, 7))
			segments(i,
				median(Y1[w]),
				y1[which.min(abs(median(Y1[w]) - x1))],
				median(Y1[w]),
				col=substring(colors[1], 1, 7))
		}
		
		w <- which(m==i & !is.na(Y2))
		if (length(w) > 1) {
			d2 <- density(Y2[w], adjust=adjust, from=min(Y2[w]), to=max(Y2[w]))
			x2 <- d2$x
			y2 <- d2$y
			y2 <- y2/sum(y2)*length(y2)
			y2 <- -y2/max(y2)/2 + i
			polygon(c(y2, i, i),
				c(x2, max(x2), min(x2)),
				col=colors[2],
				border=substring(colors[2], 1, 7))
			segments(i,
				median(Y2[w]),
				y2[which.min(abs(median(Y2[w]) - x2))],
				median(Y2[w]),
				col=substring(colors[2], 1, 7))
		}
	}
}

### Figure 1

dev.new()
layout(matrix(1:4, 2))
p <- par(mar=c(4, 4, 1, 1), mgp=c(2.3, 1, 0))

plot(as.numeric(r$Size)[w], pseudo_frame[w], pch=46, col=colors["frame"], xlab="Size (mbp)", ylab="Pseudogenes", log="x", ylim=range(pseudo_frame, pseudo_stop, na.rm=TRUE))
points(as.numeric(r$Size)[w], pseudo_stop[w], pch=46, col=colors["stop"])
for (b in c(100))
	abline(a=0, b=b, untf=TRUE, col=colors[2]) # b per mbp
legend("topleft", c("frameshift(s)", "internal stop(s)"), text.col=substring(colors[c("frame", "stop")], 1, 7), bty="n")

plot(as.numeric(r$Size)[w]*1000/as.numeric(r$Scaffolds)[w], pseudo_frame[w]/as.numeric(r$Size)[w]/1000, pch=46, col=colors["frame"], xlab="Average kbp per scaffold", ylab="Pseudogenes per kbp", log="x")
points(as.numeric(r$Size)[w]*1000/as.numeric(r$Scaffolds)[w], pseudo_stop[w]/as.numeric(r$Size)[w]/1000, pch=46, col=colors["stop"])
#legend("topleft", c("frameshift(s)", "internal stop(s)"), text.col=substring(colors[c("frame", "stop")], 1, 7), bty="n")

plot(cds_tot[w]/as.numeric(r$Size)[w]/1000, pseudo_frame[w]/as.numeric(r$Size)[w]/1000, pch=46, col=colors["frame"], xlab="Coding sequences per kbp", ylab="Pseudogenes per kbp", xlim=c(0, 1.5), ylim=c(0, 1.5), xaxs="i", yaxs="i")
points(cds_tot[w]/as.numeric(r$Size)[w]/1000, pseudo_stop[w]/as.numeric(r$Size)[w]/1000, pch=46, col=colors["stop"])
for (b in c(1))
	abline(a=0, b=b, untf=TRUE, col=colors[2], lty=2)
abline(a=0.9, b=-1, col=colors[2])
#legend("topleft", c("frameshift(s)", "internal stop(s)"), text.col=substring(colors[c("frame", "stop")], 1, 7), bty="n")

plot(as.numeric(r$GC.)[w], pseudo_frame[w]/as.numeric(r$Size)[w]/1000, pch=46, col=colors["frame"], xlab="GC-content (%)", ylab="Pseudogenes per kbp")
points(as.numeric(r$GC.)[w], pseudo_stop[w]/as.numeric(r$Size)[w]/1000, pch=46, col=colors["stop"])
#legend("topleft", c("frameshift(s)", "internal stop(s)"), text.col=substring(colors[c("frame", "stop")], 1, 7), bty="n")

### Figure 2

dev.new()
layout(matrix(1:4, 2, byrow=TRUE))
p <- par(mar=c(4, 4, 1, 1), mgp=c(2.3, 1, 0))

W <- which(as.numeric(substring(r$Release.Date, 1, 4))[w] >= 2000)
split_violin(X=as.numeric(substring(r$Release.Date, 1, 4))[w][W],
	Y1=(pseudo_frame[w]/as.numeric(r$Size)[w]/1000)[W],
	Y2=(pseudo_stop[w]/as.numeric(r$Size)[w]/1000)[W],
	colors=colors[c("frame", "stop")],
	YMAX=0.1,
	xlab="",
	ylab="Pseudogenes per kbp")

genus <- sapply(strsplit(r$X.Organism.Name[w], " "), `[`, 1L)
genus[genus=="uncultured"] <- NA_character_
genus[genus=="Candidatus"] <- NA_character_
t <- tapply((pseudo_frame + pseudo_stop)[w]/as.numeric(r$Size)[w]/1000, genus, c)
GX <- sapply(t, function(x) sum(!is.na(x)))
GY <- sapply(t, median, na.rm=TRUE)
W <- which(GX > 50) # minimum number of genomes
plot(GX[W], GY[W], xlab="Genomes per genus", ylab="Average pseudogenes per kbp", log="x", ylim=c(0, 0.2), xlim=c(min(GX[W]), max(GX[W]*10)), col=colors[5])
UQ <- sapply(t, quantile, 0.25, na.rm=TRUE)
LQ <- sapply(t, quantile, 0.75, na.rm=TRUE)
arrows(GX[W], UQ[W], GX[W], LQ[W], code=3, length=0.02, angle=90, col=colors[5])
print("Select points by clicking on the plot, and then press Esc, right-click, or click the Stop button to continue.")
I <- identify(GX, GY, names(GX), col=colors[2], font=3)
points(GX[I], GY[I], col=colors[2], pch=16)

genera <- c("Shigella", "Francisella")
for (g in genera) {
	W <- which(genus==g)
	X <- as.numeric(substring(r$Release.Date, 1, 4))[w][W] + as.numeric(substring(r$Release.Date, 6, 7))[w][W]/12 + as.numeric(substring(r$Release.Date, 9, 10))[w][W]/365
	t_frame <- tapply((pseudo_frame[w]/as.numeric(r$Size)[w]/1000)[W], X, range)
	t_stop <- tapply((pseudo_stop[w]/as.numeric(r$Size)[w]/1000)[W], X, range)
	plot(NA,, xlab="Release date", ylab="Pseudogenes per kbp", xlim=c(2015, 2021), ylim=range(unlist(t_frame), unlist(t_stop), na.rm=TRUE))
	segments(as.numeric(names(t_frame)),
		sapply(t_frame, `[`, 1L),
		as.numeric(names(t_frame)),
		sapply(t_frame, `[`, 2L),
		col=paste0(substring(colors["frame"], 1, 7), "66"),
		lwd=2)
	segments(as.numeric(names(t_stop)),
		sapply(t_stop, `[`, 1L),
		as.numeric(names(t_stop)),
		sapply(t_stop, `[`, 2L),
		col=paste0(substring(colors["stop"], 1, 7), "66"),
		lwd=2)
	legend("topleft", g, bty="n", text.font=3)
}

### Figure 3

dev.new()
layout(matrix(1:4, 2, byrow=TRUE))
p <- par(mar=c(5, 4, 4, 1), mgp=c(2.3, 1, 0))

level <- gsub(" ", "", r$Level[w])
level <- paste0(substring(level, 1, 5),
	ifelse(nchar(level) > 6,
		".",
		substring(level, 6, 6)))
split_violin(X=as.factor(level),
	Y1=pseudo_frame[w]/as.numeric(r$Size)[w]/1000,
	Y2=pseudo_stop[w]/as.numeric(r$Size)[w]/1000,
	colors=colors[c("frame", "stop")],
	YMAX=0.2,
	xlab="Release year",
	ylab="Pseudogenes per kbp",
	las=1)

coverage <- sapply(stats, grep, pattern="# Genome coverage: ", value=T)
coverage <- sapply(coverage, `[`, 1L)
coverage <- sapply(coverage, gsub, pattern="^# Genome coverage: (.*?)x$", replacement="\\1")
coverage <- as.numeric(coverage)
plot(as.numeric(coverage)[w], pseudo_frame[w]/as.numeric(r$Size)[w]/1000, pch=46, col=colors["frame"], xlab="Reported coverage", ylab="Pseudogenes per kbp", xlim=c(10, 10000), log="x")
points(as.numeric(coverage)[w], pseudo_stop[w]/as.numeric(r$Size)[w]/1000, pch=46, col=colors["stop"])

p <- par(mar=c(8, 4, 1, 1), mgp=c(2.3, 1, 0))
technology <- sapply(stats, grep, pattern="# Sequencing technology: ", value=TRUE, fixed=TRUE)
technology <- substring(sapply(technology, `[`, 1L), 25)
m <- matrix(c(grepl("Illumina", technology, ignore.case=TRUE),
		grepl("PacBio", technology, ignore.case=TRUE),
		grepl("Ion", technology, ignore.case=TRUE),
		grepl("Oxford", technology, ignore.case=TRUE),
		grepl("454", technology, ignore.case=TRUE),
		grepl("HiSeq", technology, ignore.case=TRUE),
		grepl("MiSeq", technology, ignore.case=TRUE),
		grepl("NextSeq", technology, ignore.case=TRUE)),
		ncol=8)
tech <- ifelse(m[, 1] & m[, 2] & rowSums(m[, 2:5])==1, "Illumina\n+ PacBio",
	ifelse(m[, 1]  & m[, 3] & rowSums(m[, 2:5])==1, "Illumina\n+ Nanopore",
	ifelse(m[, 1] + rowSums(m[, 2:5]) > 1, "Illumina\n+ other(s)",
	ifelse(m[, 1] & rowSums(m[, 1:5])==1, "Illumina",
	ifelse(m[, 2] & rowSums(m[, 1:5])==1, "PacBio",
	ifelse(m[, 3] & rowSums(m[, 1:5])==1, "Ion Torrent",
	ifelse(m[, 4] & rowSums(m[, 1:5])==1, "Oxford\nNanopore",
	ifelse(m[, 5] & rowSums(m[, 1:5])==1, "Roche 454",
	NA_character_))))))))
#tech <- ifelse(m[, 6], "HiSeq",
#	ifelse(m[, 7], "MiSeq",
#	ifelse(m[, 8], "NextSeq",
#	NA_character_)))
split_violin(X=as.factor(tech[w]),
	Y1=pseudo_frame[w]/as.numeric(r$Size)[w]/1000,
	Y2=pseudo_stop[w]/as.numeric(r$Size)[w]/1000,
	colors=colors[c("frame", "stop")],
	YMAX=0.2,
	xlab="",
	ylab="Pseudogenes per kbp")

assembler <- sapply(stats, grep, pattern="# Assembly method: ", value=TRUE, fixed=TRUE)
assembler <- sapply(assembler, `[`, 1L)
assembler <- sapply(assembler, gsub, pattern="^# Assembly method: (.*?) v\\..*$", replacement="\\1")
substring(assembler, 1, 1) <- toupper(substring(assembler, 1, 1))
assembler[assembler=="IDBA_UD"] <- "IDBA-UD"
assembler[assembler=="Idba_ud"] <- "IDBA-UD"
assembler[assembler=="HGAP2"] <- "HGAP"
assembler[assembler=="SOAP denovo"] <- "SOAPdenovo"
assembler[assembler=="Unknown program"] <- NA_character_
assembler[grepl(" and ", assembler, fixed=TRUE)] <- NA_character_
assembler[grepl(";", assembler, fixed=TRUE)] <- NA_character_
assembler[grepl(",", assembler, fixed=TRUE)] <- NA_character_
assembler[grepl("^CLC", assembler)] <- "CLC Workbench"
assembler[grepl("CLC Genomics Workbench", assembler)] <- "CLC Workbench"
assembler[grepl("SOAP", assembler)] <- "SOAPdenovo"
assemblers <- c("ALLPATHS", "HGAP", "MetaSPAdes", "IDBA-UD", "MEGAHIT", "Newbler", "Unicycler", "Velvet", "ABySS", "Canu", "Celera", "CLC Workbench", "SOAPdenovo", "SMRT", "A5-miseq", "MaSuRCA", "MIRA", "Platanus", "SKESA", "SPAdes")
for (a in assemblers)
	assembler[grepl(a, assembler, ignore.case=TRUE)] <- a
t <- table(assembler)
W <- which(assembler[w] %in% names(t[t > 1000]) & # at least 1000 genomes
	!grepl("Unknown", assembler[w]) &
	!grepl("in-house", assembler[w]) &
	m[w, 1] & # Illumina only
	genus[w] %in% names(GX)[GY < 0.05 & GX > 500])
split_violin(X=as.factor(assembler[w][W]),
	Y1=(pseudo_frame[w]/as.numeric(r$Size)[w]/1000)[W],
	Y2=(pseudo_stop[w]/as.numeric(r$Size)[w]/1000)[W],
	colors=colors[c("frame", "stop")],
	YMAX=0.2,
	xlab="",
	ylab="Pseudogenes per kbp")

### Figure 4

tf <- tempfile()
genera <- c("Shigella", "Francisella")
categories <- c("transposase", "transporter", "hypothetical protein", "other") # last category must be "other"

classify <- function(x) {
	if (length(x) == 0) {
		rep(NA, 4)
	} else {
		y <- rep(4, length(x))
		for (i in seq_len(length(categories) - 1))
			y[grepl(categories[i], x)] <- i
		tabulate(y, 4)
	}
}

dev.new()
layout(matrix(seq_len(length(genera)*2), ncol=2))

for (genus in genera) {
	print(genus)
	pBar <- txtProgressBar(style=3)
	
	ftps <- r$GenBank.FTP[w][grep(paste0("^", genus), r$X.Organism.Name[w])]
	ftps <- paste0(ftps,
		"/",
		sapply(strsplit(ftps, "/", fixed=TRUE), tail, n=1),
		"_genomic.gbff.gz")
	results <- vector("list", length(ftps))
	for (k in seq_along(ftps)) {
		t <- try(download.file(ftps[k], tf, quiet=TRUE),
			silent=TRUE)
		if (is(t, "try-error"))
			next
		t <- try(gbk <- readLines(tf, warn=FALSE),
			silent=TRUE)
		unlink(tf)
		if (is(t, "try-error"))
			next
		
		W <- which(substring(gbk, 1, nchar("     CDS")) == "     CDS")
		func_frame <- func_stop <- character(length(W))
		for (i in seq_along(W)) {
			g1 <- W[i]
			for (g2 in (g1 + 1L):length(gbk))
				if (substring(gbk[g2], 1, nchar("                     ")) != "                     ")
					break
			g <- grep('/pseudo', gbk[g1:g2], fixed=TRUE)[1]
			if (is.na(g))
				next
			string <- paste(substring(gbk[g1:g2], nchar("                     ")), collapse="")
			product <- gsub('.*/product="(.*?)".*', "\\1", string)
			if (product == string) # no product tag
				next
			g <- grep('/note=', gbk[g1:g2], fixed=TRUE)[1]
			if (is.na(g))
				next
			if (grepl('"frameshifted;', gbk[g1:g2][g], fixed=TRUE)) {
				func_frame[i] <- product
			} else if (grepl('"internal stop;', gbk[g1:g2][g], fixed=TRUE)) {
				func_stop[i] <- product
			}
		}
		results[[k]] <- list(frame=func_frame[func_frame != ""],
			stop=func_stop[func_stop != ""])
		
		setTxtProgressBar(pBar, k/length(ftps))
	}
	
	func_frame <- sapply(lapply(results, `[[`, "frame"), classify)
	func_stop <- sapply(lapply(results, `[[`, "stop"), classify)
	W <- which(!is.na(colSums(func_frame)) & !is.na(colSums(func_stop)))
	o <- order(as.numeric(substring(r$Release.Date, 1, 4))[w[grep(paste0("^", genus), r$X.Organism.Name[w])]][W] + as.numeric(substring(r$Release.Date, 6, 7))[w[grep(paste0("^", genus), r$X.Organism.Name[w])]][W]/12 + as.numeric(substring(r$Release.Date, 9, 10))[w[grep(paste0("^", genus), r$X.Organism.Name[w])]][W]/365)
	
	p <- par(mar=c(0, 4, 1, 0), mgp=c(2.3, 1, 0))
	barplot(func_frame[, W[o]], border=NA, col=colors[categories], ylab="Frameshifts", ylim=c(0, max(colSums(func_frame[, W]))), yaxs="i", cex.axis=1.5, cex.lab=1.5)
	p <- par(mar=c(1, 4, 0, 0), mgp=c(2.3, 1, 0))
	barplot(func_stop[, W[o]], border=NA, col=colors[categories], ylab="Internal stops", ylim=c(max(colSums(func_stop[, W])), 0), yaxs="i", cex.axis=1.5, cex.lab=1.5)
	print(length(W)) # number of annotated genomes
	legend("bottomleft", rev(categories), fill=rev(colors[categories]), bty="n", col=NA, cex=1.5)
}
