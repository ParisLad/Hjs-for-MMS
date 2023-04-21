library(BSgenome.Hsapiens.UCSC.hg19.masked)
genome <- BSgenome.Hsapiens.UCSC.hg19.masked
seqnames <- seqnames(genome)
for (seqname in seqnames)
{
	x <- genome[[seqname]]
	#x <- genome$chr1  # Kostas comment: insert chromosome here
	#myfile <- paste("C:\\Work\\GeneScreening\\hg_hj.bed", sep="") # Kostas comment: insert output filename here
	active(masks(x)) <- TRUE
	active(masks(x))["RM"] <- FALSE
	hjtable <- read.table("C:\\Work\\GeneScreening\\hj.txt", header = TRUE) # Kostas comment: insert hj.txt path here
	for (i in 1:nrow(hjtable)) # Kostas comment insert for start end for hj
	{
		#print(paste("Run HJ no. ", seqname, i, hjtable$HJ[i]))
		Lpattern <- toString(hjtable$HJ[i])
		M1 <- matchPattern(Lpattern, x, max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE, algorithm="auto")
		#print(paste("M1 = ", M1))
		myfile <- paste("C:\\Work\\GeneScreening\\hj", hjtable$HJ[i], ".bed", sep="")
		sink(myfile, append=TRUE)
		if (length(M1) != 0)
		{
			chrvalue <- seqname
			for (k in 1:length(M1))
			{
				cat(chrvalue)
				cat("\t")
				cat(start(M1[k]))
				cat("\t")
				cat(end(M1[k]))
				cat("\t")
				cat(toString(hjtable$HJ[i]))
				cat("\n")
			}
		}
		sink()
	}

}
