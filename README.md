library(BSgenome.Hsapiens.UCSC.hg19.masked)
genome <- BSgenome.Hsapiens.UCSC.hg19.masked
seqnames <- seqnames(genome)
hjtable <- read.table("C:\\Work\\GeneScreening\\hj.txt", header = TRUE)
# Kostas comment: insert hj.txt path here

for (i in 1:nrow(hjtable)) # Kostas comment insert for start end for hj
{

        for (seqname in seqnames)
        {
                x <- genome[[seqname]]
                active(masks(x)) <- TRUE
                active(masks(x))["RM"] <- FALSE

                #print(paste("Run HJ no. ", seqname, i, hjtable$HJ[i]))
                Lpattern <- toString(hjtable$HJ[i])
                M1 <- matchPattern(Lpattern, x, max.mismatch=0, min.mismatch=0,
with.indels=FALSE, fixed=TRUE, algorithm="auto")
                #print(paste("M1 = ", M1))
                myfile <- paste("C:\\Work\\GeneScreening\\hj\\", hjtable$HJ[i],
".bed", sep="")
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


        library("regioneR")
        test <- import(myfile)
        test
        Microdel.hg19_path <- system.file("tests", package = "rtracklayer")
        Microdel.hg19_bed <- file.path(Microdel.hg19_path, "Microdel.hg19.bed")
        Microdel.hg19 <- import(Microdel.hg19_bed)
        Microdel.hg19
        pt <- overlapPermTest(A=Microdel.hg19, B=test, ntimes=1000)
        pt
        lz <- localZScore(pt=pt, A=Microdel.hg19, B=test)
        plot(lz)
        plot(pt)
        output <- overlapRegions(test, Microdel.hg19)
        output




}
