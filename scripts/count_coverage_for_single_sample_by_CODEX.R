# docker run --rm -it -v /data/work/home/iwkusmirek/cnv_cov:/data/work/home/iwkusmirek/cnv_cov -w /data/work/home/iwkusmirek/cnv_cov biodatageeks/cnv-opt-codexcov Rscript count_coverage_for_single_sample_by_CODEX.R NA06984.chrom11.ILLUMINA.bwa.CEU.exome.20120522.bam codex.bed 20 11 NA06984_codex.txt

count_coverage_for_single_sample_by_CODEX <- function(bambedObj, mapqthres) {
    ref <- bambedObj$ref
    chr <- bambedObj$chr
    bamdir <- bambedObj$bamdir
    st <- start(ref)[1]
    ed <- end(ref)[length(ref)]
    Y <- matrix(NA, nrow = length(ref), ncol = 1)
    readlength <- rep(NA, 1)
    i <- 1
    bamurl <- bamdir[i]
    which <- RangesList(quack = IRanges(st - 10000, ed + 10000))
    names(which) <- as.character(chr)
    what <- c("pos", "mapq", "qwidth")
    flag <- scanBamFlag(isDuplicate = FALSE, isUnmappedQuery = FALSE,
                        isNotPassingQualityControls = FALSE, 
                        isFirstMateRead = TRUE)
    param <- ScanBamParam(which = which, what = what, flag = flag)
    bam <- scanBam(bamurl, param = param)[[1]]
    mapqfilter <- (bam[["mapq"]] >= mapqthres)
    readlength[i] <- round(mean(bam[["qwidth"]]))
    if(is.nan(readlength[i])){
      flag <- scanBamFlag(isDuplicate = FALSE, isUnmappedQuery = FALSE,
                          isNotPassingQualityControls = FALSE)
      param <- ScanBamParam(which = which, what = what, flag = flag)
      bam <- scanBam(bamurl, param = param)[[1]]
      mapqfilter <- (bam[["mapq"]] >= mapqthres)
      readlength[i] <- round(mean(bam[["qwidth"]]))
    }
    message("Getting coverage for sample ", bamurl, ": ", 
            "read length ", readlength[i], ".", sep = "")
    irang <- IRanges(bam[["pos"]][mapqfilter], width = 
                    bam[["qwidth"]][mapqfilter])
    Y[, i] <- countOverlaps(ref, irang)
    list(Y = Y, readlength = readlength)
}

library(CODEX)
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("Invalid number of arguments!!!", call.=FALSE)
}
bamFile <- args[1]
bedFile <- args[2]
mapqthres <- strtoi(args[3])
chr <- args[4]
outputFile <- args[5]

bambedObj <- getbambed(bamdir = bamFile, bedFile = bedFile, sampname = NULL, projectname = NULL, chr)
coverageObj <- count_coverage_for_single_sample_by_CODEX(bambedObj, mapqthres = mapqthres)
write.table(data.frame(coverageObj$Y), sep=",", col.names=F, row.names=F, quote=F,  file = outputFile)
