args <- commandArgs(T)
if(length(args)==0) {
  print("Usage: Rscript sangerSeq2.R list.xls")
  quit()
}


library(sangerseqR)
library(sangeranalyseR)
library(seqinr)
library(stringr)


# sangerQC <- function(inputFile, direction, outDir){
sangerQC <- function(sample, inputFile, direction, outDir){
  sangerRead <- SangerRead(inputSource           = "ABIF",
                           readFeature           = direction,
                           readFileName          = inputFile,
                           geneticCode           = GENETIC_CODE,
                           TrimmingMethod        = "M1",
                           M1TrimmingCutoff      = 0.001,
                           M2CutoffQualityScore  = NULL,
                           M2SlidingWindowSize   = NULL,
                           baseNumPerRow         = 100,
                           heightPerRow          = 200,
                           signalRatioCutoff     = 0.33,
                           showTrimmed           = TRUE)
  newSangerRead <- updateQualityParam(sangerRead,
                                      TrimmingMethod       = "M2",
                                      M1TrimmingCutoff     = NULL,
                                      M2CutoffQualityScore = 40,
                                      M2SlidingWindowSize  = 30)
  # writeFasta(newSangerRead,
  #            outputDir         = outDir,
  #            compress          = FALSE,
  #            compression_level = NA)
  name <- inputFile
  start <- newSangerRead@ QualityReport@ trimmedStartPos
  end <- newSangerRead@ QualityReport@ trimmedFinishPos
  rawseqlength<-newSangerRead@ QualityReport@ rawSeqLength
  pa <- pairwiseAlignment(primarySeq(newSangerRead)[start:end], secondarySeq(newSangerRead)[start:end],type = "global-local")
  writePairwiseAlignments(pa,file=paste0(name,".pa"))
  line <- scan(file=paste0(name,".pa"), what = character(0), sep = "\n")
  trimseqlength <- line[14]
  identity <- line[15]
  result<-rbind(rawseqlength,trimseqlength,identity)
  write.table(result, file=paste0(name,"_report"), quote=F)
  seq <- read.fasta(file = gsub("ab1$","seq",inputFile), seqtype = "DNA",
                      as.string = T, set.attributes = F)
  trimmedSeq <- getFrag(seq[[1]], begin = start, end = end)
  write.fasta(trimmedSeq, names = sample, as.string = F, nbchar = 2500,
              file.out = gsub("ab1$","fa",inputFile))
}


fastaRename <- function(ab1,sample){
  fastafile <- read.fasta(file = gsub("ab1$","fa",ab1), seqtype = "DNA",
                          as.string = T, set.attributes = F)
  # write.fasta(sequences = fastafile, names = sample, as.string = T, nbchar = 2500,
  #             file.out = paste0(s,"-QC-",str_match(ab1, "\\[(.+)\\].ab1")[2],".fa"))
  # write.fasta(sequences = fastafile, names = sample, as.string = T, nbchar = 2500,
  #             file.out = paste0(s,"-QC-",strsplit(ab1,"[-.]")[[1]][2],".fa"))
  write.fasta(sequences = fastafile, names = sample, as.string = T, nbchar = 2500,
              file.out = gsub("ab1$","fa",ab1))
  print(fastafile)
}


path <- getwd()
list <- read.delim(args[1], comment.char = "")
for(s in unique(as.character(list[,2]))){
  print(s)
  # f = list.files(pattern=paste0(s,".*F].ab1$"))
  # r = list.files(pattern=paste0(s,".*R].ab1$"))
  # f = list.files(pattern=paste0("^",s,"-.+F.*.ab1$"))
  # r = list.files(pattern=paste0("^",s,"-.+R.*.ab1$"))
  # sangerQC(f, "Forward Read", path)
  # fastaRename(f,s)
  # sangerQC(r, "Reverse Read", path)
  # fastaRename(r,s)
  files = list.files(pattern=paste0("^",s,"-\\S+.ab1$"))
  for(f in files){
    # sangerQC(f, "Forward Read", path)
    sangerQC(s, f, "Forward Read", path)
    # fastaRename(f,s)
  }
}


#生成makebasecalls之前之后的峰图以及report文件
filelist <- list.files(pattern="*_report")
m6 <- vector()
o <- length(filelist)
for(i in 1:o){
  lines <- readLines(filelist[i])
  identity <- lines[4]
  trimSeqLength <- lines[3]
  rawSeqLength <- lines[2]
  m5 <- cbind(filelist[i],identity,trimSeqLength,rawSeqLength)
  m6 <- rbind(m6, m5)
}

write.table(m6, "result.txt", quote = F)

