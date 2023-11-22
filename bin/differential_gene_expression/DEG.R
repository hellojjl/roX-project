#!/usr/bin/env Rscript
fpkmbyDESeq2 <- function(gtf_file,count_directory,count_matrix,outfilename,robust=TRUE){
  library('rtracklayer')
  gtf <- import.gff2(gtf_file)
  gtf1 <- gtf[gtf$type == 'exon']
  gtf2 <- split(gtf1,gtf1$gene_id)

  df <- read.table(count_matrix,sep='\t')
  colnames(df) <- c('name','files','cond')

  df$cond <- factor(df$cond,levels=unique(df$cond))
  library('DESeq2')
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = df,directory = count_directory,design= ~cond)
  dds <- DESeq(ddsHTSeq)
  res <- results(dds)

  rowRanges(dds) <- gtf2
  res2 <- fpkm(dds,robust=robust)
  res3 <- data.frame(res2,res)
  write.table(cbind('gene'=rownames(res3),res3),file=outfilename,quote=F,sep='\t',row.names = F)
}

args <- commandArgs(T)

if (length(args)!=5){
  message('Usage: fpkm.R <Transcriptome gtf file> <directory> <count_info_matrix> <output file name> <robust value>')
  cat('count_info_matrix should contain 3 columns: samplename, filename,groupnumber\n')
  q()
}
if(args[5]=='FALSE'){
  fpkmbyDESeq2(gtf_file=args[1],count_directory=args[2],count_matrix=args[3],outfilename=args[4],robust=FALSE)
}else{
fpkmbyDESeq2(gtf_file=args[1],count_directory=args[2],count_matrix=args[3],outfilename=args[4])
}
message('finished at ',Sys.time())
