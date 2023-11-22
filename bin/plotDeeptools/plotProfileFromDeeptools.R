#to be improved in next version: 1) sample names contain "()".

readCommentInfoFromDeeptoolsFile <- function(infile){
 
  comm_info <- readLines(infile,n = 2)
  #line1
  test <- unlist(strsplit(sub('#','',comm_info[1]),'\t'))
  info1 <- matrix(unlist(strsplit(test,':')),ncol = 2,byrow = T)
  peakinfo <- as.data.frame(info1,col.names = c('peaktype','peaknum')) #col.names not work here
  colnames(peakinfo) <- c('peaktype','peaknum')
  #line2
  test <- unlist(strsplit(sub('#','',comm_info[2]),'\t'))
  info2 <- matrix(unlist(strsplit(test,':')),ncol = 2,byrow = T)
  lengthinfo <- as.data.frame(info2)
  colnames(lengthinfo) <- c('param','length')
  #num_info <- as.data.frame(info2[,2],row.names = info2[,1])
  lengthinfo$length <- as.numeric(lengthinfo$length)
  return(list(peakinfo=peakinfo,lengthinfo=lengthinfo))
}


plotProfile <- function(infile,
                        ignore.peaktypes=NA,regionlabel=NULL,replaceSymbol=NA,
                        start_label='TSS',end_label='TTS',middle_label='center',left_label=NA,right_label=NA,nrow=NULL,
                        returnDataType='dataframe',axis.textsize=9,axis.titlesize=11,error=NA,legend.position='bottom',changeFacet=FALSE,
                        xtitle=NULL,ytitle=NA,xyscales=NULL,samplename=NULL){
  #cols: region1: sample1-bin1, ..., sample1-binN; sample2-bin1, ..., sample2-binN
  #rows: bedtype1 (autosomes) * region_num1; bedtype2 (chrX) * region_num2
  #infile <- 'H3K27me3_matrix1_dm6_chrgenes_scaleregion_chrtypes.tab'
  
  #samplenames=c('Male','roXKO','Female'),peaktypes_num,peaktypes=c('autosomes','chrX','chrY'),up=3000,down=3000,scale=5000,bs=100,
  ##version2:
  #to be improved: error plot;
  
  finfo <- readCommentInfoFromDeeptoolsFile(infile = infile) #list
  if(!is.null(regionlabel)){
    #regionlabel <- c('roX.peaks.autosomes.bed'='autosomes','roX.peaks.chrX.bed'='chrX')
    dic <- data.frame(old=names(regionlabel),new=regionlabel)
    #should check whether the old names in the dic are in the bed regions
    if(! all( names(regionlabel) %in% finfo$peakinfo$peaktype ) ){
      message('Please check whether the names for region label is in the datasets.')
      q()
    }
    for(i in seq(nrow(finfo$peakinfo))){
      finfo$peakinfo$peaktype[i] <- ifelse(finfo$peakinfo$peaktype[i] %in% dic$old,dic$new[which(dic$old==finfo$peakinfo$peaktype)],finfo$peakinfo$peaktype[i])
      
    }
  }
  
  fi <- read.delim(infile,sep = '\t',comment.char = '#') #v2
 
  #change the colnames and adjust the dataframe
  colnames(fi) <- colnames(fi)[-seq(nrow(finfo$peakinfo))] #now ncol is shorter
  df <- fi[,-seq(ncol(fi)+1-nrow(finfo$peakinfo),ncol(fi))]
  #extract samplename and determine bin numbers:
  #samplename <- unique(gsub("\\..*","",colnames(df))) #v2: https://statisticsglobe.com/r-remove-characters-before-or-after-point-in-string
  #need to improve as when the name contains "(",")", they will be transformed into "." as well.
  if(is.null(samplename)){
    a <- gsub("(.*)\\..*","\\1",colnames(df))#v3: ref:https://www.it1352.com/2372834.html
    if(is.na(replaceSymbol)){
      #message('haha')
      samplename <- unique(a)
    }else if(replaceSymbol == '()'){
      #v3
      #message('haha2')
      b <- gsub("(.*)\\.","\\1)",a)
      a <- gsub("(.*)\\.","\\1(",b)
      samplename <- unique(a)  #v3
      samplename <- samplename[-which(grepl(')',samplename) & ! grepl('\\(',samplename))]
      #message('samplename: ',samplename)
      #should consider other separators in future.
    }else{
      message('unknown replaceSymbol!')
      break() #
    }
  }
  
  
  #colnames(df) <- paste0('bin',seq(ncol(df))) #meaningless line
  df <- as.data.frame(lapply(df,as.numeric))
  peaktypes <- finfo$peakinfo$peaktype
  peaknum <- finfo$peakinfo$peaknum
  df$chrtype <- rep(peaktypes,times=peaknum)
  df$peakid <- paste0('peak',1:nrow(df))
  library(reshape2)
  ldf <- melt(df,id.vars = c('chrtype','peakid'),variable.name = 'binid',value.name = 'enrichment')
  down <- finfo$lengthinfo$length[1]
  up <- finfo$lengthinfo$length[2]
  genebody <- finfo$lengthinfo$length[3]
  binsize <- finfo$lengthinfo$length[4]
  message((finfo$lengthinfo$length[1] +finfo$lengthinfo$length[2]+finfo$lengthinfo$length[3])/finfo$lengthinfo$length[4] == (ncol(df)-2)/length(samplename))
  #message(up,genebody,down,binsize,ncol(df),length(samplename))
  if( (finfo$lengthinfo$length[1] +finfo$lengthinfo$length[2]+finfo$lengthinfo$length[3])/finfo$lengthinfo$length[4] == (ncol(df)-2)/length(samplename) ){
    #(up+down+genebody)/binsize
    binnum <- (finfo$lengthinfo$length[1] +finfo$lengthinfo$length[2]+finfo$lengthinfo$length[3])/finfo$lengthinfo$length[4]
  }
 
  ldf$sample <- factor(rep(samplename,each=binnum*nrow(df)),levels = samplename)
  ldf$binregion <- rep(seq(binnum),each=nrow(df))
  agg <- aggregate(enrichment~chrtype+sample+binregion,ldf,FUN = 'mean')
  agg$sd <- aggregate(enrichment~chrtype+sample+binregion,ldf,FUN = 'sd')$enrichment
 
  
  #whether ignore some chroms.need to imporove
  if(! is.na(ignore.peaktypes)){
    agg <- subset(agg,! agg$chrtype %in% ignore.peaktypes)
  }
  
  #agg$binid <- as.numeric(sub('bin','',agg$binid))
  if(returnDataType == 'dataframe'){
    return(agg)
  }else{
    library(ggplot2)
    
    if(!changeFacet){
      p <- ggplot(agg,aes(x=binregion,y=enrichment,color=chrtype)) +
        geom_line() +
        #labs(x=NULL) +
        facet_wrap(~sample) +
        theme(plot.background = element_blank(),panel.background = element_rect(color = 'black',fill = 'transparent'),
              panel.grid = element_blank(),strip.background = element_blank(),axis.text = element_text(color = 'black',size = axis.textsize),
              axis.title = element_text(size = axis.titlesize),
              legend.position = legend.position,legend.background = element_blank(),legend.title = element_blank(),
              legend.box.background = element_blank(),legend.key = element_rect(fill = 'transparent'))
      if(!is.null(nrow) & !is.null(xyscales)){
        p <- p +
          facet_wrap(~sample,nrow = nrow,scales = xyscales)
      }else if(!is.null(nrow) & is.null(xyscales)){
        p <- p +
          facet_wrap(~sample,scales = xyscales)
      }
    }else{
      p <- ggplot(agg,aes(x=binregion,y=enrichment,color=sample)) +
        geom_line() +
        #labs(x=NULL) +
        facet_wrap(~chrtype) +
        theme(plot.background = element_blank(),panel.background = element_rect(color = 'black',fill = 'transparent'),
              panel.grid = element_blank(),strip.background = element_blank(),axis.text = element_text(color = 'black',size = axis.textsize),
              axis.title = element_text(size = axis.titlesize),
              legend.position = legend.position,legend.background = element_blank(),legend.title = element_blank(),
              legend.box.background = element_blank(),legend.key = element_rect(fill = 'transparent'))
      if(!is.null(nrow) & !is.null(xyscales)){
        p <- p +
          facet_wrap(~chrtype,nrow = nrow,scales = xyscales)
      }else if(!is.null(nrow) & is.null(xyscales)){
        p <- p +
          facet_wrap(~chrtype,scales = xyscales)
      }
    }
    
   
    if(is.na(left_label)){
      left_label <- ifelse(up<1000,paste0('-',as.character(up),' bp'),paste0('-',as.character(up/1000),' kb'))
    }
    if(is.na(right_label)){
      right_label <- ifelse(down<1000,paste0(as.character(up),' bp'),paste0(as.character(up/1000),' kb'))
    }
    
    
    if(genebody==0){
      #reference-point mode
      p <- p +
        scale_x_continuous(breaks = c(0,up/binsize,binnum),labels = c(left_label,middle_label,right_label))

    }else{
      #scale-region mode
      p <- p +
        scale_x_continuous(breaks = c(0,up/binsize,(up+genebody)/binsize,binnum),labels = c(left_label,start_label,end_label,right_label))
      
    }
    
   
    if(is.na(ytitle)){
      p <- p +
        labs(x=xtitle)
    }else{
      p <- p +
        labs(x=xtitle,y=ytitle)
    }
    if(!is.na(error)){
      p <- p +
        geom_e
    }
    if(returnDataType == 'plot'){
      return(p)
    }else{
      print(p)
    }
  } 
  
}
