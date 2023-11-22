
plotGOFromGeneOntology <- function(file,maxNumberOfGOterms='all',plotType='bar',sheet=1,
                                   sortby='FDR',returnData=TRUE,returnDataType='plot',
                                   axisTextSize=6,axisTitleSize=8,showname='term and id'){
  if(endsWith(file,'.txt')){
    fi <- read.table(file = file,header = TRUE,sep = '\t')
  }else if(endsWith(file,'.xlsx')){
    library(openxlsx)
    fi <- read.xlsx(xlsxFile = file,sheet = sheet)
    
  }else{
    message(paste('file format of file',file,'is unkown!'))
    q()
  }
  colnames(fi) <- c('GOterms','RefList','GeneList','expected','over/under','FoldEnrichment','rawPvalue','FDR')
  

  #fi <- subset(fi,fi$'over/under' =='+')
  fi <- subset(fi,fi$`over/under`=='+' & fi$GOterms != 'Unclassified (UNCLASSIFIED)')
  if(showname=='term'){
    a <- strsplit(fi$GOterms,' \\(')
    fi$GOterms <-  matrix(unlist(a),ncol = 2,byrow = T)[,1]
  }
  ### sort the dataframe (may check whether necessary??)
  if(! sortby %in% c('FDR','FoldEnrichment')){
    message('the value of "sortby" is incorrect! Accepted values include "FDR" and "FoldEnrichment".')
    q()
  }
  if(sortby=='FDR'){
    fi <- fi[order(fi$FDR),]
  }else if(sortby=='FoldEnrichment'){
    fi <- fi[order(fi$FoldEnrichment,decreasing = TRUE),]
  }
  ##
  #if(! plotALL){
  if(maxNumberOfGOterms !='all'){
    fi <- fi[1:maxNumberOfGOterms,]
  }
  fi$GOterms <- factor(fi$GOterms,levels = rev(fi$GOterms))
  library(ggplot2)
  if(plotType=='bar'){
    p <- ggplot(fi,aes(x=GOterms,y=-log10(FDR))) +
      geom_bar(stat = 'identity',position = 'dodge') +
      labs(x=NULL,y=expression(-log[10]~FDR)) +
      coord_flip() +
      scale_y_continuous(position = 'right') +
      theme(axis.text = element_text(color = 'black',size = axisTextSize),axis.title = element_text(size = axisTitleSize))
  }else if(plotType=='bubble'){
    #need to re-order in future
    if(sortby=='FDR'){
      p <- ggplot(fi,aes(y=GOterms,x=-log10(FDR))) +
        geom_point(aes(color=-log10(FDR),size=GeneList)) +
        labs(x=expression(-log[10]~FDR),y=NULL) +
        theme(axis.text = element_text(color = 'black',size = 10)) +
        scale_color_continuous() +
        guides(color=guide_legend(title = expression(-log[10]~FDR)),size=guide_legend(title = 'gene number'))
    }else if(sortby=='FoldEnrichment'){
      p <- ggplot(fi,aes(y=GOterms,x=FoldEnrichment)) +
        geom_point(aes(color=-log10(FDR),size=GeneList)) +
        labs(x='Fold Enrichment',y=NULL) +
        theme(axis.text = element_text(color = 'black',size = 10)) +
        scale_color_continuous() +
        guides(color=guide_legend(title = expression(-log[10]~FDR)),size=guide_legend(title = 'gene number'))
    }
    #
    
  }else{
    message('Unsupported plotType!')
    q()
  }
  if(!returnData){
    print(p)
  }else{
    if(returnDataType=='dataframe'){
      return(fi)
    }else if(returnDataType=='plot'){
      return(p)
    }else{
      message('Error! Please check whether the data type to be returned is specified correctly! Accepted value is "dataframe" or "plot"')
      q()
    }
  }
  
  
}