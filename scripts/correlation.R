#!/bin/R

#this script combines all output files to create the correlation files
#this is outputted to the outfolder/correlation directory

library("ggplot2")


args <- commandArgs(trailingOnly = TRUE)
outfolder<-args[1]

#if folder doesnt exist create it
if (!dir.exists(paste0(outfolder,"/correlation"))){
  dir.create(paste0(outfolder,"/correlation"))
}

#create a correlation file for each enst_file x matrix_file (histone mark type)
for(enst_file in list.files(paste0(outfolder,"/quants/enst/"))){
  for(matrix_file in list.files(paste0(outfolder,"/matrix/"))){
    #make enst df and assign proper row/col names
    enst<-read.csv(paste0(outfolder,"/quants/enst/",enst_file), sep="\t", header=FALSE)
    row.names(enst)<-enst[,1]
    enst<-enst[-1]
    colnames(enst)<-"TPM"
    
    print(paste0(matrix_file," ",enst_file))

    #create proper matrix where if there are duplicate col names (cell types) they are averaged
    #and then saved to a singular column in matrix_uniq which is ultimately assigned proper names
    matrix<-read.csv(paste0(outfolder,"/matrix/",matrix_file), sep="\t",check.names = FALSE)
    colnames(matrix)<-gsub("\\.","-",colnames(matrix))
    duplicates<-colnames(matrix)[duplicated(colnames(matrix))]
    matrix_uniq<-matrix[,unique(colnames(matrix))]
    for(duplicate in duplicates){
      matrix_duplicate<-matrix[,colnames(matrix)==duplicate]
      duplicate_average<-rowSums(matrix_duplicate)/ncol(matrix_duplicate)
      matrix_uniq[,duplicate]<-duplicate_average
    }
    rownames(matrix_uniq)<-matrix_uniq[,1]
    matrix_uniq<-matrix_uniq[,-1]

    #make a reduced enst and matrix df where they share same cell types
    enst_red<-data.frame(enst[intersect(colnames(matrix_uniq),rownames(enst)),])
    rownames(enst_red)<-intersect(colnames(matrix_uniq),rownames(enst))
    colnames(enst_red)<-"TPM"
    matrix_red<-data.frame(matrix_uniq[,intersect(colnames(matrix_uniq),rownames(enst))])
    colnames(matrix_red)<-intersect(colnames(matrix_uniq),rownames(enst))
      
    #perform pearson/spearman correlation on this
    #there WILL be NAs. this is totally okay, just means all vectors (all cell types at this location) are 0
    #which is biologically valid but will produce warning
    pearson_corr<-as.numeric(as.character(array(cor(t(matrix_red),enst_red, use="complete.obs", method="pearson"))))
    spearman_corr<-array(cor(t(matrix_red), enst_red,use="complete.obs", method="spearman"))
    
    #create dataframe that will be printed where the cols are bins/pearson/spearman
    #and write it to outfile
    corr_matrix<-cbind("pearson"=pearson_corr,"spearman"=spearman_corr)
    rownames(corr_matrix)<-rownames(matrix_red)
    start_bin<-as.numeric(sapply(strsplit(rownames(corr_matrix),"_"),"[[",2))
    plot_matrix<-data.frame(cbind("loc"=start_bin,corr_matrix))
    plot_matrix$pearson<-as.numeric(as.character(plot_matrix$pearson))
    mark<-sapply(strsplit(sapply(strsplit(matrix_file,"_"),"[[",3),"\\."),"[[",1)
    outfile<-paste0(outfolder,"/correlation/",mark,"-",enst_file)
    print(outfile)
    write.table(plot_matrix,file=outfile,sep="\t",row.names=FALSE)

    #create a picture for this df (pearson)
    outpic<-paste0(strsplit(outfile,"\\.txt")[[1]][1],".png")
    gg<-ggplot(plot_matrix,aes(x=loc,y=spearman))+ggtitle(paste0(mark," ",enst_file))+geom_point()
    ggsave(gg,file=outpic)

  }
}


