args <- commandArgs(trailingOnly = TRUE)
outfolder<-args[1]

if (!dir.exists(paste0(outfolder,"/correlation"))){
  dir.create(paste0(outfolder,"/correlation"))
}

for(enst_file in list.files(paste0(outfolder,"/quants/enst/"))){
  for(matrix_file in list.files(paste0(outfolder,"/matrix/"))){
    enst<-read.csv(paste0(outfolder,"/quants/enst/",enst_file), sep="\t", header=FALSE)
    row.names(enst)<-enst[,1]
    enst<-enst[-1]
    colnames(enst)<-"TPM"
    
    print(paste0(matrix_file," ",enst_file))
    matrix<-read.csv(paste0(outfolder,"/matrix/",matrix_file), sep="\t",check.names = FALSE)
    colnames(matrix)<-gsub("\\.","-",colnames(matrix))
    duplicates<-colnames(matrix)[duplicated(colnames(matrix))]
    matrix_red<-matrix[,unique(colnames(matrix))]
    for(duplicate in duplicates){
      matrix_duplicate<-matrix[,colnames(matrix)==duplicate]
      duplicate_average<-rowSums(matrix_duplicate)/nrow(matrix_duplicate)
      matrix_red[,duplicate]<-duplicate_average
    }
    
    rownames(matrix_red)<-matrix_red[,1]
    matrix_red<-matrix_red[,-1]

    enst_red<-data.frame(enst[intersect(colnames(matrix),rownames(enst)),])
    rownames(enst_red)<-intersect(colnames(matrix),rownames(enst))
    colnames(enst_red)<-"TPM"
      
    pearson_corr<-as.numeric(as.character(array(cor(t(matrix_red),enst_red, use="complete.obs", method="pearson"))))
    spearman_corr<-array(cor(t(matrix_red), enst_red,use="complete.obs", method="spearman"))
    
    corr_matrix<-cbind("pearson"=pearson_corr,"spearman"=spearman_corr)
    rownames(corr_matrix)<-rownames(matrix_red)
    start_bin<-as.numeric(sapply(strsplit(rownames(corr_matrix),"_"),"[[",2))
      #spearman_corr<-cor(data.frame(cbind(enst_red,t(matrix_red[1,]))), use="complete.obs", method="spearman") 
    plot_matrix<-data.frame(cbind("loc"=start_bin,corr_matrix))
    plot_matrix$pearson<-as.numeric(as.character(plot_matrix$pearson))
    mark<-sapply(strsplit(sapply(strsplit(matrix_file,"_"),"[[",3),"\\."),"[[",1)
    outfile<-paste0(outfolder,"/correlation/",mark,"-",enst_file)
    write.table(plot_matrix,outfile,sep="\t",row.names=FALSE)
    
    #gg<-ggplot(plot_matrix,aes(x=loc,y=pearson))+ggtitle(paste0(mark," ",enst_file))+geom_point()
    #print(gg)
  }
}

