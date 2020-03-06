rm(list=ls())
runsqtl<-function(tr.f,snp.f,bed.f,out.f){
  library(sQTLseekeR)
  #index genotype file
  message("[info] sQTLseekeR")
  message("    index genotype file")
  #print(snp.f)
  genotype.indexed.f = index.genotype(snp.f)

  #read transcript count data
  message("    read transcript count data")
  te.df = read.table(tr.f, as.is=TRUE, header=TRUE, sep="\t")

  #prepare transcript expression matrix
  #remove low expression gene and transcript
  #remove single transcript in a gene
  #remove low dispersion genes
  message("    prepare transcript expression data")
  tre.df = prepare.trans.exp(te.df,min.dispersion = 0.01)

  #read and process genebed information
  gene.bed = read.table(bed.f, as.is=TRUE, sep="\t")
  colnames(gene.bed) = c("chr","start","end","geneId")

  #run the main function ,analysis sqtls(snp-gene pairs)
  message("    run sQTLseekeR")
  res.df = sqtl.seeker(tre.df, genotype.indexed.f, gene.bed, svQTL=F,verbose = T)

  #compute q-value (i.e. corrected pvalue)
  sqtls.df = sqtls(res.df, FDR=1, svQTL.removal = F)
  
  #write sqtls result with qvalue to file
  message('    write sQTLseekeR result')
  write.table(sqtls.df,file = out.f,quote = F,row.names = F,sep="\t")

  #return(sqtls.df)
}

args<-commandArgs(T)
runsqtl(args[1],args[2],args[3],args[4])