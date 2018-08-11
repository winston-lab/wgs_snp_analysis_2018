p1=read.table("D:/WINSTON LAB/Data/Whole genome sequencing/SNP files/RGC95.SNP", stringsAsFactors = F, skip=53)
p2=read.table("D:/WINSTON LAB/Data/Whole genome sequencing/SNP files/RGC98.SNP", stringsAsFactors = F, skip=53)
bedfile=read.csv("D:/Rstuff/CodingSeqAnnotation_allGenes.csv", stringsAsFactors = F)
sgd_ids=read.table("D:/WINSTON LAB/R scripts/primer_design_app/sgd_ids.tsv",stringsAsFactors = F, sep="\t", head=F)
all_combined=NULL
all_combined_2=NULL
#Read in files
sample_list=c("S1","S3","S7","S9","LF5a","LF6a","LF9b","HF5b")
for(j in sample_list)
{
  sample=read.table(paste("D:/WINSTON LAB/Data/Whole genome sequencing/SNP files/",j,".SNP",sep=""), stringsAsFactors = F, skip=53)
  
  matches_p1=NULL
  matches_p2=NULL
  
  #Find out which snps are common between the sample and parents
  for(i in 1:nrow(sample))
  {
    chr=which(p1[,1]==sample[i,1])
    
    if(!is.na(match(sample[i,2],p1[chr,2])))
    {
      matches_p1=c(matches_p1,i)
    }
    chr=which(p2[,1]==sample[i,1])
    
    if(!is.na(match(sample[i,2],p2[chr,2])))
    {
      matches_p2=c(matches_p2,i)
    }
  }
  
  #Remove snps which are common between the samples and parents
  union_matches=union(matches_p1,matches_p2)
  snps=setdiff(1:nrow(sample),union_matches)
  
  diff_snps=sample[snps,]
  
  cds_candidates= NULL
  inter_candidates=NULL
  
  #Find SNPs that are in the coding regions of genes (cds_candidates)
  #Find SNPs that are in the intergenic regions of genes (inter_candidates). In this case, they are "mapped" to the
  #nearest gene.
  for(i in 1:nrow(diff_snps))
  {
    chr=which(bedfile[,1]==diff_snps[i,1])
    var1 = head(which(diff_snps[i,2]<bedfile[chr,2]),1)
    var2 = tail(which(diff_snps[i,2]>bedfile[chr,3]),1)
    
    if(length(var1)>0 && length(var2)>0)
    {
      if(var1-var2==2)
      {
       cds_candidates=c(cds_candidates,bedfile[chr,][mean(c(var1,var2)),6])
      } else
      {
        inter_candidates=c(inter_candidates,bedfile[chr,][var1,6])
      }
      
    }
  }
  assign(paste("sys_names_",j,sep=""),cds_candidates)
  assign(paste("inter_sys_names_",j,sep=""),inter_candidates)
  assign(paste("common_names_",j,sep=""),sgd_ids[match(cds_candidates,sgd_ids[,2]),3])
  all_combined=c(all_combined,cds_candidates)
  all_combined_2=c(all_combined_2,inter_candidates)
}

#Display which genes contain the maximum number of hits
sort(table(all_combined))
sort(table(all_combined_2))
