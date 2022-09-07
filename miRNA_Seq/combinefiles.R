# readNmerge <- function(fn=c("SampleA-maturemiRNA-counts.txt","SampleB-maturemiRNA-counts.txt"), path="/analysis/maturemiRNAcounts/")


readNmerge <- function(fn, path)
  # Don't miss the last slash in path!!
{
  print(fn[2])
  nFiles <- length(fn)
  for(f in 1:nFiles)
  {
    if(f==1){
     f1 <- read.table(fn[f], sep="\t", head=T, fill=T)
    }else{
     f2 <- read.table(fn[f], sep="\t", head=T, fill=T)
     f1 <- merge(f1, f2, by="miRNA", all=T)
    }
  }	
return(f1)
}


# ## BR Mature
# outnew <- readNmerge(fn=c("70459-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/70460-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/70461-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74665-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74666-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74667-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74668-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74669-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74670-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74671-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74672-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74673-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74674-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74675-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74676-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74677-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74678-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74679-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74680-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74681-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74682-maturemiRNA-counts.txt"), path="analysis/maturemiRNAcounts/")
# write.table(outnew, file="Samples-BRMature-miRNAcounts.tsv",sep="\t",row.names=F)

# ## CR Mature
# outnew <- readNmerge(fn=c("70462-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74683-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74684-maturemiRNA-counts.txt"), path="analysis/maturemiRNAcounts/")
# write.table(outnew, file="Samples-CRMature-miRNAcounts.tsv",sep="\t",row.names=F)

# ## BR Genome
# outnew <- readNmerge(fn=c("70459-taggedBAMcounts.txt","analysis/genome-miRNAcounts/70460-taggedBAMcounts.txt","analysis/genome-miRNAcounts/70461-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74665-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74666-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74667-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74668-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74669-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74670-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74671-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74672-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74673-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74674-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74675-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74676-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74677-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74678-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74679-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74680-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74681-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74682-taggedBAMcounts.txt"), path="analysis/genome-miRNAcounts/")
# write.table(outnew, file="Samples-BRGenome-miRNAcounts.tsv",sep="\t",row.names=F)

# ## CR Genome
# outnew <- readNmerge(fn=c("70462-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74683-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74684-taggedBAMcounts.txt"), path="analysis/genome-miRNAcounts/")
# write.table(outnew, file="Samples-CRGenome-miRNAcounts.tsv",sep="\t",row.names=F)


## Genome
outnew <- readNmerge(fn=c("analysis/genome-miRNAcounts/70459-taggedBAMcounts.txt","analysis/genome-miRNAcounts/70460-taggedBAMcounts.txt","analysis/genome-miRNAcounts/70461-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74665-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74666-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74667-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74668-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74669-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74670-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74671-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74672-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74673-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74674-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74675-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74676-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74677-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74678-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74679-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74680-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74681-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74682-taggedBAMcounts.txt", "analysis/genome-miRNAcounts/70462-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74683-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74684-taggedBAMcounts.txt"), path="analysis/genome-miRNAcounts/")
write.table(outnew, file="Samples-FULL_Genome-miRNAcounts.tsv",sep="\t",row.names=F)

## Mature
outnew <- readNmerge(fn=c("analysis/maturemiRNAcounts/70459-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/70460-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/70461-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74665-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74666-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74667-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74668-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74669-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74670-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74671-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74672-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74673-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74674-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74675-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74676-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74677-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74678-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74679-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74680-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74681-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74682-maturemiRNA-counts.txt", "analysis/maturemiRNAcounts/70462-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74683-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74684-maturemiRNA-counts.txt"), path="analysis/maturemiRNAcounts/")
write.table(outnew, file="Samples-FULL_Mature-miRNAcounts.tsv",sep="\t",row.names=F)


## ALL THE RESULTS
## Mature
outnew <- readNmerge(fn=c("analysis/maturemiRNAcounts/70459-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/70460-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/70461-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74665-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74666-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74667-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74668-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74669-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74670-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74671-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74672-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74673-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74674-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74675-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74676-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74677-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74678-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74679-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74680-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74681-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74682-maturemiRNA-counts.txt", "analysis/maturemiRNAcounts/70462-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74683-maturemiRNA-counts.txt","analysis/maturemiRNAcounts/74684-maturemiRNA-counts.txt"), path="analysis/maturemiRNAcounts/")
write.table(outnew, file="Samples-FULL_Mature-miRNAcounts.tsv",sep="\t",row.names=F)

# ## CR Genome
# outnew <- readNmerge(fn=c("70462-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74683-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74684-taggedBAMcounts.txt"), path="analysis/genome-miRNAcounts/")
# write.table(outnew, file="Samples-CRGenome-miRNAcounts.tsv",sep="\t",row.names=F)


# # readNmerge <- function(fn=c("SampleA-maturemiRNA-counts.txt","SampleB-maturemiRNA-counts.txt"), path="/analysis/maturemiRNAcounts/")
# readNmerge <- function(fn=c("70459-taggedBAMcounts.txt","70460-taggedBAMcounts.txt","70461-taggedBAMcounts.txt","70462-taggedBAMcounts.txt","74665-taggedBAMcounts.txt","74666-taggedBAMcounts.txt","74667-taggedBAMcounts.txt","74668-taggedBAMcounts.txt","74669-taggedBAMcounts.txt","74670-taggedBAMcounts.txt","74671-taggedBAMcounts.txt","74672-taggedBAMcounts.txt","74673-taggedBAMcounts.txt","74674-taggedBAMcounts.txt","74675-taggedBAMcounts.txt","74676-taggedBAMcounts.txt","74677-taggedBAMcounts.txt","74678-taggedBAMcounts.txt","74679-taggedBAMcounts.txt","74680-taggedBAMcounts.txt","74681-taggedBAMcounts.txt","74682-taggedBAMcounts.txt","74683-taggedBAMcounts.txt","74684-taggedBAMcounts.txt"), path="analysis/genome-miRNAcounts/")
#   # Don't miss the last slash in path!!
# {
#   nFiles <- length(fn)
#   for(f in 1:nFiles)
#   {
#     if(f==1){
#      f1 <- read.table(paste(path, fn[f], sep=""), sep="\t", head=T, fill=T)
#     }else{
#      f2 <- read.table(fn[f], sep="\t", head=T, fill=T)
#      f1 <- merge(f1, f2, by="miRNA", all=T)
#     }
#   }	
# return(f1)
# }

# outnew <- readNmerge()
# write.table(outnew, file="Samples-genomemiRNAcounts.tsv",sep="\t",row.names=F)
