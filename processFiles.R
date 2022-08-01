processFiles = function(dt){
  res1=list()
  res2=list()
  i=1
  flag=1
  seq=''
  seqCount=1
  for(i in 1:nrow(dt)){
    k=dt[[1]][i]
    if(length(k)==0){
      print(i)
      break
    }
    if(substr(k,1,1)=='>') {
      if(flag==0) {
        res2[[seqCount]]=seq #res2 won't add when reading first line
        seqCount=seqCount+1
      }
      res1[[seqCount]]=substr(k,2,nchar(k))
      flag=1
      seq=''
    }
    else{
      if(flag==1 & substr(k,1,1)=='>') break
      seq=paste0(seq,k)
      flag=0
    }
  }
  res2[[seqCount]]=seq #add the last sequence
  print(paste0("there are in total ",seqCount," seqs!"))
  return(data.table(anno=unlist(res1),seq=unlist(res2)))
}

protein_seqs = fread('./Homo_sapiens.GRCh38.pep.all.fa', 
                     header = F, 
                     na.strings = "Char")

coding_seq = fread('./Homo_sapiens.GRCh38.cds.all.fa', 
                     header = F, 
                     na.strings = "Char")

proteinSeq = processFiles(protein_seqs)
cdsSeq = processFiles(coding_seq)
# proteinSeq need to be further formatted
txIDs = gsub(sapply(strsplit(proteinSeq$anno,split = ' '),function(x)x[5]), 
             pattern = '.*t:(.*)\\..*',replacement = '\\1')
txIDs = gsub(sapply(strsplit(cdsSeq$anno,split = ' '),function(x)x[1]), 
             pattern = '(.*)\\..*',replacement = '\\1')
rm(proteinSeq)
protein_seqs=data.table(ensembl_transcript_id=txIDs, peptide = proteinSeq$seq)
coding_seq=data.table(ensembl_transcript_id=txIDs, coding = cdsSeq$seq)

all(protein_seqs$ensembl_transcript_id == coding_seq$ensembl_transcript_id)

mRNA=fread('./Homo_sapiens.GRCh38.cdna.all.fa', 
           header = F, 
           na.strings = "Char")
cdnaSeq = processFiles(mRNA)
txIDs = gsub(sapply(strsplit(cdnaSeq$anno,split = ' '),function(x)x[1]), 
             pattern = '(.*)\\..*',replacement = '\\1')
cdna_seq=data.table(ensembl_transcript_id=txIDs, mRNA = cdnaSeq$seq)
table(coding_seq$ensembl_transcript_id %in% cdna_seq$ensembl_transcript_id)
combineSeq = data.table(TxID = protein_seqs$ensembl_transcript_id, 
                        proteinSeq =protein_seqs$peptide)
combineSeq[, coding:=coding_seq[combineSeq$TxID, coding, on='ensembl_transcript_id']]
combineSeq[, cds:=cdna_seq[combineSeq$TxID, mRNA, on='ensembl_transcript_id']]
#combineSeq$cds3utr = 
  
tmp=list()
for (x in 1:nrow(combineSeq))
  tmp[[x]]=gsub(combineSeq$cds[x], 
        pattern = paste0(".*(",combineSeq$coding[x],".*)"), 
        replacement = "\\1")

