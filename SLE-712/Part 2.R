#PART 2 
#The sequence allocated for me is 53

#to load required packages   
library("seqinr")
library("R.utils")
library("rBLAST")
library("ape")
library("ORFik")
library("Biostrings")

#Question 1. Download the whole set of E. coli gene DNA sequences and use gunzip to decompress. Use the makeblast() function to create a blast database. How many sequences are present in the E.coli  set? 
if ( ! file.exists("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa") ) {
  download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",
                destfile = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")
}
R.utils::gunzip("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz", overwrite= TRUE)
makeblastdb("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa",dbtype = "nucl","-parse_seqids")
#there are 4140 sequences present in the E.coli set. 

#Question 2
#For Query
if ( ! file.exists("sample.fa") ) {
  download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa",
                destfile = "sample.fa")
  
}
Ecoli_Sample <- read.fasta("sample.fa")
allocated_Seq <- Ecoli_Sample$`53`
allocated_Seq [1:789]
seqinr::getLength(allocated_Seq)
table(allocated_Seq)
#to calculate the proportion of GC
seqinr::GC(allocated_Seq)
#The proportion of GC is 54.24588% 

#Question 3.
source("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R")
myblastn_tab
res <- myblastn_tab(allocated_Seq, db="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
str(res)
res
hits <- as.character(res$sseqid[1:3])
hits
#We can see that the allocated sequence parfectly matches the Source E.coli sequence. 
# Percent Identity-100, E-value - 0 , bitscore - int 1517
db <- read.fasta("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
str(str(db[1:6]))

#Quistion 4 # create a copy of my allocated sequnce that has certain number of mismathces in it
seqinr::write.fasta(allocated_Seq,names="allocated_Seq",file.out = "allocated_Seq.fa")
makeblastdb("allocated_Seq.fa",dbtype="nucl", "-parse_seqids")
res <- myblastn_tab(allocated_Seq, db="allocated_Seq.fa")
res
my_allocated_mutator<-mutator(allocated_Seq,20)
res <- myblastn_tab(my_allocated_mutator, db="allocated_Seq.fa") #
res
#Question 5 # to identify blast 
cointoss < - myblastn_tab(my_allocated_mutator){
  sample(c(0,1),1,replace= TRUE)}


}
cointoss <- function(){sample(c(0,1),1,replace=TRUE}
my_allocated_mutator<-mutator(allocated_Seq,50)
res <-myblastn_tab(my_allocated_mutator, db="allocated_Seq.fa")
res       
#35 mismatches
my_allocated_mutator <-mutator(allocated_Seq,100)
res <- myblastn_tab(my_allocated_mutator, db="allocated_Seq.fa")
res
#68 mismatches
my_allocated_mutator <-mutator(allocated_Seq,200)
res <- myblastn_tab(my_allocated_mutator, db="allocated_Seq.fa")
res
#147 mismatches 
my_allocated_mutator <-mutator(allocated_Seq,150)
res <- myblastn_tab(my_allocated_mutator, db="allocated_Seq.fa")
res
#111 Mismatches 
my_allocated_mutator <-mutator(allocated_Seq,500)
res <- myblastn_tab(my_allocated_mutator, db="allocated_Seq.fa")
res
#Null
my_allocated_mutator <-mutator(allocated_Seq,210)
res <- myblastn_tab(my_allocated_mutator, db="allocated_Seq.fa")
res
#151 mismatches  
my_allocated_mutator <-mutator(allocated_Seq,215)
res <- myblastn_tab(my_allocated_mutator, db="allocated_Seq.fa")
res