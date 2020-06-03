---
title:"SLE 712 Assessment 3"
Author: "Nishat Nini Urmi"
date: 4/6/2020
output: html_document
---

##Part 1: Part 1: Importing files, data wrangling, mathematical operations, plots and saving code on GitHub 

#Question 1.Read in the file, making the gene accession numbers the row names. Show a table of values for the first six genes

x<-read.table("data/gene_expression.tsv") # 
head(x)
str(x)
x <- read.table("data/gene_expression.tsv", header = TRUE)
head(x)
str(x)
x<- read.table("data/gene_expression.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
head(x)
str(x)
x[1:6,1:2]
#Question 2. Make a new column which is the mean of the other columns. Show a table of values for the first six genes.
x$Mean_cal <- rowMeans(x)
head(x)
x[1:6,1:3]
#Question 3. List the 10 genes with the highest mean expression
order(-x$Mean_cal)
x[order(-x$Mean_cal),]
Highest_mean <-x[order(-x$Mean_cal),]
head(Highest_mean,10)
#question4  Determine the number of genes with a mean <10  
subset(x,Mean_cal<10)
#Question 5.  Make a histogram plot of the mean values in png format and paste it into your report
range(x$Mean_cal)
hist(x$Mean_cal,main="Histogram for Mean values",xlab = "mean value",border = "red",col="blue",breaks = 50,xlim=c(0,10000))
#Question6. Import this csv file into an R object. What are the column names? 
y<- read.csv("data/growth_data.csv")
head(y)
str(y) 
y <- read.csv("data/growth_data.csv", header = TRUE)
head(y)
y <- read.csv("data/growth_data.csv", header = TRUE, stringsAsFactors = FALSE)
head(y)
str(y)
subset(y,Site=="northeast")
NE <- subset(y,Site=="northeast")
SW <- subset(y,Site=="southwest")
head(NE)
head(SW)
#Name of the columns are- "site", "TreeID","Circumf_2004_cm", "Circumf_2009_cm", "Circumf_2014_cm " and "Circumf_2019_cm"

#Question7. Calculate the mean and standard deviation of tree circumference at the start and end of the study at both sites. 
mean(NE$Circumf_2004_cm)# Mean calculation of tree circumference at the start and end of the study for North East Site
mean(NE$Circumf_2019_cm)
sd(NE$Circumf_2004_cm)
sd(NE$Circumf_2019_cm)

mean(SW$Circumf_2004_cm)#Mean  and standard deviation calculation of tree circumfurence at the satrt and end of the site Southwest
mean(SW$Circumf_2019_cm)
sd(SW$Circumf_2004_cm)
sd(SW$Circumf_2019_cm)

# Question 8.Make a box plot of tree circumference at the start and end of the study at both sites.

boxplot(NE$Circumf_2004_cm,NE$Circumf_2019_cm,SW$Circumf_2004_cm,SW$Circumf_2019_cm,names=c("NE2004","NE2019","SW2004","SW2019"),ylab="Circumference(cm)",main="Growth at 2 plantation sites")
#Question 9. Calculate the mean growth over the past 10 years at each site.
NE$Growth <- (NE$Circumf_2019_cm-NE$Circumf_2009_cm ) #new column for growth over the past 10 years for NorthEast 
head(NE)
SW$Growth <- (SW$Circumf_2019_cm-SW$Circumf_2009_cm)#new column for gowth over the past 10 years for SouthEast
head(SW)
#Questoin 10..Use the t.test and wilcox.test functions to estimate the p-value that the 10 year growth is different at the two sites
t.test(SW$Growth,NE$Growth)
wilcox.test(SW$Growth,NE$Growth)
# the P-value for t-test is 1.713e-06, and for Wilcox test it is 4.626e-06.


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
#Null 



res
#Question 5
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

