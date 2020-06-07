---
title:"SLE 712 Assessment 3"
Author: "Nishat Nini Urmi"
date: 4/6/2020
output: html_document
---

##Part 1: Part 1: Importing files, data wrangling, mathematical operations, plots and saving code on GitHub 

#Question 1.Read in the file, making the gene accession numbers the row names. Show a table of values for the first six genes

#Read the file and save it as "x" and making the gene accession numbers as row names 
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/gene_expression.tsv", destfile = "gene_expression.tsv" )
x<- read.table("data/gene_expression.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
head(x) #to see the first and last part of the target table
str(x) #to see the interanl structure of "x"
#to make the table of first six genes in the table, where inside the bracket [row,colum] 
x[1:6,1:2]

#Question 2. Make a new column which is the mean of the other columns. Show a table of values for the first six genes.
#to make a new row with mean values from other rows
x$Mean_cal <- rowMeans(x)
head(x)
x[1:6,1:3] # as the new row is column three

#Question 3. List the 10 genes with the highest mean expression
 # the function to show data in a selected table from highest to lowest 
x[order(-x$Mean_cal),]
Highest_mean <-x[order(-x$Mean_cal),] # saving the list as a different object
head(Highest_mean,10) #to list the 10 genes with highest mean expression

#Question4  Determine the number of genes with a mean <10  
filteredmean <-subset(x,Mean_cal<10) #subset the genes 
nrow(filteredmean)# shows the total number of genes with a mean<10 

#Question 5.  Make a histogram plot of the mean values in png format and paste it into your report

hist(x$Mean_cal,main="Histogram for Mean values",xlab = "mean value",border = "red",col="blue",breaks = 50,xlim=c(0,10000))

#Question6. Import this csv file into an R object. What are the column names? 
#to download the file from provided link
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/growth_data.csv", destfile = "growth_data.csv" )
#saved the file as object "y"
y <- read.csv("data/growth_data.csv", header = TRUE, stringsAsFactors = FALSE)
head(y) #to determine the column names
str(y)

#Question7. Calculate the mean and standard deviation of tree circumference at the start and end of the study at both sites. 
# Subset two "sites" as two different objects for calculation
NE <- subset(y,Site=="northeast")
SW <- subset(y,Site=="southwest")
head(NE)
head(SW)
# Mean calculation of tree circumference at the start and end of the study for North East Site
mean(NE$Circumf_2004_cm)
mean(NE$Circumf_2019_cm)
sd(NE$Circumf_2004_cm)
sd(NE$Circumf_2019_cm)
#Mean  and standard deviation calculation of tree circumfurence at the satrt and end of the site Southwest
mean(SW$Circumf_2004_cm)
mean(SW$Circumf_2019_cm)
sd(SW$Circumf_2004_cm)
sd(SW$Circumf_2019_cm)

# Question 8.Make a box plot of tree circumference at the start and end of the study at both sites.
#boxplot function used with headings
boxplot(NE$Circumf_2004_cm,NE$Circumf_2019_cm,SW$Circumf_2004_cm,SW$Circumf_2019_cm,names=c("NE2004","NE2019","SW2004","SW2019"),ylab="Circumference(cm)",main="Growth at 2 plantation sites")

#Question 9. Calculate the mean growth over the past 10 years at each site.
#new column for growth over the past 10 years for NorthEast 
NE$Growth <- (NE$Circumf_2019_cm-NE$Circumf_2009_cm ) 
head(NE)
#new column for gowth over the past 10 years for SouthEast
SW$Growth <- (SW$Circumf_2019_cm-SW$Circumf_2009_cm)
head(SW)

#Questoin 10..Use the t.test and wilcox.test functions to estimate the p-value that the 10 year growth is different at the two sites
t.test(SW$Growth,NE$Growth)
wilcox.test(SW$Growth,NE$Growth)



#PART 2 
#The sequence allocated for me is 53

#to load required packages for working with sequence data   
library("seqinr")
library("R.utils")
library("rBLAST")
library("ape")
library("ORFik")
library("Biostrings")

#Question 1. Download the whole set of E. coli gene DNA sequences and use gunzip to decompress. Use the makeblast() function to create a blast database. How many sequences are present in the E.coli  set? 
#to download the file from provided link 
download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",
              destfile = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz") 
#To unzip the files
R.utils::gunzip("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz", overwrite= TRUE)
#to apply the provided BLAST function and determine the number of base pairs 
makeblastdb("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa",dbtype = "nucl","-parse_seqids")

#Question 2 Download the sample fasta sequences and read them in as above. For your allocated sequence, determine the length (in bp) and the proportion of GC bases
#To download sample file 
if ( ! file.exists("sample.fa") ) {
  download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa",
                destfile = "sample.fa")
}  
#to read the fasta file   
Ecoli_Sample <- read.fasta("sample.fa")
#renaming the allocated sequence
allocated_Seq <- Ecoli_Sample$`53`
#to calculate the length of the sequence
seqinr::getLength(allocated_Seq)
table(allocated_Seq)
#to calculate the proportion of GC
seqinr::GC(allocated_Seq)

#Question 3.You will be provided with R functions to create BLAST databases and perform blast searches. Use blast to identify what E. coli  gene your sequence matches best. Show a table of the top 3 hits including percent identity, E-value and bit scores.
#to download the provided function from the source to make BLAST databases for this experiment
source("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R")
myblastn_tab #provided function
res <- myblastn_tab(allocated_Seq, db="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
str(res)
res
#to see top three hits 
hits <- as.character(res$sseqid[1:3])
hits

#Quistion 4 You will be provided with a function that enables you to make a set number of point mutations to your sequence of interest. Run the function and write an R code to check the number of mismatches between the original and mutated sequence
# create a copy of my allocated sequnce that has certain number of mismathces in it
seqinr::write.fasta(allocated_Seq,names="allocated_Seq",file.out = "allocated_Seq.fa")
makeblastdb("allocated_Seq.fa",dbtype="nucl", "-parse_seqids")
#save it as res 
res <- myblastn_tab(allocated_Seq, db="allocated_Seq.fa")
res
#taste with a specific number to show the mismatches
my_allocated_mutator<-mutator(allocated_Seq,20)
res <- myblastn_tab(my_allocated_mutator, db="allocated_Seq.fa") 

 

#Question 5 Using the provided functions for mutating and BLASTing a sequence, determine the number and proportion of sites that need to be altered to prevent the BLAST search from matching the gene of origin. Because the mutation is random, you may need to run this test multiple times to get a reliable answer

#at first we have to create a fucntion that can mutate+blast and then summerise the result as a 0 or 1. 
myfunc <- function(myseq,nmut) { 
  mutseq <- mutator( myseq= allocated_Seq, nmut = nmut) #the sequence for mutation, nmut=nmut as it will recognize the number after the sequence as nmut 
  res <- myblastn_tab(myseq= allocated_Seq,db= "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa") #for blast
  if (is.null(res)) {myres= 0} else {myres = 1} 
  return(myres) # to summerise the result, depending on if the blast was or wasn't successful
}
myfunc(myseq,100) #applying the created function
replicate(n = 100, myfunc(myseq,100)) #to repeat the function several times that would give a vector of 100 values
n <-c(0,100,200,300)
mean(replicate(100,myfunc(myseq,100 ))) #

