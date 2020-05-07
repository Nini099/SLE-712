#Question 1.Read in the file, making the gene accession numbers the row names. Show a table of values for the first six genes
list.files()
x<-read.table("data/gene_expression.tsv")
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
#Question 3. 
order(x$Mean_cal)
order(-x$Mean_cal)
x[order(-x$Mean_cal),]
Highest_mean <-x[order(-x$Mean_cal),]
head(Highest_mean,10)
#question4 
subset(x,Mean_cal<10)

#Question6. Import this csv file into an R object. What are the column names? 
y<- read.csv("data/growth_data.csv")
head(y)
y <- read.csv("data/growth_data.csv", header = TRUE)
head(y)
y <- read.csv("data/growth_data.csv", header = TRUE, stringsAsFactors = FALSE)
head(y)
str(y)
#Name of the columns are- "site", "TreeID","Circumf_2004_cm", "Circumf_2009_cm", "Circumf_2014_cm " and "Circumf_2019_cm"

#Question7. Calculate the mean and standard deviation of tree circumference at the start and end of the study at both sites. 

mean(y$Circumf_2004_cm & y$Circumf_2019_cm)
