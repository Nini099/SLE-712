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
mean(NE$Circumf_2004_cm)#for 
mean(NE$Circumf_2019_cm)
mean(SW$Circumf_2004_cm)
mean(SW$Circumf_2019_cm)
sd(NE$Circumf_2004_cm)
sd(NE$Circumf_2019_cm)
sd(SW$Circumf_2004_cm)
sd(SW$Circumf_2019_cm)
# Question 8.Make a box plot of tree circumference at the start and end of the study at both sites.
boxplot(NE$Circumf_2004_cm,NE$Circumf_2019_cm,SW$Circumf_2004_cm,SW$Circumf_2019_cm)
boxplot(NE$Circumf_2004_cm,NE$Circumf_2019_cm,SW$Circumf_2004_cm,SW$Circumf_2019_cm,names=c("NE2004","NE2019","SW2004","SW2019"),ylab="Circumference(cm)",main="Growth at 2 plantation sites")
#Question 9. Calculate the mean growth over the past 10 years at each site.
NE$Growth <- (NE$Circumf_2019_cm-NE$Circumf_2009_cm )
head(NE$Growth)
SW$Growth <- (SW$Circumf_2019_cm-SW$Circumf_2009_cm)
head(SE$Growth)
head(y)
#Questoin 10..Use the t.test and wilcox.test functions to estimate the p-value that the 10 year growth is different at the two sites
t.test(SW$Growth,NE$Growth)
wilcox.test(SW$Growth,NE$Growth)
# the P-value for t-test is 1.713e-06, and for Wilcox test it is 4.626e-06.