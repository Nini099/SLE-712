# SLE-712
Bioinformatics Assessment part 1 and 2 
## Author  

#### Nishat Nini Urmi
#### ID: 219276546

## Motivation
For academic purposes to learn about coding language "R". By answering the question in first part assignment will let me test my ability to use R to perform file reading, oredeirng, subsetting, table generation with relevant columns and chart creation. The second part of the assessment is to learn using provided functions to determine the limit of BLAST . 

## Contents
### Data_file 
contains the "gene_exression.tsv" and "growth_data.csv" files that are needed to answer the questions on **Part 1** of the assessment. 
### Script_file 
Contains the codes for answers from both parts of the assessment with relevent comments about the functions. 

### how to use this file 
a brief description of each codes in this file with its inputs and outputs-

#### reading a (.tsv) file by making it's first column as row names 
can be done by using function read.table by  setting attributes *rownames*= 1 ,takes the value of the column , GeneID and returns into a dataframe. The function __head( )__   with n=6 returns the first six genes of the table and __str ()__ gives the structure of the data frame

function : **read.table ()** 

input : **.tsv file**

output : **data. frame**
```{r}
x<- read.table("data/gene_expression.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
head(x,6)
```

#### To make a new column that will show the mean of other two columns 
The following function assigns a new column name *Mean_cal* to the data frame  __"gene_expression.tsv"__, which is named as "x" . The symbol "$" is used in between the new name of column and the data frame ,that makes the name of the new column. The mean value can be expressed in the column by using the function __rowMeans()__ where we incorporate the dataframe to get the value. The number 6 is used in __head()__ function to see the data frame for first six columns  

function : **rowMeans()**

input : **values for new column**

output : **dataset + new column**

```{r}
x$Mean_cal <- rowMeans(x)
head(x,6)
```
#### List of 10 genes with highest mean expression 
To get the highest value in a data set, the function __order__  ("-" incorporated) is used that orders the values from the *Mean_Cal* column in a descending order. The values are saved as "Highes_mean". Then running the function with the attribute 10 give the list of 10 genes with highest mean expreassion. 
 

function : **order()**

input : **values from selected column**

outout: **data frame* 
```
{r}
Highest_mean <-x[order(-x$Mean_cal),]
head(Highest_mean,10)
```
#### Determining the number of genes that has mean value less than 10 

Filters out the mean values less than 10. the __Subset()__ is used. 

input : **data.frame, row with condition*

output : **data.frame**

Another function __nrow()__ calculates the total number of rows that satisfies the condition 

function : **nrow()**

input: **filtered data.frame**

output: **total number of rows **

```{r}
filteredmean <-subset(x,Mean_cal<10)
nrow(filteredmean)
```

#### Make a histogram plot with the selected column from data frame 

The function __hist()__ to make a histogram plot , where the values are from a selected row "Mean_cal". 

function: **hist()**

input: **data.frame**

output: **histogram plot**

```{r}
hist(x$Mean_cal,main="Histogram for Mean values",xlab = "mean value",border = "red",col="blue",breaks = 50,xlim=c(0,10000)) 
```

#### Read a .csv file and save it as a object
The function __read.csv()__ reads the .csv file, where first row is set as row name. The data frame is saved as an object "y" 

function: **read.csv()**

input: **csv.file**

output : **data.frame**
```{r}
y <- read.csv("data/growth_data.csv", header = TRUE, stringsAsFactors = FALSE)
head(y)
```
#### Calculating the mean and standard deviation 
The function __subset()__ is used to devide the two sites and save it as two seperate data.frame. Then the __mean()__ and __sd()__ functions were used to calculte seperately for the both data.frame , based on seperate columns in within the data.frame .  

function: **subset()**

input: **data.frame, site== site name in the main data frame" 

output: **seperated data.frame**
```{r}
NE <- subset(y,Site=="northeast")
SW <- subset(y,Site=="southwest")
head(NE)
head(SW)
```
function: **mean()**

input: **data.frame**

output: **calculated total mean for the column** 
```{r}
mean(NE$Circumf_2004_cm)
mean(NE$Circumf_2019_cm)
sd(NE$Circumf_2004_cm)
sd(NE$Circumf_2019_cm)
```
### Rmd file 
contains the detailed answers to the assessment questions. 
<<<<<<< HEAD
  

##Packages
```{r}
library("seqinr")
=======
## Packages
The following packages are needed to work with sequnce in R for part 2 
```{r}
library("seqinr") 
>>>>>>> 86aed0bd7f03eff6c1d7c7162aa59eb948f0afc2
library("R.utils")
library("rBLAST")
library("ape")
library("ORFik")
library("Biostrings")
<<<<<<< HEAD
 ``` 

## Installation 
=======
```
##
## Installation
>>>>>>> 86aed0bd7f03eff6c1d7c7162aa59eb948f0afc2
Rstudio 
version 1.2.1335 

LICENSE
