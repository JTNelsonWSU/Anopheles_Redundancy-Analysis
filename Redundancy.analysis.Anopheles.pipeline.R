library(vegan) #load packages
library(SoDA)
library(adespatial)

my_data <- read.table(paste("Anopheles.2L.seq.data.txt"), header = T) #load Anopheles sequence data for 2L chrommosome (respnose variables)

env <- read.table(paste("ag1000G.RDA.env.matrix.auto.input.txt"), header=T) #load environmental matrix (predictor variables; autosomal)
#env <- read.table(paste("ag1000G.RDA.env.matrix.2R.input.txt"), header=T) #load environmental matrix (predictor variables; for X chromosome only)

db <- read.table(paste("db.2R.txt"),header=T,row.names=1) #load cartesian geographic points for each sample

env.climate <- env[,3:11] #grab columns
env.climate <- as.matrix(env.climate) #df to matrix

new1 <- as.numeric(unlist(db$dbMEM.1)) #the following unlists the sample coords in the db df (use only if the cartisean points were conveted by a different source).
new2 <- as.numeric(unlist(db$dbMEM.2))
new3 <- as.numeric(unlist(db$dbMEM.3))
new4 <- as.numeric(unlist(db$dbMEM.4))
new5 <- as.numeric(unlist(db$dbMEM.5))
temp <- cbind(new1,new2,new3,new4,new5)
rownames(temp) <- rownames(db)
db <- temp


climate.rda <- rda(my_data ~ env.climate + Condition(db)) #performs RDA. predicting allele frequencies based off climate variables that are conditioned by geographic space. Conditioning to geography helps reduce false-positives in this system and is recommeneded.
RsquareAdj(climate.rda) #adjusting r2

climate.rda #prints the proportion of variance descibed by the constrained (variance explained by environment only) and unconstrained (variance explained by environment and other unknown factors) axes.

png("ag1000G.2R.RDA.scree.plot.updated.png") #prints the scree plot to png
screeplot(climate.rda)
dev.off()
###########################################

#####Post RDA analysis
#####Identifying outliers within the data

###############################################

load.rda <- scores(climate.rda, choices=c(1:5), display="species") #select RDA axes that explain ~90% of the variance (here we used the top 5 axes)
load.rda <- data.frame(load.rda)


load.rda$RDA1 <- load.rda$RDA1 * load.rda$RDA1 #we need to square the loadings to remove negative values so that we can identify outliers by the number of SD above the mean. 
load.rda$RDA2 <- load.rda$RDA2 * load.rda$RDA2
load.rda$RDA3 <- load.rda$RDA3 * load.rda$RDA3
load.rda$RDA4 <- load.rda$RDA4 * load.rda$RDA4
load.rda$RDA5 <- load.rda$RDA5 * load.rda$RDA5

write.table(load.rda, "ag1000G.RDA.2R.output.RDA12345.loadings.updated.txt",quote=F) #saves the output for the top loadings

new1 <- as.numeric(unlist(load.rda$RDA1)) #converts the loadings back to original format after squaring them
new2 <- as.numeric(unlist(load.rda$RDA2))
new3 <- as.numeric(unlist(load.rda$RDA3))
new4 <- as.numeric(unlist(load.rda$RDA4))
new5 <- as.numeric(unlist(load.rda$RDA5))

#temp <- cbind(new1,new2,new3,new4)
temp <- cbind(new1,new2,new3,new4,new5)
rownames(temp) <- rownames(load.rda)
load.rda <- temp


png("ag1000G.RDA.2R.loadings.updated.png",width=500,height=300) #prints the top loadings (top 5) to individual histograms. Better to observe the data before continuing. Make sure it looks correct.
all.hist <- par(mfrow=c(1, 5))
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 
hist(load.rda[,4], main="Loadings on RDA4") 
hist(load.rda[,5], main="Loadings on RDA5") 
# hist(load.rda[,4], main="Loadings on RDA4") 
par(all.hist)
dev.off()


outliers <- function(x,z){  #makes a function that will identify positions that have loadings that are at least 5 SD larger than the average. This cutoff is arbitrary and will vary across studies.
  lims <- mean(x) +  (z * sd(x))     # find loadings + z sd from mean loadings     
  x[x >= lims[1]]               # SNP names in these tails
}

cand1 <- outliers(load.rda[,1],5) # Applys function to the 5 RDA axes
cand2 <- outliers(load.rda[,2],5) 
cand3 <- outliers(load.rda[,3],5) 
cand4 <- outliers(load.rda[,4],5) 
cand5 <- outliers(load.rda[,5],5) 

# ncand <- length(cand1) + length(cand2) + length(cand3) + length(cand4)
ncand <- length(cand1) + length(cand2) + length(cand3) + length(cand4) + length(cand5) # concatinates the lengths and gives the total number of ourliers identified (raw list)
ncand

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1)) #combining outliers and making them a df
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
cand4 <- cbind.data.frame(rep(4,times=length(cand4)), names(cand4), unname(cand4))
cand5 <- cbind.data.frame(rep(5,times=length(cand5)), names(cand5), unname(cand5))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- colnames(cand4) <- colnames(cand5) <- c("axis","snp","loading") #changing column names

cand <- rbind(cand1, cand2, cand3, cand4, cand5) #making one df containing all outliers across RDA axes

cand$snp <- as.character(cand$snp) #converts SNPs to character (if SNPs were used as headers, this function will need to be implmented)

foo <- matrix(nrow=(ncand), ncol=9)  # 9 columns for 9 predictors
colnames(foo) <- colnames(env.climate)

my_data <- data.frame(my_data) #makes df


for (i in 1:length(cand$snp)) {  #this loop calulates the correlation coefficient for each snp against each variable. The variable with the strongest correlation is identified as the predictor. Compare predictors across RDA axes to identify which variable had the largest influence. 
  nam <- cand[i,2]
  snp.gen <- my_data[,nam]
  foo[i,] <- apply(env.climate,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo) # combine cand with predictor information 
head(cand)

length(cand$snp[duplicated(cand$snp)]) #find the number of duplicates (note that a position can be an outlier in more than 1 axis; we need to account for this)

cand <- cand[!duplicated(cand$snp),] #remove duplicates


for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,13] <- names(which.max(abs(bar[4:12]))) # gives the variable
  cand[i,14] <- max(abs(bar[4:12]))              # gives the correlation
}

colnames(cand)[13] <- "predictor" #changing columns names
colnames(cand)[14] <- "correlation"

table(cand$predictor) #prints final table
write.table(cand, "ag1000G.RDA.2R.candidates.updated.txt", row.names=F,quote=F) #prints final table as the main output for downstream analyses
