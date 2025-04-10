library(readr)
F <- read_csv("COPY Copy of Final completed results- shing(WQ + Stygo).csv")
F = as.data.frame(F)
View(F)

View(F)
count <- F[,11:ncol(F)]
View(count)
count = as.data.frame(count[,-ncol(count)])
count = count[,-ncol(count)]
View(count)

#parameters of geometric
p = 1/(colMeans(count)+1)
# allows us to find fit

d = p[1]*(1-p[1])^(0:152)
plot(d)

#running cu
t = table(count[,1])
t

sad = rep(0,152)
sad[as.numeric(names(t))+1]=t
sad
plot(cumsum(d),cumsum(sad)/sum(sad)) 
# cumsum used to calculate cumulative sum of the vector passed as argument
abline(0,1) # adds a straight line to the plot

# log series

# want to compute fishers alpha
# reference = fisher et al. 2023
#

# makes diversity values to compared to WQ
source("fisher.R")
# "fisher" used to combine the p-values from several independent tests bearing upon the 
# same overall null hypothesis
alpha = array() 
# alpha refers to the significance level, which is the probability of rejecting 
# a true null hypothesis (a Type I error), often set at 0.05 (or 5%)
# initializes alpha as an empty array
for (i in (which(rowSums(count)>0)))
  alpha[i]=fisher(count[i,count[i,]>0])
# count is assumed to be a matrix or data frame
# rowSums(count) > 0 identifies rows in count where the sum across columns is greater than zero
# which(rowSums(count) > 0) gives the row indices that satisfy this condition
# the loop iterates over these row indices (i)
#
# count[i,] extracts row i from the count matrix
# count[i,] > 0 creates a logical vector indicating which elements in row i are greater than zero
# count[i, count[i,] > 0] selects only the nonzero elements from row i
# fisher() is applied to these nonzero elements
# The result is stored in alpha[i]
sort(alpha)
# sorts alpha in ascending order and returns the result

# Summary of code above:
# Loads an external script (fisher.R).
# Initializes an empty alpha array.
# Loops over rows of count that have at least one nonzero element.
# For each of these rows:
  # Extracts nonzero values.
  # Passes them to the fisher() function.
  # Stores the result in alpha.
# Sorts alpha and returns the sorted values.

stygos = as.data.frame(F)

plot(stygos$doc,alpha,log="xy")
plot(stygos$pH,alpha,log="xy")
plot(stygos$`GW level change (m)`,alpha,log="xy")
plot(stygos$Temp,alpha,log="xy")
plot(stygos$`EC uS/cm`,alpha,log="xy")
plot(stygos$`DO mg/L`,alpha,log="xy")
plot(stygos$tn,alpha,log="xy")
# stygos is the data frame, and doc is a column or element within it

F[F=="N/A"]=NA

rownames(F) = 1:52

colnames(F) 
EV = as.data.frame(F[,4:10])
rownames(EV) = 1:nrow(EV)
EV = as.data.frame(na.omit(EV))
rownames(EV)

EV = matrix(as.numeric(unlist(EV)),nrow(EV),ncol(EV))
princomp(EV)
is.numeric(EV)
table(unlist(EV))
table(table(unlist(EV)))
princomp(as.matrix(scale(EV)))
is.numeric(EV)

# scores of PCA
p = princomp(scale(EV))$scores
princomp(scale(EV))$loadings
rownames(EV)
plot(p)

shortalpha = alpha[as.numeric(rownames(EV))]
summary(lm(log(shortalpha)~p))

F[rownames(EV),]
rownames(EV)
rownames(F)
F2 = F[rownames(EV),]
summary(lm(log(shortalpha)~p + F2$Catchment))
F

# Factor analysis
fa = factanal(EV, scores = "regression", factors = 3)
is.numeric(EV)

# point shapes to compare sites
plot(p, pch=as.integer(as.factor(F2$Catchment)))
# =="MQ" is for wanting the MQ values
w = which(F2$Catchment=="MQ")
h = chull(p[F2$Catchment=="MQ",1:2])
polygon(p[w[h],], col=hsv(h=0.6, alpha=0.2))
# !="MQ" is for not wanting the MQ values
w = which(F2$Catchment!="MQ")
h = chull(p[F2$Catchment!="MQ",1:2])
polygon(p[w[h],], col=hsv(h=0.5, alpha=0.2))

# both axes (x and y) are logarithmic

cor.test(as.numeric(stygos$doc),alpha,method="s")

# to do: poisson regression, (regression and linear modelling in class), more information on count data
# alpha is important = is significant or not
# composition differences?
# correspondence analysis = polygon 
# CEGS paper (section "alternative species richness estimators")

# describes the chances of achieving success in a series of independent trials, each having two possible outcomes

# mixed effect model?? tells something else?
# compositional difference - catchments different?
  # pca on species matrix and use those scores to do similar plot


# ask john about correspondence analysis next time we meet
# package = vegan
# function = CCA


