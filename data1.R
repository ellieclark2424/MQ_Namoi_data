library(usethis)
use_git_config(user.name = "ellieclark2424", user.email = "eloise.clark@students.mq.edu.au")
create_github_token()
install.packages("gitcreds")
library(gitcreds)
gitcreds_set() 


library(readr)
data <- read_csv("COPY Copy of Final completed results- shing(WQ + Stygo).csv")
data = as.data.frame(data)
View(data)
data$Species_richness
data = data[data$Species_richness>0,]
data = data[-4,] # got rid of row with only 2 stygos in it
count <- data[,11:ncol(data)]
count = as.data.frame(count[,-ncol(count)])
count = count[,-ncol(count)]
View(count)

data$tn=as.numeric(data$tn)
is.numeric(data$tn)
data$doc=as.numeric(data$doc)
is.numeric(data$doc)

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

data[data=="N/A"]=NA
stygos = as.data.frame(data)
plot(stygos$doc,alpha,log="xy")
plot(stygos$pH,alpha,log="xy")
plot(stygos$`GW level change (m)`,alpha,log="xy")
plot(stygos$Temp,alpha,log="xy")
plot(stygos$`EC uS/cm`,alpha,log="xy")
plot(stygos$`DO mg/L`,alpha,log="xy")
plot(stygos$tn,alpha,log="xy")
# stygos is the data frame, and doc is a column or element within it

rownames(data) = 1:52

colnames(data) 
EV = as.data.frame(data[,4:10])
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

data[rownames(EV),]
rownames(EV)
rownames(data)
F2 = data[rownames(EV),]
summary(lm(log(shortalpha)~p + F2$Catchment))
data

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
  # pca on species matrix and use those scores to do similar plot (log data set?)
   # colour the points by bore hole
   # points scale to alpha
# try correspondence analysis
# make folder "r" to put r script in
# data folder, r folder, and output folder (possible meta folder)

# ask john about correspondence analysis next time we meet
# package = vegan
# function = CCA

# Diversity - Fisher
install.packages("vegan")
library(vegan)
str(data)
sa <- data[, 11:27]
sa <- data.frame(lapply(sa, as.numeric))
data$Fishers_alpha <- fisher.alpha(sa)
head(data$Fishers_alpha)

# Trying other diversity methods
# Richness = number of species present
data$Richness <- specnumber(sa)
head(data$Richness)
# Shannon diversity
data$Shannon <- exp(diversity(sa, index = "shannon"))
head(data$Shannon)
data$Shannon
# Simpson diversity
data$Simpson <- simpson.unb(sa, inverse=T)
head(data$Simpson)
data$Simpson



# ggplots for bore ID and date
library(ggplot2)
library(dplyr)
data %>% ggplot(aes(x= pH, y= Species_richness)) +
  geom_jitter() +
  facet_grid(Bore_ID~Date)
data %>% ggplot(aes(x= doc, y= Species_richness)) +
  geom_jitter() +
  facet_grid(Bore_ID~Date)
data %>% ggplot(aes(x= DO, y= Species_richness)) +
  geom_jitter() +
  facet_grid(Bore_ID~Date)
data %>% ggplot(aes(x= EC, y= Species_richness)) +
  geom_jitter() +
  facet_grid(Bore_ID~Date)
data %>% ggplot(aes(x= Temp, y= Species_richness)) +
  geom_jitter() +
  facet_grid(Bore_ID~Date)
data %>% ggplot(aes(x= Catchment, y= pH)) +
  geom_jitter() +
  facet_grid(Bore_ID~Date)

data %>%
  ggplot(aes(x = Temp, y = Species_richness, color = Catchment)) +
  geom_jitter() +  # Adds jittered points to avoid overplotting
  facet_grid(Catchment ~ Date) +  # Facets by Catchment and Date
  theme_minimal()  # Optional: Adds a minimal theme
data %>%
  ggplot(aes(x = pH, y = Species_richness, color = Catchment)) +
  geom_jitter() +  # Adds jittered points to avoid overplotting
  facet_grid(Catchment ~ Date) +  # Facets by Catchment and Date
  theme_minimal()  # Optional: Adds a minimal theme
data %>%
  ggplot(aes(x = EC, y = Species_richness, color = Catchment)) +
  geom_jitter() +  # Adds jittered points to avoid overplotting
  facet_grid(Catchment ~ Date) +  # Facets by Catchment and Date
  theme_minimal()  # Optional: Adds a minimal theme
data %>%
  ggplot(aes(x = DO, y = Species_richness, color = Catchment)) +
  geom_jitter() +  # Adds jittered points to avoid overplotting
  facet_grid(Catchment ~ Date) +  # Facets by Catchment and Date
  theme_minimal()  # Optional: Adds a minimal theme
data %>%
  ggplot(aes(x = tn, y = Species_richness, color = Catchment)) +
  geom_jitter() +  # Adds jittered points to avoid overplotting
  facet_grid(Catchment ~ Date) +  # Facets by Catchment and Date
  theme_minimal()  # Optional: Adds a minimal theme
data %>%
  ggplot(aes(x = doc, y = Species_richness, color = Catchment)) +
  geom_jitter() +  # Adds jittered points to avoid overplotting
  facet_grid(Catchment ~ Date) +  # Facets by Catchment and Date
  theme_minimal()  # Optional: Adds a minimal theme

install.packages("lmerTest")
library(lmerTest)

# Mixed effects model
data$Date <- as.factor(data$Date)
data$Catchment <- as.factor(data$Catchment)
data$Bore_ID <- as.factor(data$Bore_ID)

library(car)
vif(lm(Species_richness ~ GW_level_change + Temp + pH + EC + DO + tn + doc + Catchment, data = data))

data_scaled <- data
numeric_vars <- c("GW_level_change", "Temp", "pH", "EC", "DO", "tn", "doc")
data_scaled[numeric_vars] <- scale(data_scaled[numeric_vars])

model <- lmer(Species_richness ~ GW_level_change + Temp + pH + EC + DO + tn + doc + Catchment + (1 | Date) + (1 | Bore_ID),
              data = data_scaled,
              control = lmerControl(optimizer = "bobyqa"))

# corrplot shows only the relationship between the variables and nothing about the causality
install.packages("corrplot")
library(corrplot)
numeric_vars <- data[, c("GW_level_change", "Temp", "pH", "EC", "DO", "tn", "doc")]
cor_matrix <- cor(numeric_vars, use = "complete.obs")
corrplot(cor_matrix, method = "color")

str(data[, c("GW_level_change", "Temp", "pH", "EC", "DO", "tn", "doc")])
numeric_vars <- data[, c("GW_level_change", "Temp", "pH", "EC", "DO", "tn", "doc")]
numeric_vars <- data.frame(lapply(numeric_vars, function(x) as.numeric(as.character(x))))

cor_matrix <- cor(numeric_vars, use = "pairwise.complete.obs")
corrplot(cor_matrix, method = "color")
install.packages("mnormt")
install.packages("psych")
library(psych)
corr_results <- corr.test(numeric_vars, use = "pairwise")
corrplot(corr_results$r, method = "color")


model <- lmer(Species_richness ~ GW_level_change + Temp + pH + EC + DO + tn + doc + Date + 
                (1 | Catchment) + (1 | Bore_ID), data = data)
summary(model)

# Syntax model = lmer(target ~ fixed1 + fixed2+ (1|random1) + (1|random2), data=data)
library(lme4)
lmer(Species_richness ~ GW_level_change + Temp + pH + EC + DO + tn + doc + (1 | Bore_ID), data = data)
modelM = lmer(log(Species_richness) ~ GW_level_change + Temp + pH + EC + DO + tn + doc + (1 | Bore_ID) + (1 | Date) + (1 | Catchment), data = data)
summary(modelM)
# positive or negative effect
ranef(modelM)
rand(modelM)

library(lme4)
modelM = lmer(log(alpha) ~ GW_level_change + Temp + pH + EC + DO + tn + doc + (1 | Bore_ID) + (1 | Date) + (1 | Catchment), data = data)
summary(modelM)
# positive or negative effect
ranef(modelM)
rand(modelM)
modelM = lmer(log(Simpson) ~ GW_level_change + Temp + pH + EC + DO + tn + doc + (1 | Bore_ID) + (1 | Date) + (1 | Catchment), data = data)
summary(modelM)
# positive or negative effect
ranef(modelM)
rand(modelM)
modelM = lmer(log(Shannon) ~ GW_level_change + Temp + pH + EC + DO + tn + doc + (1 | Bore_ID) + (1 | Date) + (1 | Catchment), data = data)
summary(modelM)
# positive or negative effect
ranef(modelM)
rand(modelM)
is.numeric(data$tn)
is.numeric(data$doc)
hist(log(alpha))
cor.test(data$pH,data$Fishers_alpha, method="s")


# Mixed effects model using catchment and bore ID as fixed effects
lmer(Species_richness ~ GW_level_change + Temp + pH + EC + DO + tn + doc + Catchment + (1 | Bore_ID), data = data)
str(data$tn)
data$tn <- as.numeric(as.character(data$tn))
# Scaling data and changing all data to numeric
vars_to_scale <- c("GW_level_change", "Temp", "pH", "EC", "DO", "tn", "doc")
for (v in vars_to_scale) {
  data[[v]] <- as.numeric(as.character(data[[v]]))
}
data[vars_to_scale] <- lapply(data[vars_to_scale], scale)
summary(data[vars_to_scale])

# Testing different random effects
model <- lmer(Species_richness ~ GW_level_change + Temp + pH + EC + DO + tn + doc + Catchment + (1 | Bore_ID), data = data)
coef(model)

lmer(Species_richness ~ GW_level_change + Temp + pH + EC + DO + tn + doc + (1 | Bore_ID), data = data)
model2 <- lmer(Species_richness ~ GW_level_change + Temp + pH + EC + DO + tn + doc + Bore_ID + (1 | Catchment), data = data)

lmer(Species_richness ~ GW_level_change + Temp + pH + EC + DO + tn + doc + (1 | Bore_ID) + (1 | Date), data = data)

model1 <- lmer(Species_richness ~ GW_level_change + Temp + pH + EC + DO + tn + doc + (1 | Catchment) + (1 | Bore_ID), data = data)
coef(model1)

lmer(Species_richness ~ GW_level_change + Temp + pH + EC + DO + tn + doc + (1 | Date), data = data)


# catchment specific affects
lmer(Species_richness ~ (GW_level_change + Temp + pH + EC + DO + tn + doc) * Catchment + (1 | Bore_ID), data = data)

# Tests:
# Main effects of each hydrochem variable
# Main effect of Catchment
# Interaction between Catchment and each hydrochem variable (i.e., whether the relationships differ between catchments)




