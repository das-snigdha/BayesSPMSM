# Code to clean the raw data and store it as an .Rdata file
rm(list = ls())

# read the raw data from .csv file
obs.full = read.csv("./data/multitooth.csv")

# add upper jaw indicator as a covariate
obs.full$Jaw = 0
obs.full$Jaw[obs.full$tooth < 16] = 1

# extract the necessary variables & rename them
obs = obs.full[,c(1,2,5,4,6,7,11,8,10)]
names(obs)[3:9] = c("x1", "x2", "x3", "x4", "x5", "C", "State")
obs = obs[order(obs$id, obs$tooth),]

# re-index the states as 0, 1, 2, 3
obs$State = obs$State - 1

# re-index the teeth as 1, 2, ..., 28
obs$tooth[obs$tooth <= 16] = obs$tooth[obs$tooth <= 16] - 1
obs$tooth[obs$tooth >= 17] = obs$tooth[obs$tooth >= 17] - 3

# save the data as a matrix
obs = as.matrix(obs)
save(obs, file = "data/PDdata.RData")
