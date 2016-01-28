
setwd("~/Data/FAR_paper/1_PrepareData/1_EBM/3_Scripts/")
source("../2_Tools/held_model_AS_algo.r")
FF_data <- read.table("../1_Data/Time_Prof_all_RF.txt", header=TRUE)

all <- rowSums(FF_data[,-1])
nat <- rowSums(subset(FF_data, select=c("Solar", "Volcanic")))

plot(FF_data$year, nat>all)
plot(FF_data$year, all, type="l", col="red")
lines(FF_data$year, nat, col="black")
plot(FF_data$year, all-nat, type="l", col="red")
