setwd("~/Data/FAR_paper/1_PrepareData/1_EBM/3_Scripts/")
source("../2_Tools/held_model_AS_algo.r")
FF_data <- read.table("../1_Data/Time_Prof_all_RF.txt", header=TRUE)

volcanic_forcing <- FF_data$Volcanic
volcanic_forcing <- c(volcanic_forcing, numeric(length(2012:2100)))
volcanic_effect <- held_model(volcanic_forcing)
volcanic_df <- data.frame("year"=1750:2100, "volcanic_effect"=volcanic_effect[-1,"T"])
write.table(volcanic_df, file="../4_Results/volcanic_effect.txt", row.names=FALSE)

nat_forcing <- apply(subset(FF_data, select=c("Volcanic", "Solar")), 1, sum)
nat_forcing <- c(nat_forcing, numeric(length(2012:2100)))
nat_effect <- held_model(nat_forcing)
nat_df <- data.frame("year"=1750:2100, "nat_effect"=nat_effect[-1,"T"])
write.table(nat_df, file="../4_Results/nat_effect.txt", row.names=FALSE)

aer_forcing <- apply(subset(FF_data, select=c("aerosolERF")), 1, sum)
aer_forcing <- c(aer_forcing, numeric(length(2012:2100)))
aer_effect <- held_model(aer_forcing)

ghg_forcing <- apply(subset(FF_data, select=c("CO2", "OtherWMGHG")), 1, sum)
ghg_forcing <- c(ghg_forcing, numeric(length(2012:2100)))
ghg_effect <- held_model(ghg_forcing)

all_forcing <- apply(FF_data[, -1], 1, sum)
all_forcing <- c(all_forcing, numeric(length(2012:2100)))
all_effect <- held_model(all_forcing)
all_df <- data.frame("year"=1750:2100, "all_effect"=all_effect[-1,"T"])

nat_aer_ghg_df <- data.frame("year"=1750:2100, "nat"=nat_effect[-1,"T"], "aer"=aer_effect[-1,"T"], "ghg"=ghg_effect[-1, "T"]) 
write.table(nat_aer_ghg_df, file="../4_Results/nat_aer_ghg_df.txt", row.names=FALSE)
