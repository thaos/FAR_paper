setwd("~/Data/FAR_paper/1_PrepareData/2_FormatData/3_Scripts/")
source("../2_Tools/ReadData_algo.R")

res_folder <- "../4_Results/"
nat <- read.table("../../1_EBM/4_Results/nat_effect.txt",
                  header=TRUE)
names(nat)[2] <- "nat"

tas_cnrm <- create_rds(model="CNRM")
tas_cnrm <- readRDS(tas_cnrm)
tas_nat_cnrm <- merge(tas_cnrm, nat, by=c("year"))
saveRDS(tas_nat_cnrm, file=paste(res_folder, "tas_nat_cnrm.rds", sep=""))

tas_ipsl <- create_rds(model="IPSL")
tas_ipsl <- readRDS(tas_ipsl)
tas_nat_ipsl <- merge(tas_ipsl, nat, by=c("year"))
saveRDS(tas_nat_ipsl, file=paste(res_folder, "tas_nat_ipsl.rds", sep=""))
