# load functions
packages <- c("scales",
              "gtable",
              "grid",
              "extRemes",
              "quantreg",
              "ncdf")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages,
                           rownames(installed.packages())),
                   repos="http://cran.us.r-project.org")
}
lapply(packages, library, character.only=TRUE)
library(devtools)
devtools::load_all("~/Data/FARpackage/FARg")

setwd("~/Data/FAR_paper/2_Analysis/3_Scripts/")
make_shortcut <- function(folder){
function(file) paste(folder, file, sep="")
}
res_path <- make_shortcut("../4_Results/")
data_path <- make_shortcut("../1_Data/")
tools_path <- make_shortcut("../2_Tools/")
source(tools_path("gam_allnat.r"), chdir=TRUE)
source(tools_path("gam_rmnat.r"), chdir=TRUE)

CO2_cc.nc  <- open.ncdf(data_path("CMIP5_CO2forcings/rcp85.nc"))
CO2 <- data.frame("year"=CO2_cc.nc$var$co2$dim[[1]]$vals, "CO2"=get.var.ncdf(CO2_cc.nc, varid="co2"))
CO2_0 <- CO2$CO2[1]
CO2$CO2 <- unlist(lapply(CO2$CO2, function(x) 5.35 * log(x/CO2$CO2[1])))

CO2_RF  <- read.table("../../1_PrepareData/1_EBM/1_Data/Time_Prof_all_RF.txt", header=TRUE)
CO2_RF$ant <- CO2_RF$aerosolERF + CO2_RF$CO2

tas_nat_cnrm <- readRDS(data_path("tas_nat_cnrm.rds"))
tas_nat_ipsl <- readRDS(data_path("tas_nat_ipsl.rds"))
l_tas_nat <- list("cnrm"=tas_nat_cnrm, "ipsl"=tas_nat_ipsl)
# l_tas_nat <- lapply(l_tas_nat, merge, y=subset(CO2, select=c("year", "CO2")), by="year")
l_tas_nat <- lapply(l_tas_nat, merge, y=subset(CO2, select=c("year", "CO2")), by="year", all.x=TRUE)
l_tas_nat <- lapply(l_tas_nat, merge, y=subset(CO2_RF, select=c("year", "ant")), by="year")



change_anomalie_period <- function(data, y="y", time="year", period=c(1850, 1879)){
  i_pre_ind <- which(data[, time] >= period[1] & data[,time ] <= period[2])
  bias <- mean(data[i_pre_ind, y])
  data[, y] %<>% `-`(., bias)
  attr(data, "bias") <- bias
  data
}
tas_pi_cnrm <- change_anomalie_period(tas_nat_cnrm, y="eur_tas")
tas_pi_ipsl <- change_anomalie_period(tas_nat_ipsl, y="eur_tas")
l_tas_pi <- lapply(l_tas_nat, change_anomalie_period, y="eur_tas")

test_methodo <- function(l_tas_nat, redef_ano=FALSE, correct_ant_bias=TRUE, reuse=FALSE){
  if(redef_ano) 
    l_tas_nat <- lapply(l_tas_nat, change_anomalie_period, y="eur_tas")
  l_samples <- lapply(l_tas_nat, boot_samples)
  l_gam_an <- lapply(l_samples, function(x){
                       print(x)
                       boot_gam_allnat(x, "eur_tas", "year", "nat", "year", correct_ant_bias=correct_ant_bias, reuse_ant_bias=reuse)
                   })
  l_gauss_fit <- lapply(l_gam_an, function(x){
                      fit_and_boot_allnat(x$l_gam_an_boot, gauss_fit)
                   })
  if(!redef_ano) z=1.6 else z=(1.6 - unlist(lapply(l_tas_nat, attr, which="bias")))
  l_far <- mapply(function(x, y, z){
                    boot_far_allnat_onperiod(x, y$l_gam_an_origin, l_time=1850:2100, xp=z)
                   },
                   x=l_gauss_fit,
                   y=l_gam_an,
                   z=z,
                   SIMPLIFY=FALSE
                   )
  boot_res_allnat_df <- lapply(l_far, function(a_far) add_param(a_far, operation=p1/p0, name="RR"))
  for (i in 1:length(boot_res_allnat_df))
    boot_res_allnat_df[[i]]["RR",,] %<>% imput_aurel(.) %T>% print()
  boot_res_allnat_far <- lapply(boot_res_allnat_df, function(bres) add_param(bres, operation=al_trans(RR), name="alRR") %>% 
                                add_param(., operation=el_trans(RR), name="elRR"))
  ic_allnat_far <- mapply(get_ic_onperiod,
                          b_onperiod=boot_res_allnat_far, 
                          method_name=names(boot_res_allnat_df),
                          ci_p=0.9,
                          SIMPLIFY=FALSE) %>% 
  do.call(rbind, .) %T>% print()
  ic_allnat_far$param <- as.character(ic_allnat_far$param)
  ic_allnat_far$param[ic_allnat_far$param == "gam_ant"] <- "x_t,ant" 
  ic_allnat_far$param[ic_allnat_far$param == "gam_nat"] <- "x_t,nat"
  ic_allnat_far$param[ic_allnat_far$param == "gam_all"] <- "x_t" 
  ic_allnat_far
}
test_010 <- test_methodo(l_tas_nat)
test_110 <- test_methodo(l_tas_nat, 1, 1, 0)
test_111 <- test_methodo(l_tas_nat, 1, 1, 1)
test_100 <- test_methodo(l_tas_nat, 1, 0, 0)
test_101 <- test_methodo(l_tas_nat, 1, 0, 1)
test_010$config <- "no center, no reuse"
test_110$config <- "center, no reuse"
test_111$config <- "center, reuse"
test_100$config <- "center, no correct"
test <- rbind(test_010, test_110, test_111, test_100)
rbind(test_00, test_11, test_11)
plot_far(test_010, axis_trans="al")
plot_far(test_110, axis_trans="al")
plot_far(test_111, axis_trans="al")
plot_far(test_100, axis_trans="al")
pdf(res_path("far_ic_app_test.pdf"))
plot_pannel_far(test, axis_trans="al")
dev.off()

ic_subset <- test[test$param %in% c("p0", "p1"), ] %T>% print()
p0001 <- ggplot(ic_subset,  aes(x=time))+
ggtitle(expression(paste(p[0], " and ", p[1], " from 1870 to 2070")))+ 
geom_point(aes(x=time,  y=Estim,  group=method,  color=method), size=1)+
geom_line(aes(x=time,  y=Estim,  group=method,  color=method), size=1)+ 
geom_ribbon(aes(x=time,  ymin=IC_inf,  ymax=IC_sup,  group=method,  fill=method),  alpha=0.2)+ 
coord_cartesian(xlim=c(1850,2100))+
facet_grid(param~config, scales="free_y") +
ylab("Estimate")+
ylim(c(0, 0.000002))+
theme(legend.position = "bottom")
pdf(file=res_path("p0p1_app_test_ylim.pdf"))
# pdf(file=res_path("p0p1_app_test.pdf"))
# pdf(file="p0p1_app_poster.pdf", height=15)
plot(p0001)
dev.off()

df_tas_pi <- rbind(cbind(l_tas_pi[[1]], method="cnrm"), cbind(l_tas_pi[[2]], method="ipsl"))
df_tas_nat <- rbind(cbind(l_tas_nat[[1]], method="cnrm"), cbind(l_tas_nat[[2]], method="ipsl"))
p0 <- ggplot(data=subset(test, param %in% c("x_t", "x_t,nat", "x_t,ant")), 
             aes(x=time))+
ggtitle(expression(paste("Covariate Estimation and Decomposition: ", x[t], " = ", x["t, ant"], " + ",  x["t, nat"]))) +
geom_point(data=df_tas_nat, aes(x=year, y=eur_tas), size=0.9, alpha=0.7) +
geom_point(data=df_tas_pi, aes(x=year, y=eur_tas), color="yellow",size=0.9, alpha=0.3) +
geom_line(aes(x=time, y=Estim, group=param, color=param), size=1) +
geom_ribbon(aes(x=time, ymin=IC_inf, ymax=IC_sup, group=param, color=param), alpha=0.4) + 
facet_grid(param~method+config) +
coord_cartesian(xlim=c(1850,2100))+
theme(legend.position="none")+
ylab("Estimate")
pdf(file=res_path("gam_decomp_app_test.pdf"), width=20)
plot(p0)
dev.off()


l_samples <- lapply(l_tas_pi, boot_samples)
l_samples <- lapply(l_tas_nat, boot_samples)
# l_samples <- lapply(l_samples, function(x){
#                       lapply(x, correct_anomalies_bias, y_name="eur_tas", time_name="year")
#                    })
l_gam_an <- lapply(l_samples, function(x){
                     print(x)
                     boot_gam_allnat(x, "eur_tas", "CO2", "nat", "year", reuse_ant_bias=TRUE)
                   })
df <- cbind(l_gam_an$cnrm$l_gam_an_origin[[1]],"clim_mod"="cnrm")
df <- rbind(df, cbind(l_gam_an$ipsl$l_gam_an_origin[[1]],"clim_mod"="ipsl"))
p_gg <- ggplot(df, aes(x=ant, y=y)) +
geom_point(alpha=0.3) +
geom_line(aes(y=gam_all), color="red", size=2)+
xlab("CO2 forcings") +
ylab("T_reg") +
ggtitle("GAM fit : T_reg=f(C02)")+
facet_grid(clim_mod~.)
p_gg

p_gg <- ggplot(df, aes(x=time, y=ant)) +
geom_point(alpha=0.3) +
geom_line(aes(y=ant), color="red", size=2)+
ylab("CO2 forcings") +
xlab("year") +
ggtitle("C02 forcings")
p_gg




# qthreshold=0.95
# l_gpd_fit <- fit_and_boot_allnat(l_gam_an$l_gam_an_boot, gpd_fit, qthreshold=qthreshold)

l_gauss_fit <- lapply(l_gam_an, function(x){
                      fit_and_boot_allnat(x$l_gam_an_boot, gauss_fit)
                   })
# pdf("qqplot_app.pdf", height=5, width=10)
# par(mfrow=c(1,2))
pdf(res_path("qqplot_app.pdf"), height=5, width=10)
source(tools_path("redefine_plot_fun.r"))
par(mfrow=c(1,2))
plot(l_gauss_fit[[1]][[1]], main="CNRM : Gauss -- QQ-Plot")
plot(l_gauss_fit[[2]][[1]], main="IPSL : Gauss -- QQ-Plot")
dev.off()

l_far <- mapply(function(x, y, z){
                  boot_far_allnat_onperiod(x, y$l_gam_an_origin, l_time=sort(unique(l_tas_nat$cnrm$year)), xp=z)
                   },
                   x=l_gauss_fit,
                   y=l_gam_an,
                   z=(1.6 - unlist(lapply(l_tas_pi, attr, which="bias"))),
                   # z=1.6,
                   SIMPLIFY=FALSE
                   )
l_far_file <- save_fast(l_far)

l_far <- readRDS("l_far.rds")
boot_res_allnat_df <- lapply(l_far, function(a_far) add_param(a_far, operation=p1/p0, name="RR"))
for (i in 1:length(boot_res_allnat_df))
  boot_res_allnat_df[[i]]["RR",,] %<>% imput_aurel(.) %T>% print()
boot_res_allnat_far <- lapply(boot_res_allnat_df, function(bres) add_param(bres, operation=al_trans(RR), name="alRR") %>% 
                              add_param(., operation=el_trans(RR), name="elRR"))
ic_allnat_far <- mapply(get_ic_onperiod,
                        b_onperiod=boot_res_allnat_far, 
                        method_name=names(boot_res_allnat_df),
                        ci_p=0.9,
                        SIMPLIFY=FALSE) %>% 
do.call(rbind, .) %T>% print()
ic_allnat_far$param <- as.character(ic_allnat_far$param)
ic_allnat_far$param[ic_allnat_far$param == "gam_ant"] <- "x_t,ant" 
ic_allnat_far$param[ic_allnat_far$param == "gam_nat"] <- "x_t,nat"
ic_allnat_far$param[ic_allnat_far$param == "gam_all"] <- "x_t" 
# ic_allnat_far <- get_ic_onperiod(b_onperiod=boot_res_allnat, method_name="gpd_fit")

l_q <- mapply(function(x, y){
                  boot_q_allnat_onperiod(x, y$l_gam_an_origin, l_time=1850:2100, p=0.01)
                   },
                   l_gauss_fit,
                   l_gam_an,
                   SIMPLIFY=FALSE
                   )
l_q_file <- save_fast(l_q)
l_q <- readRDS("l_q.rds") 
ic_allnat_q <- mapply(get_ic_onperiod,
                        b_onperiod=l_q, 
                        method_name=names(l_q),
                        ci_p=0.9,
                        SIMPLIFY=FALSE) %>% 
do.call(rbind, .) %T>% print()
ic_allnat_q$param <- as.character(ic_allnat_q$param)
ic_allnat_q$param[ic_allnat_q$param == "gam_ant"] <- "x_t,ant" 
ic_allnat_q$param[ic_allnat_q$param == "gam_nat"] <- "x_t,nat"
ic_allnat_q$param[ic_allnat_q$param == "gam_all"] <- "x_t" 
pdf(res_path("q99_app.pdf"))
plot_boot_time(ic_allnat_q, param="q", main="Quantile 99% from 1850 to 2100")
dev.off()

df_tas_nat <- rbind(cbind(l_tas_pi[[1]], method="cnrm"), cbind(l_tas_pi[[2]], method="ipsl"))
df_tas_nat <- rbind(cbind(l_tas_nat[[1]], method="cnrm"), cbind(l_tas_nat[[2]], method="ipsl"))
p0 <- ggplot(data=subset(ic_allnat_far, param %in% c("x_t", "x_t,nat", "x_t,ant")), 
             aes(x=time))+
ggtitle(expression(paste("Covariate Estimation and Decomposition: ", x[t], " = ", x["t, ant"], " + ",  x["t, nat"]))) +
geom_point(data=df_tas_nat, aes(x=year, y=eur_tas), size=0.9, alpha=0.7) +
geom_point(data=df_tas_pi, aes(x=year, y=eur_tas), color="yellow",size=0.9, alpha=0.3) +
geom_line(aes(x=time, y=Estim, group=param, color=param), size=1) +
geom_ribbon(aes(x=time, ymin=IC_inf, ymax=IC_sup, group=param, color=param), alpha=0.4) + 
facet_grid(param~method) +
coord_cartesian(xlim=c(1850,2100))+
theme(legend.position="none")+
ylab("Estimate")
pdf(file=res_path("gam_decomp_app_nat.pdf"))
plot(p0)
dev.off()

# ic_subset <- ic_allnat_far[ic_allnat_far$param %in% c("p0", "p1") & ic_allnat_far$method=="ipsl",] %T>% print()
ic_subset <- ic_allnat_far[ic_allnat_far$param %in% c("p0", "p1"), ] %T>% print()
p0001 <- ggplot(ic_subset,  aes(x=time))+
ggtitle(expression(paste(p[0], " and ", p[1], " from 1870 to 2070")))+ 
geom_point(aes(x=time,  y=Estim,  group=method,  color=method), size=1)+
geom_line(aes(x=time,  y=Estim,  group=method,  color=method), size=1)+ 
geom_ribbon(aes(x=time,  ymin=IC_inf,  ymax=IC_sup,  group=method,  fill=method),  alpha=0.2)+ 
coord_cartesian(xlim=c(1850,2100))+
# facet_grid(param~., scales="free_y") +
facet_grid(param~. ) +
ylab("Estimate")+
coord_trans(y="log10")+
# ylim(c(0, 0.000002))+
theme(legend.position = "bottom")
pdf(file=res_path("p0p1_app_nat.pdf"))
# pdf(file="p0p1_app_poster.pdf", height=15)
plot(p0001)
dev.off()

plot_far(ic_allnat_far, axis_trans="al", main=expression(paste("Relative Risk = ", p[1]/p[0], " from 1850 to 2100" )), xlim=c(1850,2100), col=c(gg_color_hue(2)))
dev.print(pdf, file=res_path("far_ic_app_nat.pdf"), height=6)
plot_far(subset(ic_allnat_far, method=="gauss_fit"), axis_trans="al", main=expression(paste("Relative Risk = ", p[1]/p[0], " from 1850 to 2100" )), xlim=c(1850,2100), col=c(gg_color_hue(2)))
dev.print(pdf, file="far_ic_app_poster.pdf", width=12, height=10)

