# args <- commandArgs(TRUE)
# print(args)
# .libPaths("/cnrm/vdr/USERS/thaos/R/x86_64-mageia-linux-gnu-library/3.0")
# 
# load functions
packages <- c("scales",
              "gtable",
              "grid")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) { install.packages(setdiff(packages,
                           rownames(installed.packages())),
                   repos="http://cran.us.r-project.org")
}
lapply(packages, library, character.only=TRUE)

source("gam_allnat.r")
source("gam_rmnat.r")
#tas_nat <- readRDS("tas_nat.rds")

s_fun <- function(time){
  (time >= 1950) * (6/150^3) * (time-1950)^3
}

simulate_data <- function(){
  tas_nat <- readRDS("../1_Data/tas_nat_cnrm.rds")
  time <- tas_nat$year
  ANT <- s_fun(time)
  natG <- tas_nat$nat
  NAT <- natG * 0.5
  v_sd <- 0.5 + 0.02 * (ANT+NAT)
  e <- unlist(lapply(v_sd, rnorm, n=1, mean=0))
  y <- ANT + NAT + e
  data.frame("time"=time, "natG"= natG, "y"=y, "ANT"=ANT, "NAT"=NAT)
}
  

compute_ic <- function(i_rep, qthreshold=0.96){
  set.seed(i_rep)
  data <- simulate_data()
  l_samples <- boot_samples(data, R=250)
  l_samples <- lapply(l_samples, correct_anomalies_bias, y_name="y", time_name="time")
  l_gam_an <- boot_gam_allnat(l_samples, "y", "time", "natG", "time", reuse_ant_bias=FALSE)
  l_gpd_fit <- fit_and_boot_allnat(l_gam_an$l_gam_an_boot, gpd_fit, qthreshold=qthreshold)
  l_gauss_fit <- fit_and_boot_allnat(l_gam_an$l_gam_an_boot, gauss_fit)
  l_fit <- list("l_gpd_fit"=l_gpd_fit, "l_gauss_fit"=l_gauss_fit)
  l_far <- lapply(l_fit, function(fit) boot_far_allnat_onperiod(fit, l_gam_an$l_gam_an_origin, 1850:2100, 1.6))
  boot_res_allnat_df <- lapply(l_far, function(a_far) add_param(a_far, operation=p1/p0, name="RR"))
  for (i in 1:length(boot_res_allnat_df))
    boot_res_allnat_df[[i]]["RR",,] %<>% imput_aurel(.) %T>% print()
  boot_res_allnat_try <- lapply(boot_res_allnat_df, function(bres) add_param(bres, operation=al_trans(RR), name="alRR") %>% 
                                add_param(., operation=el_trans(RR), name="elRR"))
  ic_allnat_try <- mapply(get_ic_onperiod,
                          b_onperiod=boot_res_allnat_try, 
                          method_name=names(boot_res_allnat_df),
                          SIMPLIFY=FALSE) %>% 
  do.call(rbind, .) %T>% print()
  # ic_allnat_try <- get_ic_onperiod(b_onperiod=boot_res_allnat, method_name="gpd_fit")
  plot_far(ic_allnat_try, axis_trans="al", main="RR = p_all / p_nat  from 1870 to  2070 \n atan log y-axis", xlim=c(1850,2100))
  file=paste("coverage_n", i_rep, ".txt", sep="")
  write.table(ic_allnat_try, file=file)
  file 
}
# compute_ic(i_rep=args[1])
# compute_ic(i_rep=84)

#a=read.table("coverage_n1.txt") 
inIC <- function(ic_df){
  tas_nat <- readRDS("../1_Data/tas_nat_cnrm.rds")
  time <- sort(unique(tas_nat$year))
  ANT <- s_fun(time)
  natG <- unique(subset(tas_nat, select=c("year", "nat")))$nat
  NAT <- natG * 0.5
  xp=1.6
  v_sd <- (0.5 + 0.02 * (ANT+NAT))
  v_sd_nat <- (0.5 + 0.02 * (NAT))
  p_nat_true <- mapply(pnorm, mean=NAT, sd=v_sd_nat, MoreArgs=list(lower.tail=FALSE, q=xp))
  p_all_true <- mapply(pnorm, mean=NAT+ANT, sd=v_sd, MoreArgs=list(lower.tail=FALSE, q=xp))
  rr_true <- p_all_true/p_nat_true
  l_ic_df <- split.data.frame(ic_df, f=ic_df$method)
  l_ic_rr <- lapply(l_ic_df, function(x) x[x$param=="RR", ])
  covered <- function(rr_true, ic_inf, ic_sup){
    ans <- (rr_true >= ic_inf & rr_true <= ic_sup)
    i_na <- which(is.na(ans) & ic_inf == 0 & is.infinite(ic_sup))
    ans[i_na] <- TRUE
    ans
  }
  l_cover <- lapply(l_ic_rr, function(x) covered(rr_true, x$IC_inf, x$IC_sup)) 
}
# b <- inIC(a)
# 
# fname <- paste("coverage/coverage_n",1:100,".txt",sep="")
# fname <- fname[file.exists(fname)]
# ic_df <- sapply(fname, function(x){
#                   read.table(x) %>% inIC()
#                           })
# inIC_df <- apply(ic_df, 1, function(x) list(simplify2array(x)))
# inIC_df <- lapply(inIC_df, `[[`, 1)
# cover_df <- lapply(inIC_df, function(x) apply(x,1,mean))
# lapply(cover_df, summary)
