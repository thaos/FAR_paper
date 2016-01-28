# load functions
packages <- c("scales",
              "gtable",
              "grid",
              "extRemes",
              "gridExtra",
              "devtools",
              "quantreg")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages,
                           rownames(installed.packages())),
                   repos="http://cran.us.r-project.org")
}
lapply(packages, library, character.only=TRUE)
devtools::load_all("~/Data/FARpackage/FARg")

setwd("~/Data/FAR_paper/2_Analysis/3_Scripts/")
make_shortcut <- function(folder){
function(file) paste(folder, file, sep="")
}
res_path <- make_shortcut("../4_Results/")
data_path <- make_shortcut("../1_Data/")
tools_path <- make_shortcut("../2_Tools/")

source(tools_path("far_gam_coverage.r"), chdir=TRUE)
source(tools_path("gam_allnat.r"), chdir=TRUE)
source(tools_path("gam_rmnat.r"), chdir=TRUE)

set.seed(1)
data <- simulate_data()

pdf(file=res_path("Ytimeseries.pdf"),height=10)
# png(file="Ytimeseries.png", height=2*480)
all_p <- ggplot(data=data, aes(x=time, y=y)) +
geom_point(alpha=0.3) +
geom_line(aes(y=ANT+NAT), size=2) + 
ggtitle(expression(paste("time series of ", y[t], " and of the covariate ", x[t]))) +
ylab(expression(y[t]))
ant_p <- ggplot(data=data, aes(x=time, y=y)) +
geom_line(aes(y=ANT), size=2) + 
ggtitle(expression(paste("time series of ", x["t, ant"])))+
ylab(expression(x["t, ant"]))
nat_p <- ggplot(data=data, aes(x=time, y=y)) +
geom_line(aes(y=NAT), size=2) + 
ggtitle(expression(paste("time series of ", x["t, nat"])))+
ylab(expression(x["t, nat"]))
err_p <- ggplot(data=data, aes(x=time, y=y)) +
geom_point(aes(y=y-NAT-ANT), alpha=0.3) +
ggtitle(expression(paste("time series of ", epsilon[t]))) +
ylab(expression(epsilon[t]))
grid.arrange(all_p, ant_p, nat_p, err_p, nrow=4)
dev.off()

l_samples <- boot_samples(data, R=250)
# Ne pas utiliser sur le cas idéaliser car on travail en absolue pas en anomalie
# l_samples <- lapply(l_samples, correct_anomalies_bias, period=c(1850,1879), y_name="y", time_name="time")
l_gam_an <- boot_gam_allnat(l_samples, "y", "time", "natG", "time", reuse_ant_bias=FALSE)
l_gpd_fit <- fit_and_boot_allnat(l_gam_an$l_gam_an_boot, gpd_fit, qthreshold=0.95)
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
                        ci_p=0.75,
                        SIMPLIFY=FALSE) %>% 
do.call(rbind, .) %T>% print()

data_t <- unique(subset(data, select=c("time", "NAT", "ANT")))
gam_all_t <- with(data_t, NAT+ANT)[1:251]
# gam_all_t2 <- predict(gam(y~s(time)+natG, data=data))[1:251]
gam_nat_t <- with(data_t, NAT)[1:251]
gam_ant_t <- with(data_t, ANT)[1:251]
xp=1.6
v_sd <- 0.5^2 + 0.02 * (gam_all_t)
v_sd_nat <- 0.5^2 + 0.02 * (gam_nat_t)
p_nat_true <- mapply(pnorm, mean=gam_nat_t, sd=v_sd_nat, MoreArgs=list(lower.tail=FALSE, q=xp))
p_all_true <- mapply(pnorm, mean=gam_all_t, sd=v_sd, MoreArgs=list(lower.tail=FALSE, q=xp))
rr_true <- p_all_true/p_nat_true
alrr_true <- al_trans(rr_true)
format_gam <- function(gam_xxx_t, name){
  ans <- as.data.frame(matrix(gam_xxx_t, nrow=length(gam_xxx_t), ncol=3))
  ans <- cbind("truth", ans)
  ans <- cbind(1850:2100, ans)
  ans <- cbind(name, ans)
  names(ans) <- c("param", "time", "method", "IC_inf", "Estim", "IC_sup")
  ans
}
gam_formated <- mapply(format_gam, 
                       gam_xxx_t=list(gam_all_t, gam_nat_t, gam_ant_t, rr_true, alrr_true),
                       name=c("x_t", "x_t,nat", "x_t,ant", "RR", "alRR"), SIMPLIFY=FALSE)
gam_formated <- do.call(rbind, gam_formated)
ic_allnat <- read.table(data_path("coverage_alpha90_linsd/coverage_n1.txt"))
# ic_allnat <- ic_allnat_try 
ic_allnat$method %<>% lapply(., function(x) substring(x,3,nchar(x))) %>% unlist()
ic_allnat %<>% rbind(., gam_formated) 
ic_allnat$param[ic_allnat$param == "gam_ant"] <- "x_t,ant" 
ic_allnat$param[ic_allnat$param == "gam_nat"] <- "x_t,nat"
ic_allnat$param[ic_allnat$param == "gam_all"] <- "x_t" 

p_gam <- ggplot(data=subset(ic_allnat, param %in% c("x_t", "x_t,nat", "x_t,ant") & method == "gpd_fit"), 
             aes(x=time))+
ggtitle(expression(paste("Covariate Estimation and Decomposition: ", x[t], " = ", x["t, ant"], " + ",  x["t, nat"]))) +
geom_line(aes(x=time, y=Estim, group=param, color=param), size=1) +
geom_line(data=subset(ic_allnat, param %in% c("x_t", "x_t,nat", "x_t,ant") & method == "truth"), aes(x=time, y=Estim, group=param), color="black", size=1) +
geom_ribbon(aes(x=time, ymin=IC_inf, ymax=IC_sup, group=param, fill=param, color=param), alpha=0.4, size=0.5) + 
facet_grid(param~., scales="free_y") +
theme(legend.position="none")+
ylab("Estimate")
pdf(file=res_path("gam_decomp.pdf"))
print(p_gam)
dev.off()

# plot_far(subset(ic_allnat, method != "gpd_fit"), axis_trans="al", main="Relative Risk = P_all / P_nat  from 1850 to  2100 \n atan log y-axis", xlim=c(1850,2100), col=c(gg_color_hue(2)[1], "black"))
# plot_far(subset(ic_allnat, method != "gauss_fit"), axis_trans="al", main="Relative Risk = P_all / P_nat  from 1850 to  2100 \n atan log y-axis", xlim=c(1850,2100), col=c(gg_color_hue(2)[2], "black"))
plot_far(ic_allnat, axis_trans="al", main=expression(paste("Relative Risk = ", p[1]/p[0], " from 1850 to 2100" )), xlim=c(1850,2100), col=c(gg_color_hue(2), "black"))
dev.print(pdf, file=res_path("far_ic.pdf"), height=6)


# fname <- paste("coverage_q95_reuse/coverage_n",1:500,".txt",sep="")
fname <- data_path(paste("coverage_alpha90_linsd/coverage_n",1:500,".txt",sep=""))
fname <- fname[file.exists(fname)]
ic_df <- sapply(fname, function(x){
                  read.table(x) %>% inIC()
             })
inIC_df <- apply(ic_df, 1, function(x) list(simplify2array(x)))
inIC_df <- lapply(inIC_df, `[[`, 1)
cover_df <- lapply(inIC_df, function(x) apply(x,1,mean))
cover_df %<>% as.data.frame() %>% 
cbind(., "time"= 1850:2100) %>% 
melt(., id.vars=c("time"), variable.name="method", value.name="coverage")
cover_df$method %<>% lapply(., function(x) substring(x,3,nchar(as.character(x)))) %>% unlist()

p_cover <- ggplot(data=cover_df, 
             aes(x=time))+
ggtitle(expression(paste("Coverage Probability of the ", CI["90%"](FAR)))) +
geom_line(aes(x=time, y=coverage, group=method, color=method), size=1) +
facet_grid(method~., scales="free_y") +
scale_y_continuous(breaks=seq(0,1,0.1), limits=c(0,1))+
theme(legend.position="none")#+
# ylim(c(0,1))
p_cover
dev.print(pdf, file=res_path("ic_cover_y01.pdf"))
# dev.print(pdf, file="ic_cover_reuse.pdf" )
# dev.print(pdf, file="ic_cover_const.pdf" )

qthreshold=0.95
l_gpd_fit <- fit_and_boot_allnat(l_gam_an$l_gam_an_boot, gpd_fit, qthreshold=qthreshold)
l_gauss_fit <- fit_and_boot_allnat(l_gam_an$l_gam_an_boot, gauss_fit)
l_fit <- list("l_gpd_fit"=l_gpd_fit, "l_gauss_fit"=l_gauss_fit)
l_far <- lapply(l_fit, function(fit) boot_far_allnat_onperiod(fit, l_gam_an$l_gam_an_origin, 1850:2100, 1.6))
boot_res_allnat_df <- lapply(l_far, function(a_far) add_param(a_far, operation=p1/p0, name="RR"))
compute_UB <- function(threshold, sigma, shape){
  ifelse(shape < 0, threshold-sigma/shape, Inf)
}

compute_UB <- function(threshold, sigma, shape){
  ifelse(shape < 0, threshold-sigma/shape, Inf)
}
boot_res_allnat_df[[1]]%>%
add_param(., operation=compute_UB(threshold0, sigma0, shape0), name="UB_nat") %>%
get_ic_onperiod(b_onperiod=., method_name="gpd_natnat", ci_p=1/3) %>%
plot_boot_time(., param="UB_nat", main="UB_nat from 1870 to 2070") %>% 
print()
