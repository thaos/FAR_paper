packages <- c("devtools", 
              "magrittr", 
              "mgcv",
              "reshape2",
              "ggplot2")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages,
                           rownames(installed.packages())),
                   repos="http://cran.us.r-project.org")
}
lapply(packages, library, character.only=TRUE)
# load_all("~/Data/FARpackage/FARg")

fit_rmnat <- function(y, ant, nat, time, fit, pre_ind=c(1850,1879), ant_pre_ind=NULL, init=NULL, knots=NULL, fixed=FALSE, ant_bias=NULL, ...){
  data=data.frame("time"=time, "y"=y, "ant"=ant, "nat"=nat)
  if(is.null(ant_pre_ind)) ant_pre_ind <- ant[time >= pre_ind[1] & time  <= pre_ind[2]]
  gam_fit <- gam_allnat(y, ant, nat, time, knots=knots, fixed=fixed)
  if(is.null(ant_bias)) ant_bias <- compute_ant_bias(gam_fit, ant_pre_ind)
  gam_an <- correct_ant_bias(gam_fit, bias=ant_bias)
  data=cbind(data, gam_an)
  data$y_rmnat <- y - data$nat
  data$gam_ant <- data$gam_all - data$gam_nat
  ans <- fit(y_rmnat, data=data, ~gam_ant, ~gam_ant, time_var="time", init=init, ...)
  ans$ant_bias <- ant_bias
  ans$gam_fit <- gam_fit
  ans$pre_ind <- pre_ind
  ans$ant_pre_ind <- ant_pre_ind
  ans$knots <- knots 
  ans$fixed <- fixed 
  ans
}

get_far_rmnat_onperiod <- function(object, xp, t0, l_time, original_data){
  p0 <- get_p_allnat_onperiod(object, xp, rep(t0, length(l_time)), original_data, ant_nat_or_all="ant") 
  p1 <- get_p_allnat_onperiod(object, xp, l_time, original_data, ant_nat_or_all="ant") 
  far <- 1 - p0[1,]/p1[1,]
  ans <- rbind(far, p0, p1)
  colnames(ans) <- l_time
  rownames(ans) <- c("FAR", paste(rownames(p0),"0", sep=""), paste(rownames(p1), "1", sep=""))
  ans
}

boot_far_rmnat_onperiod  <- function(l_fit, l_gam_an_origin, t0, l_time, xp){
  l_far <- mapply(get_far_rmnat_onperiod, object=l_fit, original_data=l_gam_an_origin, MoreArgs=list("xp"=xp, "t0"=t0, "l_time"=l_time), SIMPLIFY=FALSE)
  l_gam_an_origin <- lapply(l_gam_an_origin, function(x) unique(x[x$time %in% l_time, c("time", "gam_all", "gam_nat", "gam_ant")]))
  rbinding <- function(x, y){
    print(str(x))
    y <- t(y[,c("gam_all", "gam_nat", "gam_ant")]) 
    print(str(y))
    rbind(x, y)
  }
  l_far2 <- mapply(rbinding, l_far, l_gam_an_origin, SIMPLIFY=FALSE)
  a_far <- array(unlist(l_far2), dim=c(dim(l_far2[[1]]), length(l_far)))
  dimnames(a_far)[1:2] <- dimnames(l_far2[[1]])
  dimnames(a_far)[[3]] <- 1:length(l_far)
  a_far
}

boot_q_rmnat_onperiod  <- function(l_fit, l_gam_an_origin, t0, l_time, xp){
  l_far <- mapply(get_far_rmnat_onperiod, object=l_fit, original_data=l_gam_an_origin, MoreArgs=list("xp"=xp, "t0"=t0, "l_time"=l_time), SIMPLIFY=FALSE)
  l_gam_an_origin <- lapply(l_gam_an_origin, function(x) unique(x[x$time %in% l_time, c("time", "gam_all", "gam_nat", "gam_ant")]))
  rbinding <- function(x, y){
    print(str(x))
    y <- t(y[,c("gam_all", "gam_nat", "gam_ant")]) 
    print(str(y))
    rbind(x, y)
  }
  l_far2 <- mapply(rbinding, l_far, l_gam_an_origin, SIMPLIFY=FALSE)
  a_far <- array(unlist(l_far2), dim=c(dim(l_far2[[1]]), length(l_far)))
  dimnames(a_far)[1:2] <- dimnames(l_far2[[1]])
  dimnames(a_far)[[3]] <- 1:length(l_far)
  a_far
}


add_param <- function(b_onperiod, operation, name){
  ans <- array(dim=dim(b_onperiod)+c(1,0,0))
  ans[-dim(ans)[1],,] <- b_onperiod
  rownames_ans <- dimnames(b_onperiod)[[1]]
  dimnames(ans)[2:3] <- dimnames(b_onperiod)[2:3]
  b_onperiod %<>% melt(., varnames=c("param", "time", "bootsample")) %>%
  dcast(bootsample+time~param, data=.)
  new_param <- eval(substitute(operation), envir=b_onperiod)
  b_onperiod <- cbind(b_onperiod, new_param)
  rownames_ans %<>% append(., name)
  ans[dim(ans)[1],,] <- new_param
  dimnames(ans)[[1]] <- rownames_ans
  ans
}

get_ic_onperiod <- function(b_onperiod, method_name=NULL,  ci_p=0.95, ...){
  alpha <- 1-ci_p
  if (is.null(method_name)) method_name <- deparse(substitute(b_onperiod))
  b_onperiod %<>% melt(., varnames=c("param", "time", "bootsample"))%>% 
  cbind(., "method" = method_name) %>%
  #   tapply(b_onperiod$value, INDEX = subset(b_onperiod, select=c("param", "time", "method")), FUN=quantile, probs=c(alpha/2, 0.5, 1-alpha/2))
  aggregate(value~param+time+method, data=., FUN=quantile, probs=c(alpha/2, 0.5, 1-alpha/2), ...)
  ic <- b_onperiod$value
  colnames(ic) <- c("IC_inf", "Estim", "IC_sup")
  b_onperiod$value <- NULL
  cbind(b_onperiod, ic)
}

plot_pannel_boot_time <- function(ic_df, param="FAR", main="FAR(t1)"){ 
  env <- environment()
  ic_subset <- ic_df[ic_df$param == param, ] %T>% print()
  ic_subset %>% ggplot(.,  aes(x=time))+
  ggtitle(main)+
  geom_point(aes(x=time,  y=Estim,  group=method,  color=method), size=1)+
  geom_line(aes(x=time,  y=Estim,  group=method,  color=method), size=1)+ 
  geom_ribbon(aes(x=time,  ymin=IC_inf,  ymax=IC_sup,  group=method,  fill=method),  alpha=0.2)+ 
  facet_grid(method~config, scales="free_y") +
  # coord_cartesian( xlim = xlim, ylim = ylim)+
  theme(legend.position = "bottom")
}

plot_boot_time <- function(ic_df, param="FAR", main="FAR(t1)"){ 
  env <- environment()
  ic_subset <- ic_df[ic_df$param == param, ] %T>% print()
  ic_subset %>% ggplot(.,  aes(x=time))+
  ggtitle(main)+
  geom_point(aes(x=time,  y=Estim,  group=method,  color=method), size=1)+
  geom_line(aes(x=time,  y=Estim,  group=method,  color=method), size=1)+ 
  geom_ribbon(aes(x=time,  ymin=IC_inf,  ymax=IC_sup,  group=method,  fill=method),  alpha=0.2)+ 
  # coord_cartesian( xlim = xlim, ylim = ylim)+
  theme(legend.position = "bottom")
}

gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

source("ggplot_dual_axis.r")
plot_far <- function(ic_mat, axis_trans, main="", xlim, ylim, col=NULL){
  ticks=c(0, 1/100, 1/10, 1/5, 1/3, 1/2, 1, 2, 3, 5, 10, 100, Inf)
  if(missing(xlim)) xlim <- range(ic_mat$time)
  if(missing(ylim)) ylim <- c(-1.1, 1.1)
  inv <- get(paste(axis_trans, "_inv", sep=""))
  trans <- get(paste(axis_trans, "_trans", sep=""))
  breaks <- trans(ticks) 
  param <- paste(axis_trans,"RR", sep="")
  p1 <- plot_boot_time(ic_mat, param=param, main=main) +
  scale_y_continuous(name = "Relative Risk",
                     breaks = breaks,
                     labels = trans_format(inv, function(x) format(x, digits=2)))+
  coord_cartesian(xlim=xlim, ylim=ylim)
  if (!is.null(col)) 
    p1 <- p1 + scale_color_manual(values=col)
  p2 <- p1 + scale_y_continuous(name="FAR",
                                breaks = breaks,
                                labels = trans_format(inv, function(x) format(rrtofar(x), digits=3)))
  ggplot_dual_axis(p1,p2)
}

plot_pannel_far <- function(ic_mat, axis_trans, main="", xlim, ylim, col=NULL){
  ticks=c(0, 1/100, 1/10, 1/5, 1/3, 1/2, 1, 2, 3, 5, 10, 100, Inf)
  if(missing(xlim)) xlim <- range(ic_mat$time)
  if(missing(ylim)) ylim <- c(-1.1, 1.1)
  inv <- get(paste(axis_trans, "_inv", sep=""))
  trans <- get(paste(axis_trans, "_trans", sep=""))
  breaks <- trans(ticks) 
  param <- paste(axis_trans,"RR", sep="")
  p1 <- plot_pannel_boot_time(ic_mat, param=param, main=main) +
  scale_y_continuous(name = "Relative Risk",
                     breaks = breaks,
                     labels = trans_format(inv, function(x) format(x, digits=2)))+
  coord_cartesian(xlim=xlim, ylim=ylim)
  if (!is.null(col)) 
    p1 <- p1 + scale_color_manual(values=col)
  p1
}

al_trans <- function(x) atan(log(x))/(pi/2)
al_inv <- function(x) exp(tan((pi/2)*x))
el_trans <- function(x) ifelse(x <= 1, x-1, 1-1/x)
el_inv <- function(x) ifelse(x > 0, 1/(1-x), x + 1)
rrtofar <- function(x) 1 - 1/x
