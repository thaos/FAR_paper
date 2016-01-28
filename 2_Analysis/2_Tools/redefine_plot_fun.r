plot.gpd_fit <- function(x, ...){
  source("fitdiagnostics.R", local=TRUE)
  res <- fevd(x$y, x$data, threshold=predict(x$rq_fitted), scale.fun=x$sig_mod, type="GP", method="MLE", initial=x$results$par)
  print(all.equal(x$par,res$results$par))
  # plot(res, type="probprob", main = "Probability Plot")
  plot(res, type="qq2", xlab="(y > threshold) Empirical Quantiles", main="QQ-Plot")
  title(main="GPD -- QQ-plot")
}

plot.gauss_fit <- function(x, main="", ...){
  qqnorm <- function(y, pch = 20,
             xlab = "Standard Normal Quantiles",
             ylab = "Sample Quantiles", make.plot=TRUE, ...)
    {
      args <- list(...)
      y <- sort(na.omit(y))
      n <- length(y)
      p <- (1:length(y) - .5)/length(y)
      k <- .895 / (sqrt(n) * (1 - .01 / sqrt(n) + .85 / n))
      l <- suppressWarnings(qnorm(p - k))
      q <- qnorm(p)
      u <- suppressWarnings(qnorm(p + k))
      if(make.plot) {
        if(is.null(args$xlim)) plot(y, q, xlim = range(l, q, u, na.rm = TRUE), xlab = xlab, ylab = ylab, pch = pch, ...)
        else plot(y, q, xlab = xlab, ylab = ylab, pch = pch, ...)
        lines(y, l, lty = 2, col = "darkgray")
        lines(y, u, lty = 2, col = "darkgray")
      }
      out <- data.frame(lower=l, upper=u, qnorm=q, data=y)
      invisible(out)
    }
	param <- compute_par.gauss_fit(x, x$data)
	mu <- param[, 1]
	sig2 <- param[, 2]^2
	res <- x$y - mu
	res_std <- sort(res / sqrt(sig2))
  p_emp <- pnorm(res_std)
  p_theo <- seq_along(p_emp)/length(p_emp)
  # plot(p_emp, p_theo, xlab="Residual Empirical Probabilities", ylab="Residual Model Probabilities", main="Probability Plot")
	#abline(a=0, b=1, lty=2)
	qqn <- qqnorm(res_std, ylab="Quantiles from Model Simulated Data", xlab="Empical Quantiles", main=main)
  qqfit <- lm(qnorm~data, data=qqn)
  abline(a=coef(qqfit)[1], b=coef(qqfit)[2], col="darkgray")
	abline(a=0, b=1, col="orange", lty=2)
	legend("topleft", legend=c("1-1 line", "regression line", "95% confidence bands"), lty=c(2,1,2), col= c("orange", "darkgray", "darkgray"), lwd=1.2, bty="n")
}

