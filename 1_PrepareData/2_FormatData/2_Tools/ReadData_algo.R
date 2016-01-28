library(ncdf)

get_data  <- function(l.nc, origin="1850-01-01"){
	require(date)
	getTAS <- function(nc.file){
		nc <- open.ncdf(nc.file)
		get.var.ncdf(nc, "tas")
	}
	nc <- open.ncdf(l.nc[1])
	years  <-  Vectorize(function(x)substr(x, 1,4))(as.Date(get.var.ncdf(nc, "time"),  origin="1850-01-01"))
	years  <-  as.numeric(years)
	close.ncdf(nc)
	tas <- sapply(l.nc, getTAS)
	colnames(tas) <- paste("run", 1:length(l.nc), sep="")
	rownames(tas) <- years
	tas
}

create_rds <- function(model="CNRM", file=paste("tas_", model, ".rds", sep="")){
  folder <- list.files("../1_Data/", pattern <- model, full.names=TRUE) 
  l.gbl <- list.files(folder, pattern="gbl_sts", full.names=TRUE)
  l.eur <- list.files(folder, pattern="eur_sts", full.names=TRUE)
  d.gbl <- get_data(l.gbl)
  d.gbl.v <- as.numeric(d.gbl)
  d.eur <- get_data(l.eur)
  d.eur.v <- as.numeric(d.eur)
  year <- as.numeric(rep(rownames(d.eur), ncol(d.eur)))
  tas  <- as.data.frame(cbind(year,  d.eur.v,  d.gbl.v))
  names(tas) <- c("year", "eur_tas", "gbl_tas") 
  file <- paste("../4_Results/", file, sep="")
  saveRDS(tas, file=file)
  file
}

