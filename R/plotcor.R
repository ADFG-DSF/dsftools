#' Plot a correlation matrix
#' @description Plots a correlation matrix as shades of red (positive) or blue
#' (negative), or other colors as specified.
#'
#' The plotting component of this function is taken directly from `jagshelper::plotcor_jags()`.
#'
#' If row/column names are given with indices in square brackets
#' (e.g. `"a[1]", "a[2]"`, etc.), a single axis
#' tick will be used for all elements with a single name preceding the `[` character.
#' This was originally implemented in the `jagshelper` package, as applied to a
#' correlation matrix calculated from a data.frame
#' of MCMC output, and is left here in case it is useful.  This has the effect of
#' giving greater visual weight to single parameters, and reducing plot clutter.
#'
#' Values of correlation are overlayed, with
#' character size scaled according to the absolute correlation.
#' @param x Correlation matrix returned from `\link{cor()}`
#' @param mincor Minimum (absolute) correlation to use for text labels.  Defaults to 0 (all will be plotted)
#' @param maxn Maximum number of nodes per parameter name for text labels, to prevent plot clutter.  Defaults to 4.
#' @param maxcex Maximum character expansion factor for text labels.  Defaults to 1.
#' @param legend Whether to produce a plot legend.  Defaults to `TRUE`.
#' @param colors A vector of two colors defining the color ramps to use for negative
#' and positive colors.  Defaults to blue (negative) and red (positive).
#' @param ... Optional plotting arguments
#' @return `NULL`
#' @author Matt Tyers
#' @examples
#' ## simulating a dataset
#' n <- 100
#' xx <- rnorm(n, mean=10, sd=1)
#' dfsim <- data.frame(a = xx + rnorm(n, sd=0.5),
#'                     b = xx + rnorm(n, sd=0.5),
#'                     c = xx + rnorm(n, sd=2),
#'                     d = xx + rnorm(n, sd=4),
#'                     e = -xx + rnorm(n, sd=2),
#'                     f = -xx + rnorm(n, sd=0.5),
#'                     g = -xx + rnorm(n, sd=0.5))
#'
#' ## printing the correlation matrix
#' cor(dfsim)
#'
#' ## plotting the correlation matrix
#' plotcor(cor(dfsim))
#'
#' ## using a different color ramp
#' plotcor(cor(dfsim), colors=c("cornflowerblue","orange"))
#'
#'
#' ## heck we can plot the correlation matrix for our simulated ASL dataset
#' asldata <- sim_data$data
#' asldata$age <- as.numeric(as.factor(asldata$age))
#' plotcor(cor(asldata))  # data had NA values, leaving these in plot
#' plotcor(cor(asldata, use="na.or.complete"))  # data had NA values, removing these
#'
#' @importFrom graphics axis rect segments text
#' @export
plotcor <- function(x, mincor=0, maxn=4, maxcex=1, legend=TRUE, colors=c("#2297E6","#DF536B"), ...) {
  if(!inherits(x, "matrix")) stop("Input must be a correlation matrix returned from cor()")
  if(nrow(x) != ncol(x)) stop("Input must be a correlation matrix returned from cor()")
  if(length(colors) != 2) stop("Argument colors= must have two elements")

  dfcor <- x  # renaming this variable so the jagshelper stuff works!
  dfnames <- dimnames(dfcor)[[1]] #names(df)
  dfwhich <- sapply(strsplit(dfnames,split="[",fixed=T),FUN="[",1)
  dfhowmany <- rep(NA,length(dfnames))
  for(i in 1:length(dfnames)) dfhowmany[i] <- sum(dfwhich==dfwhich[i])
  dfdim <- cumsum(1/dfhowmany)
  dfdim1 <- c(0, dfdim[-length(dfdim)])

  xmat <- matrix(dfdim, nrow=length(dfdim), ncol=length(dfdim))
  ymat <- matrix(dfdim, nrow=length(dfdim), ncol=length(dfdim), byrow=T)
  xmat1 <- matrix(dfdim1, nrow=length(dfdim), ncol=length(dfdim))
  ymat1 <- matrix(dfdim1, nrow=length(dfdim), ncol=length(dfdim), byrow=T)

  cols <- 0*xmat
  for(i in 1:nrow(xmat)) {
    for(j in 1:i) {
      if(!is.na(dfcor[i,j])) {
        if(dfcor[i,j] > 0) {
          cols[i,j] <- cols[j,i] <- adjustcolor(colors[2],alpha.f=dfcor[i,j])
        }
        if(dfcor[i,j] < 0) {
          cols[i,j] <- cols[j,i] <- adjustcolor(colors[1],alpha.f=-dfcor[i,j])
        }
      }
      if(i==j) cols[i,j] <- 1
    }
  }

  plot(NA, xlim=(1+.1*legend)*range(0,dfdim), ylim=rev(range(0,dfdim)),
       yaxt="n", xaxt="n", ylab="", xlab="", bty='n',...=...)#,yaxs="i", xaxs="i"
  dfwhichunique <- unique(dfwhich)
  axis(side=1, at=1:length(dfwhichunique)-.5, labels=dfwhichunique, las=2)
  axis(side=2, at=1:length(dfwhichunique)-.5, labels=dfwhichunique, las=2)

  rect(xleft=xmat1, xright=xmat, ybottom=ymat1, ytop=ymat, border=cols, col=cols)
  dfcor4cex <- dfcor
  dfcor4cex[is.na(dfcor4cex)] <- 0.5
  for(i in 1:nrow(xmat)) {
    for(j in 1:i) {
      # if((dfhowmany[i]<=maxn) & (dfhowmany[j]<=maxn) & (abs(dfcor[i,j])>=mincor)) {
      if((dfhowmany[i]<=maxn) & (dfhowmany[j]<=maxn) & (abs(dfcor4cex[i,j])>=mincor)) {
        # text(x=dfdim1[i]+0.5/dfhowmany[i], y=dfdim1[j]+0.5/dfhowmany[j], labels=round(dfcor[i,j],2), cex=maxcex*abs(dfcor[i,j])^.3)
        # text(x=dfdim1[j]+0.5/dfhowmany[j], y=dfdim1[i]+0.5/dfhowmany[i], labels=round(dfcor[i,j],2), cex=maxcex*abs(dfcor[i,j])^.3)
        text(x=dfdim1[i]+0.5/dfhowmany[i], y=dfdim1[j]+0.5/dfhowmany[j], labels=round(dfcor[i,j],2), cex=maxcex*abs(dfcor4cex[i,j])^.3)
        text(x=dfdim1[j]+0.5/dfhowmany[j], y=dfdim1[i]+0.5/dfhowmany[i], labels=round(dfcor[i,j],2), cex=maxcex*abs(dfcor4cex[i,j])^.3)
        if(is.na(dfcor[i,j])) {
          text(x=dfdim1[i]+0.5/dfhowmany[i], y=dfdim1[j]+0.5/dfhowmany[j], labels="NA", cex=maxcex*abs(dfcor4cex[i,j])^.3)
          text(x=dfdim1[j]+0.5/dfhowmany[j], y=dfdim1[i]+0.5/dfhowmany[i], labels="NA", cex=maxcex*abs(dfcor4cex[i,j])^.3)
        }
      }
    }
  }
  # abline(v=0:length(dfwhichunique))
  segments(x0=rep(0, length(dfwhichunique)+1),
           x1=rep(length(dfwhichunique), length(dfwhichunique)+1),
           y0=0:length(dfwhichunique))
  segments(y0=rep(0, length(dfwhichunique)+1),
           y1=rep(length(dfwhichunique), length(dfwhichunique)+1),
           x0=0:length(dfwhichunique))

  if(legend) {
    legendby <- .25
    legendn <- 2/legendby+1

    legendl <- rep(1.05*max(dfdim),legendn)
    legendr <- rep(1.1*max(dfdim),legendn)
    legendb <- seq(from=.5*max(dfdim), to=0, length.out=legendn+1)[-legendn-1]
    legendt <- seq(from=.5*max(dfdim), to=0, length.out=legendn+1)[-1]

    legendcols <- rep(0, legendn)
    legendcors <- seq(-1,1,length.out=legendn)
    for(i in 1:(1/legendby)) {
      legendcols[i] <- adjustcolor(colors[1], alpha.f=-legendcors[i])
      legendcols[legendn+1-i] <- adjustcolor(colors[2], alpha.f=-legendcors[i])
    }

    rect(xleft=legendl, xright=legendr, ytop=legendt, ybottom=legendb, col=legendcols, border=NA)
    text(x=.5*(legendl+legendr), y=.5*(legendt+legendb), labels=legendcors, cex=.7)
  }
}

