## ID: plot.R, last updated 2025-03-17, F.Osorio

envelope.Laplace <- function(object, reps = 50, conf = 0.95, plot.it = TRUE)
{ ## simulated envelope
  envel <- function(n, center, Scatter, reps, conf) {
    conf <- 1 - conf
    # initialize progress bar
    cat(" Progress:\n")
    pb <- txtProgressBar(min = 0, max = reps, style = 3)
    elims <- matrix(0, nrow = n, ncol = reps)
    for (i in 1:reps) {
      x <- rmLaplace(n, center = center, Scatter = Scatter)
      fit <- LaplaceFit(x)
      z <- WH.Laplace(x, center = fit$center, Scatter = fit$Scatter)
      elims[,i] <- sort(z)
      # update progress bar
      setTxtProgressBar(pb, i)
    }
    close(pb)
    band <- matrix(0, nrow = n, ncol = 2)
    for (i in 1:n)
      band[i,] <- quantile(elims[i,], probs = c(conf / 2, 1 - conf / 2))
    band
  }

  n <- object$dims[1]
  z <- WH.Laplace(object$x, center = object$center, Scatter = object$Scatter)

  if (plot.it) {
    band  <- envel(n, object$center, object$Scatter, reps, conf)
    ylim <- range(z, band)
    qqnorm(z, ylim = ylim, main = "Transformed distances Q-Q plot")
    par(new = TRUE)
    qqnorm(band[,1], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
    par(new = TRUE)
    qqnorm(band[,2], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
  }
  invisible(list(transformed = z, envelope = band))
}
