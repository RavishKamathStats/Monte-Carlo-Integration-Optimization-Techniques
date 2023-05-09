# Libraries ---------------------------------------------------------------
library(Deriv); library(rootSolve); library(rgl); library(yarrr)
library(insight); library(tidyverse)
# Plotting Function -------------------------------------------------------
plotGraph   <- function(func, func2 = NULL, interval, 
                          RV = NULL,radians = F, hist = F, plotTitle = NULL){
  xlabel      <- as.expression(body(func)); 
  ylabel      <- expression(h(x));
  if (radians == T){
    plot(func, from = interval[1], to = interval[2], axes = F, 
         xlab = xlabel, ylab = ylabel, main = plotTitle)
    if(interval[2] == pi){
      vals    <- seq(0, interval[2], length.out = 5)
      labels  <- expression(0, pi/4, pi/2, 3*pi/4,pi)
    }
    else if(interval[2] == 2*pi){
      vals    <- seq(0, interval[2], length.out = 5)
      labels  <- expression(0, pi/2, pi, 3*pi/2, 2*pi)
    }
    else{
      vals    <- seq(0, interval[2], length.out = 5)
      labels  <- expression(0, pi, 2*pi, 3*pi, 4*pi)
    }
    axis(side = 1, at = vals, labels = labels)
    axis(side = 2)
    abline(h = 0)
  }
  else if (hist == T) {
    xlabel    <- as.expression(body(func2)); 
    hist(RV, 100, freq = FALSE, xlim = bound, xlab = xlabel, main = plotTitle)
    vals <- seq(interval[1], interval[2], by = 0.001)
    lines(vals, func2(vals), col = "red")
  }
  else{
      plot(func, from = interval[1], to = interval[2], 
           axes = T, lwd = 2, las = 1, ylab = ylabel, xlab = "",
           main = plotTitle)
    if(is.null(RV) == F && is.null(func2) == F){
      legend("topright", legend = c(as.expression(body(func)),
                                    as.expression(body(func2))),
             col=c("black", "red"), lty=1:2, cex=0.8)
      par(new = T)
      plot(density(RV),xlab = "", ylab = "", lty = 4, col = "red",
           lwd = 2, main= "", xaxt = "n", yaxt = "n")
    }
    else{
      title(xlab = xlabel, main = plotTitle)
    }
  }
}
plotArgMax  <- function(func, result, contour = F, steps = F, stepsCol = NULL,
                        colour = "gold", pointLines = T){
  argMax = result$ArgMax
  if(contour == T){
    if(steps == T && (is.null(start) == F)){
      argMaxSteps = as.matrix(result$iters)
      ag      = argMaxSteps[dim(argMaxSteps)[1],]
      points(argMaxSteps[1,1], argMaxSteps[1,2], col = "green", pch = 19)
      if(is.null(stepsCol)){
        lines(argMaxSteps, lwd = 2, col = "blue")
      }
      else{
        lines(argMaxSteps, lwd = 2, col = stepsCol)
      }
      points(ag[11], ag[2], col = "red", pch = 19)
      text(x = c(argMaxSteps[1,1], ag[1]),
           y = c(argMaxSteps[1,2], ag[2]),
           labels = c("Start", "ArgMax"), cex = 1.5, col = c("green", "red"))
    }else{
      if(pointLines == F){
        points(argMax[1], argMax[2], type = "p", pch = 16, cex = 2, col = colour)
      }else{
        abline(h = argMax[2], lty = 2, lwd = 1, col = "red")
        abline(v = argMax[1], lty = 2, lwd = 1, col = "red")
        points(argMax[1], argMax[2], type = "p", pch = 16, cex = 2, col = colour)
      }
    }
   
  }else if(steps == T){
    argMaxSteps = result$iters
    ag = argMaxSteps[length(argMaxSteps)]
    points(argMaxSteps[1], func(argMaxSteps[1]), col = "green", pch = 19)
    if(is.null(stepsCol)){
      lines(argMaxSteps,y = func(argMaxSteps), lwd = 2, col = "blue")
    }
    else{
      lines(argMaxSteps,y = func(argMaxSteps), lwd = 2, col = stepsCol)
    }
    points(ag, func(ag), col = "red", pch = 19)
    text(x = c(argMaxSteps[1], ag),
         y = c(func(argMaxSteps[1]), func(ag)),
         labels = c("Start", "ArgMax"), cex = 1.5, col = c("green", "red"))
  }else{
    if(pointLines == F){
      points(argMax, func(argMax), type = "p", pch = 16, cex = 2, col = colour)
    }else{
      abline(h = func(argMax), lty = 2, lwd = 1, col = "red")
      abline(v = argMax, lty = 2, lwd = 1, col = "red")
      points(argMax, func(argMax), type = "p", pch = 16, cex = 2, col = colour)
    }
  }
}
plotContour <- function(func, interval, levels = 200, colour = T,
                        plotTitle = NULL){
  x = seq(interval[1], interval[2], length = 101)
  y = x
  z = outer(x,y)
  for(i in 1:length(x)){
    for(j in 1:length(y)){
      z[i,j] = func(c(x[i], y[j]))
    }
  }
  cols        <- hcl.colors(10, "YlOrRd")
  if(colour == T){
    contour(x, y, z, nlevels = levels, col = cols, main = plotTitle)
  }else{
    contour(x, y, z, nlevels = levels, main = plotTitle)
  }
}
plotSimAnn  <- function(func, func2 = NULL, RV = NULL, simAnneal, interval, 
                       contour = F, hist = F, radians = F, colour = "blue",
                       add = F, plotTitle = NULL){
  ag        <- simAnneal$ArgMax
  props     <- simAnneal$AcceptedProps
  if (contour == TRUE){
    if(add == F){
      plotContour(func = func, interval = interval, plotTitle = plotTitle)
    }
    points(props[1,1], props[1,2], col = "green", pch = 19)
    lines(props, lwd = 2, col = colour)
    points(ag[1], ag[2], col = "red", pch = 19)
    text(x = c(props[1,1], ag[1]),
         y = c(props[1,2], ag[2]),
         labels = c("Start", "ArgMax"),
         pos = 1, col = c("green", "red"))
  }
  else if (is.null(func2) == F && is.null(RV) == F && hist == T){
    if(add == F){
      plotGraph(func, func2, interval, RV = RV, hist = T, plotTitle = plotTitle)
      }
    lines(props, func2(props), lwd = 2, col = colour)
    points(props[1], func2(props[1]), col = "green", pch = 19, cex = 2)
    points(ag[1], func2(ag[1]), col = "red", pch = 19, cex = 2)
    text(x = c(props[1], ag[1]),
         y = c(func2(props[1]), func2(ag[1])),
         labels = c("Start", "ArgMax"),
         pos = 2, col = c("green", "red"))
  }
  else{
    if(add == F){
      plotGraph(func, interval = interval, radians = radians, 
                plotTitle = plotTitle)
    }
    lines(props, func(props), lwd = 2, col = colour)
    points(props[1], func(props[1]), col = "green", pch = 19, cex = 2)
    points(ag[1], func(ag[1]), col = "red", pch = 19, cex = 2)
    text(x = c(props[1], ag[1]),
         y = c(func(props[1]), func(ag[1])),
         labels = c("Start", "ArgMax"),
         pos = 1, col = c("green", "red"))
  }
}
plotBox.Opt <- function(func, func2, interval, RV = NULL, nfuncArgsVec = 1, 
                        randomStart = F,loglik = F, B = 1000, yaxis = NULL){
  if(randomStart == T){
    res1    <- numeric()
    res3    <- numeric()
    while(length(res1) < B && length(res3) < B){
      point  <- runif(nfuncArgsVec, interval[1], interval[2])
      if(loglik == T){
        try(newt <- NewtonRaphson(func2, start = point)$ArgMax, silent = TRUE)
        res1 <- c(res1, newt)
      }
      else{
        try(newt<- NewtonRaphson(func, start = point)$ArgMax, silent = T)
        res1 = c(res1, newt)
      }
      siman <- SimAnnealing(func, nfuncArgsVec = nfuncArgsVec, 
                             start = point, interval = interval,
                             max.iter = 1000)$ArgMax
      res3 = c(res3, siman)
    }
  }else{
    if(loglik == T){
      res1    <- unlist(replicate(B,
                              NewtonRaphson(func2, interval = interval)$ArgMax))
    }else{
      res1    <- unlist(replicate(B,
                               NewtonRaphson(func, interval = interval)$ArgMax))
    }
    
    res3      <- unlist(replicate(B,
                              SimAnnealing(func, nfuncArgsVec = nfuncArgsVec,
                                               interval = bound,
                                               max.iter = 1000)$ArgMax))
  }
  if(is.null(RV)){
    res2        <- unlist(replicate(B,
                                    MC_optim(func, interval = interval)$ArgMax))
    df          <- data.frame(res1, res2, res3)
    if(is.null(yaxis)){
      boxplot(df, names = c("NewtRaph","MC Optim","Sim Anneal"))
    }
    else{
      boxplot(df, names = c("NewtRaph","MC Optim","Sim Anneal"), ylim = yaxis)
    }

  }else{
    res2        <- unlist(replicate(B,
                                    MC_optim(func, interval = interval)$ArgMax))
    res2.1      <- unlist(replicate(B,
                                    MC_optim(func, RV = RV, sample = 100,
                                             interval = interval)$ArgMax))
    df          <- data.frame(res1, res2, res2.1, res3)
    
    if(is.null(yaxis)){
      boxplot(df, names = c("NewtRaph","MC Optim", "MC Optim 2","Sim Anneal"))
    }
    
    else{
      boxplot(df, names = c("NewtRaph","MC Optim", "MC Optim 2",
                            "Sim Anneal"), ylim = yaxis)
    }
  }
  title(main = paste("Monte Carlo Optimization Boxplot \nSimulations:", B))
}
# Monte Carlo Optimization Functions --------------------------------------
NewtonRaphson <- function(func, start = NULL, interval){
  if(is.null(start)){
    sol       <- optimize(func, interval = interval, maximum = T)
    res       <- list(ArgMax = sol$maximum)
  }
  else{
    nr        <- nlm(func, p = start, gradtol = 1e-10)
    nr1       <- start
    for(i in 1:(nr$iterations)){
      nr1     <- rbind(nr1, nlm(func, p = start, iterlim = i)$estimate)
    }
    x         <- nr$estimate
    res       <- list(ArgMax = x, iters = nr1)
  }
  return(res)
}
MC_optim      <- function(func, nfuncArgsVec = 1, 
                          sample = 1000, RV = NULL, interval){
  if(nfuncArgsVec != 1){
    if(is.null(RV)){
      U       <- runif(nfuncArgsVec*sample, interval[1], interval[2])
      rv      <- matrix(U, ncol = nfuncArgsVec)
    }else{
      rv      <- RV(sample)
    }
    hh        <- func(rv)
    loc       <- which.max(hh)
    result    <- list(ArgMax = round(rv[loc,]), digits = 4)
  }else{
    if (is.null(RV)){
      rv      <- runif(sample, min = interval[1], max = interval[2])
    }else{
      rv      <- RV(sample)
    }
    hh          <- func(rv)
    loc         <- which.max(hh)
    result      <- list(ArgMax = round(rv[loc],digits = 4))
  }
  return(result)
}
SimAnnealing  <- function(func, nfuncArgsVec = 1, max.iter = NULL,
                          interval, start = NULL){
  # If no starting, choose random from interval
  if(is.null(start)){
    x       <- runif(nfuncArgsVec, min = interval[1], max = interval[2])
    x       <- matrix(x, ncol = nfuncArgsVec)
  }else{
    if(length(start) != nfuncArgsVec){
      stop("start must have same number of arguements as your function")
    }else{
      x     <- start
      x     <- matrix(x, ncol = nfuncArgsVec)
    }
  }
  # Setting the values for getting the argMax and Max
  argMax    <- xCur      <- prop    <- x
  globalMax <- maxCur    <- maxProp <- maxVal <- func(x)
  #Parameters for the While Loop  
  iter      <- diff      <- temp    <- scale  <- factor  <- 1
  if(is.null(max.iter)){
    while (diff > 1e-10){
      prop    <- x[iter,] + runif(nfuncArgsVec,-1,1)*scale[iter]
      #Accept Prop value
      if(min(prop) >= interval[1] && max(prop) <= interval[2]){
        # Probability
        rho   <- min(exp((func(prop) - maxCur[length(maxCur)])/temp[iter]), 1)
        if(runif(1) <= rho){
          xCur  <- rbind(xCur, prop)
          maxCur<- c(maxCur, func(prop)) 
        }
      }else {
        # Reject Proposed Value
        xCur  <- rbind(xCur, xCur[dim(xCur)[1],])
        maxCur<- c(maxCur, maxCur[length(maxCur)])
      }
      #Keeping tab of proposed x values 
      x       <- rbind(x, prop)
      
      # Update the best ArgMax
      if(globalMax[length(globalMax)] < maxCur[length(maxCur)]){
        argMax      <- rbind(argMax, xCur[dim(xCur)[1],])
        globalMax   <- rbind(globalMax,func(xCur[dim(xCur)[1],]))
      } 
      #Temperature update and While loop checks
      temp    <- c(temp, 1/(0.1*log(1+iter)))
      scale   <- sqrt(temp)
      if (iter > 10^4 && (length(unique(x[(iter/2):iter])) > 1)&&(diff > 1e-10)){
        diff  <- globalMax[length(globalMax)] - max(maxCur[1:(iter/2)])
      }
      iter    <- iter + 1
    }
  }else{
    while (diff > 1e-10 && iter <= max.iter ){
      prop    <- x[iter,] + rnorm(nfuncArgsVec)*scale[iter]
      #Accept Prop value
      if(min(prop) >= interval[1] && max(prop) <= interval[2]){
        # Probability
        rho   <- min(exp((func(prop) - maxCur[length(maxCur)])/temp[iter]), 1)
        if(runif(1) <= rho){
          xCur  <- rbind(xCur, prop)
          maxCur<- c(maxCur, func(prop)) 
        }
      }else {
        # Reject Proposed Value
        xCur  <- rbind(xCur, xCur[dim(xCur)[1],])
        maxCur<- c(maxCur, maxCur[length(maxCur)])
      }
      #Keeping tab of proposed x values 
      x       <- rbind(x, prop)
      
      # Update the best ArgMax
      if(globalMax[length(globalMax)] < maxCur[length(maxCur)]){
        argMax      <- rbind(argMax, xCur[dim(xCur)[1],])
        globalMax   <- rbind(globalMax,func(xCur[dim(xCur)[1],]))
      } 
      #Temperature update and While loop checks
      temp    <- c(temp, 1/(0.1*log(1+iter)))
      scale   <- sqrt(temp)
      if (iter > 10^4 && (length(unique(x[(iter/2):iter])) > 1)&&(diff > 1e-10)){
        diff  <- globalMax[length(globalMax)] - max(maxCur[1:(iter/2)])
      }
      iter    <- iter + 1
    }
  }
  res       <- list(ArgMax = argMax[dim(argMax)[1],], 
                    Max = globalMax[length(globalMax)],
                    ArgMaxvec = argMax,
                    Maxvec = globalMax,
                    AcceptedProps = xCur, 
                    iterations = iter
  )
  return(res)
}
# Example 1 ---------------------------------------------------------------
### Function
h.1     <- function(x) sin(sqrt(x))
bound   <- c(0,2*pi)
plotGraph(func = h.1, interval = bound, radians = T)

# Newton Raphson
res     <- NewtonRaphson(func = h.1, interval = bound)
plotGraph(func = h.1, interval = bound, radians = T, 
          plotTitle = "Newton-Raphson ArgMax Point")
plotArgMax(func = h.1, result = res)



# MC Optimization
res1    <- MC_optim(func = h.1, interval = bound)
plotGraph(func = h.1, interval = bound, radians = T, 
          plotTitle = "Basic MC Optimization")
plotArgMax(func = h.1, result = res1)

# Simulated Annealing
res2    <- SimAnnealing(func = h.1, nfuncArgsVec = 1, interval = bound)
plotGraph(func = h.1, interval = bound, radians = T)
plotArgMax(func = h.1, result = res2)
plotSimAnn(func = h.1, interval = bound, simAnneal = res2, radians = T,
           plotTitle = "Simmulated Annealing")

# Box Plot
plotBox.Opt(func = h.1, interval = bound, B= 500, yaxis = c(1.5,3))

# Example 3 ---------------------------------------------------------------
### Function
h.3     <- function(x) x*exp((1-x^2)/2)
bound   <- c(0,4)
plotGraph(func = h.3, interval = bound)

# Newton Raphson
res     <- NewtonRaphson(f = h.3, interval = bound)
plotGraph(func = h.3, interval = bound, 
          plotTitle = "Newton-Raphson ArgMax Point")
plotArgMax(func = h.3, result = res)

# MC Optimization
res1    <- MC_optim(func = h.3, interval = bound)
plotGraph(func = h.3, interval = bound, plotTitle = "Basic MC Optimization")
plotArgMax(func = h.3, result = res1)

# MC Optimization (using gamma RV)
Temp    <- 0.1
Z       <- integrate(function(x)exp(h.3(x)/Temp), 
                     lower = bound[1], upper = bound[2])
h.3.new <- function(x) exp(h.3(x)/Temp)/Z$value 
g       <- function (x) 1/11*dgamma(x, shape = 2, rate = 1) + 0.42
M = optimize(function(x) h.3.new(x)/g(x), 
             interval = bound, maximum = TRUE)$objective

g.rv      <- function(n){
  accept  <- numeric()
  y       <- numeric()
  while(sum(accept) < n){
    ynew  <- 1/11*rgamma(1,7,1) + 0.42
    u     <- runif(1)
    acceptnew = 1*(M*u <= h.3.new(ynew)/g(ynew))
    
    y     <- c(y, ynew)
    accept<- c(accept, acceptnew)
  }
  return(y[accept==1])
}

# Checking accept-Reject plot
rv = g.rv(10000)
hist(rv, freq = F, breaks = 100, xlim = bound, 
     main = "Gamma Density Accept-Reject \nfor MC Optimization")
plot(h.3.new, lty = 1, add = T, xlim = bound, col = "red", lwd = 2)

res1.1  <- MC_optim(func = h.3, RV = g.rv, interval = bound)
plotGraph(func = h.3, interval = bound, 
              plotTitle = "MC Optimization Using\n Gamma Density")
plotArgMax(func = h.3, result = res1.1)

# Simulated Annealing
res2    <- SimAnnealing(func = h.3, nfuncArgsVec = 1, interval = bound)
plotSimAnn(func = h.3, interval = bound, simAnneal = res2, 
           plotTitle = "Simmulated Annealing")

# BoxPlot
plotBox.Opt(func = h.3, interval = bound, RV = g.rv, nfuncArgsVec = 1,
            B = 500, yaxis = c(0.9, 1.1))

# Example 4 ---------------------------------------------------------------
# Function
h.4     <- function(x) exp(cos(3*x)*sin(4*x)^4)
h.4.neg <- function(x) -exp(-cos(x/2))*sin(5*x)
bound   <- c(0, 2*pi)

plotGraph(func = h.4, interval = bound, radians = T)

# # Newton Raphson
# res     <- NewtonRaphson(func = h.4, interval = bound)
# plotGraph(func = h.4, interval = bound, radians = T)
# plotArgMax(func = h.4, result = res)

# Using Newton Raphson with a starting point: 
start1  <- pi 
start2  <- 1/5
res.1   <- NewtonRaphson(func = h.4.neg, start = start1)
res.2   <- NewtonRaphson(func = h.4.neg, start = start2)
plotGraph(func = h.4, interval = bound, radians = T, 
          plotTitle = "Newton-Raphson Method \n(2 different starting points)")
plotArgMax(func = h.4, result = res.1, steps = T)
plotArgMax(func = h.4, result = res.2, steps = T, stepsCol = "magenta")

# MC Optimization
res1    <- MC_optim(func = h.4, interval = bound)
plotGraph(func = h.4, interval = bound, radians = T,
          plotTitle = "Basic MC Optimization")
plotArgMax(func = h.4, result = res1)

# Simulated Annealing
res2    <- SimAnnealing(func = h.4, interval = bound)
plotSimAnn(func = h.4, interval = bound, simAnneal = res2, radians = T,
           plotTitle = "Simulated Annealing")

# Box-plot
plotBox.Opt(func = h.4, interval = bound, nfuncArgsVec = 1, B = 500, 
            randomStart = T, yaxis = c(0, 3*pi))

# Example 6 ---------------------------------------------------------------
# density
dmixgamma       <- function(x, params) {
  p               <- params[1]
  alpha           <- c(params[2], params[3])
  beta            <- c(params[4], params[5])
  prob            <- c(p, 1- p)
  k               <- length(prob)
  n               <- length(x)
  den              <- function(i){
    prob[i]*dgamma(x, alpha[i], beta[i])
  }
  if( n == 1){den(1)}
  else{rowSums(vapply(1:k, den, numeric(n)))}
}
# random generation
rmixgamma       <- function(n, params) {
  p       <- params[1]
  alpha   <- c(params[2], params[3])
  beta    <- c(params[4], params[5])
  prob    <- c(p, 1- p)
  k <- sample.int(length(prob), n, replace = TRUE, prob = prob)
  rgamma(n, alpha[k], beta[k])
}
# log-likelihood
loglik.mixGamma <- function(x, params){
  pdf = dmixgamma(x, params = params)
  sum(log(pdf))
}

p               <- 1/3
alpha           <- c(38, 50)
beta            <- c(20, 11)
param           <- c(p, alpha, beta)
bound           <- c(0,8)

# Plot of the mixture model
rv = rmixgamma(1e5, params = param)
plotGraph(func = function(x) loglik.mixGamma(x, params = param),
          func2= function(x) dmixgamma(x, params = param),
          interval = bound, RV = rv, hist = T)

# Newton Raphson starting at different points
start1   <- 4
start2   <- 7
res.1    <- NewtonRaphson(func = function(x)-loglik.mixGamma(x, params = param), 
                          start = start1)
res.2    <- NewtonRaphson(func = function(x)-loglik.mixGamma(x, params = param), 
                          start = start2)
plotGraph(func = function(x) loglik.mixGamma(x, params = param),
          func2= function(x) dmixgamma(x, params = param),
          interval = bound, RV = rv, hist = T, 
          plotTitle = "Newton-Raphson Method \n(2 different starting points)")
plotArgMax(func = function(x) dmixgamma(x, params = param), result = res.1,
           steps = T)
plotArgMax(func = function(x) dmixgamma(x, params = param), result = res.2,
           steps = T, stepsCol = "magenta")


# MC Optimization
# MC Optimization (Not Good!)
plotGraph(func = function(x) loglik.mixGamma(x, params = param),
          func2= function(x) dmixgamma(x, params = param),
          interval = bound, RV = rv, hist = T,
          plotTitle = "Basic MC Optimization \n(10 Simulations)")
for (i in 1:10){
  res1    <- MC_optim(func = function(x) -loglik.mixGamma(x,params = param), 
                      interval = bound)
  plotArgMax(func = function(x) dmixgamma(x, params = param),
              result = res1, pointLines = F)
}

# Simulate Annealing with a Starting Point
start1    <- 4
start2    <- 7
res2.1    <- SimAnnealing(func = function(x) loglik.mixGamma(x, params = param),
                   nfuncArgsVec = 1, start = start1, interval = bound)
res2.2    <- SimAnnealing(func = function(x) loglik.mixGamma(x, params = param),
                        nfuncArgsVec = 1, start = start2, interval = bound)
plotSimAnn(func = function(x) loglik.mixGamma(x, params = param),
           func2 = function(x) dmixgamma(x, params = param),
           RV = rv, simAnneal = res2.1, interval = bound, hist = T,
           plotTitle = "Simmulated Anneal at 2 Different Points")

plotSimAnn(func = function(x) loglik.mixGamma(x, params = param),
           func2 = function(x) dmixgamma(x, params = param),
           RV = rv, simAnneal = res2.2, interval = bound, hist = T, add = T, 
           colour = "orange")

plotBox.Opt(func = function(x) loglik.mixGamma(x, params = param),
            interval = bound, nfuncArgsVec = 1, loglik = T,
            func2 = function(x) -loglik.mixGamma(x, params = param), B = 10,
            randomStart = T
            )

# Example 7 ---------------------------------------------------------------
dmixlogis  = function(x, params){
  p        = params[1]
  location = c(params[2], params[3])
  scale    = c(params[4], params[5])
  p*dlogis(x) + (1-p)*dlogis(x)
}
rmixlogis  = function(n, params){
  if (length(params) != 5){
    stop("params must have 1 value for prob,
       and 2 each for location and scale.")
  }
  p        = params[1]
  location = c(params[2], params[3])
  scale    = c(params[4], params[5])
  if (length(n) == 1){
    rv = c(location[1] + scale[1]*rlogis(n),
           location[2] + scale[2]*rlogis(n)
    )
  }else{
    rv = c(location[1] + scale[1]*rlogis(n[1]),
           location[2] + scale[2]*rlogis(n[2]) 
    )
  }
  return(rv)
}
loglik.mixlogis = function(X, rv){
  sum(log(prob*dlogis(rv-X[1]) + (1-prob)*dlogis(rv-X[2])))
}

prob     = 1/4
loc      = c(2, 5)
scal     = c(1, 3)
numOfrv  = c(100, 300)
param = c(prob, loc, scal)
bound = c(-5,15)
mixlogis.rv = rmixlogis(n = numOfrv, params = param)

# Plot of the Mixture model
plotContour(func = function(X) loglik.mixlogis(X, rv = mixlogis.rv), 
            interval = bound, colour = F)

# Newton Raphson starting at different points
start1   <- c(10,-1)
start2   <- c(-1,10)
res.1    <- NewtonRaphson(
                    func = function(X) -loglik.mixlogis(X,rv = mixlogis.rv), 
                    start = start1
                    )
res.2    <- NewtonRaphson(
  func = function(X) -loglik.mixlogis(X,rv = mixlogis.rv), 
  start = start2)

plotContour(func = function(X) loglik.mixlogis(X, rv = mixlogis.rv),
            interval = bound, colour = T,
            plotTitle = "Newton-Raphson Method at 2 different Points")
plotArgMax(result = res.1, contour = T, steps = T, stepsCol = "purple")
plotArgMax(result = res.2, contour = T, steps = T, stepsCol = "blue")


# MC Optimization
# MC Optimization (Not Good!)
plotContour(func = function(X) loglik.mixlogis(X, rv = mixlogis.rv),
            interval = bound, colour = T, 
            plotTitle = "Basic MC Optimization \n(20 Simulations)")
for (i in 1:20){
  res1    <- MC_optim(
            func = function(X) -loglik.mixlogis(X, rv = mixlogis.rv),
            interval = bound, nfuncArgsVec = 2 )
  plotArgMax (result = res1, contour = T, colour = "blue", pointLines = F)
}




# Simulate Annealing with a Starting Point
start1   <- c(10,-1)
start2   <- c(-1,10)
set.seed(21)
res2.1  <- SimAnnealing(
                      func = function(X) loglik.mixlogis(X, rv = mixlogis.rv),
                      nfuncArgsVec = 2, 
                      start = start1, 
                      interval = bound
                      )
set.seed(31)
res2.2  <- SimAnnealing(
                        func = function(X) loglik.mixlogis(X, rv = mixlogis.rv),
                        nfuncArgsVec = 2, 
                        start = start2, 
                        interval = bound
                        )
plotSimAnn(func =function(X) loglik.mixlogis(X,rv = mixlogis.rv), 
           simAnneal = res2.1, interval = bound, contour = T, 
           plotTitle = "Simulated Annealing at 2 Different Points")
plotSimAnn(func =function(X) loglik.mixlogis(X,rv = mixlogis.rv), 
           simAnneal = res2.2, interval = bound, contour = T, add = T,
           colour = "purple")    
              

# Example 2 Not Using -----------------------------------------------------
### Function 
h.2     <- function(x) (-x^3-x^5-1)/(1+x^2)
bound   <- c(0,1)
plotGraph(func = h.2, interval = bound)


# Newton Raphson
res     <- NewtonRaphson(f = h.2, interval = bound)
plotGraph(func = h.2, interval = bound)
plotArgMax(func = h.2, result = res)

# MC Optimization
res1    <- MC_optim(func = h.2, interval = bound)
plotGraph(func = h.2, interval = bound)
plotArgMax(func = h.2, result = res1)

# Simulated Annealing
res2    <- SimAnnealing(func = h.2, nfuncArgsVec = 1, interval = bound)
plotGraph(func = h.2, interval = bound)
plotSimAnn(func = h.2, interval = bound, simAnneal = res2)

# Box Plot
plotBox.Opt(func = h.2, interval = bound, B= 500)
# Example 5 Not Using -----------------------------------------------------
### Function
schwefel <- function (X) {
  sum(X*sin(sqrt(abs(X))))
}

bound    <- c(-3*pi, 3*pi)
plotContour(func = schwefel, interval = bound)

# Newton Raphson
res      <- NewtonRaphson(func = schwefel, interval = bound)
plotContour(func = schwefel, interval = bound)
plotArgMax(func = schwefel, result = res)

# Newton Raphson At a different point
start    <- c(-8,-8)
res.1    <- NewtonRaphson(func= function (x) -schwefel(x), start, bound)
plotContour(func = schwefel, interval = bound)
plotArgMax (func = schwefel, result = res.1, contour = T, steps = T)

start    <- c(5,-5)
res.1    <- NewtonRaphson(func= function (x) -schwefel(x), start, bound)
plotContour(func = schwefel, interval = bound)
plotArgMax (func = schwefel, result = res.1, contour = T, steps = T)

# MC Optimization (Not Good!)
res1     <- MC_optim(func = schwefel, interval = bound, nfuncArgsVec = 2)
plotContour(func = schwefel, interval = bound)
plotArgMax (result = res1, contour = T)

#Simulated Annealing
res2     <- SimAnnealing(func = schwefel, nfuncArgsVec = 2, interval = bound)
plotSimAnn(func = schwefel, simAnneal = res2, interval = bound, contour = T)

plotBox.Opt(func = schwefel, interval = bound, nfuncArgsVec = 1, B = 50, 
            randomStart = T, yaxis = c(-3*pi, 3*pi))
