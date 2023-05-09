# Libraries ---------------------------------------------------------------
library(tictoc); library(ggplot2)
# Plotting Function -------------------------------------------------------
plotGraph     <- function(func, func2 = NULL, interval, 
                       RV = NULL,radians = F, hist = F){
  xlabel      <- as.expression(body(func)); 
  ylabel      <- expression(h(x));
  if (radians == T){
    plot(func, from = interval[1], to = interval[2], axes = F, 
         xlab = xlabel, ylab = ylabel)
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
    tit       <- "Distribution"
    xlabel    <- as.expression(body(func2)); 
    hist(RV, 100, freq = FALSE, xlim = bound, xlab = xlabel, main = tit)
    vals <- seq(interval[1], interval[2], by = 0.001)
    lines(vals, func2(vals), col = "red")
  }
  else{
    plot(func, from = interval[1], to = interval[2], 
         axes = T, lwd = 2, las = 1, ylab = ylabel, xlab = "")
    if(is.null(RV) == F && is.null(func2) == F){
      legend("topright", legend = c(as.expression(body(func)),
                                    as.expression(body(func2))),
             col=c("black", "red"), lty=1:2, cex=0.8)
      par(new = T)
      plot(density(RV),xlab = "", ylab = "", lty = 4, col = "red",
           lwd = 2, main= "", xaxt = "n", yaxt = "n")
    }
    else if (is.null(func2) == F){
      curve(func2, add = T, lty = 2, lwd = 4, col = "red")
      legend("bottomright", legend = c(as.expression(body(func)),
                                    as.expression(body(func2))),
             col=c("black", "red"), lty=1:2, cex=0.8)
    }
    else{
      title(xlab = xlabel)
    }
  }
}
plotBox.Integ <- function(func, ImpFunc, ImpRV, Stratas, 
                       B = 1000, n = 100, yaxis = NULL,interval){
  true        <- integrate(func, interval[1],  interval[2])$value
  res1        <- unlist(replicate(B,
                          MC(func = func, interval = interval, sim = n)$Ihat))
  res2        <- unlist(replicate(B, 
                               MC_imp(func, func2 = ImpFunc, RV = ImpRV,
                                      interval = interval, sim = n)$Ihat))
  res3        <- unlist(replicate(B, 
                               MC_strat(func, interval = interval,
                                        stratas = Stratas, sim = n)$Ihat))
  
  df          <- data.frame(res1, res2, res3)
  if(is.null(yaxis)){
    boxplot(df, names = c("Naive MC","Importance MC","Stratifed MC"))
  }else{
    boxplot(df, names = c("Naive MC","Importance MC","Stratifed MC")
            ,ylim = yaxis)
  }
  title(main = paste("Monte Carlo IntegrationBoxplot \nSample Size:", n, 
                     "  Simulations:", B))
  abline(h = true, col = "red", lwd = 2, lty = 2)
}
# Monte Carlo Functions ---------------------------------------------------
#Monte Carlo Estimate Functions (Naive, Important, Stratified)
MC        <- function(func, interval, sim = 10^5, plots = F){
  rv      <- runif(sim, interval[1], interval[2])
  x       <- func(rv)*(interval[2]-interval[1])
  Ihat    <- mean(x)
  SE2     <- var(x)/sim
  true    <- integrate(func, interval[1], interval[2])$value
  if (plots == T){
    X     <- cumsum(x)/(1:sim)
    SE    <- sqrt(cumsum(x - X)^2)/(1:sim)
    alpha <- qnorm(.05/2,lower.tail = F)
    title <- "Naive Monte Carlo Integration\nMean & CI Range"
    plot(X, xlab = "Iterations", type = "l", lwd = 2,
           ylim = mean(x)+20*c(-SE[sim], SE[sim]), ylab = "", main = title)
      lines(X+alpha*SE, col = "gold", lwd = 2)
      lines(X-alpha*SE, col ="gold", lwd = 2)
      abline(h = true, lty = 3, lwd = 2)
      legend("topright", legend = c("Estimated Value",
                                    "95% C.I.",
                                    "Expected Value"),
             col=c("black", "gold"), lty=1:1:3, lwd = 2)
  }
  return(list(Ihat = Ihat, SE2 = SE2, I = true))
}
MC_imp    <- function(func, func2, RV, interval, sim = 10^5, plots = F){
  rv      <- RV(sim)
  weight  <- (1/(interval[2] - interval[1]))/func2(rv)
  x       <- (interval[2] - interval[1])*func(rv)*weight
  Ihat    <- mean(x)
  SE2     <- var(x)/sim
  true    <- integrate(func, interval[1], interval[2])$value
  if(plots == T){
    X     <- cumsum(x)/(1:sim)
    SE    <- sqrt(cumsum(x - X)^2)/(1:sim); alpha = qnorm(.05/2,lower.tail = F)
    title <- "Importance Sampling Monte Carlo Integration\nMean & CI Range"
    plot(X, xlab = "Iterations", type = "l", lwd = 2,
          ylim = mean(x)+20*c(-SE[sim], SE[sim]), ylab = "", main = title)
     lines(X+alpha*SE, col = "gold", lwd = 2)
     lines(X-alpha*SE, col ="gold", lwd = 2)
     abline(h = true, lty = 3, lwd = 2)
     legend("topright", legend = c("Estimated Value",
                                   "95% C.I.",
                                  "Expected Value"),
            col=c("black", "gold"), lty=1:1:3, lwd = 2)
  }
  return(list(Ihat = Ihat, SE2 = SE2, I = true))
}
MC_strat  <- function(func, stratas, interval, sim = 10^5){
  if(stratas > sim) stop("sim must be larger than stratas")
  stratLen<- diff(interval)/stratas
  strats  <- seq(from = interval[1], interval[2], length.out = stratas + 1)
  hh      <- vector(mode = "list", length = stratas)
  for(j in 1:stratas){
      hh[[j]]  <- func(runif(sim/stratas, strats[j], strats[j+1]))
  }
  t       <- lapply(lapply(hh, mean), FUN = "*", stratLen)
  tvar    <- lapply(lapply(hh, sd), FUN = "*", stratLen)
  Ihat    <- sum(unlist(t))
  SE2     <- (sum(stratLen^2*unlist(tvar)^2/sim))^2
  true    <- integrate(func, interval[1], interval[2])$value
  return(list(Ihat = Ihat, SE2 = SE2, I = true))
}


# Example 1 ---------------------------------------------------------------
h.1     <- function(x) sqrt(x^3 + sqrt(x))
f       <- function(x, m) 1
bound   <- c(0, 1)

# Plot 
plotGraph(func = h.1, interval = bound)

# Naive MC
set.seed(2)
res     <- MC(func = h.1, interval = bound, sim = 10^5, plots = T)

# Importance Sampling
g     <- function(x) 2*x
g.rv  <- function(n){
  u = runif(n)
  sqrt(u)
}
plotGraph(func = h.1, interval = bound)
curve(g, add = T, lty = 3, col = "blue", lwd = 3)

func_test1 = function(x){
  res <- h.1(x)*f(x)/g(x)
}

curve(expr = func_test1, from = 0, to  = pi, lty = 4, col = "red", lwd = 4,
      ylab = "|h|f/g", ylim = c(0,10))

set.seed(10)
res1    <- MC_imp(func = h.1, func2 = g, RV = g.rv, 
                  interval = bound, sim = 10^5, plot = F)

# Stratified Sampling
res2    <- MC_strat(func = h.1, stratas = 4, 
                    interval = bound, sim = 10^5, plot = T)


# Box Plots
plotBox.Integ(func = h.1, ImpFunc = g, ImpRV = g.rv, 
              interval = bound, Stratas = 4, n = 10^3, yaxis = c(0.8,1.2))

# Table Results
estimate     <- c(res$Ihat, res1$Ihat, res$Ihat)
error        <- c(res$SE2, res1$SE2, res1$SE2)
df           <- data.frame(estimate, error)
rownames(df) <- c("Naive Monte Carlo", "Importance Sampling", "Stratified Sampling")


# Example 2 ---------------------------------------------------------------
h.2     <- function(x) exp(-5*(x-3)^4)
bound   <- c(1, 5)
plotGraph(func = h.2, interval = bound)

# Naive MC
res     <- MC(func = h.2, interval = bound, sim = 10^5, plots = T)

# Importance Sampling

# Plot of g and h.2
g       <- function(x) dnorm(x, mean = 3, sd = 1)
g.rv    <- function(x) rnorm(x, mean = 3, sd = 1)
set.seed(2)
rv      <- g.rv(10^6)
plotGraph(func = h.2, func2 = g, interval = bound, RV = rv)

set.seed(13)
res1    <- MC_imp(func = h.2,func2 = g,RV = g.rv, interval = bound, 
                  sim = 10^5, plots = T)

# Stratified Sampling
res2    <- MC_strat(func = h.2, stratas = 8, interval = bound, sim = 10^5)




plotBox.Integ(func = h.2, ImpFunc = g, ImpRV = g.rv, 
              interval = bound, Stratas = 8, n = 10^5)

# Table Results
estimate     <- c(res$Ihat, res1$Ihat, res$Ihat)
error        <- c(res$SE2, res1$SE2, res1$SE2)
df           <- data.frame(estimate, error)
rownames(df) <- c("Naive Monte Carlo", "Importance Sampling", "Stratified Sampling")




# Example 3 ---------------------------------------------------------------
h.3     <- function(x) exp(-x)/(1+x^2)
bound   <- c(0, 20)
plotGraph(func = h.3, interval =bound)

# Naive MC
res     <- MC(func = h.3, interval = bound, sim = 10^5, plots = F)

# Importance Sampling
f       <- function(x) 1 ## Uniform(0,1)
g.0     <- function(x) 1 ## Uniform(0,1)
g.1     <- function(x) exp(-x) ##Exp(1)
g.2     <- function(x) 1/((1 + x^2)*pi) ## Cauchy(0,1)
g.3     <- function(x) exp(-x)/(1 - exp(-1)) 
g.4     <- function(x) 4/((1 + x^2)*pi) ## Cauchy(0,1) x 4

func0   <- function(x){
  res   <- h.3(x)*f(x)/g.0(x)
  return(res)
}
func1   <- function(x){
  res   <- h.3(x)*f(x)/g.1(x)
  return(res)
}
func2   <- function(x){
  res   <- h.3(x)*f(x)/g.2(x)
  return(res)
}
func3   <- function(x){
  res   <- h.3(x)*f(x)/g.3(x)
  return(res)
}
func4   <- function(x){
  res = h.3(x)*f(x)/g.4(x)
  return(res)
}

# my.df <-data.frame(x = c(lower,upper))
# ggplot(my.df, aes(x=x)) + xlim(c(lower,upper)) + ylim(c(0,3)) +
#   stat_function(fun = func4, geom = "point", aes(colour = "4/((1+x^2)*pi")) +
#   stat_function(fun = func3, geom = "point", aes(colour = "exp(-x)/(1-e^(-1)"))+ 
#   stat_function(fun = func2, geom = "point", aes(colour = "1/(1+x^2)*pi")) +
#   stat_function(fun = func1, geom = "point", aes(colour = "exp(-1)")) +
#   stat_function(fun = func0, geom = "point", aes(colour = "1"))   + 
#   scale_colour_manual("Different g Functions", 
#                       values = c("red", "blue", "green", "orange", "gold")) + 
#   ggtitle("|h|f/g Plot for each g function")

par(bg = "lightgrey")
curve(expr = func0, from = 0, to  = 3, lty = 4, col = "darkred", lwd = 4,
      ylab = "|h|f/g")
curve(expr = func1, add = T, lty = 4, col = "darkgreen", lwd = 4)
curve(expr = func2, add = T, lty = 4, col = "darkorange", lwd = 4)
curve(expr = func3, add = T, lty = 4, col = "darkmagenta", lwd = 4)
curve(expr = func4, add = T, lty = 4, col = "navy", lwd = 4)
legend("topright", legend = c(as.expression(body(g.0)),
                              as.expression(body(g.1)),
                              as.expression(body(g.2)),
                              as.expression(body(g.3)),
                              as.expression(body(g.4))
                              ),col=c("darkred",
                                      "darkgreen",
                                      "darkorange",
                                      "darkmagenta",
                                      "navy"
                                      ),lty=4, lwd = 4)



# Choosing g.3 as our importance function
g.3.rv  <- function(x, interval = bound){
  u     <- runif(x, interval[1], interval[2])
  x     <- -log(1 - (1 - exp(-1))*u)
}

# h.3 and g.3 Graph
RV      <- g.3.rv(10^6, interval = bound)
plotGraph(func = h.3, func2 = g.3, interval = bound, RV = RV)

# Importance Sampling vs. Naive MC
res1    <- MC_imp(func = h.3, func2 = g.3, RV = g.3.rv,
                  interval = bound, sim = 10^5, plots = T)

# Stratified Sampling
res2    <- MC_strat(func = h.3, interval = bound, stratas = 4, 
                    sim = 10^5, plots = T)

# Box plot Comparison
plotBox.Integ(func = h.3, ImpFunc = g.3, ImpRV = g.3.rv, 
              interval = bound, Stratas = 4, yaxis = c(0.4,0.6))

# Comparing Sys Time ------------------------------------------------------
func = h.3; interval = c(0,1); n = 10^5; B = 1000;ImpFunc = g.3
ImpRv = g.3.rv; Statas = 4
tic("Naive Monte Carlo time")
test1   <- unlist(replicate(B,
                  MC(func = func, interval = interval, sim = n)$Ihat))
toc()
tic("Importance Sampling Monte Carlo  time")
test2   <- unlist(replicate(B, 
                    MC_imp(func, func2 = ImpFunc, RV = ImpRV,
                                   interval = interval, sim = n)$Ihat))
toc()
tic("Stratified Sampling Monte Carlo time")
test2   <- unlist(replicate(B, 
                            MC_strat(func, interval = interval,
                                     stratas = Stratas, sim = n)$Ihat))
toc()

h = function (x) 1/(x^2 + 5*x + 6)
bound = c(0,1)

curve(h, from = 0, to = 1)

MC_strat(func = h, strata = 16, interval = bound, sim = 10^4, plot = T)

dd = function (x) 1/( x + 2.5 )*(x + 2.8)


curve(dd, col = "red", lty = 3, ylim = c(0,2))
curve(ii, add = T)


