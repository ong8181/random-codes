####
#### Twin surrogate function (Rcpp version)
####

# Load packages
library(Rcpp); packageVersion("Rcpp") # 1.0.7
library(tseriesChaos); packageVersion("tseriesChaos") # 0.1.13.1
library(rEDM); packageVersion("rEDM") # 0.7.4

# Load Cpp function
sourceCpp("TwinSurrogate.cpp")

# Load data and visialize
ds <- read.csv("data/logistic.csv")
y <- ds$Y
plot(y, type = "l")

# Determine the best embedding dimension
y_std <- as.numeric(scale(y)) # Normalize
E_test <- simplex(y_std, E = 1:20)
E_y <- E_test[which.min(E_test$mae),'E'] # Use MAE as an index of forecasting skill; E = 2

# Generate twin surrogate (iter = 100)
twin_y <- TwinSurrogateRcpp(y_std,
                            dim = E_y,
                            num.iter = 100,
                            surrogate.option = "random",
                            initial.point = "twin",
                            distance.method = "norm", # "norm" is the default option
                            s.update = "on",
                            n.twin.threshold = 10,
                            output.message = T) # "output.message = T" shows the number of twins

# Visualize an example
par(las = 1)
plot(y_std, type = "l", lwd = 2, ylab = "Value", ylim=c(-2,2.5))
lines(twin_y[1], col = 2)
legend(110, 2.3, c("Original", "Twin surrogate"), col=c(1,2), lty=1, lwd=c(2,1))
