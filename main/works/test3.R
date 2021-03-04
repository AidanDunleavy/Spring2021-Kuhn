# 6.2.a pp 252-256
# 2^(4-1) design;  2^3 = 8-run, k = 4-factors, p = 2 
library(FrF2)
design <- FrF2( 8, 4, generators = "ABC", randomize = FALSE)
design
y <- runif(8, 0, 1) # "dummy" values for response
aliases( # (.)^3 gives saturated model up to 3-way factors
  lm( y ~ (.)^3, data = design)
) 

# 6.2.c.1 pp 252-256
# estimated effects, for 2^(4-1) design, melo study
melo <- FrF2(8, 4, 
             generators = "ABC",
             randomize = FALSE)
weight <- c(4.70, 14.67, 1.71, 3.73, 9.47, 7.61, 0.15, 4.78); 
library(DoE.base)
melo <- add.response( melo , weight ); melo
m <- lm( weight ~ 1 + A + B + C + D
         + A:B + A:C + A:D, data = melo) 
summary(m)

# 6.2.c.2 half-normal plot of 2^(4-1) design, for melo study
library(daewr)
LGB(coef(m)[-1], alpha = .1) # plot and report
halfnorm(
  coef(m)[-1], 
  names, 
  alpha = 0.25, 
  refline = TRUE
) 

# 6.2.d fit model after removing 3 |smallest| effects
A <- (as.numeric( melo$A) - 1.5 ) / .5; A
B <- (as.numeric( melo$B) - 1.5 ) / .5; B
C <- (as.numeric( melo$C) - 1.5 ) / .5; C
D <- (as.numeric( melo$D) - 1.5 ) / .5; D
AC <- A*C; AC # create AC interaction by elementwise mult
melo.AC <- cbind( melo, AC); melo.AC
library(DoE.base)
m <- lm( weight ~ 1 + A + B + D + AC, data = melo.AC); summary(m)  

# 6.2.f fit full factorial model A + B + D + BD
BD <- B*D; BD # create BD interaction by elementwise mult = BD interaction
melo.BD <- as.data.frame(cbind( A, B, D, BD, weight)); melo.BD
library(DoE.base)
m <- lm( weight ~ A*B*D, data = melo.BD); summary(m)  
m <- lm( weight ~ A + B + D + BD, data = melo.BD); summary(m)  
m <- lm( weight ~ B + D + BD, data = melo.BD); summary(m)  

# 6.2.g BD interaction plot
melo <- FrF2(8, 4, 
             generators = "ABC",
             factor.names = list(
               sucrose = c(150, 250), 
               temperature = c(20, 30), 
               yeast.conc = c(2.0, 5.0),
               agitation = c(50, 100)),
             randomize = FALSE)
weight <- c(4.70, 14.67, 1.71, 3.73, 9.47, 7.61, 0.15, 4.78); 
library(DoE.base)
melo <- add.response( melo , weight ); melo
interaction.plot( # weight best when temp low, agitation high
  melo$temperature, melo$agitation, melo$weight, 
  type="b",pch=c(24,18,22), leg.bty="o",
  main="Interaction Plot of Levan Weight for Temperature x Agitation",
  xlab="Temperature (degrees C)", 
  ylab="Average Levan Weight"
)
# average levan weight for two different sucrose levels
cellmeans <- tapply( melo$weight, melo$sucrose, mean); cellmeans
