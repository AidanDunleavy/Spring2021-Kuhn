# Lecture Notes 1-B
# 2 Completely Randomized Designs with One Factor

library(later)
library(htmltools)
library(digest)
library(promises)
library(classInt)
library(units)
library(agricolae)
library(AlgDesign)
library(BsMD)
library(car)
library(daewr)
library(DoE.base)
library(FrF2)
library(GAD)
library(gmodels)
library(leaps)
library(lme4)
library(lsmeans)
library(mixexp)
library(multcomp)
library(Vdgraph)
library(pracma)
library(zip)

# 2.1 Introduction 
# 2.2 Replication and Randomization 

# drug types
set.seed(7638)
f <- factor( rep ( c("A", "B", "C" ), each = 5)); f # 5 x 3 = 15 drug types
fac <- sample ( f, 15 ) # SRS 15 from 15 patients
eu <- 1:15 # eu means experimental unit
plan <- data.frame ( patient = eu, drug = fac ); plan

# heart rates
set.seed(38)
f <- factor( rep ( c("NM-S", "MT-S", "H-S", 
                      "NM-M", "MT-M", "H-M",
                      "NM-L", "MT-L", "H-L"), each = 2)); f 
fac <- sample ( f, 18 ) # SRS from artists x volumes
eu <- 1:18 # eu means experimental unit
plan <- data.frame ( student = eu, artist.volume = fac ); plan

# 2.3 A Historical Example 
# 2.4 Linear Model for CRD 

# 2.4.1 drug doses A
y <- drug.response.A <- c(5.90, 5.92, 5.91, 5.89, 5.88,
                        5.51, 5.50, 5.50, 5.49, 5.50,
                        5.01, 5.00, 4.99, 4.98, 5.02)
dose <- factor( rep ( c(10, 20, 30 ), each = 5)); dose
patient <- as.numeric(1:15); patient; class(patient)
drug <- data.frame(patient, dose, drug.response.A); drug

m1 <- aov( drug.response.A ~ dose, data = drug )
summary( m1 )

# 2.4.2 drug doses B
y <- drug.response.B <- c(5.90, 4.42, 7.51, 7.89, 3.78,
                          6.31, 3.54, 4.73, 7.20, 5.72,
                          4.52, 6.93, 4.48, 5.55, 3.52)
dose <- factor( rep ( c(10, 20, 30 ), each = 5)); dose
patient <- as.numeric(1:15); patient; class(patient)
drug <- data.frame(patient, dose, drug.response.B); drug

m2 <- aov( drug.response.B ~ dose, data = drug )
summary( m2 )

# 2.4.4 matrix representation, drug doses B
X1 <- c(1,1,1,1,1,
        1,1,1,1,1,
        1,1,1,1,1)
X2 <- c(0,0,0,0,0,
        1,1,1,1,1,
        0,0,0,0,0)
X3 <- c(0,0,0,0,0,
        0,0,0,0,0,
        1,1,1,1,1)
X <- as.matrix(data.frame(cbind(X1,X2,X3))); X

ssE <- t(y) %*% (diag(15) - X %*% solve(t(X) %*% X) %*% t(X)) %*% y; ssE # alternative formula
sigma2.est <- ssE/(15-3); sigma2.est

m1 <- aov( drug.response.A ~ dose, data = drug )
summary( m1 )

beta.est <- solve(t(X)%*%X) %*% t(X) %*% y; beta.est

m2 <- aov( drug.response.B ~ dose, data = drug )
summary( m2 )

m3 <- lm(drug.response.B ~ dose, data = drug)
summary(m3)

ssE <- t(y) %*% y - t(beta.est) %*% t(X) %*% y; ssE 

# 2.5 Verifying Assumptions of the Linear Model 
# constant variance, normality of residuals

m1 <- aov( drug.response.A ~ dose, data = drug )
summary( m1 )

par( mfrow = c(2,2) )
plot( m1, which = 5 ) # should be constant stnd resid for different factor levels
plot( m1, which = 1 ) # should be constant for different fitted values
plot( m1, which = 2 ) # normal prob plot for normality should be linear
plot( residuals(m1) ~ patient, 
      main = "Residuals vs Exp. Unit", 
      font.main = 1, 
      data = drug)
abline( h=0, lty = 2 ) # should be constant
par(mfrow = c(1,1))

# 2.6 Analysis Strategies When Assumptions Are Violated 

# 2.6.1 box-box transformations
y <- drug.response.A <- c(5.90, 5.92, 5.91, 5.89, 5.88,
                          5.51, 5.50, 5.50, 5.49, 5.50,
                          5.01, 5.00, 4.99, 4.98, 5.02)
x <- dose <- factor( rep ( c(10, 20, 30 ), each = 5)); dose

library(MASS) # find lambda for transformation
bc <- boxcox(m1); bc
lambda <- bc$x[ which.max( bc$y ) ]; lambda

# 2.6.2 run ANOVA on box-cox transformed data
library(MASS) 
tdrug <- transform(drug, tdrug.response = drug.response.A^(1.232323) ); tdrug
m1t <- aov( tdrug.response ~ dose, data = tdrug )
summary(m1t) 

par( mfrow = c(2,2) ) # plots slightly better with transformed data
plot( m1t, which = 5 ) # should be constant stnd resid for different factor levels
plot( m1t, which = 1 ) # should be constant for different fitted values
plot( m1t, which = 2 ) # normal prob plot for normality should be linear
plot( residuals(m1t) ~ patient, 
      main = "Residuals vs Exp. Unit", 
      font.main = 1, 
      data = drug)
abline( h=0, lty = 2 ) # should be constant
par(mfrow = c(1,1))

# 2.6.2 distribution-based transformations

# 2.6.3 weighted least squares, p 35
 
with( drug, { # "with" means use all statements in {} with data "drug"
  std <- sqrt( tapply( drug.response.A, dose, var) )
  weights <- rep( 1/std, each = 5 ) 
  m1w <- lm( drug.response.A ~ dose, weights = weights, data = drug )
  anova( m1w )
})

# 2.6.4 generalized linear model: error distribution not normal

# 2.7 Determining the Number of Replicates 

m1 <- aov( drug.response.A ~ dose, data = drug )
summary( m1 ) # drug response, case A

library(daewr)
rmin <- 2 # smallest number of replicates
rmax <- 6 # largest number of replicates
alpha <- rep(0.05, rmax - rmin +1)
sigma <- sqrt(0.0002)
nlev <- 3
nreps <- rmin:rmax
Delta <- 5
power <- Fpower1(alpha, nlev, nreps, Delta, sigma); power  

m2 <- aov( drug.response.B ~ dose, data = drug )
summary( m2 ) # drug response, case B

library(daewr)
rmin <- 2 # smallest number of replicates
rmax <- 6 # largest number of replicates
alpha <- rep(0.05, rmax - rmin +1)
sigma <- sqrt(2.3332)
nlev <- 3
nreps <- rmin:rmax
Delta <- 5
power <- Fpower1(alpha, nlev, nreps, Delta, sigma); power  


# 2.8 Comparison of Treatments after the F-test 

# 2.8.1 pre-planned comparisons, qualitative factors
library(daewr)
m1 <- aov( drug.response.A ~ dose, data = drug )
summary( m1 )
con1 <- matrix( c(1, -1/2, -1/2, 
                 0, 1, -1 ), 3, 2 )
L1 <- t(con1)
rownames(L1) <- c("10 vs 1/2(20 + 30)", "-20 vs 30"); L1

options( digits = 3)
library(gmodels)
fit.contrast(m1, "dose", L1)

# 2.8.2 linear, quadratic trends, quantitative data
contrasts(drug$dose) <- contr.poly(3)
contrasts(drug$dose)

m1 <- aov( drug.response.A ~ dose, data = drug )
summary.lm(m1)

# Tukey's HSD method of comparison

m1 <- aov( drug.response.A ~ dose, data = drug )
m1.tukey <- TukeyHSD( m1, ordered = T )
m1.tukey

# SNK method of pairwise comparison

library(agricolae)
compare <- SNK.test( m1, "dose",  alpha = 0.05 )
print(compare)

# Dunnett method, compare to control (best), defaults to first level

summary(drug)
library(multcomp)
drug.dun <- glht(m1, linfct = mcp( dose = "Dunnett"),
                  alternative = "less")
summary(drug.dun)


# 2.9 Review of Important Concepts

### Test 1 questions

# 2.2.f pp 51-53
y <- flight.time <- c(
  5.2, 5.1, 5.3, 5.3, 5.2, 5.2, 5.0, 5.2,
  5.2, 5.1, 5.1, 5.0, 5.2, 5.1, 5.1, 5.2,
  5.1, 5.1, 5.0, 5.0, 5.1, 5.0, 4.9, 5.1,
  4.8, 4.9, 4.8, 5.0, 4.8, 4.9, 4.9, 4.9
); y; class(y)
x <- wing.length <- factor(c(rep(4,8), rep(4.75, 8), rep(5.5,8), rep(6,8))); x; class(x)
helicopter <- as.numeric(1:32); helicopter; class(helicopter)
paper.heli <- data.frame(helicopter, wing.length, flight.time); paper.heli

set.seed(7638) # randomization
f <- factor( rep ( c(4, 4.75, 5.5, 6 ), each = 8)); f
fac <- sample ( f, 32 )
eu <- 1:32
plan <- data.frame ( helicopter = eu, wing.length = fac ); plan

# 2.2.h pp 51-53
# see 2.2.f for paper.heli data
m.heli <- aov( flight.time ~ wing.length, data = paper.heli )
summary( m.heli )

# 2.2.i pp 51-53
# see 2.2.f for paper.heli data
par( mfrow = c(1,2) )
plot( m.heli, which = 5 ) # should be constant for different fitted values
plot( m.heli, which = 2 ) # normal prob plot for normality should be linear
par(mfrow = c(1,1))

# 2.2.j pp 51-53
# see 2.2.f for paper.heli data
contrasts(paper.heli$wing.length) <- contr.poly(4)
contrasts(paper.heli$wing.length)

m.heli1 <- aov( flight.time ~ wing.length, data = paper.heli )
summary.lm(m.heli1)

# 2.4.b pp 51-53
y <- rise.height <- c(
  11.4, 11.0, 11.3, 9.5,
  27.8, 29.2, 26.8, 26.0,
  47.6, 47.0, 47.3, 45.5,
  61.6, 62.4, 63.0, 63.9
); y; class(y)
x <- bak.powder <- factor(c(rep(0.25,4), rep(0.5, 4), rep(0.75,4), rep(1,4))); x; class(x)
biscuit <- as.numeric(1:16); biscuit; class(biscuit)
dough.biscuit <- data.frame(biscuit, bak.powder, rise.height); dough.biscuit

m.biscuit <- aov( rise.height ~ bak.powder, data = dough.biscuit )
summary( m.biscuit )

# 2.4.c pp 51-53
# see 2.4.b for dough.biscuit data
contrasts(dough.biscuit$bak.powder) <- contr.poly(4)
contrasts(dough.biscuit$bak.powder)

m.biscuit1 <- aov( rise.height ~ bak.powder, data = dough.biscuit )
summary.lm(m.biscuit1)

# 2.4.d pp 51-53
# see 2.4.b for dough.biscuit data
m.biscuit <- aov( rise.height ~ bak.powder, data = dough.biscuit )
summary( m.biscuit ) # variance estimate equals MSE

# 2.4.d pp 51-53
# see 2.4.b for dough.biscuit data
par( mfrow = c(1,2) )
plot( m.biscuit, which = 5 ) # should be constant for different fitted values
plot( m.biscuit, which = 2 ) # normal prob plot for normality should be linear
par(mfrow = c(1,1))

# 2.6.a and 2.6.b pp 51-53
library(daewr)
rmin <- 2 # smallest number of replicates
rmax <- 10 # largest number of replicates
alpha <- rep(0.05, rmax - rmin + 1)
sigma <- sqrt(2.1)
nlev <- 4
nreps <- rmin:rmax
Delta <- 2*sigma
power <- Fpower1(alpha, nlev, nreps, Delta, sigma)  
power

# 2.6.c pp 51-53
library(daewr)
rmin <- 6 # smallest number of replicates
rmax <- 11 # largest number of replicates
alpha <- rep(0.05, rmax - rmin + 1)
sigma <- sqrt(2.1)
nlev <- 4 # try 2, 4 and 8
nreps <- rmin:rmax
Delta <- 2*sigma
power <- Fpower1(alpha, nlev, nreps, Delta, sigma)  
power

