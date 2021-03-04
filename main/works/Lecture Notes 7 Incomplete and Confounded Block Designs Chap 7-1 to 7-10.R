# Lecture Notes 7
# Chapter 7.1-1.10 Incomplete and Confounded Block Designs

library(EMSaov)
library(AlgDesign)
library(BsMD)
library(car)
library(partitions)
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
library(later)
library(htmltools)
library(digest)
library(promises)
library(classInt)
library(units)
library(agricolae)
library(daewr)

# 7.2 Balanced Incomplete Block (BIB) Designs

# minimum number of blocks for BIB design
library(daewr)
BIBsize(4,2)

# D-optimal design for BIB
# b = 6 blocks of size k = 2, t = 4 levels for each factor
library(AlgDesign)
BIB <- optBlock( ~ ., withinData = factor(1:4), blocksizes = rep(2, 6))
des <- BIB$rows # save 
dim(des) <- NULL
des <- matrix(des, nrow = 6, ncol = 2, byrow = TRUE,
dimnames = list(
  c( "Block1", "Block2", "Block3", "Block4", "Block5", "Block6"), 
  c("unit1", "unit2")))
des

# 7.3 analysis of BIB mous data.frame
mous <- c(
  1, 2, 3, 7, 8, 9,
  1, 4, 5, 7, 10, 11,
  2, 4, 6, 8, 10, 12,
  3, 5, 6, 9, 11, 12
)
temp <- as.character(c(
  0, 0, 0, 0, 0, 0,
  10, 10, 10, 10, 10, 10,
  20, 20, 20, 20, 20, 20,
  30, 30, 30, 30, 30, 30
))
roc <- c(
  10.3, 9.1, 6.1, 9.3, 9.4, 5.1,
  1.8, 9.8, 4.2, 2.8, 8.8, 5.2,
  6.5, 8.1, 8.1, 5.5, 9.1, 7.1,
  18.1, 17.2, 14.1, 12.1, 15.2, 14.4
)
rocd <- data.frame(cbind(mous,temp,roc)); rocd

mod1 <- aov( roc ~ mous + temp, data = rocd)
summary(mod1)

library(lsmeans)
lsmeans(mod1, pairwise ~ temp, adjust = ("tukey"))

# number of replicates for treatment levels 
# r = bk/t = 12(2)/4 = 6
library(daewr) 
rmin <- 3 # smallest number of replicates
rmax <- 12 # largest number of replicates
sigma <- sqrt(4.67) # sqrt of msE
alpha <- .05
Delta <- 5 # ROC differs by 5
nlev <-  6 # four ROC levels
nreps <- c(rmin:rmax)
power <- Fpower1(alpha, nlev, nreps, Delta, sigma)
options(digits = 5)
power

# number of replicates for blocks 
# b = tr/k = 4(6)/2 = 12
library(daewr) 
rmin <- 3 # smallest number of replicates
rmax <- 12 # largest number of replicates
sigma <- sqrt(4.67) # sqrt of msE
alpha <- .05
Delta <- 5 # ROC differs by 5
nlev <-  12 # 12 mice (blocks)
nreps <- c(rmin:rmax)
power <- Fpower1(alpha, nlev, nreps, Delta, sigma)
options(digits = 5)
power

# 7.4 BTIB and PBIB designs
mous <- c(
  1, 2, 3, 4, 5, 6,
  1, 2, 3, 4, 5, 6
)
temp <- as.character(c(
  0, 0, 0, 0, 0, 0,
  10, 10, 20, 20, 30, 30
))
roc <- c(
  10.3, 9.1, 6.1, 9.8, 4.2, 8.1,
  1.8, 6.5, 18.1, 8.1, 17.2, 14.1
)
rocm <- data.frame(cbind(mous,temp,roc)); rocm

modm <- aov( roc ~ mous + temp, data = rocd)
library(lsmeans)
lsmeans(modm,pairwise~temp,adjust=("tukey"))

# 7.4 generalized cyclic incomplete block design
# b = t = 6, k = 2
library(agricolae)
treat <- c(1, 2, 3, 4, 5, 6)
des <- design.cyclic(treat, k = 2, r = 2)
des$book

# 7.5 Row and Column designs
library(agricolae)
treat <- c(1, 2, 3, 4, 5, 6)
RCD <- design.cyclic(treat, k = 4, r = 4, rowcol = TRUE, seed = 1)
RCD$book

# 7.6 confounded 2^k designs
# CCBF design: 2^4 = 16-run, k = 4 factors, 1 confounding block ABD 
library(FrF2)
d.ccbf <- FrF2( 16, 4, blocks = c("ABCD"),
               alias.block.2fis = TRUE, randomize = FALSE)
d.ccbf
# lm fit
y <- c(60, 64, 60, 51, 69, 61, 41, 57, 61, 56, 52, 44, 61, 66, 45, 63)
d <- add.response(d.ccbf, response = y) 
fit <- lm( y ~ Blocks + A * B * C * D, data = d)
fit
# half-normal plot of effects
effects <- coef(fit)
effects <- effects[3:17] # strip off intercept and block, leave factors
effects <- effects[ !is.na(effects) ] # remove NAs
library(daewr)
halfnorm(effects, names(effects), alpha=.25)

# 7.6 confounded 2^k-p designs
# CCBFF design: 2^6-2 = 16-run, k = 6 factors, 2 confounding blocks ABC, ABD 
library(FrF2)
d.ccbff <- FrF2(16, 6, generators = c("ABC", "ABD"), # recall,  gives minimum aberration design
            blocks = c("AB", "AC"), alias.block.2fis = TRUE, 
            randomize = FALSE)
d.ccbff

# 7.6 no block, for aliases
d.ccbff.noblock <- FrF2(16, 6, generators = c("ABC", "ABD"),  
                randomize = FALSE)
d.ccbff.noblock
library(DoE.base)
generators(d.ccbff.noblock) 
y <- runif(16, 0, 1) 
aliases( lm( y ~ (.)^3, data = d.ccbff.noblock)) # up to 3-way factors

# 7.6 d.ccbff design + responses
y <- c(60, 64, 60, 51, 69, 61, 41, 57, 61, 56, 52, 44, 61, 66, 45, 63)
d <- add.response(d.ccbff, response = y) 
fit <- lm( y ~ Blocks + A * B * C * D * E * F, data = d)
fit
# 7.6 half-normal plot of effects
effects <- coef(fit)
effects <- effects[5:17] # strip off intercept and block, leave factors
effects <- effects[ !is.na(effects) ]; effects # remove NAs
library(daewr)
halfnorm(effects, names(effects), alpha=.25)

# 7.6 choosing block defining contrasts
sun.design <- FrF2(16, 6, blocks = 4,
      alias.block.2fis = TRUE, randomize = FALSE)
summary(sun.design)

# 7.7 confounding 3-level and p-level factorial designs
library(AlgDesign)
Blockdes <- gen.factorial(3, nVars = 3, factors = "all",
  varNames = c("A", "B", "C"))
Blockdes <- optBlock( ~ A + B + C + A:B + A:C + B:C, 
  withinData = Blockdes,
  blocksizes = c(rep(5, 5)), 
  criterion = "D")
Blockdes

# 7.7 confounding p-level factorial designs
library(AlgDesign)
Blockdes <- gen.factorial(5, nVars = 4, factors = "all",
                          varNames = c("A", "B", "C", "D"))
Blockdes <- optBlock( ~ A + B + C + D + A:B + A:C +
                        A:D + B:C + B:D + C:D, 
                      withinData = Blockdes,
                      blocksizes = c(rep(20, 20)), 
                      criterion = "D")
Blockdes

# 7.8 cook and nachtsheim's D_s efficient method
# for blocking mixed level designs
library(AlgDesign)
des <- gen.factorial(levels = c(4, 2, 5), # three factors levels 4, 2, 5 
        factors = 'all',
        varNames = c("A", "B", "C"))
bdes <- optBlock( ~ A + B + C + A:B + A:C + B:C, # one/two-factor estimable
        withinData = des, 
        blocksizes = c(rep(12, 6)), # create 6 blocks of size 12 each
        criterion = "D") # use D efficient method
bdes
bdes$Blocks

# 7.8 blocking for orthogonal arrays
library("DoE.base") # determine number of OA runs, if possible
show.oas(factors = list(nlevels = c(2, 3, 4), number = c(3, 1, 3)))
library("DoE.base") # main effects, block orthogonal wrt one another
fnames =c ("A", "B", "C", "D", "E", "F", "Block")
BlockOA <- oa.design(nlevels = c(2, 2, 2, 3, 4, 4, 4),
    factor.names = fnames, seed=104, nruns = 48) # 4 blocks of size 12
BlockOA <- BlockOA[order(BlockOA$Block), ]; BlockOA

library(DoE.base) # 4 blocks of size 6, but sacrifice orthogonality
library(AlgDesign)
fnames <- c("A", "B", "C", "D", "E", "F")
cand <- oa.design(nlevels = c(2, 2, 2, 3, 4, 4, 4),
                  factor.names = fnames, 
                  randomize = TRUE, 
                  seed = 104,
                  nruns = 48)
bdes <- optBlock( ~ A + B + C + D + E + F,
                  withinData = cand, 
                  blocksizes = c(rep(6, 4)), # 4 blocks of size 6
                  criterion = "D")
bdes

# 7.9 partially confounded blocked factorial design
library(AlgDesign)
Blockdes <- gen.factorial(2, nVars = 2, factors = "all",
                          varNames = c("A","B"))
Blockdes <- optBlock( ~ A + B + A:B, withinData = Blockdes,
                      blocksizes = c(rep(2,6)), criterion = "D")
Blockdes

# 7.9 partially confounded mixed level factorial
library(AlgDesign)
des <- gen.factorial(levels = c(4, 2, 5), # three factors levels 4, 2, 5 
                     factors = 'all',
                     varNames = c("A", "B", "C"))
bdes <- optBlock( ~ A + B + C + A:B + A:C + B:C, # one/two-factor estimable
                  withinData = des, 
                  blocksizes = c(rep(12, 6)), # create 6 blocks of size 12 each
                  criterion = "D") # use D efficient method
bdes
bdes$Blocks
Block <- as.factor(c(rep(1:6, each = 12))) # blocks
roc <- c(8.9, 3.4, 5.6, 5.9, 2.8, 9.0, 7.7, 8.1, 8.9, 5.8, 6.7, 7.1,
         3.9, 5.4, 4.6, 4.9, 4.8, 8.0, 6.7, 5.1, 7.9, 4.8, 5.7, 8.3,
         8.2, 4.4, 6.6, 5.7, 4.8, 7.1, 5.7, 9.1, 7.9, 4.3, 8.7, 8.6,
         7.1, 6.4, 7.6, 7.9, 8.8, 9.3, 7.2, 8.3, 8.6, 5.6, 4.7, 8.5,
         4.4, 6.4, 8.6, 8.2, 7.8, 6.0, 8.7, 7.1, 7.9, 7.8, 3.7, 7.9,
         8.2, 8.4, 9.6, 7.1, 2.3, 5.7, 5.7, 6.7, 5.9, 8.8, 7.7, 8.0
         )
bdesign <- cbind(Block, bdes$design,roc); bdesign # store block design for lm  routine
library(car)
modf <- lm(roc ~ Block + A + B + C + A:B + A:C + B:C + A:B:C, 
  data = bdesign,
  contrasts = list(A = contr.sum, 
        B = contr.sum, 
        C = contr.sum,
        Block = contr.sum))
summary.aov(modf, type="III")
library(lsmeans)
lsmeans(modf, pairwise ~ A, adjust = ("tukey"))

###################################################
# Homework 7 and Test 4 questions
###################################################

# 7.2.b pp 303--305
# BIB design: minimum number of blocks for BIB design
library(daewr)
BIBsize(10,3); BIBsize(10,4); BIBsize(10,5); BIBsize(10,6)

# 7.2.c.1 pp 303--305
# D-optimal BIB design, for b = 30, k = 3
library(AlgDesign)
BIB <- optBlock( ~ ., withinData = factor(1:10), 
    blocksizes = rep(3, 10)) # size 3, from 10 treatments
des <- BIB$rows # save 
dim(des) <- NULL
des <- matrix(des, 
    nrow = 30, # 30 blocks
    ncol = 3, byrow = TRUE, # each block size 3
    dimnames = list(
    c("B1", "B2", "B3", "B4", "B5", 
      "B6", "B7", "B8", "B9", "B10", 
    "B11", "B12", "B13", "B14", "B15", 
    "B16", "B17", "B18", "B19", "B20", 
    "B21", "B22", "B23", "B24", "B25", 
    "B26", "B27", "B28", "B29", "B30"), 
    c("treat1", "treat2", "treat3")))
des

# 7.2.c.2 pp 303--305
# D-optimal BIB design, for b = 15, k = 4
library(AlgDesign)
BIB <- optBlock( ~ ., withinData = factor(1:10), 
                 blocksizes = rep(4, 10)) # size 4, from 10 treatments
des <- BIB$rows # save 
dim(des) <- NULL
des <- matrix(des, 
    nrow = 15, # 15 blocks
    ncol = 4, byrow = TRUE, # each block size 4
    dimnames = list(
    c("B1", "B2", "B3", "B4", "B5", 
      "B6", "B7", "B8", "B9", "B10", 
      "B11", "B12", "B13", "B14", "B15"), 
    c("treat1", "treat2", "treat3", "treat4")))
des

# 7.2.d pp 303--305
# generalized cyclic incomplete block design
# b = t = 10, k = 3
library(agricolae)
treat <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
des <- design.cyclic(treat, k = 3, r = 3)
des$book

# 7.6.b pp 303--305
# alias pattern for 2^(6-1) design
library(FrF2)
design <- FrF2( 32, 6, generators = "ABCDE", randomize = FALSE)
design
y <- runif(32, 0, 1)
aliases( lm( y ~ (.)^3, 
    data = design)) # (.)^3 gives saturated model up to 3-way factors

# 7.10 partially confounded mixed level factorial
library(AlgDesign)
des <- gen.factorial(levels = c(2, 3, 3), # three factors levels 2, 3, 3 
                     factors = 'all',
                     varNames = c("A", "B", "C"))
bdes <- optBlock( ~ A*B*C, # one/two/three-factor estimable
                  withinData = des, 
                  blocksizes = c(rep(6, 6)), # create 6 blocks of size 6 each
                  criterion = "D") # use D efficient method
bdes
bdes$Blocks
Block <- as.factor(c(rep(1:6, each = 6))) # blocks
roc <- c(8.9, 3.4, 5.6, 5.9, 2.8, 9.0, 
         7.7, 8.1, 8.9, 5.8, 6.7, 7.1,
         3.9, 5.4, 4.6, 4.9, 4.8, 8.0, 
         6.7, 5.1, 7.9, 4.8, 5.7, 8.3,
         8.2, 4.4, 6.6, 5.7, 4.8, 7.1, 
         5.7, 9.1, 7.9, 4.3, 8.7, 8.6
)
bdesign <- cbind(Block, bdes$design,roc); bdesign # store block design for lm  routine
library(car)
modf <- lm(roc ~ Block + A + B + C + A:B + A:C + B:C + A:B:C, 
           data = bdesign,
           contrasts = list(A = contr.sum, 
                            B = contr.sum, 
                            C = contr.sum,
                            Block = contr.sum))
summary.aov(modf, type="III")
library(lsmeans)
lsmeans(modf, pairwise ~ A, adjust = ("tukey"))

