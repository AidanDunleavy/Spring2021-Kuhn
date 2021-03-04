# Lecture Notes 5
# Chapter 6.1-6.5 Fractional Factorial Designs

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

# 6.1 Introduction
# 6.2 Half-Fractions of 2^k Designs

# 6.2 full 2^4 design, 1 replicate, for sump-pump example
Des <- expand.grid( A = c(-1, 1), # 2^4 design
                    B = c(-1, 1), # 2 x 2 x 2 x 2 = 16 configurations
                    C = c(-1, 1),
                    D = c(-1, 1)) # 16 x 1 replicates = 16
Des[] <- lapply(Des, factor) # convert to factors

defects <- c( # 
  61, 69, 60, 45, 64, 66, 44, 57, 
  60, 61, 52, 41, 56, 61, 51, 63); 
twok1rdata <- data.frame(cbind(Des,defects)); twok1rdata # data twok1r

# 6.2 half-fractional 8-run 2^(4-1) design where D = ABC
library(FrF2)
design <- FrF2( 8, 4, generators = "ABC", randomize = FALSE)
design
design.info(design) # info related to 2^(4-1) design

# 6.2 color map of correlations
library(FrF2)
design <- FrF2(8, 4, randomize = FALSE)
library(daewr)
colormap(design, mod=3)

# 6.2 alias pattern for 2^(4-1) design
library(FrF2)
y <- runif(8, 0, 1)
aliases( lm( y ~ (.)^3, 
          data = design)) # (.)^3 gives saturated model up to 3-way factors

# 6.2 fracional 2^(4-1) design, 1 replicate, for sump-pump example
library(FrF2)
sump.pump <- FrF2(8, 4, 
      generators = "ABC",
      factor.names = list(
        temp = c(0, 30), 
        click = c(60, 120), 
        hard = c(0.5, 0.8),
        surge = c(50, 150)),
      randomize = FALSE)
defects <- c(61, 61, 52, 45, 56, 66, 44, 63); 
library(DoE.base)
sump.pump <- add.response( sump.pump , defects )
m <- lm( defects ~ 1 + temp + click + hard + surge 
         + temp:click + temp:hard + temp:surge, data = sump.pump)
summary(m)

# 6.2 half-normal plot of 2^(4-1) design, for sump-pump example
library(daewr)
LGB(coef(m)[-1], rpt = TRUE) # plot and report

# 6.2 interaction plot temp x clicks
interaction.plot(
  sump.pump$temp, sump.pump$click, sump.pump$defects, 
  type="b",pch=c(24,18,22), leg.bty="o",
  main="Interaction Plot for On-Off Clicks by Temperature",
  xlab="Temperature (degrees C)", 
  ylab="Average Number of Defects")

# 6.3 Quarter and Higher Fractions of $2^k$ Designs

# 6.3  2^(5-2) design: 2^3 = 8-run, k = 5-factors, p = 2 generators D = AB and E = AC
library(FrF2)
design <- FrF2( 8, 5, generators = c("AB", "AC"), randomize = FALSE)
design

# 6.3 aliasing pattern
library(FrF2)
y <- runif(8, 0, 1)
aliases( lm( y ~ (.)^5, data = design)) # up to 5-way factors

# 6.3  2^(6-3) design: 2^3 = 8-run, k = 6-factors, 
# 6.3                  p = 3 generators D = AB, E = AC, F = BC
library(FrF2)
design <- FrF2( 8, 6, generators = c("AB", "AC", "BC"), randomize = FALSE)
design

# 6.3 aliasing pattern
library(FrF2)
y <- runif(8, 0, 1)
aliases( lm( y ~ (.)^5, data = design)) # up to 5-way factors

# 6.4 Criteria for Choosing Generators for $2^{k-p}$ Designs

# 6.4  2^(5-2) design: 2^3 = 8-run, k = 5-factors, p = 2 
library(FrF2)
#      generators D = AB and E = AC
design <- FrF2( 8, 5, generators = c("AB", "AC"), randomize = FALSE)
y <- runif(8, 0, 1) # select 8 numbers between 0 to 1, to add to factorial design
aliases( lm( y ~ (.)^5, data = design)) # up to 5-way factors
#      generators D = AB and E = BC
design <- FrF2( 8, 5, generators = c("AB", "BC"), randomize = FALSE)
y <- runif(8, 0, 1) 
aliases( lm( y ~ (.)^5, data = design)) # up to 5-way factors

# 6.4  minimum aberration: 2^(5-2) design: 2^3 = 8-run, k = 5-factors, p = 2 
library(FrF2)
design <- FrF2( 8, 5, randomize = FALSE) # gives min.aberr when generators NOT specified
library(DoE.base)
generators(design) 
y <- runif(8, 0, 1) 
aliases( lm( y ~ (.)^5, data = design)) # up to 5-way factors

# 6.4 maximum clear: 2^(5-2) design: 2^3 = 8-run, k = 5-factors, p = 2 
design <- FrF2( 8, 5, MaxC2= TRUE, randomize = FALSE) # gives maximum clear design
library(DoE.base)
generators(design) 
y <- runif(8, 0, 1) 
aliases( lm( y ~ (.)^5, data = design)) # up to 5-way factors

# 6.4  minimum aberration: 2^(6-2) design: 2^4 = 16-run, k = 6-factors, p = 2 
library(FrF2)
design <- FrF2( 16, 6, randomize = FALSE) # gives min.aberr when generators NOT specified
library(DoE.base)
generators(design) 
y <- runif(16, 0, 1) 
aliases( lm( y ~ (.)^3, data = design)) # up to 3-way factors
design <- FrF2( 16, 6, MaxC2= TRUE) # gives maximum clear design
library(DoE.base)
generators(design) 
y <- runif(16, 0, 1) 
aliases( lm( y ~ (.)^3, data = design)) # up to 3-way factors

# 6.4 sump-pump: 2^(6-2) maximum clear, minimum aberration design
# 2^(6-2) design: 2^4 = 16-run, k = 6-factors, p = 2
library(FrF2)
sump.pump <- FrF2( 16, 6, randomize = FALSE) # gives min.aberr 
sump.pump
defects <- c(61, 69, 60, 45, 64, 66, 44, 57,
             60, 61, 52, 41, 56, 61, 51, 63); defects
library(DoE.base)
sump.pump <- add.response( sump.pump , defects ); sump.pump
m.sump.pump.62 <- lm( defects ~ (.)^3 , data = sump.pump) # up to 3-ways
summary(m.sump.pump.62)
library(daewr)
cfs <- coef(m.sump.pump.62)[c(2:12,17,19,27,29)]
names <- names(cfs)
LGB(cfs, rpt = FALSE)
halfnorm(cfs, names, alpha = .25, refline = FALSE)
# 2^(5-1) design: 2^4 = 16-run, k = 5-factors, p = 1
library(FrF2)
sump.pump <- FrF2( 16, 5, 
    factor.names = c("A","B","C","D","E"),
    randomize = FALSE) # gives min.aberr 
generators(sump.pump)
sump.pump
defects <- c(61, 69, 60, 45, 64, 66, 44, 57,
             60, 61, 52, 41, 56, 61, 51, 63); defects
library(DoE.base)
sump.pump <- add.response( sump.pump , defects ); sump.pump
m.sump.pump.51 <- lm( defects ~ (.)^2 , data = sump.pump) # up to 2-ways
summary(m.sump.pump.51)
cfs <- coef(m.sump.pump.51)[c(2:16)]
names <- names(cfs)
halfnorm(cfs, names, alpha = .25, refline = FALSE)
# interaction plot of DE
Temp <- 15*((as.numeric(sump.pump$A)- 1.5 ) / 0.5 ) + 15
Hard <- 0.15*((as.numeric(culture2$H) - 1.5) / 0.5) + 0.65
interaction.plot(Temp, Hard, 
    sump.pump$defects, 
    type = "b", pch=c(18,24,22), 
    leg.bty="o",
    main="Interaction Plot for Water Temperature and Hardness",
    xlab = "Temperature (Celcius)", ylab = "Hardness")

# 6.5 Augmeting Fractional Factorials

# 6.5.1 foldover design for specified main factors in 2^5-2
# 2^3 = 8-run, k = 5-factors, p = 2 generators D = AB and E = AC
library(FrF2)
design <- FrF2( 8, 5, 
    generators = c("AB", "AC"),
    randomize = FALSE)
design
y <- runif(8, 0, 1) # dummary y values, to allow aliases to work 
aliases( lm( y ~ (.)^5, data = design)) # up to 5-way factors
designf <- fold.design(design, columns = c(1)) # foldover A
designf 
designf$fold <- NULL # drop "fold" factor
designf
y <- runif(16, 0, 1) # dummary y values, to allow aliases to work 
aliases( lm( y ~ (.)^5, data = designf)) # up to 5-way factors

design.mirror <- fold.design(design, columns = 'full') # mirror design
design.mirror 
design.mirror$fold <- NULL # drop "fold" factor
design.mirror
y <- runif(16, 0, 1) # dummary y values, to allow aliases to work 
aliases( lm( y ~ (.)^5, data = design.mirror)) # up to 5-way factors

# 6.5.1 2^(5-2) design: 2^3 = 8-run, k = 5-factors, p = 2
# 6.5.1 start with this design, then clear main factors
library(FrF2)
sump.pump <- FrF2( 8, 5, 
                   factor.names = c("A","B","C","D","E"),
                   randomize = FALSE) # gives min.aberr 
sump.pump # design only, no data
defects <- c(61, 69, 60, 45, 64, 66, 44, 57); defects
library(DoE.base)
sump.pump.data <- add.response( sump.pump , defects ); sump.pump.data
m.sump.pump.52 <- lm( defects ~ (.)^2 , data = sump.pump.data) # up to 2-ways
summary(m.sump.pump.52)
cfs <- coef(m.sump.pump.52)[c(2:6,11,13)]
names <- names(cfs)
library(daewr)
halfnorm(cfs, names, alpha = .25, refline = FALSE) # B aliased with BE
# clear main factors with mirror design
sump.pump
sump.pump.mirror <- fold.design(design, columns = 'full') # mirror design
sump.pump.mirror 
sump.pump.mirror$fold <- NULL # drop "fold" factor
sump.pump.mirror
defects <- c(61, 69, 60, 45, 64, 66, 44, 57,
             60, 61, 52, 41, 56, 61, 51, 63); defects
sump.pump.data <- add.response( sump.pump.mirror , defects ); sump.pump.data
m.sump.pump.mirror <- lm( defects ~ (.)^3 , data = sump.pump.mirror) # up to 2-ways
summary(m.sump.pump.mirror)
cfs <- coef(m.sump.pump.mirror)[c(2:13,17:19)]
names <- names(cfs)
library(daewr)
halfnorm(cfs, names, alpha = .25, refline = FALSE) # main effects cleared

# 6.5.2 2^(5-2) design: 2^3 = 8-run, k = 5-factors, p = 2
# 6.5.2 start with this design and D-optimize
library(FrF2)
sump.pump <- FrF2( 
  8, 5, 
  factor.names = c("A","B","C","D","E"),
  randomize = FALSE) # gives min.aberr 
defects <- c(61, 69, 60, 45, 64, 66, 44, 57); defects
sump.pump # design only, no data

augm <- fold.design(sump.pump); augm # augmented mirror design

# prepare data as numeric for optFederov function
# for D-optimal or A-optimal design
A <- (as.numeric( augm$A) - 1.5 ) / .5
B <- (as.numeric( augm$B) - 1.5 ) / .5
C <- (as.numeric( augm$C) - 1.5 ) / .5
D <- (as.numeric( augm$D) - 1.5 ) / .5
E <- (as.numeric( augm$E) - 1.5 ) / .5
Block <- augm$fold
augmn <- data.frame(A, B ,C, D, E, Block); augmn
# list of 2^5 = 32 candidate points of full factorial
library(AlgDesign)
cand <- gen.factorial( 
  levels = 2, 
    nVar = 5, 
  varNames = c("A","B", "C", "D", "E")
); cand
# combine original, mirror and candidate into data.frame
Block <- rep('cand', 32)
cand <- data.frame( 
  A=cand$A, B=cand$B, C=cand$C,
  D=cand$D, E=cand$E, Block
)
all <- rbind( augmn, cand); all
fr<-1:16 # "freeze" augmented orginal + mirror design
# choose from 32 candidates for D-optimal design
optim <- optFederov( # optimizes for A + AC + AD + BE 
  ~ B + I(A*C) + I(A*D) + I(B*E), 
  data = all, 
  nTrials = 24, # add 8 runs to 16
  criterion = "D", # D-optimal
  nRepeats = 10, 
  augment = TRUE, 
  rows=fr # frozen augmented
); optim
optim <- as.data.frame(optim); dim(optim)
defects <- c(61, 69, 60, 45, 64, 66, 44, 57,
             60, 61, 52, 41, 56, 61, 51, 63,
             60, 60, 52, 40, 51, 59, 54, 61); defects; dim(defects)
# sump.pump.data <- add.response( optim , defects ); sump.pump.data
sump.pump.data <- cbind(optim, defects); 
sump.pump.data <- sump.pump.data[c(5:10,12)] # extract factors, blocks, defects
names(sump.pump.data)[1:6] <- c("A","B","C","D","E","block") # rename
sump.pump.data

# fit model B + AC + AD + BE + block                      
m.sump.pump.D.optimal <- lm( 
  defects ~ B + A:C + A:D + B:E  + block , 
  data = sump.pump.data) 
summary(m.sump.pump.D.optimal)

# A-optimal choose from 32 candidates for A-optimal design
optim <- optFederov( # optimizes for A + AC + AD + BE 
  ~ B + I(A*C) + I(A*D) + I(B*E), 
  data = all, 
  nTrials = 24, # add 8 runs to 16
  criterion = "A", # A-optimal
  nRepeats = 10, 
  augment = TRUE, 
  rows=fr # frozen augmented
); optim
optim <- as.data.frame(optim); dim(optim)
defects <- c(61, 69, 60, 45, 64, 66, 44, 57,
             60, 61, 52, 41, 56, 61, 51, 63,
             60, 60, 52, 40, 51, 59, 54, 61); defects; dim(defects)
# sump.pump.data <- add.response( optim , defects ); sump.pump.data
sump.pump.data <- cbind(optim, defects); 
sump.pump.data <- sump.pump.data[c(5:10,12)] # extract factors, blocks, defects
names(sump.pump.data)[1:6] <- c("A","B","C","D","E","block") # rename
sump.pump.data

# A-optimal fit model B + AC + AD BE + block                      
m.sump.pump.A.optimal <- lm( 
  defects ~ B + A:C + A:D + B:E  + block , 
  data = sump.pump.data) 
summary(m.sump.pump.A.optimal)

###################################################
# Homework 5 and Test 3 questions
###################################################

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

# 6.4.a.1 2^(8-4) design: 2^4 = 16-run, k = 8-factors, p = 4 
# generators
library(FrF2)
design <- FrF2( 16, 8, randomize = FALSE) # min abbr design
library(DoE.base)
generators(design) 
design

# 6.4.a.2 2^(8-4) design: 2^4 = 16-run, k = 8-factors, p = 4 
# generators
library(FrF2)
design <- FrF2( 16, 8, MaxC2= TRUE) # gives maximum clear design
library(DoE.base)
generators(design) 
design

# 6.4.c 2^(8-4) design: 2^4 = 16-run, k = 8-factors, p = 4 
# aliases
library(FrF2)
design <- FrF2( 16, 8, randomize = FALSE) # min abbr design
library(DoE.base)
generators(design) 
design
y <- runif(16, 0, 1) # "dummy" values for response
aliases( # (.)^3 gives saturated model up to 3-way factors
  lm( y ~ (.)^3, data = design)
) 

# 6.4.d 2^(8-4) design: 2^4 = 16-run, k = 8-factors, p = 4 
# estimated effects
library(FrF2)
design <- FrF2( 16, 8, randomize = FALSE, seed = 5) # min abbr design
library(DoE.base)
generators(design) 
design
flight.time <- c(
  5.2, 4.9, 5.0, 5.0, 4.9, 5.0, 5.1, 4.8,
  4.9, 5.1, 5.1, 5.1, 5.2, 5.0, 5.1, 5.0); flight.time
library(DoE.base)
helicopter <- add.response( design , flight.time ); helicopter
m.helicopter.84 <- lm( flight.time ~ (.)^3 , data = helicopter) # up to 3-ways
summary(m.helicopter.84)

# 6.4.e,f,g 2^(8-4) design: 2^4 = 16-run, k = 8-factors, p = 4 
# normal plot
library(FrF2)
design <- FrF2( 16, 8, randomize = FALSE, seed = 5) # min abbr design
library(DoE.base)
generators(design) 
design
flight.time <- c(
  5.2, 4.9, 5.0, 5.0, 4.9, 5.0, 5.1, 4.8,
  4.9, 5.1, 5.1, 5.1, 5.2, 5.0, 5.1, 5.0); flight.time
library(DoE.base)
helicopter <- add.response( design , flight.time ); helicopter
m.helicopter.84 <- lm( 
  flight.time ~ 1 + A + B + C + D + E + F + G + H + A:B + A:C + 
    A:D + A:E + A:F + A:G + A:H, data = helicopter) # up to 3-ways
summary(m.helicopter.84)
library(daewr)
LGB(coef(m.helicopter.84)[-1], alpha = .1) # plot and report
halfnorm(
  coef(m.helicopter.84)[-1], 
  names, 
  alpha = 0.4, 
  refline = TRUE
) 

# 6.6.a 2^(7-3)_III design: 2^3 = 8-run, k = 7-factors, p = 3 
library(FrF2)
sand.filter <- FrF2( 
  8, 7, 
  randomize = FALSE) # gives min.aberr 
sand.filter # design only, no data
augm <- fold.design(sand.filter); augm # augmented mirror design
# data frame for 2^(7-3)_III design with mirror and percentage.arsenic
A <- (as.numeric( augm$A) - 1.5 ) / .5
B <- (as.numeric( augm$B) - 1.5 ) / .5
C <- (as.numeric( augm$C) - 1.5 ) / .5
D <- (as.numeric( augm$D) - 1.5 ) / .5
E <- (as.numeric( augm$E) - 1.5 ) / .5
F <- (as.numeric( augm$F) - 1.5 ) / .5
G <- (as.numeric( augm$G) - 1.5 ) / .5
Block <- augm$fold
percent.arsenic <- c(69.95, 58.65, 56.25, 53.25, 
    94.40, 73.45, 10.00, 2.11,
    16.20, 52.85, 9.05, 31.1, 
    7.40, 9.90, 10.85, 48.75); percent.arsenic
augmn <- data.frame(A, B ,C, D, E, F, G, Block, percent.arsenic); augmn

# 6.6.b 2^(7-3)_III design: 2^3 = 8-run, k = 7-factors, p = 3 
library(FrF2)
sand.filter <- FrF2( 
  8, 7, 
  randomize = FALSE) # gives min.aberr 
sand.filter # design only, no data
library(DoE.base)
generators(sand.filter) 
y <- runif(8, 0, 1) # "dummy" values for response
aliases( # (.)^5 gives saturated model up to 5-way factors
  lm( y ~ (.)^5, data = sand.filter)
) 

# 6.6.c 2^(7-3)_III design: 2^3 = 8-run, k = 7-factors, p = 3 
# estimated effects and half-normal plot original data
library(FrF2)
sand.filter <- FrF2( 
  8, 7, 
  randomize = FALSE) # gives min.aberr 
sand.filter # design only, no data
library(DoE.base)
percent.arsenic <- c(69.95, 58.65, 56.25, 53.25, 
    94.40, 73.45, 10.00, 2.11); percent.arsenic
library(DoE.base)
sand.filter <- add.response( sand.filter , percent.arsenic ); sand.filter
m.sand.filter.87 <- lm( 
  percent.arsenic ~ 1 + A + B + C + D + E + F + G, 
  data = sand.filter) # up to 1-ways
summary(m.sand.filter.87)
library(daewr)
LGB(coef(m.sand.filter.87)[-1], alpha = .1) # plot and report
halfnorm(
  coef(m.sand.filter.87)[-1], 
  names, 
  alpha = 0.10, 
  refline = TRUE
) 

# 6.6.d 2^(7-3)_III design: 2^3 = 8-run, k = 7-factors, p = 3 
# defining relation of original + mirror
library(FrF2)
sand.filter <- FrF2( 
  8, 7, 
  randomize = FALSE) # gives min.aberr 
sand.filter # design only, no data
augm <- fold.design(sand.filter); augm # augmented mirror design
augm$fold <- NULL # remove fold variable 
library(DoE.base)
generators(sand.filter) 
y <- runif(16, 0, 1) # "dummy" values for response
aliases( # (.)^7 gives saturated model up to 7-way factors
  lm( y ~ (.)^7, data = augm)
) 

# 6.6.e 2^(7-3)_III design: 2^3 = 8-run, k = 7-factors, p = 3 
# estimated effects half-normal original + mirror + block
library(FrF2)
sand.filter <- FrF2( 
  8, 7, 
  randomize = FALSE) # gives min.aberr 
sand.filter # design only, no data
augm <- fold.design(sand.filter); augm # augmented mirror design
A <- (as.numeric( augm$A) - 1.5 ) / .5
B <- (as.numeric( augm$B) - 1.5 ) / .5
C <- (as.numeric( augm$C) - 1.5 ) / .5
D <- (as.numeric( augm$D) - 1.5 ) / .5
E <- (as.numeric( augm$E) - 1.5 ) / .5
F <- (as.numeric( augm$F) - 1.5 ) / .5
G <- (as.numeric( augm$G) - 1.5 ) / .5
Block <- augm$fold
# block <- c(rep(-1,8),rep(1,8)); block
percent.arsenic <- c(69.95, 58.65, 56.25, 53.25, 
                     94.40, 73.45, 10.00, 2.11,
                     16.20, 52.85, 9.05, 31.1, 
                     7.40, 9.90, 10.85, 48.75); percent.arsenic
augmn <- data.frame(A, B ,C, D, E, F, G, Block, percent.arsenic); augmn
m.sand.filter.87.mir <- lm( 
  percent.arsenic ~ 1 + Block + A + B + C + D + E + F + G +  
    A:B + A:C + A:D + A:E + A:F + A:G, 
  data = augmn) 
summary(m.sand.filter.87.mir)
library(daewr)
LGB(coef(m.sand.filter.87.mir)[-1], alpha = .1) # plot and report
halfnorm(
  coef(m.sand.filter.87.mir)[-1], 
  names, 
  alpha = 0.1, 
  refline = TRUE
) 

# 6.8.a,c 2^(5-2) design: 2^3 = 8-run, k = 5-factors, p = 2 
# estimated effects and half-normal plot original data
library(FrF2)
prince <- FrF2( 
  8, 5, 
  generators = c("AC", "BC"),
  randomize = FALSE) # generators D = AC, E = BC
prince # design only, no data
library(DoE.base)
assay.response <- c(1.31091, 1.43201, 1.29951, 1.37199, 1.33566, 1.46820, 1.39023, 1.41531); assay.response
library(DoE.base)
prince <- add.response( prince , assay.response ); prince
# first check for where estimated effects are:
# summary(lm( assay.response ~ (.)^5, data = prince))
m.prince <- lm( 
  assay.response ~ 1 + A + B + C + D + E + A:B + A:E, 
  data = prince) 
summary(m.prince)
library(daewr)
LGB(coef(m.prince)[-1], alpha = .1) # plot and report
halfnorm(
  coef(m.prince)[-1], 
  names, 
  alpha = 0.3, 
  refline = TRUE
) 

# 6.8.b,c 2^(5-2) design: 2^3 = 8-run, k = 5-factors, p = 2 
# estimated effects original data
library(FrF2)
prince <- FrF2( 
  8, 5, 
  generators = c("AC", "BC"),
  randomize = FALSE) # based on generators AC, BC 
prince # design only, no data
library(DoE.base)
generators(prince) 
y <- runif(8, 0, 1) # "dummy" values for response
aliases( # (.)^5 gives saturated model up to 5-way factors
  lm( y ~ (.)^5, data = prince)
) 

# 6.8.d 2^(5-2) design: 2^3 = 8-run, k = 5-factors, p = 2 
# estimated effects half-normal original + mirror + block
library(FrF2)
prince <- FrF2( 
  8, 5, 
  generators = c("AC", "BC"),
  block.name = c("Block"),
  blocks = c("original", "mirror"),
  randomize = FALSE) # based on generators AC, BC 
prince # design only, no data
augm <- fold.design(prince); augm # augmented mirror design
Block <- augmn$fold
# block <- c(rep(-1,8),rep(1,8)); block
assay.response <- c(1.31091, 1.43201, 1.29951, 1.37199, 
    1.33566, 1.46820, 1.39023, 1.41531,
    1.31702, 1.38881, 1.32222, 1.36248, 
    1.33826, 1.32654, 1.32654, 1.34635); assay.response
augmn <- data.frame(A, B , C, D, E, Block, assay.response); augmn
# first check for where estimated effects are:
summary(lm( assay.response ~ 1 + Block + A*B*C*D*E, data = augmn))
# summary(lm(assay.response ~ (.)^5, data = augmn))
m.prince.mir <- lm( 
  assay.response ~ 1 + Block + A + B + C + D + E + 
    A:B + A:C + B:C + A:D + B:D + C:D + C:E + 
    A:B:C + B:C:D,
  data = augmn) 
summary(m.prince.mir)

# 6.8.e 2^(5-2) design: 2^3 = 8-run, k = 5-factors, p = 2 
# defining relationship for original + mirror 
library(FrF2)
prince <- FrF2(
  8, 5,
  generators = c("AC", "BC"),
  randomize = FALSE) # based on generators AC, BC
prince # design only, no data
augm <- fold.design(prince); augm # augmented mirror design
Block <- augm$fold; augm; 
colnames(augm) <- c("A", "B","C","Block","D","E")
library(DoE.base)
y <- runif(16, 0, 1) # "dummy" values for response
aliases( # (.)^5 gives saturated model up to 5-way factors
  lm( y ~ (.)^5, data = augm)
)

# 6.8.f 2^(5-2) design: 2^3 = 8-run, k = 5-factors, p = 2 
# half-normal plot original + mirror 
library(daewr)
LGB(coef(m.prince.mir)[-1], alpha = .1) # plot and report
halfnorm(
  coef(m.prince.mir)[-1], 
  names, 
  alpha = 0.1, 
  refline = TRUE
) 
