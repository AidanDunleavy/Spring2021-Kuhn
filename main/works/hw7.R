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
Blockdes <- gen.factorial(2, nVars = 3, factors = "all", # 2^3 design
                          varNames = c("A","B","C"))
Blockdes <- optBlock( ~ A + B + C + A:B + A:C + B:C + A:B:C, withinData = Blockdes,
                      blocksizes = c(rep(4,14)), criterion = "D")
Blockdes

###########
des <- gen.factorial(levels = c(2, 2, 2), # three factors levels 2, 2, 2 
                     factors = 'all',
                     varNames = c("A", "B", "C"))
bdes <- optBlock( ~ A + B + C + A:B + A:C + B:C, # one/two-factor estimable
                  withinData = des, 
                  blocksizes = c(rep(4, 14)), # create 14 blocks of size 4 each
                  criterion = "D") # use D efficient method
bdes
bdes$Blocks
Block <- as.factor(c(rep(1:14, each = 4))) # blocks
set.seed(1234)
roc <- runif(56, min=5, max=9)

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

#########################
library(daewr)
BIBsize(10,4)

##########################################
# 7.1.b

library(AlgDesign)
BIB <- optBlock( ~ ., withinData = factor(1:10), blocksizes = rep(4, 15))
des <- BIB$rows
dim(des) <- NULL
des <- matrix(des, nrow = 15, ncol = 4, byrow = TRUE,
      dimnames = list(c( "Block1", "Block2", "Block3", "Block4",
                         "Block5", "Block6", "Block7", "Block8", "Block9",
                         "Block10", "Block11", "Block12", "Block13", "Block14", "Block15"), c("unit1", "unit2", "unit3", "unit4")))
des



##############################
# 7.5.a / make a 2^6 design with 

library(DoE.base)
show.oas(factors=list(nlevels=c(2),number=c(6)))

des <- oa.design(nlevels = c(2, 2, 2, 2, 2, 2),nruns = 8,
                 columns = "min3", randomize = TRUE, seed = 104)
des


library(DoE.base)
library(AlgDesign)
fnames <- c("A", "B", "C", "D", "E", "F")
cand <- oa.design( nlevels = c(2, 2, 2, 2, 3, 3),
                   factor.names = fnames, 
                     nruns = 24, seed = 104)
cand
bdes <- optBlock( ~ A + B + C + D + E + F + A:B + A:C + A:D + A:E + A:F + B:C + B:D + B:E + B:F + C:D + C:E + C:F + D:E + D:F + E:F ,
                  withinData = cand, 
                  blocksizes = c(rep(8, 3)),
                  criterion = "D")
bdes




library(FrF2)
design <- FrF2( 16, 6, blocks = c("AB","BD","EF"),
                alias.block.2fis = TRUE, randomize = FALSE)
design



library(AlgDesign)
Blockdes <- gen.factorial(3, nVars = 3, factors = "all",
                          varNames = c("A", "B", "C"))
Blockdes <- optBlock( ~ A + B + C + A:B + A:C + B:C, 
                      withinData = Blockdes,
                      blocksizes = c(rep(5, 5)), 
                      criterion = "D")
Blockdes


library(DoE.base) # 8 blocks of size 3, but sacrifice orthogonality
library(AlgDesign)
fnames <- c("A", "B", "C", "D", "E", "F")
cand <- oa.design(nlevels = c(2, 2, 2, 2, 3, 3),
                  factor.names = fnames, 
                  randomize = TRUE, 
                  seed = 104,
                  nruns = 24)
bdes <- optBlock( ~ A + B + C + D + E + F,
                  withinData = cand, 
                  blocksizes = c(rep(3, 8)), # 8 blocks of size 3
                  criterion = "D")
bdes
