# 7.2.b pp 303--305
# BIB design: minimum number of blocks for BIB design
library(daewr)
BIBsize(12,3); BIBsize(12,4); BIBsize(12,5); BIBsize(12,6)

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