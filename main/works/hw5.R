######################
#Problem 1.a
#####################

#2^(3-1) design
library(FrF2)
design <- FrF2( 16, 7, generators = c("ABC", "ABD", "ACD"), randomize = FALSE)
design


######################
#Problem 1.b
#####################

library(FrF2)
y <- runif(16, 0, 1)
ali = aliases( lm( y~ (.)^4, data = design)) # the 2 is equal to 3-1
ali


# 2^3 = 8-run, k = 5-factors, p = 2 generators D = AB and E = AC
library(FrF2)
design <- FrF2( 16, 7, 
                generators = c("ABC", "ABD", "ACD"),
                randomize = FALSE)
design
y <- runif(16, 0, 1) # dummary y values, to allow aliases to work 
aliases( lm( y ~ (.)^7, data = design)) # up to 5-way factors
designf <- fold.design(design, columns = c(2)) # foldover A
designf 
designf$fold <- NULL # drop "fold" factor
designf
y <- runif(32, 0, 1) # dummary y values, to allow aliases to work 
aliases( lm( y ~ (.)^7, data = designf)) # up to 5-way factors

design.mirror <- fold.design(design, columns = 'full') # mirror design
design.mirror 
design.mirror$fold <- NULL # drop "fold" factor
design.mirror
y <- runif(32, 0, 1) # dummary y values, to allow aliases to work 
aliases( lm( y ~ (.)^4, data = design.mirror)) # up to 5-way factors
