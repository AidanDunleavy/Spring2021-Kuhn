library(rsm) 
ccd.des <- ccd(defects~x1+x2, n0=c(3,2), alpha="rotatable", randomize=FALSE)
ccd.des

des <- as.data.frame(ccd.des); des # convert to data.frame
set.seed(1234)
des$defects <- 30*rgamma(13, shape = 2, rate = 1.1); des # add defects response to data.frame
des
library(Vdgraph)
#+ fig.width=5, fig.height=5
Vdgraph(des[ , 3:4])



library(rsm)
ccd.des <- ccd(defects~x1+x2, n0=c(3,2), alpha="faces", randomize=FALSE)
ccd.des

des <- as.data.frame(ccd.des); des # convert to data.frame
set.seed(1234)
des$defects <- 30*rgamma(13, shape = 2, rate = 1.1); des # add defects response to data.frame


library(Vdgraph)
#+ fig.width=5, fig.height=5
Vdgraph(des[ , 3:4])

####################
# PRoblem 5
#####################


library(rsm) 
# n0 = c(#center points, #axial points)
ccd.des <- ccd(defects~x1+x2+x3, n0 = c(7,6), alpha="faces", randomize=FALSE)
ccd.des


# 
# SOdes2 <- ccd (3, n0 = c(4,6), alpha = "rotatable", inscribed = TRUE, coding = list (
#   x1 ~ (Temp - 150)/10, x2 ~ (Pres - 50)/5, x3 ~ Feedrate - 4))
# SOdes2




library(rsm)
bbd.des <- bbd(defects ~ x1 + x2 + x3, n0 = 15, randomize = FALSE)
bbd.des


bbd.act <- bbd(defects ~ x1 + x2 + x3,
               n0 = 5, # 5 center points
               coding=list(x1~(temp-15)/15, # actual levels
                           x2~(click-90)/30,
                           x3~(hard-0.65)/0.15),randomize = FALSE)
bbd.act
act <- as.coded.data(bbd.act); act # save as coded data
act$defects <- c(61, 69, 60, 45, 64, 66, 44, 57,
                 60, 61, 52, 41, 56, 61, 51, 63,
                 60); act # add defects response to coded data



# 
# library(AlgDesign)
# nonlin.desD <- optFederov(~ -1 + dfdT + dfdk + dfdc0,
#                           data=grid,
#                           nTrials=3,
#                           center=TRUE,
#                           criterion="D",
#                           nRepeats=20)
# nonlin.desD$design



library(AlgDesign)
d <- data.frame(var=paste("x",1:3,sep=""),  
                low=-1,
                high=1,
                center=0,
                nLevels=10,
                round=1,
                factor=FALSE); d
constFcn <- function(x){x[2] >= -2*x[1] - 2 && x[2] <= -2*x[1] + 1}
desCon <- optMonteCarlo(~-1+x1+x2+x3, 
                        d, 
                        constraint=constFcn,
                        criterion = "D"); desCon

desCon

# x1=c(-1,1,-1,1,-1,1,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
# x2=c(-1,-1,1,1,0,0,0,0,-1,1,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
# x3=c(0,0,0,0,-1,-1,1,1-1,-1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
x1 = c(ccd.des$x1,0)
x2 = c(ccd.des$x2,0)
x3 = c(ccd.des$x3,0)
library(AlgDesign)
data<-cbind(x1,x2,x3)
nonlin.desD <- optFederov(~ -1 + x1 + x2 + x3,
                          data=data,
                          nTrials=10,
                          center=TRUE, 
                          criterion="D",
                          nRepeats=1)
nonlin.desD$design


desCon
descon <- as.data.frame(desCon)
descon
# (d)

library(Vdgraph)
#+ fig.width=5, fig.height=5
Vdgraph(nonlin.desD$design)

Vdgraph(descon[,5:7])

Vdgraph(bbd.des[, 3:5])

Vdgraph(ccd.des[, 3:5])




# part e

cat <- expand.grid(c(1,2,3), c(1,2,3), c(1,2,3))
cat <- expand.grid(c(-1,0,1), c(-1,0,1), c(-1,0,1))
set.seed(1234)
y <- distance <- 4+3*runif(27); y

catapult <- data.frame(cat, y)
names(catapult) <- c("start", "stop", "pivot", "dist")
catapult

# part f

library(rsm)
catPQ <- rsm(dist ~ FO(start, stop, pivot)+PQ(start, stop, pivot), data=catapult)
summary(catPQ)

contour(catPQ, ~ start + stop + pivot)

max(y)
which(y==max(y))



