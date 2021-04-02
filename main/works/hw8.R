data = read.table("https://raw.githubusercontent.com/AidanDunleavy/Spring2021-Kuhn/main/main/lib/8-8-5-c.txt", header = TRUE)

R <- factor(c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), 
              rep(5,2), rep(6,2), rep(7,2), rep(8,2))); R
Block <- factor(rep(R,3)); Block
FilmType <- factor(rep(c(1,2), times = 8)); FilmType
P <- factor(rep(c(1,2,3), each = 16)); P
# Pressure subplot factor, 3 levels
film <- c(data$P1, data$P2, data$P3)
film.data <- data.frame(Block,FilmType,P,film); film.data
library(lme4)
spmodel <- lmer( film ~ (1|Block) + FilmType + (1|Block:FilmType) + P + FilmType:P, data = film.data)
summary(spmodel) # summary random effects
anova(spmodel) # summary fixed effects
spmodel <- aov(film ~ Block*FilmType + P + FilmType:P + Error(Block:FilmType/P), data = film.data ) 
summary(spmodel, type = 3) # error degrees of freedom for FP and site factors
1 - pf(3.3634,1,7) # p-value for whole plot FilmType factor
1 - pf(4.2261,2,28) # p-value for subplot site factor


par( mfrow = c(1,1) )
with(film.data,
     (interaction.plot(
       FilmType,
       P,
       film,
       type = "b",
       pch = c(18,22,6),
       leg.bty = "o",
       main = "Film Type by Pressure Interaction Plot",
       xlab = "Film Type", 
       ylab = "Quality"))
)
' this is my answer'


spmodel2 <- lmer( film ~ (1|Block) + FilmType + (1|Block:FilmType) + P + FilmType:P, data = film.data)
res <- resid(spmodel2)
plot(fitted(spmodel2), res,
     ylab="Residual", 
     xlab="Predicted Value", 
     main="Residuals by Predicted for Y")
abline(0,0)

# Normal Plot (Look for Linear Relationship)
qqnorm(res, 
       ylab="Standardized Residuals", 
       xlab="Normal Scores", 
       main="QQ Plot of Residuals for Y")
abline(0,1/5)



spmodel2 <- lmer( film ~ (1|Block) + FilmType + (1|Block:FilmType) + P + FilmType:P, data = film.data)
res <- resid(spmodel2)
plot(fitted(spmodel2), res,
     ylab="Residual", 
     xlab="Predicted Value", 
     main="Residuals by Predicted for Y")
abline(0,0)