#' ---
#' title: HW 2
#' author: Aidan Dunleavy
#' date: 01/26/2021
#' ---
#' PROBLEM 1
#' 
#' (a)
D <- expand.grid(Beat = X3_10_1$Beat[1:3],
                 Course = X3_10_1$Course[1:3])
D[] <- lapply(D, factor)
Score = X3_10_1$Score
ScoreData = data.frame(cbind(D,Score))
ScoreData
with(ScoreData, 
     (interaction.plot(
       Course, 
       Beat, 
       Score, 
       type = "b", 
       pch = c(18,24,22),
       leg.bty = "o",
       main = "Interaction Plot of Course length and Beat",
       xlab = "Course Length",
       ylab = "Score"))
)
help("interaction.plot")
