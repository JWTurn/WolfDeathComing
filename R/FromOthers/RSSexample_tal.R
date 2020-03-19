
#here is a rather comples iSSF. 'habitat' is a factor with multiple levels.
mymodel <- clogit(case ~ step_length + I(log(step_length)) + I(log(step_length)):time2death +
                    I(cos(turn_angl)) + I(cos(turn_angl)):time2death +
                    I(cos(turn_angl)):I(log(step_length)) +
                    I(log(dist2road+1)) + I(log(dist2road+1)):time2death +
                    habitat + habitat:time2death +
                    habitat:I(log(step_length)) + habitat:I(cos(turn_angl)) +
                    strata(cluster_ID), data = mydata)

#now let's set up 100 potentiale step 1's. The choice of values is up to you - these are examples
#in this example I left everything constant but distance to road
x_1 <- data.frame(step_length = mean(mydata$step_length),
                  turn_angl = acos(mean(cos(mydata$turn_angl))),
                  time2death = mean(mydata$time2death),
                  dist2road = seq(from = 0, to = max(mydata$dist2road), length.out = 100),
                  habitat = factor('forest', levels = levels(mydata$habitat)))
#if you cannot apply predict to your model, you would have to spell it out here,
#but only the linear part (the addative combination of covariates and parametrs, without exp())
log.predicted_1 <- predict(mymodel, newdata = x_1) #type = "link" should be the default
#now let's set up the alternative (reference) step. There is only one alternative step
x_2 <- data.frame(step_length = mean(mydata$step_length),
                  turn_angl = acos(mean(cos(mydata$turn_angl))),
                  time2death = mean(mydata$time2death),
                  dist2road = mean(mydata$dist2road),
                  habitat = factor('forest', levels = levels(mydata$habitat)))
log.predicted_2 <- predict(mymodel, newdata = x_2)
#and here are 100 logRSS values as fucntion of dist2road
logRSS <- log.predicted_1 - log.predicted_2
plot(logRSS ~ x_1$dist2road, col = 'red')

#you can now repeat this proccess for differnt conditions. for example:
x_1 <- data.frame(step_length = mean(mydata$step_length),
                  turn_angl = acos(mean(cos(mydata$turn_angl))),
                  time2death = 0,
                  dist2road = seq(from = 0, to = max(mydata$dist2road), length.out = 100),
                  habitat = factor('forest', levels = levels(mydata$habitat)))
log.predicted_1 <- predict(mymodel, newdata = x_1) #type = "link" should be the default
x_2 <- data.frame(step_length = mean(mydata$step_length),
                  turn_angl = acos(mean(cos(mydata$turn_angl))),
                  time2death = 0,
                  dist2road = mean(mydata$dist2road),
                  habitat = factor('forest', levels = levels(mydata$habitat)))
log.predicted_2 <- predict(mymodel, newdata = x_2)
logRSS <- log.predicted_1 - log.predicted_2
points(logRSS ~ x_1$dist2road, col = 'blue')
