simu <- data.frame(D = rep(NA,0),E = rep(NA,0),C=rep(NA,0))
simu <- rbind(simu,data.frame(D = rep(1,646),E = rep(1,646),C = rep(1,646))) #1-1-1
simu <- rbind(simu,data.frame(D = rep(1,30),E = rep(0,30),C = rep(1,30))) #1-0-1
simu <- rbind(simu,data.frame(D = rep(0,4410),E = rep(1,4410),C = rep(1,4410))) #0-1-1  
simu <- rbind(simu,data.frame(D = rep(0,59),E = rep(0,59),C = rep(1,59))) #0-0-1
simu <- rbind(simu,data.frame(D = rep(1,404),E = rep(1,404),C = rep(0,404))) #1-1-0
simu <- rbind(simu,data.frame(D = rep(1,820),E = rep(0,820),C = rep(0,820))) #1-0-0
simu <- rbind(simu,data.frame(D = rep(0,860),E = rep(1,860),C = rep(0,860))) #0-1-0
simu <- rbind(simu,data.frame(D = rep(0,871),E = rep(0,871),C = rep(0,871))) #0-0-0

#crude RR
RR_crude <- with(simu,
                 (sum(D == 1 & E == 1)/sum(E == 1))/(sum(D == 1 & E == 0)/sum(E == 0)))
RR_C0_crude <- with(simu[simu$C==0,],
                    (sum(D == 1 & E == 1)/sum(E == 1))/(sum(D == 1 & E == 0)/sum(E == 0)))
RR_C1_crude <-with(simu[simu$C==1,],
                    (sum(D == 1 & E == 1)/sum(E == 1))/(sum(D == 1 & E == 0)/sum(E == 0)))

#IP weighting
#calc Pr(A)
swmodel <- glm(E~C,family = binomial(),data=simu)
E.w <- predict(swmodel,type = "response")

num <- sum(simu$E==1)/nrow(simu)

simu$sw <- ifelse(simu$E==1,num/E.w,(1-num)/(1-E.w))

#binomial + log link (then exp(beta) ~ RR)
ipmodel <- glm(D~E,family = binomial(link = "log"),data = simu,weights = sw)
summary(ipmodel)

#poisson + log link (then exp(beta) ~ RR)
ipmodelps <- glm(D~E,family = poisson(link = "log"),data = simu,weights = sw)
summary(ipmodelps)

#NB
ipmodelnb <- MASS::glm.nb(D~E,data = simu,weights = sw)
summary(ipmodelnb)

RR_adj_ipw <- c(exp(coefficients(ipmodel)[2]),
            exp(coefficients(ipmodelps)[2]),
            exp(coefficients(ipmodelnb)[2]))

##standardization
simu_std <- rbind(simu,simu,simu)
simu_std$E[8101:(8100*2)] <- 0
simu_std$E[(8100*2+1):(8100*3)] <- 1
simu_std$D[8101:(8100*2)] <- NA
simu_std$D[(8100*2+1):(8100*3)] <- NA
simu_std$status <- c(rep(-1,8100),rep(0,8100),rep(1,8100))

stdmodel <- glm(D~ E*C,family = binomial(link = "log"),data = simu_std,na.action = na.omit)
simu_std$pred1 <- predict(stdmodel,newdata = simu_std,type = "response")

stdmodelps <- glm(D ~ E*C,family = poisson(link = "log"),data = simu_std)
simu_std$pred2 <- predict(stdmodelps,newdata = simu_std,type = "response")

stdmodelnb <- MASS::glm.nb(D ~ E*C,data = simu_std)
simu_std$pred3 <- predict(stdmodelnb,newdata = simu_std,type = "response")

RR_adj1 <- sum(simu_std$pred1[simu_std$status==1])/sum(simu_std$pred1[simu_std$status==0])
RR_adj2 <- sum(simu_std$pred2[simu_std$status==1])/sum(simu_std$pred2[simu_std$status==0])
RR_adj3 <- sum(simu_std$pred3[simu_std$status==1])/sum(simu_std$pred3[simu_std$status==0])
RR_adj_std <- c(RR_adj1,RR_adj2,RR_adj3)

#estimate crude RR using regression

crudemodel <- glm(D~E,family = binomial(link = "log"),data = simu)
crudemodelps <- glm(D~E,family = poisson(link = "log"),data = simu)
crudemodelnb <- MASS::glm.nb(D~E,data = simu)

RR_crude_reg <- c(
  exp(coefficients(crudemodel)[2]),
  exp(coefficients(crudemodelps)[2]),
  exp(coefficients(crudemodelnb)[2])
)
RR_crude_reg

#Estimate "Fox"-adjusted RR. Conditioning on distribution of C within E=1. pseudo-Standardization used.
#Construct two pseudo datasets with approx. corrected Pr(C|E=1)
Pr_C_E1 <- nrow(simu[simu$E == 1 & simu$C == 1,])/nrow(simu[simu$E==1,])

simu_pseudo_treated <- data.frame(
  D = rep(NA,8100),
  E = rep(1,8100),
  C = c(rep(1,round(8100*Pr_C_E1)),rep(0,(8100-round(8100*Pr_C_E1))))
)
simu_pseudo_untreated <- data.frame(
  D = rep(NA,8100),
  E = rep(0,8100),
  C = c(rep(1,round(8100*Pr_C_E1)),rep(0,(8100-round(8100*Pr_C_E1))))
)

#Using previous regression model to calculate "conditional standardized" RR
simu_pseudo_treated$pred1 <- predict(stdmodel,newdata = simu_pseudo_treated,type = "response")
simu_pseudo_treated$pred2 <- predict(stdmodelps,newdata = simu_pseudo_treated,type = "response")
simu_pseudo_treated$pred3 <- predict(stdmodelnb,newdata = simu_pseudo_treated,type = "response")

simu_pseudo_untreated$pred1 <- predict(stdmodel,newdata = simu_pseudo_untreated,type = "response")
simu_pseudo_untreated$pred2 <- predict(stdmodelps,newdata = simu_pseudo_untreated,type = "response")
simu_pseudo_untreated$pred3 <- predict(stdmodelnb,newdata = simu_pseudo_untreated,type = "response")

RR_condstd <- c(
  sum(simu_pseudo_treated$pred1)/sum(simu_pseudo_untreated$pred1),
  sum(simu_pseudo_treated$pred2)/sum(simu_pseudo_untreated$pred2),
  sum(simu_pseudo_treated$pred3)/sum(simu_pseudo_untreated$pred3)
)
RR_condstd
## Here RR_condstd (standardized conditioning on the distribution of C|E=1) precisely equals to 
## What reported in the book.
