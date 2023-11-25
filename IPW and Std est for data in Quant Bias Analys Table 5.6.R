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


#IP weightening
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

##standarization
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
