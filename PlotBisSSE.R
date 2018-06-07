pars <- c(0.1, 0.2, 0.04, 0.03, 0.02, 0.01)
set.seed(5)
phy <- tree.bisse(pars, max.t = 45, x0 = 0)

brate2 <- seq(0.05, 0.5, 0.01)

lnlikely <- sapply(brate2, function(x) biSSE.likelyhood(pars[1], x, pars[3], 
                                pars[4], pars[5], pars[6], phy, phy$tip.state))

plot(brate2, lnlikely, xlab = "Birth rate while in state 1", 
     ylab = "Log Likelihood")
title("BiSSE: Log likelihood vs lambda1", 
      sub = paste0("lambda0 = ", pars[1], "; mu0 =", pars[3], "; mu1 = ", 
                   pars[4], "; q01 = ", pars[5], "; q10 = ", pars[6]))

likely <- exp(lnlikely)
plot(brate2, likely, xlab = "Birth rate while in state 1", 
     ylab = "Likelihood")
title("BiSSE: Likelihood vs lambda1", 
      sub = paste0("lambda0 = ", pars[1], "; mu0 =", pars[3], "; mu1 = ", 
                   pars[4], "; q01 = ", pars[5], "; q10 = ", pars[6]))
