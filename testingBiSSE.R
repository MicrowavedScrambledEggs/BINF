trstr = "(((chimp:1,human:1):1,gorilla:2):1,orang:3);"
tr = read.tree(file="", text=trstr)
endStates <- c(1,1,1,1)
names(endStates) <- tr$tip.label
bisse_func <- make.bisse(tr, endStates, sampling.f = c(1,1), strict = FALSE)

compareBiSSEFunctions <- function(b0, b1, u0, u1, q01, q10, phy, tipStates)
{
  my.bisse.out <- biSSE.likelyhood(b0,b1,u0,u1,q01,q10,phy,tipStates, 
                                   output.branches = TRUE)
  # Compare to diversitree's function
  parms <- c(b0,b1,u0,u1,q01,q10)
  names(parms) = c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")
  diversitree.out <- bisse_func(pars=parms, root=ROOT.FLAT, root.p=NULL, 
                                intermediates=TRUE, condition.surv=FALSE)
  cat(paste("log likelyhood:\n\tMy function:", my.bisse.out[1], 
            "\tdiversitree function:", diversitree.out[1], "\n\n"))
  cat("Relative probabilities at tip end of branches:\n\tMy function\n")
  print(attributes(my.bisse.out)$init)
  cat("\tdiversitree function:\n")
  print(attributes(diversitree.out)$intermediates$init)
  cat("\nRelative probabilities at root end of branches:\n\tMy function\n")
  print(attributes(my.bisse.out)$base)
  cat("\tdiversitree function:\n")
  print(attributes(diversitree.out)$intermediates$base)
  cat("\n\n")
  return(list(my.bisse.out, diversitree.out))
}

brate1 <- 0.35
brate2 <- 0.4
drate1 <- 0.1
drate2 <- 0.2
qrate1 <- 0.3
qrate2 <- 0.2

# Single birthrate, no death rate, no state transition 
yuleBiSSE <- compareBiSSEFunctions(brate1, brate1,0,0,0,0,tr,endStates)

# Single birthrate, single death rate, no state transition 
birthdeath <- compareBiSSEFunctions(brate1, brate1, drate1, drate1, 0, 0, 
                                    tr, endStates)

# Single birthrate, single death rate, single state transition 
stateIndependent <- compareBiSSEFunctions(brate1, brate1, drate1, drate1, 
                                          qrate1, qrate1, tr, endStates)

# Single birthrate, no death rate, Single state transition 
yuleStates <- compareBiSSEFunctions(brate1, brate1, 0, 0, 
                                    qrate1, qrate1, tr, endStates)

# varying everything
stateDependent <- compareBiSSEFunctions(brate1, brate2, drate1, drate2, 
                                        qrate1, qrate2, tr, endStates)


# Testing with a random tree
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(2)
phy <- tree.bisse(pars, max.t = 50, x0 = 0)

lik <- make.bisse(phy, phy$tip.state, strict = FALSE)
diversitree.out <- lik(pars, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, 
                       condition.surv=FALSE)
my.bisse.out <- biSSE.likelyhood(pars[1], pars[2], pars[3], pars[4], pars[5], 
                                 pars[6], phy, phy$tip.state)
cat(paste("log likelyhood:\n\tMy function:", my.bisse.out[1], 
          "\tdiversitree function:", diversitree.out[1], "\n\n"))
