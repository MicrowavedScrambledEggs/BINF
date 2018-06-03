library(deSolve)
library(diversitree)

# BiSSE likelyhoodfunction
# Takes birth rate and death rate for each character state (b0, b1, u0, u1), a q matrix
# for the characer transition rates, the tree, and a vector of states for the tree tips
# The vector of states has to have a names attribue matching the names in phy$tip.label to
# states
biSSE.likelyhood <- function(b0, b1, u0, u1, q01, q10, phy, endstates)
{
  # Uses the differentials of the probabilities each event could happen between a small 
  # time period
  
  # D_n1(t) : The likelyhood that an extant species has the character state 1 at time t
  # Starting state: if tip has state 1 D_n1(0) = 1; D_n0(0) = 0
  #                 if tip has state 0 D_n1(0) = 0; D_n0(0) = 1
  # E_0(t) : The likelyhood that a species goes extinct at time t while in state 0
  # The likelyhood that something stays at state 0 on a branch between time t and
  # time t + small_t and does not go extinct
  # D_n0(t+small_t) = p(did not go extinct)
  #                   * ( 
  #                         p(No state change and no speciation)
  #                       + p(State change and no speciation)
  #                       + p(No state change and speication then extinction)
  #                       + p(No state change and speication then extinction)
  #                     )
  # We don't include p(state change and speciation then extinction) because small_t is too
  # small for the probability of both occuring to be significant.
  # Any speciation events have to be followed by extinction because it would mean a branch
  # that is not observed in the given tree. We include speciation then extinction twice as
  # extinction could have ouccured in either of the two decendants from a speciation.
  # D_n0(t+small_t) = (1 - u0*small_t)
  #                   * ( 
  #                         (1 -  q[0,1]*small_t) * (1 - b0*small_t) * D_n0(t)
  #                       + q[0,1] * small_t * (1 - b0*small_t) * D_n1(t)
  #                       + (1 -  q[0,1]*small_t) * b0 * small_t * E_0(t) * D_n0(t)
  #                       + (1 -  q[0,1]*small_t) * b0 * small_t * E_0(t) * D_n0(t) 
  #                     )
  # Because terms of order small_t^2 are negligably small we can drop them, giving:
  # D_n0(t+small_t) = (1 - (b0 + u0 + q[0,1])*small_t) * D_n0(t) + q[0,1] * small_t 
  #                   * D_n1(t) + 2 * b0 * small_t * E_0(t) * D_n0(t) 
  #
  # Likewise:
  # D_n1(t+small_t) = (1 - (b1 + u1 + q[1,0])*small_t) * D_n1(t) + q[1,0] * small_t 
  #                   * D_n0(t) + 2 * b1 * small_t * E_1(t) * D_n1(t) 
  #
  # As said in Madison et al 2007 "Estimating a Binary Character's Effect on Speciation 
  # and Extinction": "Dividing [D_N0 (t+small_t)-D_n0(t)] and [D_n1(t+small_t)-D_n1(t)] by 
  # the time interval, small_t, and taking the limit as small_t goes to zero, we can derive 
  # two coupled differential equations: "
  # dD_n0/dt = -(b0 + u0 + q[0,1]) * D_n0(t) + q[0,1] * D_n1(t) + 2 * b0 * E_0(t) * D_n0(t)
  # dD_n1/dt = -(b1 + u1 + q[1,0]) * D_n1(t) + q[1,0] * D_n0(t) + 2 * b0 * E_1(t) * D_n1(t)
  # 
  # With solutions for E_0(t) and E_1(t), those equations can be numerically intergrated 
  # along a branch. Therefore with the probabilites at the leaf end of a branch, the 
  # D_n0(t) and D_n1(t) at the root end of the branch can be derived.
  
  parms = c(b0, b1, u0, u1, q01, q10)
  names(parms) = c("b0", "b1", "u0", "u1", "q01", "q10")
  
  # Code from Nick Matze========================================
  
  # t = current time point of integration
  # y = state variable we are tracking (named)  MUST HAVE NAMES!!!
  # parms = model parameters (named)            MUST HAVE NAMES!!!
  
  define_BiSSE_eqns_in_R <- function(t, y, parms)
  {
    with(data=
           as.list(c(y, parms)), 
         { # expr
           # When the limit is taken as deltaT goes to 0, the
           # change in E0 in dt is:
           # probs of: 
           # - extinction
           # - no change but later extinction
           # - state change, then extinction
           # - no change, speciation, extinction of both
           dE0t <- u0 - (u0 + q01 + b0)*E0t + q01*E1t + b0*(E0t)^2
           dE1t <- u1 - (u1 + q10 + b1)*E1t + q10*E0t + b1*(E1t)^2
           
           # probs of:
           # - no change
           # - character change, no speciation
           # - speciation followed by extinction of 1
           # - extinction (prob of observed clade = 0, since the clade is extant)
           dD0t <- -1*(b0 + u0 + q01)*D0t + q01*D1t + 2*b0*E0t*D0t + 0
           dD1t <- -1*(b1 + u1 + q10)*D1t + q10*D0t + 2*b1*E1t*D1t + 0
           
           # Return the list of coupled differential equations
           res <- c(dE0t, dE1t, dD0t, dD1t)
           return(list(res))
           #return(list_of_diff_eqns)
         }
    )
  }
  
  # ============================================================================
  
  # Reorder phy's edge matrix so that when itterating through calculating likelyhoods 
  # for nodes we don't start from a node we haven't already calculated likelyhoods for 
  phy <- reorder(phy, order = "pruningwise")
  
  edges <- phy$edge
  node_times <- branching.times(phy)
  num_tips <- length(phy$tip.label)
  num_nodes <- phy$Nnode + num_tips
  
  # Table of state likelyhoods (E_0, E_1, D_n0, D_n1) for each node
  # Initialise everything to 0 except the D's for the tips, which we know from endStates
  E_D_likelyhood <- matrix(data = 0, nrow = num_nodes, ncol = 4)
  E_D_likelyhood[1:num_tips, 4] <- sapply(phy$tip.label, function(x) endstates[which(names(endstates) == x)])
  E_D_likelyhood[which(E_D_likelyhood[,4] == 0), 3] <- 1
  
  numSteps = 1000 #number of steps to calc E's and D's down each branch
  # because we have to calculate likelyhoods at nodes using both the node's descending
  # branches, our loop will traverse down two branches at a time
  loop_edge_indexes <- seq(1,length(phy$edge.length), 2)
  
  for(i in loop_edge_indexes)
  {
    j <- i + 1 # the right branch edge number
    
    # node numbers for left tip, right tip and ancestor
    left_tip_node <- phy$edge[i,2]
    right_tip_node <- phy$edge[j,2]
    ancestor <- phy$edge[i,1]
    
    # time values for lsoda 
    left_times <- seq(from=0, to=phy$edge.length[i], by= phy$edge.length[i]/numSteps)
    right_times <- seq(from=0, to=phy$edge.length[j], by= phy$edge.length[j]/numSteps)
    
    # initial likelyhoods for lsoda on left branch
    y = E_D_likelyhood[left_tip_node,]
    names(y) = c("E0t", "E1t", "D0t", "D1t")
    # likelyhoods down left branch 
    left_likely <- lsoda(y, left_times, define_BiSSE_eqns_in_R, parms)
    
    # initial likelyhoods for lsoda on right branch
    y = E_D_likelyhood[right_tip_node,]
    names(y) = c("E0t", "E1t", "D0t", "D1t")
    # likelyhoods down right branch
    right_likely <- lsoda(y, right_times, define_BiSSE_eqns_in_R, parms)
    
    # Likelyhood of the brach states at the root-facing end
    # (Col indexes for state likelyhoods are +1 in lsoda output because time column is 
    # at col index 1)
    D_Left_asc0 <- left_likely[nrow(left_likely), 4]
    D_Left_asc1 <- left_likely[nrow(left_likely), 5]
    D_Right_asc0 <- right_likely[nrow(right_likely), 4]
    D_Right_asc1 <- right_likely[nrow(right_likely), 5]
    # Probability ancestor is in state 0
    D_anc0 <- D_Left_asc0 * D_Right_asc0 * b0
    # Probability ancestor is in state 0
    D_anc1 <- D_Left_asc1 * D_Right_asc1 * b1
    E_D_likelyhood[ancestor,3] <- D_anc0
    E_D_likelyhood[ancestor,4] <- D_anc1
  }
  
  # likelyhood of the tree given the parameters and assuming that it's a 50%/50% prob
  # that the root was either state
  lazy_likely <- E_D_likelyhood[ancestor, 3] + E_D_likelyhood[ancestor, 4]
  print(E_D_likelyhood)
  print(lazy_likely)
  return(log(lazy_likely))
}

trstr = "(((chimp:1,human:1):1,gorilla:2):1,orang:3);"
tr = read.tree(file="", text=trstr)
plot(tr)
endStates <- c(1,1,1,1)
names(endStates) <- tr$tip.label

# Single birthrate, no death rate, no state transition 
biSSE.likelyhood(0.222222222,0.222222222,0,0,0,0,tr,endStates)
# Compare to diversitree's function
bisse_func <- make.bisse(tr, endStates, sampling.f=c(1,1), strict=FALSE)
parms <- c(0.222222222,0.222222222,0,0,0,0)
names(parms) = c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")
bisse_func(pars=parms, root=ROOT.GIVEN, root.p=c(0,1), intermediates=TRUE, condition.surv=FALSE)
bisse_func(pars=parms, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)

# Single birthrate, single death rate, no state transition 
birthrate <- 0.222222222
deathrate <- 0.1

b0 <- birthrate
b1 <- birthrate
u0 <- deathrate
u1 <- deathrate
q01 <- 0
q10 <- 0

biSSE.likelyhood(b0,b1,u0,u1,q01,q10,tr,endStates)
# Compare to diversitree's function
parms <- c(b0,b1,u0,u1,q01,q10)
names(parms) = c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")
bisse_func(pars=parms, root=ROOT.GIVEN, root.p=c(0,1), intermediates=TRUE, condition.surv=FALSE)
bisse_func(pars=parms, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)

# Single birthrate, single death rate, single state transition 
birthrate <- 0.222222222
deathrate <- 0.1

b0 <- birthrate
b1 <- birthrate
u0 <- deathrate
u1 <- deathrate
q01 <- 0.05
q10 <- 0.05

biSSE.likelyhood(b0,b1,u0,u1,q01,q10,tr,endStates)
# Compare to diversitree's function
parms <- c(b0,b1,u0,u1,q01,q10)
names(parms) = c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")
bisse_func(pars=parms, root=ROOT.GIVEN, root.p=c(0,1), intermediates=TRUE, condition.surv=FALSE)
bisse_func(pars=parms, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)

# Testing with a random tree
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(2)
phy <- tree.bisse(pars, max.t = 50, x0 = 0)

lik <- make.bisse(phy, phy$tip.state)
lik(pars)
biSSE.likelyhood(pars[1], pars[2], pars[3], pars[4], pars[5], pars[6], phy, phy$tip.state)



