library(deSolve)
library(diversitree)

# BiSSE likelihood function
# Takes birth rate and death rate for each character state (b0, b1, u0, u1), a q values 
# for the characer transition rates, the tree, and a vector of states for the tree tips
# The vector of states has to have a names attribue matching the names in phy$tip.label to
# match the tips to their states.
# Equations taken from Madison et al 2007 "Estimating a Binary Character's Effect on 
# Speciation and Extinction"
# Outputs the log likelihood of the parameters given the tree and states. This is 
# calculated from the mean of the likelihoods for the root being in each state, as we
# assume a 0.5 probability of the root being in either state (which may not be the case
# for all phylogenies)
# If output.branches = TRUE then matricies storing the relative probabilites at either
# side of each branch will be included in the output as attributes
biSSE.likelyhood <- function(b0, b1, u0, u1, q01, q10, phy, tipstates, 
                             output.branches = FALSE)
{
  # Uses the differentials of the probabilities each event could happen between a small 
  # time period
  
  # D_n1(t) : The likelyhood that an extant species has the character state 1 at time t
  # Starting state: if tip has state 1 D_n1(0) = 1; D_n0(0) = 0
  #                 if tip has state 0 D_n1(0) = 0; D_n0(0) = 1
  # E_0(t) : The likelyhood that a species extinct by time 0 is in state 0 at time t 
  # E_0 and E_1 at time 0 for all tips is 0 as they are extant at time 0
  
  # The likelyhood that a species stays in state 0 on a branch between time t and
  # time t + small_t and is extant at time 0 is
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
  #
  # The likelyhood that a species extinct by time 0 is in state 0 on a branch between time t 
  # and time t + small_t equals:
  # E_0(t+small_t) = u0 * small_t   # Goes extinct in small_t
  #                  # + No state change and no speciation
  #                  + (1 - u0 * small_t) * (1 - q01 * small_t) * (1 - b0 * small_t) * E_o(t)
  #                  # + State change and no speciation
  #                  + (1 - u0 * small_t) * q01 * small_t * (1 - b0 * small_t) * E_1(t)
  #                  # + No state change and speciation
  #                  + (1 - u0 * small_t) * (1 - q01 * small_t) * (b0 * small_t) * E_o(t)^2
  #
  # Because terms of order small_t^2 are negligably small we can drop them, giving:
  # E_0(t+small_t) = u0 * small_t + (1 - (u0 + q01 + b0) * small_t) * E_0(t)
  #                  + q01 * small_t * E_1(t) + b0 * small_t * E_0(t)^2
  #
  # Once again, divding the change in E0  by the time interval and taking the limit as
  # small_t goes to zero gives us the differential equations:
  # dE_0/dt = u0 - (u0 + q01 + b0) * E_0(t) _ q01 * E_1(t) + b0 * E_0(t)^2
  # 
  # The same can be done for E_1. With the four coupled differential equations solutions can be
  # found from numberical integration
  
  parms = c(b0, b1, u0, u1, q01, q10)
  names(parms) = c("b0", "b1", "u0", "u1", "q01", "q10")
  
  # Create a function for use by the lsoda() function from the deSolve package for numerical 
  # intergration. Returns the values got from applying the differntial equations for E0, E1, D0 
  # and D1 described above
  # Code from Nick Matze========================================
  # t = current time point of integration
  # y = state variables we are tracking (named "E0t", "E1t", "D0t", "D1t" respectively)
  # parms = model parameters (named "b0", "b1", "u0", "u1", "q01", "q10" respectively)
  
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
  # Initialise everything to 0 except the D's for the tips, which we know from tipstates
  E_D_likelyhood <- matrix(data = 0, nrow = num_nodes, ncol = 4)
  E_D_likelyhood[1:num_tips, 4] <- sapply(phy$tip.label, function(x) tipstates[which(names(tipstates) == x)])
  E_D_likelyhood[which(E_D_likelyhood[,4] == 0), 3] <- 1
  
  if(output.branches){
    E_D_init <- E_D_likelyhood 
    E_D_base <- matrix(data = NA, nrow = num_nodes, ncol = 4)
  }
  
  numSteps = 1000 #number of steps to calc E's and D's down each branch
  # because we have to calculate likelyhoods at nodes using both the node's daughter
  # branches, our loop will traverse down two branches at a time
  loop_edge_indexes <- seq(1,length(phy$edge.length), 2)
  
  for(i in loop_edge_indexes)
  {
    j <- i + 1 # the right branch edge number
    
    # node numbers for node at tip end of left branch, right branch and ancestor
    left_tip_node <- phy$edge[i,2]
    right_tip_node <- phy$edge[j,2]
    ancestor <- phy$edge[i,1]
    
    # time values for lsoda 
    left_times <- seq(from=0, to=phy$edge.length[i], by= phy$edge.length[i]/numSteps)
    right_times <- seq(from=0, to=phy$edge.length[j], by= phy$edge.length[j]/numSteps)
    
    # initial likelyhoods for lsoda on left branch
    y <- E_D_likelyhood[left_tip_node,]
    names(y) = c("E0t", "E1t", "D0t", "D1t")
    # likelyhoods down left branch 
    left_likely <- lsoda(y, left_times, define_BiSSE_eqns_in_R, parms)
    
    # initial likelyhoods for lsoda on right branch
    y <- E_D_likelyhood[right_tip_node,]
    names(y) = c("E0t", "E1t", "D0t", "D1t")
    # likelyhoods down right branch
    right_likely <- lsoda(y, right_times, define_BiSSE_eqns_in_R, parms)
    
    # Likelyhood of the brach states at the root-facing end
    # (Col indexes for state likelyhoods are +1 in lsoda output because time column is 
    # at col index 1)
    if(output.branches){
      
    }
    D_Left_root <- left_likely[nrow(left_likely), 4:5]
    D_Right_root <- right_likely[nrow(right_likely), 4:5]
    
    # Likelihood ancestor is in state 0
    D_anc0 <- D_Left_root[1] * D_Right_root[1] * b0
    # Likelyhood ancestor is in state 1
    D_anc1 <- D_Left_root[2] * D_Right_root[2] * b1
    E_D_likelyhood[ancestor,3] <- D_anc0
    E_D_likelyhood[ancestor,4] <- D_anc1
    
    # Extinction likelyhoods at node is just the mean between the two branches
    E_anc0 <- (left_likely[nrow(left_likely), 2] + right_likely[nrow(right_likely), 2]) / 2
    E_anc1 <- (left_likely[nrow(left_likely), 3] + right_likely[nrow(right_likely), 3]) / 2
    E_D_likelyhood[ancestor,1] <- E_anc0
    E_D_likelyhood[ancestor,2] <- E_anc1

    if(output.branches)
    {
      E_D_base[left_tip_node,] <- left_likely[nrow(left_likely), 2:5]
      E_D_base[right_tip_node,] <- left_likely[nrow(left_likely), 2:5]
      
      # relative probabilities of each state on each branch at the root facing end
      rel_Left_root <- D_Left_root / sum(D_Left_root)
      rel_Right_root <- D_Right_root / sum(D_Right_root)
      E_D_base[left_tip_node,3:4] <- rel_Left_root
      E_D_base[right_tip_node,3:4] <- rel_Right_root
      
      # Relative Probability ancestor is in state 0
      relD_anc0 <- rel_Left_root[1] * rel_Right_root[1] * b0
      # Relative Probability ancestor is in state 1
      relD_anc1 <- rel_Left_root[2] * rel_Right_root[2] * b1
      E_D_init[ancestor,] <- c(E_anc0, E_anc1, relD_anc0, relD_anc1)
    }
  }
  
  # likelyhood of the parameters given the tree and states assuming that it's a 50%/50% prob
  # that the root was either state
  lazy_likely <- mean(E_D_likelyhood[ancestor, 3:4])
  log_lazy_likely <- log(lazy_likely)
  if(output.branches){
    attributes(log_lazy_likely) <- list(init = t(E_D_init), base = t(E_D_base))
  }
  return(log_lazy_likely)
}




