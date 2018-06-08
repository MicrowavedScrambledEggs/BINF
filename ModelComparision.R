# Needs the bd.log.likelyhood.function() from BirthDeathLikelyhood.R
# and biSSE.likelyhood() from BiSSELikelyhood.R
library(ape)

# Takes a phylogeny and a birth rate and outputs the log likelihood of
# the birth rate given the tree under the yule model
yule.log.likely <- function(brate, phy)
{
  # The probability that a the ith birth event happened is proportional
  # to i*brate (except for the 1st birth event which is 1 as otherwise
  # would not have a tree) so for the n birth events in the tree it is:
  # n!*brate^(n-1). The probability no extra birth events happened 
  # along branch j is proportional to exp(-brate*branch.length[j]) so
  # for all branches it is exp(-brate*sum(branch.length))
  # Combining then converting to log likelihood we get:
  return(-brate * sum(phy$edge.length) + lfactorial(phy$Nnode) 
         + (phy$Nnode - 1) * log(brate))
}

finches <- read.nexus("http://www.r-phylo.org/w/images/0/02/Geospiza.nex")

optim_yule <- optimise(yule.log.likely, phy = finches, interval = c(0,100), 
                       maximum = TRUE)
optim_bd <- optim(c(2.5,2.25), control = list(fnscale=-1),
                  function(p) bd.log.likelyhood.function(p[1],p[2],finches))

# BiSSE ==================================================================

# load tree
shore_birds <- read.nexus("http://www.zoology.ubc.ca/prog/diversitree/doc/files/Thomas-tree.nex")
# shore_birds <- read.nexus("Thomas-tree.nex")
# It's a bit big for how slow my BiSSE runs so will take a subtree
sub.shore.bird <- subtrees(shore_birds)
# Find the biggest tree where all nodes have 2 daughters as that's the
# only types of trees that work with my function
bin_trees <- c()
n_tips <- c()
for(i in 1:length(sub.shore.bird)){
  if(is.binary.tree(sub.shore.bird[[i]])) {
    tip_n <- sub.shore.bird[[i]]$Ntip
    cat("The", i, "th tree is binary. No of tips:", tip_n, "\n")
    bin_trees <- c(bin_trees, i)
    n_tips <- c(n_tips, tip_n)
  }
}
shore.bird <- sub.shore.bird[[bin_trees[which.max(n_tips)]]]

# load state data
sex_dimorph <- read.csv("http://www.zoology.ubc.ca/prog/diversitree/doc/files/Lislevand-states.csv",
                        as.is = TRUE)
states <- sex_dimorph$dimorph
names(states) <- sex_dimorph$species
states <- states[shore.bird$tip.label]
names(states) <- shore.bird$tip.label

# Convert to binary states. Count absolute difference in body length 
# > 7.5% as sexually dimorphic. Count NA's as state 0 for now
states.075 <- rep(0, length(states))
states.075[which(abs(states) > 0.075)] <- 1
names(states.075) <- names(states) 

# Test
testBisse <- biSSE.likelyhood(0.4,0.3,0.2,0.1,0.07,0.04, shore.bird, states.075, 
                 output.branches = TRUE)
prefunc <- make.bisse(shore.bird, states.075)

# every parameter different
optim_ve <- optim(c(0.04,0.03,0.02,0.01,0.07,0.04), control = list(fnscale=-1),
                  method=c("L-BFGS-B"), lower = c(0.000001,0.000001,0,0,0,0),
                  upper = c(10,10,10,10,10,10),
                  fn = function(x) biSSE.likelyhood(x[1],x[2],x[3],x[4],x[5],x[6],
                                                    shore.bird, states.075))
# Single death and transition rates
optim_dt <- optim(c(0.04,0.03,0.02,0.07), control = list(fnscale=-1),
                  method=c("L-BFGS-B"), lower = c(0.000001,0.000001,0,0),
                  upper = c(1,1,1,1),
                  fn = function(x) biSSE.likelyhood(x[1],x[2],x[3],x[3],x[4],x[4],
                                                    shore.bird, states.075))
# No death, single birth rate
optim_y <- optim(c(0.04,0.07,0.04), control = list(fnscale=-1), method=c("L-BFGS-B"), 
                 lower = c(0.000001,0,0), upper = c(1,1,1),
                  fn = function(x) biSSE.likelyhood(x[1],x[1],0,0,x[2],x[3],
                                                    shore.bird, states.075))
# Single transition rate
optim_t <- optim(c(0.04,0.03,0.02,0.01,0.07), control = list(fnscale=-1),
                 method=c("L-BFGS-B"), lower = c(0.000001,0.000001,0,0,0),
                 upper = c(1,1,1,1,1),
                  fn = function(x) biSSE.likelyhood(x[1],x[2],x[3],x[4],x[5],x[5],
                                                    shore.bird, states.075))
