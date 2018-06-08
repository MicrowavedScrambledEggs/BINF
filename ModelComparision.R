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
shore.bird <- sub.shore.bird[[5]]

# load state data
sex_dimorph <- read.csv("http://www.zoology.ubc.ca/prog/diversitree/doc/files/Lislevand-states.csv",
                        as.is = TRUE)
states <- sex_dimorph$dimorph
names(states) <- sex_dimorph$species
states <- states[shore.bird$tip.label]
names(states) <- shore.bird$tip.label

# Convert to binary states. Count absolute difference in body length 
# > 15% as sexually dimorphic
states.15 <- (abs(states) > 0.15) + 0

# Test
testBisse <- biSSE.likelyhood(0.4,0.3,0.2,0.1,0.07,0.04, shore.bird, states.15, 
                 output.branches = TRUE)

# every parameter different
optim_ve <- optim()
