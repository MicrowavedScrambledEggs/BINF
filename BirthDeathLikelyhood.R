library(ape)
finches <- read.nexus("http://www.r-phylo.org/w/images/0/02/Geospiza.nex")

## Birth death model likelyhood function
# Takes a birth rate, a death rate, and an object of class phylo 
bd.likelyhood.function = function(brate, drate, phy)
{
  N = length(phy$tip.label) # Number of linearges
  t = branching.times(phy) # time since linearge birth
  
  a = drate / brate # extinction rate
  r = brate - drate # speciation rate
  
  # P(t, T) Likelyhood of a lineage alive at time t is still alive at later time T
  # r /(r*exp(-r*(T-t)))
  
  # Probability that the ith extant lineage to be born was born at time t_i
  # (i - 1)*brate*P(t_i, T)
  # The prob that the 1st and 2nd lineages were born is 1 as otherwise we wouldn't have a tree
  
  # Probability that a lineage does not give birth to another lineage in time period t
  # (1 - u_t) where u_t is the probability that a lineage DOES give birth to another lineage
  # in time t
  # u_t = (brate*(1-exp(-r*t)))/(brate - drate*exp(-r*t))
  
  cat("b:", brate, "d:", drate, "N:", N, "a", a, "r", r, "\n")
  if(r < 0 || a > 1)
    {likelyhood = 1e-100}
  else
    {likelyhood = gamma(N)*r**(N-2)*exp(r*sum(t[2:length(t)]))*(1-a)**N*prod(1/(exp(r*t)-a)**2)}
  if (likelyhood == 0) cat("Shock Horror!", brate, drate, "\n")
  print(log(likelyhood))
  return(log(likelyhood))
}

# optim(c(2.5,2.25), function(p) bd.likelyhood.function(p[1],p[2],finches), 
#       lower=c(0,0), control = list(fnscale=-1))