library(ape)
finches <- read.nexus("http://www.r-phylo.org/w/images/0/02/Geospiza.nex")

## Birth death model likelyhood function
# Takes a birth rate, a death rate, and an object of class phylo 
bd.likelyhood.function = function(brate, drate, phy)
{
  N = length(phy$tip.label) # Number of linearges
  t = branching.times(phy) # time since linearge birth
  # will have to edit if nodes have > 2 daughers
  
  a = drate / brate # extinction rate
  r = brate - drate # speciation rate
  cat("b:", brate, "d:", drate, "N:", N, "a", a, "r", r, "\n")
  if(r < 0 || a > 1)
    {likelyhood = 1e-100}
  else
    {likelyhood = gamma(N)*r**(N-2)*exp(r*sum(t[2:length(t)]))*(1-a)**N*prod(1/(exp(r*t)-a)**2)}
  if (likelyhood == 0) cat("Shock Horror!", brate, drate, "\n")
  print(log(likelyhood))
  return(log(likelyhood))
}

optim(c(2.5,2.25), function(p) bd.likelyhood.function(p[1],p[2],finches), 
      lower=c(0,0), control = list(fnscale=-1))