## Birth death model likelyhood function
# Takes a birth rate, a death rate, and an object of class phylo 
bd.likelyhood.function = function(brate, drate, phy)
{
  N = length(phy$tip.label) # Number of linearges
  t = branching.times(phy) # time since linearge birth
  # will have to edit if nodes have > 2 daughers
  
  a = drate / brate
  r = brate - drate
  
  likelyhood = gamma(N)*r**(N-2)*exp(r*sum(t[2:length(t)]))*(1-a)*prod(1/(exp(r*t)-a)**2)
  return(likelyhood)
}