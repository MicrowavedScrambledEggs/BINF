# BiSSE likelyhoodfunction
# Takes birth rate and death rate for each character state (b0, b1, d0, d1), a q matrix
# for the characer transition rates, and the tree

biSSE.likelyhood <- function(b0, b1, d0, d1, q)
{
  # Uses the differentials of the probabilities each event could happen between a small 
  # time period
  
  # D_n1(t) : The likelyhood that an extant species has the character state 1 at time t
  # Starting state: if tip has state 1 D_n1(0) = 1 D_n0(0) = 0
  #                 if tip has state 0 D_n1(0) = 0 D_n0(0) = 1
  # The likelyhood that something stays at state 0 on a branch between time t and
  # time t + small_t and does not go extinct
  # D_n0(t+small_t) = p(did not go extinct)
  #                   * ( p(No state change and no speciation)
  #                   + p(State change and no speciation)
  #                   + p(No state change and speication then extinction)
  #                   + p(No state change and speication then extinction))
  # We don't include p(state change and speciation then extinction) because small_t is too
  # small for the probability of both occuring to be significant.
  # Any speciation events have to be followed by extinction because it would mean a branch
  # that is not observed in the given tree. We include speciation then extinction twice as
  # extinction could have ouccured in either of the two decendants from a speciation.
  # D_n0(t+small_t) = (1 - d0*small_t)
}