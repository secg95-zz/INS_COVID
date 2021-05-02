get_expected = function(beta, N0, tau1, tau2) {
  "
  Returns
  -------
  expected_NN : numeric vector
    Daily expected number of new infectious cases
  expected_N : numeric vector
    Daily expected number of infectious cases.
  "
  daysUntilFirstSintoms=0
  daysUntilFirstSintoms=ceiling(1/tau1)-1
  beta=c(rep(beta[1],daysUntilFirstSintoms),beta)
  moments = length(beta)
  expected_NN = rep(0, moments)
  expected_N = rep(0, moments)
  expected_XN = rep(0, moments)
  expected_X = rep(0, moments)
  X_previous = N0
  N_previous = 0
  
  for(t in 1:moments){
    expected_N[t] = N_previous * pexp(1, rate=tau2, lower.tail=FALSE) + X_previous * pexp(1, rate=tau1, lower.tail=TRUE)
    expected_X[t] = N_previous * beta[t] + X_previous * pexp(1, rate=tau1, lower.tail=FALSE)
    expected_NN[t] = + X_previous * pexp(1, rate=tau1, lower.tail=TRUE)
    expected_XN[t] = N_previous * beta[t]
    N_previous = expected_N[t]
    X_previous = expected_X[t]
  }
  expected = list("NN"=expected_NN[(daysUntilFirstSintoms+1):moments], "N"=expected_N[(daysUntilFirstSintoms+1):moments],"XN"=expected_XN[(daysUntilFirstSintoms+1):moments], "X"=expected_X[(daysUntilFirstSintoms+1):moments])
  return(expected)
}
