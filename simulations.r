separable_beta = function(steps, phi, I0) {
  discrete_pchisq = diff(pchisq(0:steps, df=4))
  beta = function(t, tau) {
      return(phi * discrete_pchisq[tau])
  }
  # simulate a sequence of incident cases
  I = I0
  expected_I = I0
  for (t in 2:steps) {
    I = c(I, 0) # so I[t] = 0
    expected_I = c(expected_I, 0)
    for (tau in 1:(t - 1)) {
      I[t] = I[t] + (rpois(1, lambda=beta(t, tau)) * I[t - tau])
      expected_I[t] = expected_I[t] + beta(t, tau) * expected_I[t - tau]
    }
  }
  return(list("I"=I, "expected_I"=expected_I, "residual"=I - expected_I))
}

isolating_measures = function(steps, phi, tau0, I0) {
  t0 = floor(steps / 2)
  discrete_pchisq = diff(pchisq(0:steps, df=4))
  beta = function(t, tau) {
    if (t <= t0) {
      return(phi * discrete_pchisq[tau])
    } else {
      if (tau > tau0) {
        return(0)
      } else {
        return ((phi / 2) * (discrete_pchisq[tau] / pchisq(tau0, df=4)))
      }
    }
  }
  # simulate a sequence of incident cases
  I = I0
  expected_I = I0
  for (t in 2:steps) {
    I = c(I, 0) # so I[t] = 0
    expected_I = c(expected_I, 0)
    for (tau in 1:(t - 1)) {
      I[t] = I[t] + (rpois(1, lambda=beta(t, tau)) * I[t - tau])
      expected_I[t] = expected_I[t] + beta(t, tau) * expected_I[t - tau]
    }
  }
  return(list("I"=I, "expected_I"=expected_I))
}