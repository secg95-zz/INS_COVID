simulate_bayesian = function(steps, phi, I0) {
  discrete_pchisq = diff(pchisq(0:steps, df=6))
  beta = function(t, tau) {
      return(phi * discrete_pchisq[tau])
  }
  # simulate a sequence of incident cases
  I = I0
  expected = I0
  expected_conditional = I0
  for (t in 2:steps) {
    I = c(I, 0) # so I[t] = 0
    expected = c(expected, 0)
    expected_conditional = c(expected_conditional, 0)
    for (tau in 1:(t - 1)) {
      I[t] = I[t] + (rpois(1, lambda=beta(t, tau)) * I[t - tau])
      expected[t] = expected[t] + beta(t, tau) * expected[t - tau]
      expected_conditional[t] = expected_conditional[t] + beta(t, tau) * I[t - tau]
    }
  }
  return(list("I"=I, "expected"=expected, "expected_conditional"=expected_conditional,
              "omega"=discrete_pchisq))
}

source("betaStateSpace/serial_interval_construction.R")
simulate_bayesian2 = function(steps, phi, I0) {
  discrete_pchisq = weightConstructionSerialInterval(lags=steps)
  beta = function(t, tau) {
    return(phi * discrete_pchisq[tau])
  }
  # simulate a sequence of incident cases
  I = I0
  expected = I0
  expected_conditional = I0
  for (t in 2:steps) {
    I = c(I, 0) # so I[t] = 0
    expected = c(expected, 0)
    expected_conditional = c(expected_conditional, 0)
    for (tau in 1:(t - 1)) {
      I[t] = I[t] + (rpois(1, lambda=beta(t, tau)) * I[t - tau])
      expected[t] = expected[t] + beta(t, tau) * expected[t - tau]
      expected_conditional[t] = expected_conditional[t] + beta(t, tau) * I[t - tau]
    }
  }
  return(list("I"=I, "expected"=expected, "expected_conditional"=expected_conditional,
              "omega"=discrete_pchisq))
}

simulate_inseparable_beta = function(steps, phi, tau0, I0) {
  t0 = floor(steps / 2)
  discrete_pchisq = diff(pchisq(0:steps, df=6))
  beta = function(t, tau) {
    if (t <= t0) {
      return(phi * discrete_pchisq[tau])
    } else {
      if (tau > tau0) {
        return(0)
      } else {
        return (phi * discrete_pchisq[tau] / pchisq(tau0, df=6))
      }
    }
  }
  omega = function(t, tau) {
    if (t <= t0) {
      return(discrete_pchisq[tau])
    } else {
      if (tau > tau0) {
        return(0)
      } else {
        return (discrete_pchisq[tau] / pchisq(tau0, df=4))
      }
    }
  }
  # simulate a sequence of incident cases
  I = I0
  expected = I0
  expected_conditional = I0
  for (t in 2:steps) {
    I = c(I, 0) # so I[t] starts at 0
    expected = c(expected, 0)
    expected_conditional = c(expected_conditional, 0)
    for (tau in 1:(t - 1)) {
      I[t] = I[t] + (rpois(1, lambda=beta(t, tau)) * I[t - tau])
      expected[t] = expected[t] + beta(t, tau) * expected[t - tau]
      expected_conditional[t] = expected_conditional[t] + beta(t, tau) * I[t - tau]
    }
  }
  return(list("I"=I, "expected"=expected, "expected_conditional"=expected_conditional,
              "omega"=omega))
}