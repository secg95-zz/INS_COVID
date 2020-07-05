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

simulate_impartial = function(steps, beta, tau1, tau2, I0) {
  "
  Simulate Cerón's impartial model.
  
  Parameters
  ----------
  steps : numeric
    Number of steps in the simulation.
  beta : numeric vector
    Series where beta[t] is the expected number of people an infectious patient
    will infect in step t.
  tau1 : numeric
    Expected number of days since exposure to infection onset.
  tau2 : numeric
    Expected number of days a patient remains infectious.
  I0 : numeric
    Number of intial just-infectious cases.
  "
  # initialize simulated vectors to be returned
  E = 0
  N = I0
  I = I0
  expected_I = I0
  expected_E = 0
  expected_N = I0
  # carry out simulation and append to those vectors
  for (t in 1:steps) {
    I = c(I, rbinom(1, size=E[t], prob=1 / tau1))
    E = c(E, E[t] - I[t + 1] + rpois(1, lambda=beta[t] * N[t]))
    N = c(N, N[t] + I[t + 1] - rbinom(1, size=N[t], prob=1 / tau2))
    # track expected values
    expected_I = c(expected_I, expected_E[t] / tau1)
    expected_E = c(expected_E, expected_E[t] - expected_I[t + 1] + expected_N[t] * beta[t])
    expected_N = c(expected_N, expected_I[t + 1] + expected_N[t] * (1 - (1 / tau2)))
  }
  # calculate theoretical values
  R = beta * tau2
  f_inc = dgeom(1:steps, p=1 / tau1)
  f_inf = dgeom(1:steps, p=1 / tau2)
  omega = NULL
  for (tau in 1:steps) {
    omega = c(omega, 0) # so omega[tau] = 0
    for (tau_prime in 0:(tau - 1)) {
      omega[tau] = omega[tau] +
                   pgeom(tau_prime - 1, p=1 / tau2, lower.tail=FALSE) *
                   dgeom(tau - tau_prime - 1, p=1 / tau1)
    }
  }
  omega = omega / tau2
  return(list(
    "I"=I, "E"=E, "N"=N, "expected_I"=expected_I, "R"=R, "f_inc"=f_inc,
    "f_inf"=f_inf, "omega"=omega, "steps"=steps
  ))
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