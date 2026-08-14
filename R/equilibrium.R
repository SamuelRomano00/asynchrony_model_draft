#' Endemic equilibria and calibration
#'
#' Simulations are started at the endemic equilibrium of the model without
#' intervention, so that any change in incidence during the run is attributable
#' to the intervention schedule rather than to a transient from an arbitrary
#' initial condition.

# ---------------------------------------------------------------------------
# SIS model
# ---------------------------------------------------------------------------

#' Endemic equilibrium prevalence of the two-patch SIS model
#'
#' Solves the steady state of the SIS model in the absence of intervention
#' (`omega = 0`).
#'
#' @param R0 Numeric vector `c(R0_1, R0_2)`.
#' @param N Numeric vector `c(N_1, N_2)` of population sizes.
#' @param r Recovery rate, shared across patches.
#' @param p_12,p_21 Mobility fractions in \[0, 0.5\], see [mobility_matrix()].
#' @return Numeric vector `c(x_1, x_2)` of equilibrium prevalences, or
#'   `c(0, 0)` when the system-level reproduction number is at most 1, in which
#'   case the disease-free state is the only equilibrium.
compute_equilibrium_prevalence <- function(R0, N, r, p_12, p_21) {

  N_1 <- N[1]
  N_2 <- N[2]
  # Mobility matrix entries, named m11..m22 rather than a..d: `c` and `F` would
  # shadow the concatenation function and FALSE inside this function.
  m11 <- 1 - p_12   # patch-1 residents staying in patch 1
  m12 <- p_12       # patch-1 residents visiting patch 2
  m21 <- p_21       # patch-2 residents visiting patch 1
  m22 <- 1 - p_21   # patch-2 residents staying in patch 2

  R0_1 <- R0[1]
  R0_2 <- R0[2]

  # Effective transmission coefficients of the reduced two-variable system,
  # i.e. the algebraic form of the model's K[]/pwlk[] force of infection for
  # exactly two patches at steady state. Each denominator is the total
  # population PRESENT where transmission happens, mixing residents and
  # visitors: patch 1 sees a*N_1 + c*N_2, patch 2 sees b*N_1 + d*N_2.
  # k1_x and k2_x are built from patch-1 residents and so carry the x terms;
  # k1_y and k2_y are built from patch-2 residents and carry the y terms.
  k1_x <- r * R0_1 * m11 * N_1 / (m11 * N_1 + m21 * N_2)
  k1_y <- r * R0_1 * m21 * N_2 / (m11 * N_1 + m21 * N_2)
  k2_x <- r * R0_2 * m12 * N_1 / (m12 * N_1 + m22 * N_2)
  k2_y <- r * R0_2 * m22 * N_2 / (m12 * N_1 + m22 * N_2)
  
  # Quadratic coefficients of the two equilibrium equations.
  A <- (m11 * k1_x + m12 * k2_x - r)
  B <- (m11 * k1_x + m12 * k2_x)   # eq1: coefficient of x
  C <- (m11 * k1_y + m12 * k2_y)   # eq1: coefficient of y
  Fyy <- (m21 * k1_y + m22 * k2_y) # eq2: coefficient of y
  G <- (m21 * k1_x + m22 * k2_x)   # eq2: coefficient of x
  E <- Fyy - r

  eqs <- function(vars) {
    x <- vars[1]
    y <- vars[2]
    c(A * x - B * x^2 - C * x * y + C * y,
      E * y - Fyy * y^2 - G * x * y + G * x)
  }

  jacobian <- function(vars) {
    x <- vars[1]
    y <- vars[2]
    matrix(c(A - 2 * B * x - C * y, -C * x + C,
             -G * y + G,            E - 2 * Fyy * y - G * x),
           nrow = 2, byrow = TRUE)
  }

  # x = y = 0 is always a root, since every term of both equations carries a
  # factor x or y. Three safeguards keep the solver away from it when an
  # endemic equilibrium exists.

  # 1. Threshold test. J(0) reads off the Jacobian above at x = y = 0, so the
  #    system-level R0 is the spectral radius of [[B, C], [G, F]] / r. When it
  #    is at most 1 the disease-free state is the only equilibrium and zeros
  #    are the correct answer, not a solver failure.
  M <- matrix(c(B, C, G, Fyy), nrow = 2, byrow = TRUE)
  R0_system <- max(Mod(eigen(M, only.values = TRUE)$values)) / r
  if (R0_system <= 1) {
    return(c(0, 0))
  }

  # 2. Solve in logit space. Substituting x = plogis(u) puts the trivial root
  #    at u = -Inf, out of reach of any finite Newton step, and enforces the
  #    (0, 1) box constraint automatically.
  eqs_logit <- function(u) eqs(plogis(u))
  jacobian_logit <- function(u) {
    x <- plogis(u)
    jacobian(x) %*% diag(x * (1 - x))   # chain rule: d(plogis)/du = x*(1-x)
  }

  # 3. Multi-start with validation. The first guess is the isolated-patch SIS
  #    equilibrium; a candidate is accepted only if its residual is small and
  #    it is interior.
  start_isolated <- pmin(pmax(1 - 1 / R0, 0.01), 0.99)
  starts <- list(start_isolated, c(0.5, 0.5), c(0.1, 0.1), c(0.9, 0.9), c(0.01, 0.5))

  for (s in starts) {
    sol <- try(nleqslv(qlogis(s), eqs_logit, jac = jacobian_logit,
                       method = "Newton",
                       control = list(ftol = 1e-20, xtol = 1e-20, maxit = 200)),
               silent = TRUE)
    if (inherits(sol, "try-error")) next

    x_sol <- plogis(sol$x)
    if (max(abs(eqs(x_sol))) < 1e-10 && max(x_sol) > 1e-8 && max(x_sol) < 1 - 1e-9) {
      return(x_sol)
    }
  }

  # Fallback: forward integration. The endemic equilibrium is globally
  # asymptotically stable on [0,1]^n minus the origin, so integrating from any
  # positive point cannot converge anywhere else. Slow, hence the last resort.
  warning("compute_equilibrium_prevalence(): Newton failed from every start ",
          "despite R0_system = ", signif(R0_system, 4), " > 1; ",
          "falling back to forward integration.")
  x_int <- c(0.5, 0.5)
  dt <- 0.1 / r
  for (k in 1:200000) {
    k1 <- eqs(x_int)
    if (max(abs(k1)) < 1e-14) break
    k2 <- eqs(x_int + dt / 2 * k1)
    k3 <- eqs(x_int + dt / 2 * k2)
    k4 <- eqs(x_int + dt * k3)
    x_int <- pmin(pmax(x_int + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4), 0), 1)
  }
  x_int
}

#' Reproduction numbers implied by a target equilibrium prevalence
#'
#' The inverse of [compute_equilibrium_prevalence()]: given an observed or
#' targeted equilibrium prevalence, recover the `R0` values that produce it.
#'
#' @param x0 Numeric vector `c(x0_1, x0_2)` of target equilibrium prevalences.
#' @param N Numeric vector `c(N_1, N_2)` of population sizes.
#' @param r Recovery rate.
#' @param p_12,p_21 Mobility fractions in \[0, 0.5\].
#' @return Numeric vector `c(R0_1, R0_2)`.
compute_R0 <- function(x0, N, r, p_12, p_21) {

  N_1 <- N[1]
  N_2 <- N[2]
  m11 <- 1 - p_12   # patch-1 residents staying in patch 1
  m12 <- p_12       # patch-1 residents visiting patch 2
  m21 <- p_21       # patch-2 residents visiting patch 1
  m22 <- 1 - p_21   # patch-2 residents staying in patch 2
  X_1 <- x0[1]
  X_2 <- x0[2]

  # Prevalence met in each patch, same mixing of residents and visitors as in
  # the model, treated as known here since x0 is the target.
  K_1 <- (m11 * N_1 * X_1 + m21 * N_2 * X_2) / (m11 * N_1 + m21 * N_2)
  K_2 <- (m12 * N_1 * X_1 + m22 * N_2 * X_2) / (m12 * N_1 + m22 * N_2)

  # Linear in the unknowns (R0_1, R0_2), obtained by rearranging the SIS
  # equilibrium condition into odds form, X / (1 - X), for each patch. No
  # analytic Jacobian is needed for a linear system.
  eqs <- function(vars) {
    c(m11 * K_1 * vars[1] + m12 * K_2 * vars[2] - X_1 / (1 - X_1),
      m21 * K_1 * vars[1] + m22 * K_2 * vars[2] - X_2 / (1 - X_2))
  }

  sol <- nleqslv(c(1, 1), eqs)
  c(sol$x[1], sol$x[2])
}

# ---------------------------------------------------------------------------
# Ross-Macdonald model
# ---------------------------------------------------------------------------

#' Mosquito-to-human ratio implied by a target R0
#'
#' Inverts `R0 = a^2 * b * c * m * exp(-mu * tau) / (r * mu)`, so that the
#' Ross-Macdonald model can be driven by the same `R0` lever as the SIS model.
#'
#' @param R0 Target reproduction number, scalar or per patch.
#' @param a Biting rate.
#' @param b Transmission probability, mosquito to human.
#' @param c_rate Transmission probability, human to mosquito.
#' @param r Human recovery rate.
#' @param mu Mosquito mortality.
#' @param tau Extrinsic incubation period, in days.
#' @return The mosquito-to-human ratio `m`, matching the shape of `R0`.
R0_to_m <- function(R0, a, b, c_rate, r, mu, tau) {
  R0 * r * mu / (a^2 * b * c_rate * exp(-mu * tau))
}

#' Endemic equilibrium of the two-patch Ross-Macdonald model
#'
#' Solves the human and mosquito steady states jointly, without intervention.
#'
#' @inheritParams R0_to_m
#' @param N Numeric vector `c(N_1, N_2)` of population sizes.
#' @param p_12,p_21 Mobility fractions in \[0, 0.5\].
#' @return A list with `X`, the human equilibrium prevalence per patch, and
#'   `Y`, the mosquito sporozoite prevalence per patch.
compute_equilibrium_prevalence_RM <- function(N, a, b, c_rate, r, mu, tau, R0, p_12, p_21) {

  m <- R0_to_m(R0, a, b, c_rate, r, mu, tau)

  N_1 <- N[1]; N_2 <- N[2]
  p <- mobility_matrix(p_12, p_21)
  eip_survival <- exp(-mu * tau)

  D1 <- p[1, 1] * N_1 + p[2, 1] * N_2
  D2 <- p[1, 2] * N_1 + p[2, 2] * N_2

  system <- function(vars) {
    X1 <- vars[1]; X2 <- vars[2]
    Y1 <- vars[3]; Y2 <- vars[4]

    foi1 <- a * b * (p[1, 1] * m[1] * Y1 + p[1, 2] * m[2] * Y2)
    foi2 <- a * b * (p[2, 1] * m[1] * Y1 + p[2, 2] * m[2] * Y2)

    K1 <- (p[1, 1] * N_1 * X1 + p[2, 1] * N_2 * X2) / D1
    K2 <- (p[1, 2] * N_1 * X1 + p[2, 2] * N_2 * X2) / D2

    c(foi1 * (1 - X1) - r * X1,
      foi2 * (1 - X2) - r * X2,
      a * c_rate * K1 * (eip_survival - Y1) - mu * Y1,
      a * c_rate * K2 * (eip_survival - Y2) - mu * Y2)
  }

  # Start from the exact isolated-patch equilibrium rather than its linear
  # approximation, which can send the solver to the trivial X = 0 root when
  # a*c/mu leaves the linear regime.
  R0_mean <- mean(R0)
  X_start <- pmin(0.99, pmax(0.001, mu * (R0_mean - 1) / (a * c_rate + R0_mean * mu)))
  Y_start <- pmin(0.99, pmax(0.0001, a * c_rate * X_start * eip_survival / (mu + a * c_rate * X_start)))

  sol <- nleqslv(rep(c(X_start, Y_start), each = 2), system,
                 method = "Broyden",
                 control = list(ftol = 1e-14, xtol = 1e-14, maxit = 500))

  if (sol$termcd != 1) {
    warning("compute_equilibrium_prevalence_RM(): solver may not have converged ",
            "(termcd = ", sol$termcd, ").")
  }

  list(X = c(sol$x[1], sol$x[2]), Y = c(sol$x[3], sol$x[4]))
}

#' Reproduction number of the Ross-Macdonald model from its next-generation matrix
#'
#' Independent check on the `R0` fed into [R0_to_m()]. The state ordering is
#' `(X_1, X_2, Y_1, Y_2)`, linearised around the disease-free equilibrium.
#'
#' @inheritParams compute_equilibrium_prevalence_RM
#' @return The spectral radius of the next-generation matrix, squared.
compute_R0_RM <- function(N, a, b, c_rate, r, mu, tau, R0, p_12, p_21) {

  m <- R0_to_m(R0, a, b, c_rate, r, mu, tau)

  p <- mobility_matrix(p_12, p_21)
  eip_survival <- exp(-mu * tau)
  D <- c(p[1, 1] * N[1] + p[2, 1] * N[2],
         p[1, 2] * N[1] + p[2, 2] * N[2])

  F_mat <- matrix(0, 4, 4)
  for (i in 1:2) {
    for (j in 1:2) {
      F_mat[i, 2 + j] <- a * b * p[i, j] * m[j]
      F_mat[2 + i, j] <- a * c_rate * eip_survival * p[j, i] * N[j] / D[i]
    }
  }

  V_mat <- diag(c(r, r, mu, mu))
  NGM <- F_mat %*% solve(V_mat)

  max(abs(eigen(NGM)$values))^2
}
