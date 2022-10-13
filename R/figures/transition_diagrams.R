library(tidyverse)
library(splines)
library(patchwork)

simulate_rate <- function(Tn, t_star, Omega, P_tilde, f, epsilon = rep(0, Tn)) {
  y <- numeric(Tn)
  
  y[t_star] <- Omega
  for(i in (t_star + 1):Tn) {
    y[i] <- boot::inv.logit(boot::logit(y[i - 1]) + f(y[i - 1] / P_tilde) + epsilon[i])
  }
  
  for(i in (t_star - 1):1) {
    y[i] <- boot::inv.logit(boot::logit(y[i + 1]) - f(y[i + 1] / P_tilde) - epsilon[i + 1])
  }
  
  tibble(year = 1:Tn, y0 = y)
}

spline_transition_function <- function(coefs, num_knots = 7, spline_degree = 2, constraints = TRUE) {
  if(constraints == TRUE) {
    knots <- c(seq(0, 1, length.out = num_knots), 1000)
    grid <- c(seq(from = 0, to = 1, by = .01/2), 1000) # generating inputs
  }
  else {
    knots <- seq(0, 1, length.out = num_knots)
    grid <- seq(0, 1, 0.01 / 2)
  }
  
  B <- t(bs(grid, knots = knots, degree = spline_degree, intercept = FALSE))
  B <- B[1:(nrow(B) - 1), ]
  num_grid <- length(grid)
  num_basis <- nrow(B)
  
  if(constraints == TRUE) {
    y <- c(coefs, rep(0, spline_degree + 1)) %*% B
  }
  else {
    y <- c(coefs %*% B)
  }
  approxfun(grid, y)
}


simulate_rate_tfr <- function(Tn, f, Omega) {
  y <- numeric(Tn)
  
  y[1] <- Omega
  for(i in 2:Tn) {
    y[i] <- y[i - 1] + f(y[i - 1])
  }
  tibble(year = 1:Tn, y0 = y)
}

plot_tfr_transition_function_example <- function(d, delta1, delta2, delta3, delta4) {
  grid <- seq(0, 8, 0.01)
  f <- Vectorize(function(eta) {
    ifelse(
      eta > 1, 
      d / (1 + exp(-2*log(9) / delta1 * (eta - (delta1 + delta2 + delta3 + delta4) + 0.5 * delta1))) - 
        d / (1 + exp(-2*log(9) / delta3 * (eta - delta4 - 0.5 * delta3))),
      0
    )
  })
  
  transition_function <- tibble(
    x = grid,
    y = f(grid)
  )
  
  data <- simulate_rate_tfr(Tn = 20, f = f, 6.5)
  
  arrows <- arrow(length = unit(0.1, "inches"), ends = "both")
  arrow_y <- -0.02
  arrow_y2 <- 0.06
  segment_size <- 0.2
  segment_lwd <- 1
  adj <- 0.07
  
  p1 <- ggplot(transition_function) +
    # f curve
    geom_line(aes(x, y)) +
    
    geom_segment(x = -8, xend = 0, y = d, yend = d, lty = 2, color = "blue", size = segment_size, lwd = segment_lwd) +
    annotate(x = 0, y = -d - 0.03, label = expression(-d[c]), geom = "text", color = "blue") +
    
    geom_segment(xend = 0, x = -(delta4 - adj), y = arrow_y, yend = arrow_y, lty = 1, color = "blue", arrow = arrows, size = segment_size) +
    annotate(x = delta4 / 2, y = arrow_y2, label = expression(Delta[c4]), geom = "text", color = "blue") +
    
    geom_segment(xend = -(delta4 + adj), x = -(delta4 + delta3 - adj), y = arrow_y, yend = arrow_y, lty = 1, color = "blue", arrow = arrows, size = segment_size) +
    annotate(x = delta4 + delta3 / 2, y = arrow_y2, label = expression(Delta[c3]), geom = "text", color = "blue") +
    
    geom_segment(xend = -(delta4 + delta3 + adj), x = -(delta4 + delta3 + delta2 - adj), y = arrow_y, yend = arrow_y, lty = 1, color = "blue", arrow = arrows, size = segment_size) +
    annotate(x = delta4 + delta3 + delta2 / 2, y = arrow_y2, label = expression(Delta[c2]), geom = "text", color = "blue") +
    
    geom_segment(xend = -(delta4 + delta3 + delta2 + adj), x = -(delta4 + delta3 + delta2 + delta1 - adj), y = arrow_y, yend = arrow_y, lty = 1, color = "blue", arrow = arrows, size = segment_size) +
    annotate(x = delta4 + delta3 + delta2 + delta1 / 2, y = arrow_y2, label = expression(Delta[c1]), geom = "text", color = "blue") +
    
    labs(x = expression(eta[ct]~(reversed)), y = expression(f[dl]~(reversed)), title = "Transition Function") +
    scale_y_reverse() +
    scale_x_reverse(limits = c(8.5, 0)) +
    cowplot::theme_cowplot() +
    theme(legend.position = "none", plot.title = element_text(face = "italic"))
  
  p2 <- ggplot(data, aes(year, y0)) +
    geom_line() +
    #geom_segment(x = t_star, xend = t_star, y = -1, yend = Omega, lty = 2, color = "blue") +
    #annotate(x = t_star + 5, y = Omega, label = expression(Omega[c]), geom = "text", color = "blue") +
    cowplot::theme_cowplot() +
    labs(x = expression(t), y = expression(eta[ct]), title = "TFR") +
    theme(legend.position = "none", plot.title = element_text(face = "italic"))
  
  
  list(
    plot = ggpubr::ggarrange(p1, p2, ncol = 2),
    data = data,
    transition_function = transition_function
  )
}


plot_fpem_transition_function_example <- function(omega, t_star = 40, P = 0.8, Omega = 0.2, show_P = TRUE) {
  grid <- seq(0.001, 0.999, 0.01)
  f <- Vectorize(function(eta) {
    ifelse(
      eta < P, 
      (eta - P) * omega / (P * (eta - 1)),
      #boot::logit(P * boot::inv.logit(boot::logit(eta / P) + omega)) - boot::logit(eta), 
      0
    )
  })
  
  transition_function <- tibble(
    x = grid,
    y = f(grid)
  )
  
  data <- simulate_rate(Tn = 80, t_star = t_star, P = 1, Omega = Omega, f = f)
  
  p1 <- ggplot(transition_function, aes(x, y)) +
    # f curve
    geom_line() +
    
    labs(x = expression(eta[ct]), y = expression(f[FPEM]), title = "Transition Function") +
    cowplot::theme_cowplot() +
    theme(legend.position = "none", plot.title = element_text(face = "italic"))
  
  p2 <- ggplot(data, aes(year, y0)) +
    geom_line() +
    geom_segment(x = t_star, xend = t_star, y = -1, yend = Omega, lty = 2, color = "blue") +
    annotate(x = t_star + 5, y = Omega, label = expression(Omega[c]), geom = "text", color = "blue") +
    cowplot::theme_cowplot() +
    labs(x = expression(t), y = expression(eta[ct]), title = "mCPR") +
    theme(legend.position = "none", plot.title = element_text(face = "italic"))
  
  
  if(show_P == TRUE) {
    p1 <- p1 + geom_vline(xintercept = P, lty = 2, color = "blue") +
      annotate(x = P + 0.07, y = 0.075, label = expression(P[c]^u), geom = "text", color = "blue")
    
    p2 <- p2 + geom_hline(yintercept = P, lty = 2, color = "blue") +
      annotate(x = 80, y = 0.85, label = expression(P[c]^u), geom = "text", color = "blue")
  }
  
  list(
    plot = ggpubr::ggarrange(p1, p2, ncol = 2),
    data = data,
    transition_function = transition_function
  )
}

plot_transition_function_example <- function(beta, knots = 7, spline_degree = 2, t_star = 40, P = 0.8, Omega = 0.2, constraints = TRUE, show_P = TRUE) {
  grid <- seq(0, 1, 0.01)
  f <- spline_transition_function(beta, num_knots = knots, spline_degree = spline_degree, constraints = constraints)
  
  transition_function <- tibble(
    x = grid,
    y = f(grid / P)
  )
  
  basis_functions <- map(1:length(beta), function(i) {
    b <- rep(0, length(beta))
    #b[i] <- 0.02
    b[i] <- beta[i]
    tibble(
      x = grid,
      y = spline_transition_function(b, num_knots = knots, spline_degree = spline_degree, constraints = constraints)(grid/P)
    )
  })  %>%
    bind_rows(.id = "j")
  
  data <- simulate_rate(Tn = 80, t_star = t_star, P = P, Omega = Omega, f = f)
  
  p1 <- ggplot(transition_function, aes(x, y)) +
    # Basis functions
    geom_line(data = basis_functions, aes(x = x, y = y, group = j), color = "gray") +
    # Spline curve
    geom_line() +
    labs(x = expression(eta[ct]), y = expression(f[b]), title = "Transition Function") +
    cowplot::theme_cowplot() +
    theme(legend.position = "none", plot.title = element_text(face = "italic"))
  
  p2 <- ggplot(data, aes(year, y0)) +
    geom_line() +
    geom_segment(x = t_star, xend = t_star, y = -1, yend = Omega, lty = 2, color = "blue") +
    annotate(x = t_star + 5, y = Omega, label = expression(Omega[c]), geom = "text", color = "blue") +
    cowplot::theme_cowplot() +
    labs(x = expression(t), y = expression(eta[ct]), title = "CPR") +
    theme(legend.position = "none", plot.title = element_text(face = "italic"))
  
  
  if(show_P == TRUE) {
    p1 <- p1 + geom_vline(xintercept = P, lty = 2, color = "blue") +
      annotate(x = P + 0.07, y = 0.075, label = expression(P[c]^u), geom = "text")
    
    p2 <- p2 + geom_hline(yintercept = P, lty = 2, color = "blue") +
      annotate(x = 80, y = 0.85, label = expression(P[c]^u), geom = "text", color = "blue")
  }
  
  #ggpubr::ggarrange(p1, p2, ncol = 2)
  list(p1, p2)
}

plot_splines <- function(knots = 7, spline_degree = 2) {
  grid <- seq(0, 1, 0.01)
  beta <- rep(1, knots - 1)
  constraints <- TRUE
  f <- spline_transition_function(beta, num_knots = knots, spline_degree = spline_degree, constraints = constraints)
  
  basis_functions <- map(1:length(beta), function(i) {
    b <- rep(0, length(beta))
    b[i] <- 1
    tibble(
      x = grid,
      y = spline_transition_function(b, num_knots = knots, spline_degree = spline_degree, constraints = constraints)(grid)
    )
  })  %>%
    bind_rows(.id = "j")
  
  ggplot(basis_functions, aes(x, y)) +
    geom_line(aes(group = j)) +
    cowplot::theme_cowplot() +
    labs(x = expression(eta[ct] / P[c]^u), y = expression(f[b])) +
    theme(plot.title = element_text(face = "plain", hjust = 0.5))
}

((plot_splines(5, 2) + ggtitle("K = 5, d = 2") + labs(x = "")) +
    (plot_splines(7, 2) + ggtitle("K = 7, d = 2")) + labs(x = "", y = "")) /
  ((plot_splines(5, 3) + ggtitle("K = 5, d = 3")) + 
     (plot_splines(7, 3) + ggtitle("K = 7, d = 3")) + labs(y = "")) +
  plot_layout()
ggsave("plots/spline_setups.pdf", height = 6, width = 8)

plot_transition_function_example(beta = c(0.15, 0.15, 0.1, 0.12, 0.08, 0.05, 0.1, 0.1), Omega = 0.3, P = 1, show_P = FALSE, constraints = FALSE)
ggsave("plots/spline_transition_function_example_unconstrained.pdf", height = 4)

plots <- plot_transition_function_example(beta = c(0.05, 0.07, 0.10, 0.06, 0.05, 0.06), P = 1)
ggsave(plot = plots[[1]], "plots/spline_transition_function_example_1.pdf", height = 4, width = 5)
ggsave(plot = plots[[2]], "plots/spline_transition_function_example_2.pdf", height = 4, width = 5)

plot_transition_function_example(beta = rep(0.08, 6), knots = 7)
ggsave("plots/spline_logistic_transition_function_example.pdf", height = 4)

fpem <- plot_fpem_transition_function_example(omega = 0.15, P = 0.8, Omega = 0.2)
fpem$plot
ggsave("plots/fpem_transition_function_example.pdf", height = 4)

P <- 0.8
omega <- 0.15

final_logistic <- readRDS("server-output/final_logistic")
x <- final_logistic$stan_data$spline_maxima[1:6]

rates <- Vectorize(function(eta, P, omega) (eta - P) * omega / (P * (eta - 1)), vectorize.args = "eta")

b <- rates(x, P, omega)

p <- plot_transition_function_example(beta = b, knots = 7, t_star = 40, Omega = 0.2)
p[[1]] <- p[[1]] + 
  geom_line(data = fpem$transition_function, aes(x = x, y = y), lty = 2) +
  annotate("text", x = 0.61, y = 0.11, label = expression(f[FPEM]), size = 4) +
  annotate("text", x = 0.49, y = 0.08, label = "Approximation", hjust = 1, size = 4)
p[[1]] + p[[2]]

ggsave("plots/spline_logistic_transition_function_example.pdf", height = 4)


plot_tfr_transition_function_example(1.2, 1, 2, 2, 2)$plot
ggsave("plots/tfr_transition_function_example.pdf", height = 5)
