library(targets)

analyze_fit <- function(fit, only_most_recent = FALSE) {
	data_tibble <- tibble::tibble(i = 1:length(fit$held_out), held_out = fit$held_out, y = fit$stan_data$y, c = fit$stan_data$country, t = fit$stan_data$time)

  data <- fit$posteriors$generated_quantities %>%
    ungroup() %>%
    left_join(data_tibble, by = "i") %>%
    left_join(fit$country_index, by = "c") %>%
    left_join(fit$time_index, by = "t")

	n_not_left_out <- fit$data %>%
	  bind_cols(tibble(held_out = fit$held_out)) %>%
		group_by(name_country) %>%
		summarize(n_left_out = sum(held_out), n_not_left_out = sum(!held_out))

	countries_with_all_left_out <- n_not_left_out %>%
		filter(n_not_left_out == 0)

	data <- data %>%
		filter(.data$held_out == 1)

	data <- data %>%
		filter(!(name_country %in% countries_with_all_left_out$name_country))

  if(only_most_recent == TRUE) {
		data <- data %>%
		  ungroup() %>%
			group_by(name_country) %>%
			filter(year == max(year)) %>%
		  ungroup()
  }
	
  results <- data %>%
    group_by(.data$i, .data$y) %>%
    summarize(
      q2.5  = quantile(y_pred, 0.025),
      q5    = quantile(y_pred, 0.05),
      q10   = quantile(y_pred, 0.10),
      q25   = quantile(y_pred, 0.25),
      q50   = quantile(y_pred, 0.50),
      q75   = quantile(y_pred, 0.75),
      q90   = quantile(y_pred, 0.90),
      q95   = quantile(y_pred, 0.95),
      q97.5 = quantile(y_pred, 0.975),
    ) %>%
    mutate(
			below = y < q2.5,
			above = y > q97.5,
			within = .data$y >= q2.5 & .data$y <= q97.5,
			ci_width = q97.5 - q2.5,
			error = y - q50,
			below_q5 = y < q5,
			below_q10 = y < q10,
			below_q25 = y < q25,
			below_q50 = y < q50,
			above_q50 = y > q50,
			above_q75 = y > q75,
			above_q90 = y > q90,
			above_q95 = y > q95,
		) %>%
		ungroup() %>%
		summarize(
		  below_q5  = mean(below_q5),
		  below_q10 = mean(below_q10),
		  below_q25 = mean(below_q25),
		  below_q50 = mean(below_q50),
		  above_q50 = mean(above_q50),
		  above_q75 = mean(above_q75),
		  above_q90 = mean(above_q90),
		  above_q95 = mean(above_q95),
		  
			below = mean(below), 
			above = mean(above), 
			within = mean(within), 
			ci_width = mean(ci_width),
			mean_error = mean(error),
			mean_abs_error = mean(abs(error)),
			median_error = median(error),
			median_abs_error = median(abs(error)),
			n = n()
		)

	results
}

