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

tar_load(starts_with("cv_fit_cutoff2010"))
tab2010 <- tribble(
	~title, ~fit,
	"B-spline ($d=2$, $K=5$)", cv_fit_cutoff2010_spline_2_5,
	"B-spline ($d=2$, $K=7$)", cv_fit_cutoff2010_spline_2_7,
	"B-spline ($d=3$, $K=5$)", cv_fit_cutoff2010_spline_3_5,
	"B-spline ($d=2$, $K=7$)", cv_fit_cutoff2010_spline_3_7,
	"Logistic", cv_fit_cutoff2010_logistic_2_7
) %>%
	mutate(results = map(fit, analyze_fit, only_most_recent = TRUE)) %>%
	select(-fit) %>%
	unnest(results)

tab2010_a <- tab2010 %>%
  select(title, below, within, above, ci_width, median_error, median_abs_error) %>%
	mutate_if(is.numeric, ~signif(. * 100, 3)) %>%
	mutate_at(vars(below, within, above), ~paste0(., "\\%"))

tab2010_b <- tab2010 %>%
  select(title, below_q5, below_q10, below_q25, below_q50, above_q50, above_q75, above_q90, above_q95) %>%
	mutate_if(is.numeric, ~paste0(signif(. * 100, 3), "\\%"))

#tar_load(starts_with("cv_fit_cutoff2015"))
tab2015 <- tribble(
	~title, ~fit,
	"B-spline ($d=2$, $K=5$)", cv_fit_cutoff2015_spline_2_5,
	"B-spline ($d=2$, $K=7$)", cv_fit_cutoff2015_spline_2_7,
	"B-spline ($d=3$, $K=5$)", cv_fit_cutoff2015_spline_3_5,
	"B-spline ($d=2$, $K=7$)", cv_fit_cutoff2015_spline_3_7,
	"Logistic", cv_fit_cutoff2015_logistic_2_7
) %>%
	mutate(results = map(fit, analyze_fit, only_most_recent = TRUE)) %>%
	select(-fit) %>%
	unnest(results)

tab2015_a <- tab2015 %>%
  select(title, below, within, above, ci_width, median_error, median_abs_error) %>%
	mutate_if(is.numeric, ~signif(. * 100, 3)) %>%
	mutate_at(vars(below, within, above), ~paste0(., "\\%"))

tab2015_b <- tab2015 %>%
  select(title, below_q5, below_q10, below_q25, below_q50, above_q50, above_q75, above_q90, above_q95) %>%
	mutate_if(is.numeric, ~paste0(signif(. * 100, 3), "\\%"))


print("Loading random hold-out runs...")
tar_load(starts_with("cv_fit_holdout_spline_2_5"))
tar_load(starts_with("cv_fit_holdout_spline_2_7"))
tar_load(starts_with("cv_fit_holdout_spline_3_5"))
tar_load(starts_with("cv_fit_holdout_spline_3_7"))
tar_load(starts_with("cv_fit_holdout_logistic_2_7"))
print("Finished loading...")

tab_random <- tribble(
	~index, ~title, ~fit,
	1, "B-spline ($d=2$, $K=5$)", cv_fit_holdout_spline_2_5_1,
	1, "B-spline ($d=2$, $K=7$)", cv_fit_holdout_spline_2_7_1,
	1, "B-spline ($d=3$, $K=5$)", cv_fit_holdout_spline_3_5_1,
	1, "B-spline ($d=3$, $K=7$)", cv_fit_holdout_spline_3_7_1,
	1, "Logistic", cv_fit_holdout_logistic_2_7_1,

	2, "B-spline ($d=2$, $K=5$)", cv_fit_holdout_spline_2_5_2,
	2, "B-spline ($d=2$, $K=7$)", cv_fit_holdout_spline_2_7_2,
	2, "B-spline ($d=3$, $K=5$)", cv_fit_holdout_spline_3_5_2,
	2, "B-spline ($d=2$, $K=7$)", cv_fit_holdout_spline_3_7_2,
	2, "Logistic", cv_fit_holdout_logistic_2_7_2,

	3, "B-spline ($d=2$, $K=5$)", cv_fit_holdout_spline_2_5_3,
	3, "B-spline ($d=2$, $K=7$)", cv_fit_holdout_spline_2_7_3,
	3, "B-spline ($d=3$, $K=5$)", cv_fit_holdout_spline_3_5_3,
	3, "B-spline ($d=2$, $K=7$)", cv_fit_holdout_spline_3_7_3,
	3, "Logistic", cv_fit_holdout_logistic_2_7_3,

	4, "B-spline ($d=2$, $K=5$)", cv_fit_holdout_spline_2_5_4,
	4, "B-spline ($d=2$, $K=7$)", cv_fit_holdout_spline_2_7_4,
	4, "B-spline ($d=3$, $K=5$)", cv_fit_holdout_spline_3_5_4,
	4, "B-spline ($d=2$, $K=7$)", cv_fit_holdout_spline_3_7_4,
	4, "Logistic", cv_fit_holdout_logistic_2_7_4,

	5, "B-spline ($d=2$, $K=5$)", cv_fit_holdout_spline_2_5_5,
	5, "B-spline ($d=2$, $K=7$)", cv_fit_holdout_spline_2_7_5,
	5, "B-spline ($d=3$, $K=5$)", cv_fit_holdout_spline_3_5_5,
	5, "B-spline ($d=2$, $K=7$)", cv_fit_holdout_spline_3_7_5,
	5, "Logistic", cv_fit_holdout_logistic_2_7_5
) %>%
	mutate(results = map(fit, analyze_fit, only_most_recent = TRUE)) %>%
	select(-fit) %>%
	unnest(results) %>%
	ungroup() %>%
	group_by(title) %>%
	summarize_if(is.numeric, mean)

tab_random_a <- tab_random %>%
  select(title, below, within, above, ci_width, median_error, median_abs_error) %>%
	mutate_if(is.numeric, ~signif(. * 100, 3)) %>%
	mutate_at(vars(below, within, above), ~paste0(., "\\%"))

tab_random_b <- tab_random %>%
  select(title, below_q5, below_q10, below_q25, below_q50, above_q50, above_q75, above_q90, above_q95) %>%
	mutate_if(is.numeric, ~paste0(signif(. * 100, 3), "\\%"))
