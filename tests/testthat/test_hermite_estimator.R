context("hermite_estimator")
library(hermiter)
library(magrittr)

test_that("hermite_estimator constructor returns correct class", {
  hermite_est <- hermite_estimator(N = 10, standardize = TRUE)
  expect_is(hermite_est, "hermite_estimator")
})

test_that("batch updates of hermite_estimator work as expected", {
  test_observations <- c(
    0.3719336,
    0.8484731,
    -2.663014,
    -2.028416,
    -2.429764,
    -1.677079,
    1.434098,
    -1.453405,
    -1.133242,
    0.963844,
    1.46609,
    0.371181,
    2.135272,
    -0.7328963,
    0.8465673,
    -2.168731,
    -0.269106,
    -1.885169,
    -0.07963116,
    0.1244462,
    0.1165929,
    1.460038,
    -0.06261991,
    0.07363522,
    0.03415375,
    -1.65761,
    2.058115,
    0.9413341,
    -1.759675,
    2.214421
  )
  target_coeff_vec_standardized <-
    c(
      0.5076928,
      0.03212483,
      0.08900023,
      -0.0701138,
      -0.07802765,
      0.08559936,
      -0.009535837,
      -0.0861154,
      0.06687268,
      0.07748299,
      -0.08368591
    )
  
  target_coeff_vec_unstandardized <-
    c(
      0.3918339,
      0.04782401,
      0.1507684,
      -0.1122261,
      0.1331726,
      0.02097852,
      -0.1256438,
      0.009395677,
      0.05290998,
      -0.01767049,
      -0.02434694
    )
  hermite_est <- hermite_estimator(N = 10, standardize = TRUE)
  hermite_est <- hermite_est %>% update_batch(test_observations)
  expect_equal(target_coeff_vec_standardized,
               get_coefficients(hermite_est),
               tolerance = 1e-07)
  
  hermite_est <- hermite_estimator(N = 10, standardize = FALSE)
  hermite_est <- hermite_est %>% update_batch(test_observations)
  expect_equal(target_coeff_vec_unstandardized,
               get_coefficients(hermite_est),
               tolerance = 1e-07)
})

test_that("sequential updates of hermite_estimator work as expected", {
  test_observations <- c(
    0.3719336,
    0.8484731,
    -2.663014,
    -2.028416,
    -2.429764,
    -1.677079,
    1.434098,
    -1.453405,
    -1.133242,
    0.963844,
    1.46609,
    0.371181,
    2.135272,
    -0.7328963,
    0.8465673,
    -2.168731,
    -0.269106,
    -1.885169,
    -0.07963116,
    0.1244462,
    0.1165929,
    1.460038,
    -0.06261991,
    0.07363522,
    0.03415375,
    -1.65761,
    2.058115,
    0.9413341,
    -1.759675,
    2.214421
  )
  target_coeff_vec_standardized <-
    c(
      0.52421532,
      0.08216913,
      0.02253855,
      0.01674095,
      -0.06203412,
      -0.01720743,
      0.01073874,
      -0.02520607,
      0.03651330,
      0.06876952,
      -0.05918178
    )
  target_coeff_vec_unstandardized <-
    c(
      0.3918339,
      0.04782401,
      0.1507684,
      -0.1122261,
      0.1331726,
      0.02097852,
      -0.1256438,
      0.009395677,
      0.05290998,
      -0.01767049,
      -0.02434694
    )
  hermite_est <- hermite_estimator(N = 10, standardize = TRUE)
  for (idx in seq_along(test_observations)) {
    hermite_est <-
      hermite_est %>% update_sequential(test_observations[idx])
  }
  expect_equal(target_coeff_vec_standardized,
               get_coefficients(hermite_est),
               tolerance = 1e-07)
  
  hermite_est <- hermite_estimator(N = 10, standardize = FALSE)
  for (idx in seq_along(test_observations)) {
    hermite_est <-
      hermite_est %>% update_sequential(test_observations[idx])
  }
  expect_equal(target_coeff_vec_unstandardized,
               get_coefficients(hermite_est),
               tolerance = 1e-07)
})

test_that("sequential updates of exponentially weighted hermite_estimator 
          work as expected",
          {
            test_observations <- c(
              0.3719336,
              0.8484731,
              -2.663014,
              -2.028416,
              -2.429764,
              -1.677079,
              1.434098,
              -1.453405,
              -1.133242,
              0.963844,
              1.46609,
              0.371181,
              2.135272,
              -0.7328963,
              0.8465673,
              -2.168731,
              -0.269106,
              -1.885169,
              -0.07963116,
              0.1244462,
              0.1165929,
              1.460038,
              -0.06261991,
              0.07363522,
              0.03415375,
              -1.65761,
              2.058115,
              0.9413341,
              -1.759675,
              2.214421
            )
            target_coeff_vec_standardized <-
              c(
                0.3720845,
                0.0329111,
                0.07880838,
                -0.0480489,
                -0.001239091,
                0.03989803,
                -0.08362694,
                -0.04300433,
                0.1199047,
                0.04287386,
                -0.1147223
              )
            target_coeff_vec_unstandardized <-
              c(
                0.3140483,
                0.03968127,
                0.1031914,
                -0.0542355,
                0.1318214,
                0.03501422,
                -0.1379236,
                -0.005167008,
                0.06536224,
                -0.02967193,
                -0.02234358
              )
            hermite_est <-
              hermite_estimator(N = 10,
                                standardize = TRUE,
                                exp_weight_lambda = 0.05)
            for (idx in seq_along(test_observations)) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations[idx])
            }
            expect_equal(target_coeff_vec_standardized,
                         get_coefficients(hermite_est),
                         tolerance = 1e-07)
            
            hermite_est <-
              hermite_estimator(N = 10,
                                standardize = FALSE,
                                exp_weight_lambda = 0.05)
            for (idx in seq_along(test_observations)) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations[idx])
            }
            expect_equal(target_coeff_vec_unstandardized,
                         get_coefficients(hermite_est),
                         tolerance = 1e-07)
          })

test_that("hermite_estimators combine consistently", {
  test_observations <- c(
    0.3719336,
    0.8484731,
    -2.663014,
    -2.028416,
    -2.429764,
    -1.677079,
    1.434098,
    -1.453405,
    -1.133242,
    0.963844,
    1.46609,
    0.371181,
    2.135272,
    -0.7328963,
    0.8465673,
    -2.168731,
    -0.269106,
    -1.885169,
    -0.07963116,
    0.1244462,
    0.1165929,
    1.460038,
    -0.06261991,
    0.07363522,
    0.03415375,
    -1.65761,
    2.058115,
    0.9413341,
    -1.759675,
    2.214421
  )
  hermite_est <-
    hermite_estimator(N = 10, standardize = FALSE) %>% 
    update_batch(test_observations)
  hermite_est_1 <-
    hermite_estimator(N = 10, standardize = FALSE) %>% 
    update_batch(test_observations[1:10])
  hermite_est_2 <-
    hermite_estimator(N = 10, standardize = FALSE) %>% 
    update_batch(test_observations[11:20])
  hermite_est_3 <-
    hermite_estimator(N = 10, standardize = FALSE) %>% 
    update_batch(test_observations[21:30])
  hermite_comb <-
    combine_hermite(list(hermite_est_1, hermite_est_2, hermite_est_3))
  expect_equal(hermite_est, hermite_comb, tolerance = 1e-07)
})

test_that("probability density estimation works as expected", {
  test_observations <- c(
    0.3719336,
    0.8484731,
    -2.663014,
    -2.028416,
    -2.429764,
    -1.677079,
    1.434098,
    -1.453405,
    -1.133242,
    0.963844,
    1.46609,
    0.371181,
    2.135272,
    -0.7328963,
    0.8465673,
    -2.168731,
    -0.269106,
    -1.885169,
    -0.07963116,
    0.1244462,
    0.1165929,
    1.460038,
    -0.06261991,
    0.07363522,
    0.03415375,
    -1.65761,
    2.058115,
    0.9413341,
    -1.759675,
    2.214421
  )
  x <- seq(-2, 2, by = 0.5)
  hermite_est <-
    hermite_estimator(N = 10, standardize = FALSE) %>% 
    update_batch(test_observations)
  pdf_vals <- hermite_est %>% dens(x)
  target_pdf_vals_unstandardized <-
    c(
      0.2700498,
      0.2387219,
      0.03433206,
      0.124966,
      0.3581052,
      0.2772852,
      0.170135,
      0.2005586,
      0.1451746
    )
  expect_equal(pdf_vals, target_pdf_vals_unstandardized, tolerance = 1e-07)
  
  hermite_est <-
    hermite_estimator(N = 10, standardize = TRUE) %>%
    update_batch(test_observations)
  pdf_vals <- hermite_est %>% dens(x)
  target_pdf_vals_standardized <-
    c(
      0.2495028,
      0.1759633,
      0.09314036,
      0.1551446,
      0.2806341,
      0.2997986,
      0.2198349,
      0.1576933,
      0.1306844
    )
  expect_equal(pdf_vals, target_pdf_vals_standardized, tolerance = 1e-07)
  
  pdf_vals <- hermite_est %>% dens(x, clipped=T)
  expect_equal(pdf_vals, target_pdf_vals_standardized, tolerance = 1e-07)
  
  hermite_est <-
    hermite_estimator(N = 10,
                      standardize = FALSE,
                      exp_weight_lambda = 0.1)
  for (idx in seq_along(test_observations)) {
    hermite_est <-
      hermite_est %>% update_sequential(test_observations[idx])
  }
  pdf_vals <- hermite_est %>% dens(x)
  target_pdf_vals_unstandardized <-
    c(
      0.2632709,
      0.2498718,
      -0.0264928,
      0.1122043,
      0.4209524,
      0.2511466,
      0.07670513,
      0.1983007,
      0.2352683
    )
  expect_equal(pdf_vals, target_pdf_vals_unstandardized, tolerance = 1e-07)
  
  hermite_est <-
    hermite_estimator(N = 10,
                      standardize = FALSE,
                      exp_weight_lambda = 0.1)
  for (idx in seq_along(test_observations)) {
    hermite_est <-
      hermite_est %>% update_sequential(test_observations[idx])
  }
  pdf_vals <- hermite_est %>% dens(x, clipped = TRUE)
  target_pdf_vals_unstandardized_clipped <-
    c(
      0.2632709,
      0.2498718,
      1e-8,
      0.1122043,
      0.4209524,
      0.2511466,
      0.07670513,
      0.1983007,
      0.2352683
    )
  expect_equal(pdf_vals, target_pdf_vals_unstandardized_clipped,
               tolerance = 1e-07)
  
  hermite_est <-
    hermite_estimator(N = 10,
                      standardize = TRUE,
                      exp_weight_lambda = 0.1)
  for (idx in seq_along(test_observations)) {
    hermite_est <-
      hermite_est %>% update_sequential(test_observations[idx])
  }
  pdf_vals <- hermite_est %>% dens(x)
  target_pdf_vals_standardized <-
    c(
      0.224348,
      0.2392328,
      0.04693996,
      0.03530093,
      0.2728397,
      0.3533988,
      0.1876257,
      0.1247126,
      0.2253963
    )
  expect_equal(pdf_vals, target_pdf_vals_standardized, tolerance = 1e-07)
  hermite_est <-
    hermite_estimator(N = 10)
  expect_equal(hermite_est %>% dens(x),NA)
})

test_that("cumulative distribution function estimation works as expected",
          {
            test_observations <- c(
              0.3719336,
              0.8484731,
              -2.663014,
              -2.028416,
              -2.429764,
              -1.677079,
              1.434098,
              -1.453405,
              -1.133242,
              0.963844,
              1.46609,
              0.371181,
              2.135272,
              -0.7328963,
              0.8465673,
              -2.168731,
              -0.269106,
              -1.885169,
              -0.07963116,
              0.1244462,
              0.1165929,
              1.460038,
              -0.06261991,
              0.07363522,
              0.03415375,
              -1.65761,
              2.058115,
              0.9413341,
              -1.759675,
              2.214421
            )
            hermite_est <-
              hermite_estimator(N = 10, standardize = FALSE) %>%
              update_batch(test_observations)
            cdf_from_pdf <- stats::integrate(
              f = function(x) {
                hermite_est %>% dens(x)
              },
              lower = -Inf,
              upper = 0.5
            )$value
            cdf_est <- hermite_est %>% cum_prob(0.5)
            expect_equal(cdf_est, 0.6549575, tolerance = 1e-07)
            expect_equal(cdf_from_pdf, cdf_est, tolerance = 1e-07)
            
            cdf_est <- hermite_est %>% cum_prob(0.5, clipped=TRUE)
            expect_equal(cdf_est, 0.6549575, tolerance = 1e-07)
            
            cdf_est <- hermite_est %>% cum_prob(3, clipped=TRUE)
            expect_equal(cdf_est, 1, tolerance = 1e-07)
            
            hermite_est <-
              hermite_estimator(N = 10, standardize = TRUE) %>%
              update_batch(test_observations)
            cdf_from_pdf <- stats::integrate(
              f = function(x) {
                hermite_est %>% dens(x)
              },
              lower = -Inf,
              upper = 0.5
            )$value
            cdf_est <- hermite_est %>% cum_prob(0.5)
            expect_equal(cdf_est, 0.6013645, tolerance = 1e-07)
            expect_equal(cdf_from_pdf, cdf_est, tolerance = 1e-07)
            
            hermite_est <-
              hermite_estimator(N = 10,
                                standardize = FALSE,
                                exp_weight_lambda = 0.1)
            for (idx in seq_along(test_observations)) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations[idx])
            }
            cdf_from_pdf <- stats::integrate(
              f = function(x) {
                hermite_est %>% dens(x)
              },
              lower = -Inf,
              upper = 0.5
            )$value
            cdf_est <- hermite_est %>% cum_prob(0.5)
            expect_equal(cdf_est, 0.6132811, tolerance = 1e-07)
            expect_equal(cdf_from_pdf, cdf_est, tolerance = 1e-07)
            
            hermite_est <-
              hermite_estimator(N = 10,
                                standardize = TRUE,
                                exp_weight_lambda = 0.1)
            for (idx in seq_along(test_observations)) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations[idx])
            }
            cdf_from_pdf <- stats::integrate(
              f = function(x) {
                hermite_est %>% dens(x)
              },
              lower = -Inf,
              upper = 0.5
            )$value
            cdf_est <- hermite_est %>% cum_prob(0.5)
            expect_equal(cdf_est, 0.4344541, tolerance = 1e-07)
            expect_equal(cdf_from_pdf, cdf_est, tolerance = 1e-07)
            hermite_est <-
              hermite_estimator(N = 10)
            expect_equal(hermite_est %>% cum_prob(0.5),NA)
            
            hermite_est <-
              hermite_estimator(N = 10, standardize = FALSE) %>%
              update_batch(test_observations)
            cdf_from_pdf <- 1- stats::integrate(
              f = function(x) {
                hermite_est %>% dens(x)
              },
              lower = 1,
              upper = Inf
            )$value
            cdf_est <- hermite_est %>% cum_prob_quantile_helper(1)
            expect_equal(cdf_est, 0.769892, tolerance = 1e-06)
            expect_equal(cdf_from_pdf, cdf_est, tolerance = 1e-06)
          })

test_that("quantile estimation works as expected", {
  test_observations <- c(
    0.3719336,
    0.8484731,
    -2.663014,
    -2.028416,
    -2.429764,
    -1.677079,
    1.434098,
    -1.453405,
    -1.133242,
    0.963844,
    1.46609,
    0.371181,
    2.135272,
    -0.7328963,
    0.8465673,
    -2.168731,
    -0.269106,
    -1.885169,
    -0.07963116,
    0.1244462,
    0.1165929,
    1.460038,
    -0.06261991,
    0.07363522,
    0.03415375,
    -1.65761,
    2.058115,
    0.9413341,
    -1.759675,
    2.214421
  )
  hermite_est <-
    hermite_estimator(N = 10, standardize = TRUE) %>%
    update_batch(test_observations)
  quantiles_est <- hermite_est %>% quant(c(0.25, 0.5, 0.75))
  expect_equal(quantiles_est,
               c(-1.54823599,  0.04145506,  0.90889172),
               tolerance = 1e-07)
  expect_equal(hermite_estimator(N = 10, standardize = TRUE) %>%
                 quant(c(0.25, 0.5, 0.75)), NA)
  cum_prob_check <-
    hermite_est %>% cum_prob_quantile_helper((hermite_est %>% quant(0.75)
                                              - hermite_est$running_mean) /
                                              sqrt(hermite_est$running_variance 
                                                    / (hermite_est$num_obs -1))
    )
  expect_equal(cum_prob_check, 0.75, tolerance = 0.001)
  
  hermite_est <- hermite_estimator(N = 10, standardize = TRUE)
  for (idx in seq_along(test_observations)) {
    hermite_est <-
      hermite_est %>% update_sequential(test_observations[idx])
  }
  quantiles_est <- hermite_est %>% quant(c(0.25, 0.5, 0.75))
  expect_equal(quantiles_est,
               c(-0.9205331,  0.1404577,  1.1469203),
               tolerance = 1e-07)
  cum_prob_check <-
    hermite_est %>% cum_prob_quantile_helper((hermite_est %>% quant(0.75)
                                              - hermite_est$running_mean) /
                                               sqrt(hermite_est$running_variance
                                                    / (hermite_est$num_obs -1))
    )
  expect_equal(cum_prob_check, 0.75, tolerance = 0.001)
  
  hermite_est <-
    hermite_estimator(N = 10,
                      standardize = TRUE,
                      exp_weight_lambda = 0.1)
  for (idx in seq_along(test_observations)) {
    hermite_est <-
      hermite_est %>% update_sequential(test_observations[idx])
  }
  quantiles_est <- hermite_est %>% quant(c(0.25, 0.5, 0.75))
  expect_equal(quantiles_est, c(-1.396857,  0.360104,  1.437417), 
               tolerance = 1e-06)
  cum_prob_check <-
    hermite_est %>% cum_prob_quantile_helper((hermite_est %>% quant(0.75) -
                                            hermite_est$running_mean) /
                                              sqrt(hermite_est$running_variance)
    )
  expect_equal(cum_prob_check, 0.75, tolerance = 0.001)
})
