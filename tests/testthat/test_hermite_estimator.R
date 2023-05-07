## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

context("hermite_estimator_univar")
library(hermiter)
library(magrittr)

get_eps <- function(){
  return(1e-3)
}

test_that("hermite_estimator constructor returns correct class", {
  hermite_est <- hermite_estimator(N = 10, standardize = TRUE)
  expect_is(hermite_est, "hermite_estimator_univar")
})

test_that("error trapping for hermite_estimator work as expected", {
  expect_error(hermite_estimator(N = 10, est_type="other"))
  expect_error(hermite_estimator(N = "a", standardize = TRUE))
  expect_error(hermite_estimator(N = 80, standardize = TRUE))
  expect_error(hermite_estimator(N = -1, standardize = TRUE))
  expect_error(hermite_estimator(N = 10, standardize = 8))
  expect_error(hermite_estimator(N = 10, exp_weight_lambda = -0.5))
  expect_error(hermite_estimator(N = 10, exp_weight_lambda = 1.5))
  expect_error(hermite_estimator(N = 10, exp_weight_lambda = "a"))
  hermite_est <- hermite_estimator(N = 10, standardize = TRUE)
  expect_true(is.na(cum_prob(hermite_est, x=2)))
  expect_true(is.na(quant(hermite_est, p=0.5)))
  expect_error(update_sequential(hermite_est, "a"))
  expect_error(hermite_estimator(N = 10, standardize = TRUE,
                                 observations = c("a","b")))
  expect_error(dens(hermite_est, "a"))
  expect_error(dens(hermite_est, numeric(0)))
  expect_error(cum_prob(hermite_est, "a"))
  expect_error(cum_prob(hermite_est, numeric(0)))
  expect_error(quant(hermite_est, "a"))
  expect_error(quant(hermite_est, numeric(0)))
  expect_error(quant(hermite_est, 0.9, algorithm = "a"))
  expect_error(dens(hermite_est, c()))
  expect_error(cum_prob(hermite_est, c()))
  expect_error(quant(hermite_est, c()))
  expect_error(quant(hermite_est, p=2))
  expect_error(quant(hermite_est, p=-0.5))
  expect_error(merge_hermite(list()))
  hermite_est <- hermite_estimator(N = 10, standardize = FALSE, 
                                   observations = c(1,2,3))
  expect_error(quant(hermite_est, p=0.5))
  expect_error(spearmans(hermite_est))
  # Univariate estimator specific error trapping tests
  expect_error(hermite_estimator(N = 10, est_type="other"))
  expect_error(hermite_estimator(N = "a", standardize = TRUE))
  expect_error(hermite_estimator(N = 80, standardize = TRUE))
  expect_error(hermite_estimator(N = -1, standardize = TRUE))
  expect_error(hermite_estimator(N = 10, standardize = 8))
  expect_error(hermite_estimator(N = 10, exp_weight_lambda = -0.5))
  expect_error(hermite_estimator(N = 10, exp_weight_lambda = 1.5))
  expect_error(hermite_estimator(N = 10, exp_weight_lambda = "a"))
  hermite_est <- hermite_estimator(N = 10, standardize = TRUE)
  expect_true(is.na(cum_prob(hermite_est, x=2)))
  expect_true(is.na(quant(hermite_est, p=0.5)))
  expect_error(update_sequential(hermite_est, "a"))
  expect_error(hermite_estimator(N = 10, standardize = TRUE, 
                                        observations = c("a","b")))
  hermite_est <- hermite_estimator(N = 10, standardize = TRUE, 
                                          observations = c(1,2,3))
  expect_error(dens(hermite_est, "a"))
  expect_error(cum_prob(hermite_est, "a"))
  expect_error(quant(hermite_est, "a"))
  expect_error(dens(hermite_est, c()))
  expect_error(cum_prob(hermite_est, c()))
  expect_error(quant(hermite_est, c()))
  expect_error(merge_hermite_univar(list()))
  hermite_est <- hermite_estimator(N = 10, standardize = FALSE, 
                                          observations = c(1,2,3))
  expect_error(quant(hermite_est, p=0.5))
  expect_error(spearmans(hermite_est))
  expect_error(kendall(hermite_est))
  hermite_est_1 <- hermite_estimator(N = 10, standardize = TRUE)
  hermite_est_2 <- hermite_estimator_bivar(N = 10, standardize = TRUE)
  expect_error(merge_pair(hermite_est_1,hermite_est_2))
  expect_error(merge_hermite(list(hermite_est_1,hermite_est_2)))
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
  hermite_est <- hermite_estimator(N = 10, standardize = TRUE, 
                                   observations = test_observations)
  expect_equal(target_coeff_vec_standardized,
               hermite_est$coeff_vec,
               tolerance = get_eps())
  expect_equal(mean(test_observations), 
               hermite_est$running_mean,tolerance = get_eps())
  expect_equal(sd(test_observations),sqrt(hermite_est$running_variance / 
(hermite_est$num_obs-1)),tolerance = get_eps())
  hermite_est <- hermite_estimator(N = 10, standardize = FALSE, 
                                   observations = test_observations)
  expect_equal(target_coeff_vec_unstandardized,
               hermite_est$coeff_vec,
               tolerance = get_eps())
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
               hermite_est$coeff_vec,
               tolerance = get_eps())
  hermite_est <- hermite_estimator(N = 10, standardize = TRUE)
  hermite_est <-
      hermite_est %>% update_sequential(test_observations)
  expect_equal(target_coeff_vec_standardized,
               hermite_est$coeff_vec,
               tolerance = get_eps())
  hermite_est <- hermite_estimator(N = 10, standardize = FALSE)
  for (idx in seq_along(test_observations)) {
    hermite_est <-
      hermite_est %>% update_sequential(test_observations[idx])
  }
  expect_equal(target_coeff_vec_unstandardized,
               hermite_est$coeff_vec,
               tolerance = get_eps())
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
                         hermite_est$coeff_vec,
                         tolerance = get_eps())
            
            hermite_est <-
              hermite_estimator(N = 10,
                                standardize = FALSE,
                                exp_weight_lambda = 0.05)
            for (idx in seq_along(test_observations)) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations[idx])
            }
            expect_equal(target_coeff_vec_unstandardized,
                         hermite_est$coeff_vec,
                         tolerance = get_eps())
          })

test_that("hermite_estimators merge consistently", {
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
    hermite_estimator(N = 10, standardize = FALSE, 
                      observations = test_observations)
  hermite_est_1 <-
    hermite_estimator(N = 10, standardize = FALSE, 
                      observations = test_observations[1:10])
  hermite_est_2 <-
    hermite_estimator(N = 10, standardize = FALSE, 
                      observations = test_observations[11:20])
  hermite_est_3 <-
    hermite_estimator(N = 10, standardize = FALSE, 
                      observations = test_observations[21:30])
  hermite_merged <-
    merge_hermite_univar(list(hermite_est_1, hermite_est_2, hermite_est_3))
  expect_equal(hermite_est, hermite_merged, tolerance = get_eps())
  hermite_merged <-
    merge_hermite(list(hermite_est_1, hermite_est_2, hermite_est_3))
  expect_equal(hermite_est, hermite_merged, tolerance = get_eps())
  hermite_est <-
    hermite_estimator(N = 10, standardize = TRUE, 
                      observations = test_observations)
  hermite_est_1 <-
    hermite_estimator(N = 10, standardize = TRUE, 
                      observations = test_observations[1:10])
  hermite_est_2 <-
    hermite_estimator(N = 10, standardize = TRUE, 
                      observations = test_observations[11:20])
  hermite_est_3 <-
    hermite_estimator(N = 10, standardize = TRUE, 
                      observations = test_observations[21:30])
  hermite_merged <-
    merge_hermite(list(hermite_est_1))
  expect_equal(hermite_merged,hermite_est_1)
  
  hermite_merged <-
    merge_hermite_univar(list(hermite_est_1))
  expect_equal(hermite_merged,hermite_est_1)
  
  hermite_merged <-
    merge_hermite_univar(list(hermite_est_1, hermite_est_2, hermite_est_3))
  hermite_merged_gen <-
    merge_hermite(list(hermite_est_1, hermite_est_2, hermite_est_3))
  expect_equal(hermite_merged,hermite_merged_gen)
  expect_equal(hermite_est$running_mean, hermite_merged$running_mean, 
               tolerance = get_eps())
  expect_equal(hermite_est$running_variance, hermite_merged$running_variance, 
               tolerance = get_eps())
  expect_equal(hermite_merged$num_obs,30, tolerance = get_eps())
  target_coeffs <-
    c(
      0.507692830456503,
      0.0321243433159265,
      0.0890037870794549,
      -0.0701338033329833,
      -0.0779386178176305,
      0.0851268944086588,
      -0.00850373459838627,
      -0.0885734911174612,
      0.0729495282288348,
      0.0615199678827195,
      -0.117430497748064
    )
  expect_equal(hermite_merged$coeff_vec, target_coeffs, 
               tolerance = get_eps())
  hermite_est <-
    hermite_estimator(N = 10, standardize = TRUE, 
                      observations = test_observations)
  hermite_est_1 <-
    hermite_estimator(N = 10, standardize = TRUE, 
                      observations = test_observations[1:15])
  hermite_est_2 <-
    hermite_estimator(N = 10, standardize = TRUE, 
                      observations = test_observations[16:30])
  hermite_merged <-
    merge_pair(hermite_est_1, hermite_est_2)
  expect_equal(hermite_est$running_mean, hermite_merged$running_mean, 
               tolerance = get_eps())
  expect_equal(hermite_est$running_variance, hermite_merged$running_variance, 
               tolerance = get_eps())
  expect_equal(hermite_merged$num_obs,30, tolerance = get_eps())
  target_coeffs <-
    c(
      0.507692794353622,
      0.0321248726193445,
      0.0889995081092394,
      -0.070111256131885,
      -0.0780342084589457,
      0.0856973484243477,
      -0.0101767682224489,
      -0.0848642344123028,
      0.0639666686621888,
      0.0889173735738952,
      -0.127949647659276
    )
  expect_equal(hermite_merged$coeff_vec, target_coeffs, 
               tolerance = get_eps())
  hermite_est_1 <-
    hermite_estimator(N = 10, standardize = FALSE, 
                      observations = test_observations[1:10])
  hermite_est_2 <-
    hermite_estimator(N = 20, standardize = FALSE, 
                      observations = test_observations[11:20])
  hermite_est_3 <-
    hermite_estimator(N = 10, standardize = TRUE, 
                      observations = test_observations[21:30])
  hermite_est_4 <-
    hermite_estimator(N = 10, standardize = FALSE, exp_weight_lambda = 0.01)
  hermite_est_5 <-
    hermite_estimator(N = 20, standardize = TRUE)
  expect_error(merge_hermite(list(hermite_est_1, hermite_est_2)))
  expect_error(merge_hermite(list(hermite_est_1, hermite_est_3)))
  expect_error(merge_hermite(list(hermite_est_1, hermite_est_4)))
  expect_error(merge_hermite(list(hermite_est_3, hermite_est_5)))
  expect_error(merge_hermite_univar(list(hermite_est_1, hermite_est_2)))
  expect_error(merge_hermite_univar(list(hermite_est_1, hermite_est_3)))
  expect_error(merge_hermite_univar(list(hermite_est_1, hermite_est_4)))
  expect_error(merge_hermite_univar(list(hermite_est_3, hermite_est_5)))
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
    hermite_estimator(N = 10, standardize = FALSE, 
                      observations = test_observations)
  pdf_vals <- hermite_est %>% dens(x, accelerate_series = FALSE)
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
  expect_equal(pdf_vals, target_pdf_vals_unstandardized, tolerance = get_eps())
  
  hermite_est <-
    hermite_estimator(N = 10, standardize = TRUE, 
                      observations = test_observations)
  pdf_vals <- hermite_est %>% dens(x, accelerate_series = FALSE)
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
  expect_equal(pdf_vals, target_pdf_vals_standardized, tolerance = get_eps())
  
  pdf_vals <- hermite_est %>% dens(x, clipped=T, accelerate_series = FALSE)
  expect_equal(pdf_vals, target_pdf_vals_standardized, tolerance = get_eps())
  
  hermite_est <-
    hermite_estimator(N = 10,
                      standardize = FALSE,
                      exp_weight_lambda = 0.1)
  for (idx in seq_along(test_observations)) {
    hermite_est <-
      hermite_est %>% update_sequential(test_observations[idx])
  }
  pdf_vals <- hermite_est %>% dens(x, accelerate_series = FALSE)
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
  expect_equal(pdf_vals, target_pdf_vals_unstandardized, tolerance = get_eps())
hermite_est <-
    hermite_estimator(N = 10,
                      standardize = FALSE,
                      exp_weight_lambda = 0.1)
  for (idx in seq_along(test_observations)) {
    hermite_est <-
      hermite_est %>% update_sequential(test_observations[idx])
  }
  pdf_vals <- hermite_est %>% dens(x, accelerate_series = TRUE)
  target_pdf_vals_unstandardized <-
    c(
      0.25810091,
      0.23093187,
      -0.00677207,
      0.11099507,
      0.40381941,
      0.27076818,
      0.07213663,
      0.18380665,
      0.25292388
    )
  expect_equal(pdf_vals, target_pdf_vals_unstandardized, tolerance = get_eps())

  hermite_est <-
    hermite_estimator(N = 10,
                      standardize = FALSE,
                      exp_weight_lambda = 0.1)
  for (idx in seq_along(test_observations)) {
    hermite_est <-
      hermite_est %>% update_sequential(test_observations[idx])
  }
  pdf_vals <- hermite_est %>% dens(x, clipped = TRUE, accelerate_series = FALSE)
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
               tolerance = get_eps())
  
  hermite_est <-
    hermite_estimator(N = 10,
                      standardize = TRUE,
                      exp_weight_lambda = 0.1)
  for (idx in seq_along(test_observations)) {
    hermite_est <-
      hermite_est %>% update_sequential(test_observations[idx])
  }
  pdf_vals <- hermite_est %>% dens(x, accelerate_series = FALSE)
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
  expect_equal(pdf_vals, target_pdf_vals_standardized, tolerance = get_eps())
  hermite_est <-
    hermite_estimator(N = 10)
  pdf_vals <- hermite_est %>% dens(x, accelerate_series = FALSE)
  expect_equal(length(pdf_vals), length(x))
  expect_true(all(is.na(pdf_vals)))
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
              hermite_estimator(N = 10, standardize = FALSE, 
                                observations = test_observations)
            cdf_from_pdf <- stats::integrate(
              f = function(x) {
                hermite_est %>% dens(x, accelerate_series = FALSE)
              },
              lower = -Inf,
              upper = 0.5
            )$value
            cdf_est <- hermite_est %>% cum_prob(0.5, accelerate_series = FALSE)
            expect_equal(cdf_est, 0.6549575, tolerance = get_eps())
            expect_equal(cdf_from_pdf, cdf_est, tolerance = get_eps())
            
            cdf_est <- hermite_est %>% cum_prob(0.5, clipped=TRUE, 
                                                accelerate_series = FALSE)
            expect_equal(cdf_est, 0.6549575, tolerance = get_eps())
            
            cdf_est <- hermite_est %>% cum_prob(3, clipped=TRUE,
                                                accelerate_series = FALSE)
            expect_equal(cdf_est, 1, tolerance = get_eps())
            
            hermite_est <-
              hermite_estimator(N = 10, standardize = TRUE, 
                                observations = test_observations)
            cdf_from_pdf <- stats::integrate(
              f = function(x) {
                hermite_est %>% dens(x, accelerate_series = FALSE)
              },
              lower = -Inf,
              upper = 0.5
            )$value
            cdf_est <- hermite_est %>% cum_prob(0.5, accelerate_series = FALSE)
            expect_equal(cdf_est, 0.6013645, tolerance = get_eps())
            expect_equal(cdf_from_pdf, cdf_est, tolerance = get_eps())
            
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
                hermite_est %>% dens(x, accelerate_series = FALSE)
              },
              lower = -Inf,
              upper = 0.5
            )$value
            cdf_est <- hermite_est %>% cum_prob(0.5, accelerate_series = FALSE)
            expect_equal(cdf_est, 0.6132811, tolerance = get_eps())
            expect_equal(cdf_from_pdf, cdf_est, tolerance = get_eps())
            
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
                hermite_est %>% dens(x, accelerate_series = FALSE)
              },
              lower = -Inf,
              upper = 0.5
            )$value
            cdf_est <- hermite_est %>% cum_prob(0.5, accelerate_series = FALSE)
            expect_equal(cdf_est, 0.4344541, tolerance = get_eps())
            expect_equal(cdf_from_pdf, cdf_est, tolerance = get_eps())
            
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
                hermite_est %>% dens(x, accelerate_series = TRUE)
              },
              lower = -Inf,
              upper = 0.5
            )$value
            cdf_est <- hermite_est %>% cum_prob(0.5, accelerate_series = TRUE)
            expect_equal(cdf_est, 0.4955906, tolerance = get_eps())
            expect_equal(cdf_from_pdf, cdf_est, tolerance = get_eps())
            hermite_est <-
              hermite_estimator(N = 10)
            expect_equal(hermite_est %>% cum_prob(0.5, 
                                                accelerate_series = FALSE),NA)
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
    hermite_estimator(N = 0, standardize = TRUE, 
                      observations = test_observations)
  quantiles_est <- hermite_est %>% quant(c(0.25, 0.5, 0.75))
  expect_equal(quantiles_est,
               c(-1.0771798, -0.1580109,  0.7745023),
               tolerance = get_eps())
  quantiles_est <- hermite_est %>% quant(c(0.25, 0.5, 0.75), 
                                         algorithm = "bisection")
  expect_equal(quantiles_est,
               c(-1.0771734, -0.1513388,  0.7744958),
               tolerance = get_eps())
  quantiles_est <- hermite_est %>% quant(c(0.25, 0.5, 0.75), 
                                         accelerate_series = FALSE)
  expect_equal(quantiles_est,
               c(-1.0771798, -0.1580109,  0.7745023),
               tolerance = get_eps())
  quantiles_est <- hermite_est %>% quant(c(0.25, 0.5, 0.75), 
                                         algorithm = "bisection",
                                         accelerate_series = FALSE)
  expect_equal(quantiles_est,
               c(-1.0771734, -0.1513388,  0.7744958),
               tolerance = get_eps())
  hermite_est <-
    hermite_estimator(N = 2, standardize = TRUE, 
                      observations = test_observations)
  quantiles_est <- hermite_est %>% quant(c(0.25, 0.5, 0.75))
  expect_equal(quantiles_est,
               c(-1.1394204, 0.0169664,  1.0930044),
               tolerance = get_eps())
  quantiles_est <- hermite_est %>% quant(c(0.25, 0.5, 0.75), 
                                         algorithm = "bisection")
  expect_equal(quantiles_est,
               c(-1.13935192, 0.01706313,  1.09273375),
               tolerance = get_eps())
  quantiles_est <- hermite_est %>% quant(c(0.25, 0.5, 0.75), 
                                         accelerate_series = FALSE)
  expect_equal(quantiles_est,
               c(-1.353599,  0.154451,  1.290673),
               tolerance = get_eps())
  quantiles_est <- hermite_est %>% quant(c(0.25, 0.5, 0.75), 
                                         algorithm = "bisection",
                                         accelerate_series = FALSE)
  expect_equal(quantiles_est,
               c(-1.3533800,  0.1420511,  1.2907235),
               tolerance = get_eps())
  hermite_est <-
    hermite_estimator(N = 10, standardize = TRUE, 
                      observations = test_observations)
  quantiles_est <- hermite_est %>% quant(c(0.25, 0.5, 0.75))
  expect_equal(quantiles_est,
               c(-1.55096352,  0.05649595,  0.93885643),
               tolerance = get_eps())
  quantiles_est <- hermite_est %>% quant(c(0.25, 0.5, 0.75), 
                                         algorithm = "bisection")
  expect_equal(quantiles_est,
               c(-1.55081667,  0.05632926,  0.93898751),
               tolerance = get_eps())
  quantiles_est <- hermite_est %>% quant(c(0.25, 0.5, 0.75), 
                                         accelerate_series = FALSE)
  expect_equal(quantiles_est,
               c(-1.31098344,  0.04145163,  0.90887905),
               tolerance = get_eps())
  quantiles_est <- hermite_est %>% quant(c(0.25, 0.5, 0.75), 
                                         algorithm = "bisection",
                                         accelerate_series = FALSE)
  expect_equal(quantiles_est,
               c(-1.31079556,  0.04139707,  0.90912314),
               tolerance = get_eps())
  hermite_est <- hermite_estimator(N = 10, standardize = TRUE)
  for (idx in seq_along(test_observations)) {
    hermite_est <-
      hermite_est %>% update_sequential(test_observations[idx])
  }
  quantiles_est <- hermite_est %>% quant(c(0.25, 0.5, 0.75))
  expect_equal(quantiles_est,
               c(-0.8253848,  0.1408785,  1.1352633),
               tolerance = get_eps())
  hermite_est <-
    hermite_estimator(N = 10,
                      standardize = TRUE,
                      exp_weight_lambda = 0.1)
  for (idx in seq_along(test_observations)) {
    hermite_est <-
      hermite_est %>% update_sequential(test_observations[idx])
  }
  quantiles_est <- hermite_est %>% quant(c(0.25, 0.5, 0.75))
  expect_equal(quantiles_est, c(-0.8790995,  0.4444830,  1.4835615), 
               tolerance = get_eps())
  quantiles_est <- hermite_est %>% quant(c(0.25, 0.5, 0.75), 
                                         algorithm = "bisection")
  expect_equal(quantiles_est, c(-0.8793468 , 0.4445265,  1.4835805),
               tolerance = get_eps())
  hermite_est <-
    hermite_estimator(N = 20,
                      standardize = TRUE, observations = c(1:4))
  quantiles_est <- hermite_est %>% quant(c(0.5))
  expect_equal(quantiles_est,2.496089, 
               tolerance = get_eps())
  quantiles_est <- hermite_est %>% quant(c(0.5), algorithm="bisection")
  expect_equal(quantiles_est,2.5, 
               tolerance = get_eps())
  hermite_est <- hermite_estimator(N = 10, standardize = TRUE, 
                                   observations = c(1,1))
  expect_equal(quant(hermite_est, p=0.5), 1)
})

test_that("convenience and utility functions work as expected", {
  hermite_poly_vals <- as.vector(hermite_polynomial_N(N=0,x=c(2)))
  target_values <-
    c(
      1
    )
  expect_equal(hermite_poly_vals,target_values,tolerance=get_eps())
  hermite_poly_vals <- as.vector(hermite_polynomial_N(N=1,x=c(2)))
  target_values <-
    c(
      1,
      4
    )
  expect_equal(hermite_poly_vals,target_values,tolerance=get_eps())
  hermite_poly_vals <- as.vector(hermite_polynomial_N(N=6,x=c(2)))
  target_values <-
    c(
      1,
      4,
      14,
      40,
      76,
      -16,
      -824
    )
  expect_equal(hermite_poly_vals,target_values,tolerance=get_eps())
  hermite_function_vals <- as.vector(hermite_function_N(N=0,x=c(2)))
  target_values <-
    c(
      0.101653788306418
    )
  expect_equal(hermite_function_vals,target_values,tolerance=get_eps())
  hermite_function_vals <- as.vector(hermite_function_N(N=1,x=c(2)))
  target_values <-
    c(
      0.101653788306418,
      0.287520332179079
    )
  expect_equal(hermite_function_vals,target_values,tolerance=get_eps())
  hermite_function_vals <- as.vector(hermite_function_N(N=6,x=c(2)))
  target_values <-
    c(
      0.101653788306418,
      0.287520332179079,
      0.503160581313389,
      0.586898420428556,
      0.39424986030507,
      -0.0262468952793101,
      -0.390206540413716
    )
  expect_equal(hermite_function_vals,target_values,tolerance=get_eps())
  target_integral <- stats::integrate(f=function(t)
    {hermite_function_N(N=6,t)[7,]}, lower=-Inf, upper=5)$value
  hermite_int_lower_val <- hermite_int_lower(N=6,x=5)[7]
  expect_equal(target_integral,hermite_int_lower_val,tolerance=get_eps())
  target_integral <- stats::integrate(f=function(t){
    hermite_function_N(N=6,t)[7,]}, lower=2, upper=Inf)$value
  hermite_int_upper_val <- hermite_int_upper(N=6,x=2)[7]
  expect_equal(target_integral,hermite_int_upper_val,tolerance=get_eps())
  target_integral <- stats::integrate(f=function(t){
    hermite_function_N(N=6,t)[7,]}, lower=-Inf, upper=Inf)$value
  hermite_int_full <- hermite_int_full(N=6)[7]
  expect_equal(target_integral,hermite_int_full,tolerance=get_eps())
  target_integral <- stats::integrate(function(x){x*exp(-x^2)}, 
                                      lower=-Inf,upper=Inf)$value
  quad_val <- gauss_hermite_quad_100(function(x){x})
  expect_equal(target_integral,quad_val,tolerance=get_eps())
  hermite_function_sum_vals <- as.vector(
    rowSums(hermite_function_N(N=6,x=c(1))))
  expect_equal(hermite_function_sum_vals,
                         hermite_function_sum_N(N=6,x=c(1)),tol=get_eps())
  hermite_function_sum_vals <- as.vector(
    rowSums(hermite_function_N(N=6,x=c(1,2))))
  expect_equal(hermite_function_sum_vals,
               hermite_function_sum_N(N=6,x=c(1,2)),tol=get_eps())
  hermite_function_sum_vals <- as.vector(
    rowSums(hermite_function_N(N=6,x=c(1,2,3))))
  expect_equal(hermite_function_sum_vals,
               hermite_function_sum_N(N=6,x=c(1,2,3)),tol=get_eps())
  hermite_function_sum_vals <- as.vector(
    rowSums(hermite_function_N(N=6,x=c(1,2,3,4))))
  expect_equal(hermite_function_sum_vals,
               hermite_function_sum_N(N=6,x=c(1,2,3,4)),tol=get_eps())
  hermite_function_sum_vals <- as.vector(
    rowSums(hermite_function_N(N=6,x=seq(-4,4,by=0.5))))
  expect_equal(hermite_function_sum_vals,
               hermite_function_sum_N(N=6,x=seq(-4,4,by=0.5)),tol=get_eps())
  h_input_mat <- matrix(data=rep(1,1*20),nrow=1,ncol=20)
  output_length <- length(series_calculate(h_input_mat,1, 
                                           accelerate_series = F))
  expect_equal(output_length,20)
  output_length <- length(series_calculate(h_input_mat,1, 
                                           accelerate_series = T))
  expect_equal(output_length,20)
  h_input_mat <- matrix(data=rep(1,3*20),nrow=3,ncol=20)
  output_length <- length(series_calculate(h_input_mat,1:3, 
                                           accelerate_series = F))
  expect_equal(output_length,20)
  output_length <- length(series_calculate(h_input_mat,1:3, 
                                           accelerate_series = T))
  expect_equal(output_length,20)
  h_input_mat <- matrix(data=rep(1,7*20),nrow=7,ncol=20)
  output_length <- length(series_calculate(h_input_mat,1:7, 
                                           accelerate_series = F))
  expect_equal(output_length,20)
  output_length <- length(series_calculate(h_input_mat,1:7, 
                                           accelerate_series = T))
  expect_equal(output_length,20)
  h_input_mat <- matrix(data=rep(1,13*20),nrow=13,ncol=20)
  output_length <- length(series_calculate(h_input_mat,1:13, 
                                           accelerate_series = F))
  expect_equal(output_length,20)
  output_length <- length(series_calculate(h_input_mat,1:13,
                                          accelerate_series = T))
  expect_equal(output_length,20)
})

test_that("Print and Summary work as expected", {
  expect_equal(capture.output(print(hermite_estimator())), 
               c("Univariate Hermite Estimator:",
                 "N = 50", "Standardize observations = TRUE",
                 "Exponential weighting for coefficents = FALSE",
                 "Number of observations = 0"))
  h_est <- hermite_estimator(observations = c(1,2,3))
  expect_equal(capture.output(summary(h_est)), 
               c('Univariate Hermite Estimator:','N = 50',
                 'Standardize observations = TRUE',
                 'Exponential weighting for coefficents = FALSE',
                 'Number of observations = 3','','Mean = 2',
                 'Standard Deviation = 1','Estimated Quantiles:',
      paste0('    10%    20%    30%    40%    50%    60%   70%    ',
                 '80%    90%'), 
      paste0(' 0.9375 1.0375 1.1449 1.8919 2.0007 2.1081 2.295 2.9625 3.0625')))
})
