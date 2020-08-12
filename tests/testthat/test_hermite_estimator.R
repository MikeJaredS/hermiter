context("Hermite Estimator")
library(hermiter)
library(magrittr)

test_that("hermite_estimator constructor returns correct class", {
  hermite_est <- hermite_estimator(N = 10, normalize = TRUE)
  expect_is(hermite_est, "hermite_estimator")
})

test_that("batch updates of hermite_estimator work as expected",
          {
            test_observations <-
              c(
                0.3719336,
                0.8484731,-2.663014,-2.028416,-2.429764,-1.677079,
                1.434098,-1.453405,-1.133242,
                0.963844,
                1.46609,
                0.371181,
                2.135272,-0.7328963,
                0.8465673,-2.168731,-0.269106,-1.885169,-0.07963116,
                0.1244462,
                0.1165929,
                1.460038,-0.06261991,
                0.07363522,
                0.03415375,-1.65761,
                2.058115,
                0.9413341,-1.759675,
                2.214421
              )
            target_coeff_vec_normalized <-
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
            
            target_coeff_vec_unnormalized <-
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
            hermite_est <- hermite_estimator(N = 10, normalize = T)
            hermite_est <-
              hermite_est %>% update_batch(test_observations)
            expect_equal(target_coeff_vec_normalized, hermite_est$coeff_vec, tolerance =
                           1e-7)
            
            hermite_est <- hermite_estimator(N = 10, normalize = F)
            hermite_est <-
              hermite_est %>% update_batch(test_observations)
            expect_equal(target_coeff_vec_unnormalized,
                         hermite_est$coeff_vec,
                         tolerance = 1e-7)
          })

test_that("sequential updates of hermite_estimator work as expected",
          {
            test_observations <-
              c(
                0.3719336,
                0.8484731,-2.663014,-2.028416,-2.429764,-1.677079,
                1.434098,-1.453405,-1.133242,
                0.963844,
                1.46609,
                0.371181,
                2.135272,-0.7328963,
                0.8465673,-2.168731,-0.269106,-1.885169,-0.07963116,
                0.1244462,
                0.1165929,
                1.460038,-0.06261991,
                0.07363522,
                0.03415375,-1.65761,
                2.058115,
                0.9413341,-1.759675,
                2.214421
              )
            target_coeff_vec_normalized <-
              c(
                0.5492528,
                0.08216913,
                0.004834347,
                0.01674095,-0.04670184,-0.01720743,-0.003257658,-0.02520607,
                0.04960573,
                0.06876952,-0.07160235
              )
            target_coeff_vec_unnormalized <-
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
            hermite_est <- hermite_estimator(N = 10, normalize = T)
            for (idx in c(1:length(test_observations))) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations[idx])
            }
            expect_equal(target_coeff_vec_normalized, hermite_est$coeff_vec, tolerance =
                           1e-7)
            
            hermite_est <- hermite_estimator(N = 10, normalize = F)
            for (idx in c(1:length(test_observations))) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations[idx])
            }
            expect_equal(target_coeff_vec_unnormalized,
                         hermite_est$coeff_vec,
                         tolerance = 1e-7)
          })

test_that("sequential updates of exponentially weighted hermite_estimator work as expected",
          {
            test_observations <-
              c(
                0.3719336,
                0.8484731,-2.663014,-2.028416,-2.429764,-1.677079,
                1.434098,-1.453405,-1.133242,
                0.963844,
                1.46609,
                0.371181,
                2.135272,-0.7328963,
                0.8465673,-2.168731,-0.269106,-1.885169,-0.07963116,
                0.1244462,
                0.1165929,
                1.460038,-0.06261991,
                0.07363522,
                0.03415375,-1.65761,
                2.058115,
                0.9413341,-1.759675,
                2.214421
              )
            target_coeff_vec_normalized <-
              c(
                0.3720845, 0.0329111, 0.07880838, -0.0480489, -0.001239091, 0.03989803, -0.08362694, -0.04300433, 0.1199047, 0.04287386, -0.1147223
              )
            target_coeff_vec_unnormalized <-
              c(0.3140483, 0.03968127, 0.1031914, -0.0542355, 0.1318214, 0.03501422, -0.1379236, -0.005167008, 0.06536224, -0.02967193, -0.02234358
              )
            hermite_est <- hermite_estimator(N = 10, normalize = TRUE, exp_weight_lambda = 0.05)
            for (idx in c(1:length(test_observations))) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations[idx])
            }
            expect_equal(target_coeff_vec_normalized, hermite_est$coeff_vec, tolerance =
                           1e-7)
            
            hermite_est <- hermite_estimator(N = 10, normalize = FALSE, exp_weight_lambda = 0.05)
            for (idx in c(1:length(test_observations))) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations[idx])
            }
            expect_equal(target_coeff_vec_unnormalized,
                         hermite_est$coeff_vec,
                         tolerance = 1e-7)
          })

test_that("hermite_estimators combine consistently", {
})

test_that("probability density estimation works as expected", {
})

test_that("cumulative distribution function estimation works as expected", {
})

test_that("quantile estimation works as expected", {
})