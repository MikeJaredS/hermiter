## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

context("utils_general")
library(hermiter)
library(magrittr)

get_eps <- function() {
  return(1e-3)
}

test_that("density generics work correctly", {
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
  # Standardize = FALSE
  hermite_est <- hermite_estimator(observations = test_observations,
                                   standardize = FALSE)
  expect_error(density(hermite_est))
  dens_output <-
    density(
      hermite_est,
      x_lower = -3,
      x_upper = 3,
      clipped = T
    )
  expect_equal(
    capture.output(print(dens_output, digits = 5)),
    c(
      " x[1:100] =     -3, -2.9394, -2.8788,  ..., 2.9394,      3",
      " d[1:100] =  1e-08,  1e-08, 0.012479,  ..., 0.061606, 0.062612"
    )
  )
  # Standardize = TRUE
  hermite_est <- hermite_estimator(observations = test_observations,
                                   standardize = TRUE)
  dens_output <- density(hermite_est)
  expect_equal(
    capture.output(print(dens_output, digits = 5)),
    c(
      " x[1:99] = -2.795, -2.669, -2.5501,  ..., 2.2876, 2.3505",
      " d[1:99] = 0.076129, 0.081658, 0.087815,  ..., 0.17415, 0.14304"
    )
  )
  test_observations <- c(
    -0.37826482129403,
    -1.47945641842633,
    0.586716971868732,
    -0.665669486496455,
    -1.43166984568682,
    0.474563723449495,
    -1.11341217646342,
    1.36541799968729,
    -0.24770649613921,
    1.4481757755572,
    -0.529665159959186,
    0.686121642323148,
    1.01383912025988,
    -0.780484609763183,
    -0.545088136349466,
    0.135846098413723,
    -0.240679160926683,
    1.77147252004363,
    -1.45865271681631,
    0.0628601340760367,
    1.07810464037276,
    -0.17475267390839,
    -0.66888270234312,
    -0.408256095712664,
    0.172717270523465,
    -0.741493500626988,
    1.01132229766528,
    -0.959035155862129,
    -0.482739678200718,
    -0.753321211065154,
    0.503938515418142,
    -1.78785564896096,
    2.47357639665998,
    -0.489738640544134,
    -0.714875932607426,
    0.806157676535886,
    -1.00656011023483,
    -0.984503617692169,
    1.30774013514267,
    0.440505965420727,
    -0.650310967710673,
    -1.66222913387161,
    1.25306046766581,
    0.0171057800672902,
    -0.563511566403471,
    0.388015625901842,
    0.66092470022605,
    1.65884426783205,
    -0.123975093954792,
    -0.552324416383275,
    1.18682631574925,
    0.435917095119776,
    -0.732475285222769,
    0.0837044467479744,
    0.0011521929124057,
    0.0224946049862673,
    1.64913440622687,
    -2.35583419045975,
    -0.350200566468806,
    0.578709836500825
  )
  test_observations_mat <-
    matrix(test_observations,
           nrow = 30,
           ncol = 2,
           byrow = F)
  # Standardize = FALSE
  hermite_est <-
    hermite_estimator(observations = test_observations_mat,
                      standardize = FALSE,
                      est_type = "bivariate")
  expect_error(density(hermite_est))
  dens_output <-
    density(
      hermite_est,
      x_lower = c(-3, -3),
      x_upper = c(3, 3),
      clipped = T
    )
  expect_equal(
    capture.output(print(dens_output, digits = 5)),
    c(
      "  x1      x2          d",
      "1 -3 -3.0000 0.00235936",
      "2 -3 -2.6842 0.00000001",
      "3 -3 -2.3684 0.00467963",
      "...",
      "                       ",
      "399 3 2.6842 0.00000001",
      "400 3 3.0000 0.00185688"
    )
  )
  # Standardize = TRUE
  hermite_est <-
    hermite_estimator(observations = test_observations_mat,
                      standardize = TRUE,
                      est_type = "bivariate")
  dens_output <- density(hermite_est)
  expect_equal(
    capture.output(print(dens_output, digits = 5)),
    c(
      "       x1       x2         d",
      "1 -1.4578 -1.87119  0.219245",
      "2 -1.4578 -1.26122 -0.065646",
      "3 -1.4578 -0.89678  0.178301",
      "...",
      "                            ",
      "360 1.4854 1.4755 -0.0042704",
      "361 1.4854 1.7168  0.0294621"
    )
  )
  
})

test_that("density generics work correctly", {
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
  # Standardize = FALSE
  hermite_est <- hermite_estimator(observations = test_observations,
                                   standardize = FALSE)
  expect_error(hcdf(hermite_est))
  cdf_output <-
    hcdf(
      hermite_est,
      x_lower = -3,
      x_upper = 3,
      clipped = T
    )
  expect_equal(
    capture.output(print(cdf_output, digits = 5)),
    c(
      " x[1:100] =     -3, -2.9394, -2.8788,  ..., 2.9394,      3",
      " p[1:100] =  1e-08,  1e-08,  1e-08,  ..., 0.99601, 0.99988"
    )
  )
  expect_equal(
    capture.output(summary(cdf_output, digits = 5)),
    c(
      "Hermite CDF estimates:",
      "30 observations",
      "100 evaluation points",
      " x[1:100] =     -3, -2.9394, -2.8788,  ..., 2.9394,      3",
      " p[1:100] =  1e-08,  1e-08,  1e-08,  ..., 0.99601, 0.99988"
    )
  )
  # Standardize = TRUE
  hermite_est <- hermite_estimator(observations = test_observations,
                                   standardize = TRUE)
  cdf_output <- hcdf(hermite_est)
  expect_equal(
    capture.output(print(cdf_output, digits = 5)),
    c(
      " x[1:99] = -2.795, -2.669, -2.5501,  ..., 2.2876, 2.3505",
      " p[1:99] = 0.0099982, 0.019999, 0.029998,  ..., 0.98027, 0.99028"
    )
  )
  expect_equal(
    capture.output(summary(cdf_output, digits = 5)),
    c(
      "Hermite CDF estimates:",
      "30 observations",
      "99 evaluation points",
      " x[1:99] = -2.795, -2.669, -2.5501,  ..., 2.2876, 2.3505",
      " p[1:99] = 0.0099982, 0.019999, 0.029998,  ..., 0.98027, 0.99028"
    )
  )
  test_observations <- c(
    -0.37826482129403,
    -1.47945641842633,
    0.586716971868732,
    -0.665669486496455,
    -1.43166984568682,
    0.474563723449495,
    -1.11341217646342,
    1.36541799968729,
    -0.24770649613921,
    1.4481757755572,
    -0.529665159959186,
    0.686121642323148,
    1.01383912025988,
    -0.780484609763183,
    -0.545088136349466,
    0.135846098413723,
    -0.240679160926683,
    1.77147252004363,
    -1.45865271681631,
    0.0628601340760367,
    1.07810464037276,
    -0.17475267390839,
    -0.66888270234312,
    -0.408256095712664,
    0.172717270523465,
    -0.741493500626988,
    1.01132229766528,
    -0.959035155862129,
    -0.482739678200718,
    -0.753321211065154,
    0.503938515418142,
    -1.78785564896096,
    2.47357639665998,
    -0.489738640544134,
    -0.714875932607426,
    0.806157676535886,
    -1.00656011023483,
    -0.984503617692169,
    1.30774013514267,
    0.440505965420727,
    -0.650310967710673,
    -1.66222913387161,
    1.25306046766581,
    0.0171057800672902,
    -0.563511566403471,
    0.388015625901842,
    0.66092470022605,
    1.65884426783205,
    -0.123975093954792,
    -0.552324416383275,
    1.18682631574925,
    0.435917095119776,
    -0.732475285222769,
    0.0837044467479744,
    0.0011521929124057,
    0.0224946049862673,
    1.64913440622687,
    -2.35583419045975,
    -0.350200566468806,
    0.578709836500825
  )
  test_observations_mat <-
    matrix(test_observations,
           nrow = 30,
           ncol = 2,
           byrow = F)
  # Standardize = FALSE
  hermite_est <-
    hermite_estimator(observations = test_observations_mat,
                      standardize = FALSE,
                      est_type = "bivariate")
  expect_error(hcdf(hermite_est))
  cdf_output <-
    hcdf(
      hermite_est,
      x_lower = c(-3, -3),
      x_upper = c(3, 3),
      clipped = T
    )
  expect_equal(
    capture.output(print(cdf_output, digits = 5)),
    c(
      "  x1      x2          p",
      "1 -3 -3.0000 0.00071279",
      "2 -3 -2.6842 0.00092093",
      "3 -3 -2.3684 0.00081187",
      "...",
      "                    ",
      "399 3 2.6842 0.98625",
      "400 3 3.0000 0.98795"
    )
  )
  expect_equal(
    capture.output(summary(cdf_output, digits = 5)),
    c(
      "Hermite CDF estimates:",
      "30 observations",
      "400 evaluation points",
      "  x1      x2          p",
      "1 -3 -3.0000 0.00071279",
      "2 -3 -2.6842 0.00092093",
      "3 -3 -2.3684 0.00081187",
      "...",
      "                    ",
      "399 3 2.6842 0.98625",
      "400 3 3.0000 0.98795"
    )
  )
  # Standardize = TRUE
  hermite_est <-
    hermite_estimator(observations = test_observations_mat,
                      standardize = TRUE,
                      est_type = "bivariate")
  cdf_output <- hcdf(hermite_est)
  expect_equal(
    capture.output(print(cdf_output, digits = 5)),
    c(
      "       x1       x2         p",
      "1 -1.4578 -1.87119 0.0040857",
      "2 -1.4578 -1.26122 0.0137779",
      "3 -1.4578 -0.89678 0.0133331",
      "...",
      "                         ",
      "360 1.4854 1.4755 0.86669",
      "361 1.4854 1.7168 0.90147"
    )
  )
  expect_equal(
    capture.output(summary(cdf_output, digits = 5)),
    c(
      "Hermite CDF estimates:",
      "30 observations",
      "361 evaluation points",
      "       x1       x2         p",
      "1 -1.4578 -1.87119 0.0040857",
      "2 -1.4578 -1.26122 0.0137779",
      "3 -1.4578 -0.89678 0.0133331",
      "...",
      "                         ",
      "360 1.4854 1.4755 0.86669",
      "361 1.4854 1.7168 0.90147"
    )
  )
  
})

test_that("quantile, median and IQR generics work correctly", {
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
  # Standardize = FALSE
  hermite_est <- hermite_estimator(observations = test_observations,
                                   standardize = FALSE)
  expect_error(quantile(hermite_est))
  # Standardize = TRUE
  hermite_est <- hermite_estimator(observations = test_observations,
                                   standardize = TRUE)
  quantiles_est <- quantile(hermite_est)
  target_quantiles <- c(-14.64905596,
                       -1.57302308,
                       0.04009478,
                       1.01244707,
                       2.43500549)
  expect_equal(quantiles_est, target_quantiles, get_eps())
  quantiles_est_original <- quant(hermite_est, p = seq(0, 1, 0.25))
  expect_equal(quantiles_est, quantiles_est_original, get_eps())
  
  expect_equal(median(hermite_est), 0.04009478, get_eps())
  expect_equal(median(hermite_est), quant(hermite_est, p = 0.5), get_eps())
  
  expect_equal(IQR(hermite_est), 2.58547, get_eps())
  quartiles <- quant(hermite_est, p = c(0.25, 0.75))
  expect_equal(IQR(hermite_est), quartiles[2] - quartiles[1] , get_eps())
  
  expect_equal(IQR(test_observations), 2.524678, get_eps())
  
})


test_that("cor wrapper works correctly", {
  test_observations <- c(
    -0.37826482129403,
    -1.47945641842633,
    0.586716971868732,
    -0.665669486496455,
    -1.43166984568682,
    0.474563723449495,
    -1.11341217646342,
    1.36541799968729,
    -0.24770649613921,
    1.4481757755572,
    -0.529665159959186,
    0.686121642323148,
    1.01383912025988,
    -0.780484609763183,
    -0.545088136349466,
    0.135846098413723,
    -0.240679160926683,
    1.77147252004363,
    -1.45865271681631,
    0.0628601340760367,
    1.07810464037276,
    -0.17475267390839,
    -0.66888270234312,
    -0.408256095712664,
    0.172717270523465,
    -0.741493500626988,
    1.01132229766528,
    -0.959035155862129,
    -0.482739678200718,
    -0.753321211065154,
    0.503938515418142,
    -1.78785564896096,
    2.47357639665998,
    -0.489738640544134,
    -0.714875932607426,
    0.806157676535886,
    -1.00656011023483,
    -0.984503617692169,
    1.30774013514267,
    0.440505965420727,
    -0.650310967710673,
    -1.66222913387161,
    1.25306046766581,
    0.0171057800672902,
    -0.563511566403471,
    0.388015625901842,
    0.66092470022605,
    1.65884426783205,
    -0.123975093954792,
    -0.552324416383275,
    1.18682631574925,
    0.435917095119776,
    -0.732475285222769,
    0.0837044467479744,
    0.0011521929124057,
    0.0224946049862673,
    1.64913440622687,
    -2.35583419045975,
    -0.350200566468806,
    0.578709836500825
  )
  test_observations_mat <-
    matrix(test_observations,
           nrow = 30,
           ncol = 2,
           byrow = F)
  cor_matrix <-
    cor(test_observations_mat, method = "hermite.spearman")
  target_matrix <-
    matrix(c(1.0000000, 0.4991461, 0.4991461, 1.0000000),
           nrow = 2,
           byrow = T)
  expect_equal(cor_matrix, target_matrix, get_eps())
  
  cor_matrix <-
    cor(test_observations_mat, method = "hermite.kendall")
  target_matrix <-
    matrix(c(1.0000000, 0.3478115, 0.3478115, 1.0000000),
           nrow = 2,
           byrow = T)
  expect_equal(cor_matrix, target_matrix, get_eps())
  
  expect_equal(
    stats::cor(test_observations_mat),
    cor(test_observations_mat, method = "pearson")
  )
  expect_equal(
    stats::cor(test_observations_mat, method = "kendall"),
    cor(test_observations_mat, method = "kendall")
  )
  expect_equal(
    stats::cor(test_observations_mat, method = "spearman"),
    cor(test_observations_mat, method = "spearman")
  )
  
  cor_value <- cor(x = 1:100,
                   y = 201:300,
                   method = "hermite.spearman")
  expect_equal(cor_value, 0.9825907, get_eps())
  
  cor_value <- cor(x = 1:100,
                   y = 201:300,
                   method = "hermite.kendall")
  expect_equal(cor_value, 0.9223132, get_eps())
  
  cor_value <- cor(x = 1:100,
                   y = 201:300,
                   method = "pearson")
  expect_equal(cor_value, 1, get_eps())
  
  cor_value <- cor(x = 1:100,
                   y = 201:300,
                   method = "kendall")
  expect_equal(cor_value, 1, get_eps())
  
  cor_value <- cor(x = 1:100,
                   y = 201:300,
                   method = "spearman")
  expect_equal(cor_value, 1, get_eps())
  
  test_observations_mat_1 <-
    matrix(1:90,
           ncol = 3)
  test_observations_mat_2 <-
    matrix(91:180,
           ncol = 3)
  cor_matrix <-
    cor(x = test_observations_mat_1, y = test_observations_mat_2, 
        method = "hermite.spearman")
  target_matrix <- matrix(rep(0.9758105, 9), ncol = 3)
  expect_equal(cor_matrix, target_matrix, get_eps())
  
  cor_matrix <-
    cor(x = test_observations_mat_1, y = test_observations_mat_2, 
        method = "hermite.kendall")
  target_matrix <- matrix(rep(0.9189, 9), ncol = 3)
  expect_equal(cor_matrix, target_matrix, get_eps())
  
  cor_matrix <-
    cor(x = test_observations_mat_1, y = test_observations_mat_2, 
        method = "pearson")
  target_matrix <- matrix(rep(1, 9), ncol = 3)
  expect_equal(cor_matrix, target_matrix, get_eps())
  
  cor_matrix <-
    cor(x = test_observations_mat_1, y = test_observations_mat_2, 
        method = "kendall")
  target_matrix <- matrix(rep(1, 9), ncol = 3)
  expect_equal(cor_matrix, target_matrix, get_eps())
  
  cor_matrix <-
    cor(x = test_observations_mat_1, y = test_observations_mat_2, 
        method = "spearman")
  target_matrix <- matrix(rep(1, 9), ncol = 3)
  expect_equal(cor_matrix, target_matrix, get_eps())
  
  expect_error(cor(x = test_observations_mat_1, method = "unknown"))
  expect_error(cor(y = test_observations_mat_1, method = "hermite.spearman"))
  expect_error(cor(x = c(1, 2, 3, 4), method = "hermite.spearman"))
})
