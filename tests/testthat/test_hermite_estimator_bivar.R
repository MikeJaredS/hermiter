## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

context("hermite_estimator_bivar")
library(hermiter)
library(magrittr)

get_eps <- function(){
  return(1e-3)
}

test_that("hermite_estimator_bivar constructor returns correct class", {
  hermite_est <-
    hermite_estimator(N = 10,
                      standardize = TRUE,
                      est_type = "bivariate")
  expect_is(hermite_est, "hermite_estimator_bivar")
})

test_that("error trapping for hermite_estimator_bivar work as expected", {
  expect_error(hermite_estimator(N = "a", standardize = TRUE, 
                                 est_type="bivariate"))
  expect_error(hermite_estimator(N = 80, standardize = TRUE, 
                                 est_type="bivariate"))
  expect_error(hermite_estimator(N = -1, standardize = TRUE, 
                                 est_type="bivariate"))
  expect_error(hermite_estimator(N = 10, standardize = 8, 
                                 est_type="bivariate"))
  expect_error(hermite_estimator(N = 10, exp_weight_lambda = -0.5, 
                                 est_type="bivariate"))
  expect_error(hermite_estimator(N = 10, exp_weight_lambda = 1.5, 
                                 est_type="bivariate"))
  expect_error(hermite_estimator(N = 10, exp_weight_lambda = "a", 
                                 est_type="bivariate"))
  hermite_est <- hermite_estimator(N = 10, standardize = TRUE, 
                                   est_type="bivariate")
  expect_error(update_sequential(hermite_est, "a"))
  expect_error(update_sequential(hermite_est, c(1,2,3)))
  expect_error(hermite_estimator(N = 10, standardize = TRUE, 
                                       observations = c("a","b"), 
                                 est_type="bivariate"))
  expect_error(hermite_estimator(N = 10, standardize = TRUE, 
                                       observations = c(1,2,3), 
                                 est_type="bivariate"))
  expect_error(hermite_estimator(N = 10, standardize = TRUE, 
                                       observations = matrix(c(1,2,3), 
                                        nrow=1,ncol=3), est_type="bivariate"))
  expect_error(dens(hermite_est,x=c()))
  expect_error(dens(hermite_est,x=c("a","b")))
  expect_error(dens(hermite_est,x=c(1,2,3)))
  expect_error(dens(hermite_est,x=matrix(c(1,2,3), nrow=1,ncol=3)))
  expect_error(cum_prob(hermite_est,x=c()))
  expect_error(cum_prob(hermite_est,x=c("a","b")))
  expect_error(cum_prob(hermite_est,x=c(1,2,3)))
  expect_error(cum_prob(hermite_est,x=matrix(c(1,2,3), nrow=1,ncol=3)))
  expect_error(quant(hermite_est))
  expect_error(merge_hermite_bivar(list()))
  hermite_est_1 <- hermite_estimator(N = 10, standardize = TRUE, 
                                     est_type="bivariate")
  hermite_est_2 <- hermite_estimator_univar(N = 10, standardize = TRUE)
  expect_error(merge_pair(hermite_est_1,hermite_est_2))
})

test_that("batch updates of hermite_estimator_bivar work as expected", {
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
  target_coeff_vec_standardized_x <-
    c(
      0.51961614764726,
      -0.0840330894386698,
      0.0437816793075618,
      0.131708361054797,
      -0.0333025722055132,
      -0.0930519251739091,
      -0.00707318784009528,
      0.058082502006704,
      0.013178999363968,
      -0.04474864378677,
      0.00282139633198586
    )
  target_coeff_vec_standardized_y <-
    c(
      0.536646509839032,
      -0.00251088825721887,
      -0.0172921504640147,
      0.0145284939379268,
      0.0220432422329034,
      -0.0293833771230316,
      -0.00873537430706128,
      0.0321419554321903,
      -0.0174613292270531,
      -0.0191365163924595,
      0.0293993051247836
    )
  
  target_coeff_vec_unstandardized_x <-
    c(
      0.542468993917852,
      -0.112156489591074,
      0.0240427460518777,
      0.114702892134999,
      -0.0745566716013653,
      -0.0915076243058221,
      0.0423861774862675,
      0.0682739544693091,
      -0.0171370670453854,
      -0.0519546168194523,
      0.0128991995695314
    )
  
  target_coeff_vec_unstandardized_y <-
    c(
      0.514269409817415,
      0.0134314987076482,
      0.0103702696517751,
      0.0184650339583687,
      0.0251952541085223,
      -0.0312867094537839,
      -0.00130336810122759,
      0.0288137305459172,
      -0.0266496927060038,
      -0.00612587622958993,
      0.0340885404266732
    )
  hermite_est <-
    hermite_estimator(N = 10,
                      standardize = TRUE,
                      est_type = "bivariate", 
                      observations = test_observations_mat)
  expect_equal(target_coeff_vec_standardized_x,
               hermite_est$coeff_vec_x,
               tolerance = get_eps())
  expect_equal(target_coeff_vec_standardized_y,
               hermite_est$coeff_vec_y,
               tolerance = get_eps())
  expect_equal(0.3512403, sum(hermite_est$coeff_mat_bivar), tolerance = 
                 get_eps())
  expect_equal(mean(test_observations_mat[, 1]),
               hermite_est$running_mean_x,
               tolerance = get_eps())
  expect_equal(mean(test_observations_mat[, 2]),
               hermite_est$running_mean_y,
               tolerance = get_eps())
  expect_equal(sd(test_observations_mat[, 1]),
               sqrt(hermite_est$running_variance_x / (hermite_est$num_obs - 1)),
               tolerance = get_eps())
  expect_equal(sd(test_observations_mat[, 2]),
               sqrt(hermite_est$running_variance_y / (hermite_est$num_obs - 1)),
               tolerance = get_eps())
  hermite_est <-
    hermite_estimator(N = 10,
                      standardize = FALSE,
                      est_type = "bivariate", 
                      observations = test_observations_mat)
  expect_equal(target_coeff_vec_unstandardized_x,
               hermite_est$coeff_vec_x,
               tolerance = get_eps())
  expect_equal(target_coeff_vec_unstandardized_y,
               hermite_est$coeff_vec_y,
               tolerance = get_eps())
  expect_equal(0.3421721, sum(hermite_est$coeff_mat_bivar), tolerance = 
                 get_eps())
  hermite_est <-
    hermite_estimator(N = 10,
                      standardize = FALSE,
                      est_type = "bivariate", 
                      observations = c(0.1,0.4))
  expect_equal(0.3864793, sum(hermite_est$coeff_mat_bivar), tolerance = 
                 get_eps())
})

test_that("sequential updates of hermite_estimator_bivar work as expected",
          {
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
            target_coeff_vec_standardized_x <-
              c(
                0.535556916970579,
                -0.0507365295710208,
                -0.00156670275185479,
                0.157674444810238,
                -0.0606628649219171,
                -0.158694468305431,
                0.0441839726891637,
                0.117784644482994,
                -0.0251888871496467,
                -0.072754839522615,
                0.0151690464895391
              )
            target_coeff_vec_standardized_y <-
              c(
                0.572391675202002,
                0.00405311761349712,
                -0.100314203427144,
                0.0373548871733545,
                0.0168349948864179,
                -0.0884315184187305,
                0.0177333970783137,
                0.0951726612600252,
                -0.0346444220769052,
                -0.0687423875926508,
                0.0363964386028525
              )
            
            target_coeff_vec_unstandardized_x <-
              c(
                0.542468993917852,
                -0.112156489591074,
                0.0240427460518777,
                0.114702892134999,
                -0.0745566716013653,
                -0.0915076243058221,
                0.0423861774862675,
                0.0682739544693091,
                -0.0171370670453854,
                -0.0519546168194523,
                0.0128991995695314
              )
            
            target_coeff_vec_unstandardized_y <-
              c(
                0.514269409817415,
                0.0134314987076482,
                0.0103702696517751,
                0.0184650339583687,
                0.0251952541085223,
                -0.0312867094537839,
                -0.00130336810122759,
                0.0288137305459172,
                -0.0266496927060038,
                -0.00612587622958993,
                0.0340885404266732
              )
            hermite_est <-
              hermite_estimator(N = 10,
                                standardize = TRUE,
                                est_type = "bivariate")
            for (idx in seq_len(nrow(test_observations_mat))) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations_mat[idx, ])
            }
            expect_equal(target_coeff_vec_standardized_x,
                         hermite_est$coeff_vec_x,
                         tolerance = get_eps())
            expect_equal(target_coeff_vec_standardized_y,
                         hermite_est$coeff_vec_y,
                         tolerance = get_eps())
            hermite_est <-
              hermite_estimator(N = 10,
                                standardize = TRUE,
                                est_type = "bivariate")
            hermite_est <-
                hermite_est %>% update_sequential(test_observations_mat)
            expect_equal(target_coeff_vec_standardized_x,
                         hermite_est$coeff_vec_x,
                         tolerance = get_eps())
            expect_equal(target_coeff_vec_standardized_y,
                         hermite_est$coeff_vec_y,
                         tolerance = get_eps())
            hermite_est <-
              hermite_estimator(N = 10,
                                standardize = FALSE,
                                est_type = "bivariate")
            for (idx in seq_len(nrow(test_observations_mat))) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations_mat[idx, ])
            }
            expect_equal(target_coeff_vec_unstandardized_x,
                         hermite_est$coeff_vec_x,
                         tolerance = get_eps())
            expect_equal(target_coeff_vec_unstandardized_y,
                         hermite_est$coeff_vec_y,
                         tolerance = get_eps())
          })

test_that(
  "sequential updates of exponentially weighted hermite_estimator_bivar work as 
  expected",
  {
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
    target_coeff_vec_standardized_x <-
      c(
        0.530208153780545,
        -0.165941454379716,
        0.00607735015668864,
        0.194361633570863,
        -0.0908992345500389,
        -0.153683009259818,
        0.0821158510894216,
        0.0941532107052644,
        -0.0577436136870325,
        -0.0405289503799201,
        0.0372587571492455
      )
    target_coeff_vec_standardized_y <-
      c(
        0.5200324474035,
        -0.00827836295930781,
        -0.0557573666657626,
        0.0164207797661815,
        0.0833294366215888,
        -0.0812202951629325,
        -0.0671812180262155,
        0.0768030098379967,
        0.0267241063753009,
        -0.0277947231660045,
        -0.00253844809930239
      )
    
    target_coeff_vec_unstandardized_x <-
      c(
        0.542879997302748,
        -0.183440688713213,
        -0.0139071822407013,
        0.186309323167361,
        -0.0877891862176087,
        -0.147334159928653,
        0.0962429827300817,
        0.100966777115889,
        -0.0829918346143299,
        -0.0601496377415389,
        0.0664110392747781
      )
    
    target_coeff_vec_unstandardized_y <-
      c(
        0.513290877803333,
        0.0587833207693598,
        -0.0639598126648297,
        -0.00544884001564581,
        0.114061515378859,
        -0.0749995930467861,
        -0.0673110185132834,
        0.045620915276054,
        0.0141991476040826,
        0.0220769591254195,
        0.00551671755582751
      )
    hermite_est <-
      hermite_estimator(
        N = 10,
        standardize = TRUE,
        exp_weight_lambda = 0.1,
        est_type = "bivariate"
      )
    for (idx in seq_len(nrow(test_observations_mat))) {
      hermite_est <-
        hermite_est %>% update_sequential(test_observations_mat[idx, ])
    }
    expect_equal(target_coeff_vec_standardized_x,
                 hermite_est$coeff_vec_x,
                 tolerance = get_eps())
    expect_equal(target_coeff_vec_standardized_y,
                 hermite_est$coeff_vec_y,
                 tolerance = get_eps())
    hermite_est <-
      hermite_estimator(
        N = 10,
        standardize = FALSE,
        exp_weight_lambda = 0.1,
        est_type = "bivariate"
      )
    for (idx in seq_len(nrow(test_observations_mat))) {
      hermite_est <-
        hermite_est %>% update_sequential(test_observations_mat[idx, ])
    }
    expect_equal(target_coeff_vec_unstandardized_x,
                 hermite_est$coeff_vec_x,
                 tolerance = get_eps())
    expect_equal(target_coeff_vec_unstandardized_y,
                 hermite_est$coeff_vec_y,
                 tolerance = get_eps())
  }
)

test_that("bivariate hermite estimators merge consistently", {
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
  hermite_est <-
    hermite_estimator(N = 10,
                      standardize = FALSE,
                      est_type = "bivariate", 
                      observations = test_observations_mat)
  hermite_est_1 <-
    hermite_estimator(N = 10,
                      standardize = FALSE,
                      est_type = "bivariate", 
                      observations = test_observations_mat[1:10, ])
  hermite_est_2 <-
    hermite_estimator(N = 10,
                      standardize = FALSE,
                      est_type = "bivariate", 
                      observations = test_observations_mat[11:20, ])
  hermite_est_3 <-
    hermite_estimator(N = 10,
                      standardize = FALSE,
                      est_type = "bivariate", 
                      observations = test_observations_mat[21:30, ])
  hermite_merged <-
    merge_hermite_bivar(list(hermite_est_1))
  expect_equal(hermite_est_1, hermite_merged, tolerance = get_eps())
  hermite_merged <-
    merge_hermite(list(hermite_est_1))
  expect_equal(hermite_est_1, hermite_merged, tolerance = get_eps())
  hermite_merged <-
    merge_hermite_bivar(list(hermite_est_1, hermite_est_2, hermite_est_3))
  expect_equal(hermite_est, hermite_merged, tolerance = get_eps())
  hermite_merged <-
    merge_hermite(list(hermite_est_1, hermite_est_2, hermite_est_3))
  expect_equal(hermite_est, hermite_merged, tolerance = get_eps())
  hermite_est <-
    hermite_estimator(N = 10,
                      standardize = TRUE,
                      est_type = "bivariate", 
                      observations = test_observations_mat)
  hermite_est_1 <-
    hermite_estimator(N = 10,
                      standardize = TRUE,
                      est_type = "bivariate", 
                      observations = test_observations_mat[1:10, ])
  hermite_est_2 <-
    hermite_estimator(N = 10,
                      standardize = TRUE,
                      est_type = "bivariate", 
                      observations = test_observations_mat[11:20, ])
  hermite_est_3 <-
    hermite_estimator(N = 10,
                      standardize = TRUE,
                      est_type = "bivariate", 
                      observations = test_observations_mat[21:30, ]) 
  hermite_merged <-
    merge_pair(hermite_est_1, hermite_est_2)
  expect_equal(0.3236014, sum(hermite_merged$coeff_mat_bivar), tolerance =
                 get_eps())
  hermite_merged <-
    merge_hermite_bivar(list(hermite_est_1, hermite_est_2, hermite_est_3))
  hermite_merged_gen <-
    merge_hermite(list(hermite_est_1, hermite_est_2, hermite_est_3))
  expect_equal(hermite_merged,hermite_merged_gen)
  expect_equal(hermite_est$running_mean_x,
               hermite_merged$running_mean_x,
               tolerance = get_eps())
  expect_equal(hermite_est$running_variance_x,
               hermite_merged$running_variance_x,
               tolerance = get_eps())
  expect_equal(hermite_est$running_mean_y,
               hermite_merged$running_mean_y,
               tolerance = get_eps())
  expect_equal(hermite_est$running_variance_y,
               hermite_merged$running_variance_y,
               tolerance = get_eps())
  expect_equal(hermite_merged$num_obs, 30, tolerance = get_eps())
  
  hermite_est_univar_1_x <-
    hermite_estimator(N = 10, standardize = T, 
                      observations = test_observations_mat[1:10, 1])
  hermite_est_univar_1_y <-
    hermite_estimator(N = 10, standardize = T, 
                      observations = test_observations_mat[1:10, 2])
  hermite_est_univar_2_x <-
    hermite_estimator(N = 10, standardize = T, 
                      observations = test_observations_mat[11:20, 1])
  hermite_est_univar_2_y <-
    hermite_estimator(N = 10, standardize = T, 
                      observations = test_observations_mat[11:20, 2])
  hermite_est_univar_3_x <-
    hermite_estimator(N = 10, standardize = T, 
                      observations = test_observations_mat[21:30, 1])
  hermite_est_univar_3_y <-
    hermite_estimator(N = 10, standardize = T, 
                      observations = test_observations_mat[21:30, 2])
  hermite_merged_univar_x <-
    merge_hermite(list(
      hermite_est_univar_1_x,
      hermite_est_univar_2_x,
      hermite_est_univar_3_x
    ))
  hermite_merged_univar_y <-
    merge_hermite(list(
      hermite_est_univar_1_y,
      hermite_est_univar_2_y,
      hermite_est_univar_3_y
    ))
  expect_equal(hermite_merged$coeff_vec_x,
               hermite_merged_univar_x$coeff_vec,
               tolerance = get_eps())
  expect_equal(hermite_merged$coeff_vec_y,
               hermite_merged_univar_y$coeff_vec,
               tolerance = get_eps())
  expect_equal(sum(hermite_merged$coeff_mat_bivar), 0.3444092,
               tolerance = get_eps())
  hermite_est_1 <-
    hermite_estimator(N = 10, standardize = FALSE, est_type = "bivariate", 
                      observations = test_observations_mat[1:10,])
  hermite_est_2 <-
    hermite_estimator(N = 20, standardize = FALSE, est_type = "bivariate",
                      observations = test_observations_mat[11:20,])
  hermite_est_3 <-
    hermite_estimator(N = 10, standardize = TRUE, est_type = "bivariate",
                      observations = test_observations_mat[21:30,])
  hermite_est_4 <-
    hermite_estimator(N = 10, standardize = FALSE,exp_weight_lambda = 0.01,
                      est_type = "bivariate")
  expect_error(merge_hermite(list(hermite_est_1, hermite_est_2)))
  expect_error(merge_hermite(list(hermite_est_1, hermite_est_3)))
  expect_error(merge_hermite(list(hermite_est_1, hermite_est_4)))
  expect_error(merge_hermite_bivar(list(hermite_est_1, hermite_est_2)))
  expect_error(merge_hermite_bivar(list(hermite_est_1, hermite_est_3)))
  expect_error(merge_hermite_bivar(list(hermite_est_1, hermite_est_4)))
})

test_that("bivariate probability density estimation works as expected", {
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
  x_vals <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
  x_mat <- matrix(x_vals,
                  nrow = 5,
                  ncol = 2,
                  byrow = FALSE)
  
  hermite_est <-
    hermite_estimator(N = 10,
                      standardize = FALSE,
                      est_type = "bivariate", 
                      observations = test_observations_mat)
  pdf_vals <- hermite_est %>% dens(x_mat)
  target_pdf_vals_unstandardized <-
    c(
      0.0228867713027167,
      -0.0138748970102507,
      0.00945499275336495,
      0.00219206360051186,
      -0.00045534467461236
    )
  expect_equal(pdf_vals, target_pdf_vals_unstandardized, tolerance = get_eps())
  
  hermite_est <-
    hermite_estimator(N = 10,
                      standardize = TRUE,
                      est_type = "bivariate", 
                      observations = test_observations_mat)
  pdf_vals <- hermite_est %>% dens(x_mat)
  target_pdf_vals_standardized <-
    c(
      0.0143447460506261,
      -0.00682650501682584,
      0.00423940498866387,
      0.0028586637344242,
      -0.000421779530883988
    )
  expect_equal(pdf_vals, target_pdf_vals_standardized, tolerance = get_eps())
  
  target_pdf_vals_standardized <-
    c(0.0143447460506261,
      1e-08,
      0.00423940498866387,
      0.0028586637344242,
      1e-08)
  pdf_vals <- hermite_est %>% dens(x_mat, clipped = T)
  expect_equal(pdf_vals, target_pdf_vals_standardized, tolerance = get_eps())
  
  hermite_est <-
    hermite_estimator(
      N = 10,
      standardize = FALSE,
      exp_weight_lambda = 0.1,
      est_type = "bivariate"
    )
  for (idx in seq_len(nrow(test_observations_mat))) {
    hermite_est <-
      hermite_est %>% update_sequential(test_observations_mat[idx, ])
  }
  pdf_vals <- hermite_est %>% dens(x_mat)
  target_pdf_vals_unstandardized <-
    c(
      -0.00191362976526509,
      0.000315150156501173,
      0.00635052422650019,
      0.00163105161328314,
      -0.00106262644848195
    )
  expect_equal(pdf_vals, target_pdf_vals_unstandardized, tolerance = get_eps())
  
  pdf_vals <- hermite_est %>% dens(x_mat, clipped = TRUE)
  target_pdf_vals_unstandardized_clipped <-
    c(1e-08,
      0.000315150156501173,
      0.00635052422650019,
      0.00163105161328314,
      1e-08)
  expect_equal(pdf_vals,
               target_pdf_vals_unstandardized_clipped,
               tolerance = get_eps())
  
  hermite_est <-
    hermite_estimator(
      N = 10,
      standardize = TRUE,
      exp_weight_lambda = 0.1,
      est_type = "bivariate"
    )
  for (idx in seq_len(nrow(test_observations_mat))) {
    hermite_est <-
      hermite_est %>% update_sequential(test_observations_mat[idx, ])
  }
  pdf_vals <- hermite_est %>% dens(x_mat)
  target_pdf_vals_standardized <-
    c(
      -0.0113625236490677,
      0.00496291270778249,
      -0.00119028464006515,
      -0.000447725087235188,
      0.000498681697647531
    )
  expect_equal(pdf_vals, target_pdf_vals_standardized, tolerance = get_eps())
  hermite_est <-
    hermite_estimator(N = 10, est_type = "bivariate")
  pdf_vals <- hermite_est %>% dens(x_mat)
  expect_equal(length(pdf_vals), nrow(x_mat))
  expect_true(all(is.na(pdf_vals)))
})

test_that("bivariate cumulative probability estimation works as expected",
          {
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
            hermite_est <-
              hermite_estimator(N = 10,
                                standardize = FALSE,
                                est_type = "bivariate", 
                                observations = test_observations_mat)
            integrand_pdf <- function(h_est, x, y) {
              dens(h_est, c(x, y))
            }
            integrand_pdf_vec <-
              Vectorize(integrand_pdf, vectorize.args = c("x", "y"))
            inner_func <- function(x, y) {
              integrand_pdf_vec(hermite_est, x, y)
            }
            inner_integral <- Vectorize(function(y) {
              integrate(function(x) {
                inner_func(x, y)
              }, lower = -Inf, upper = 0.5)$value
            })
            target_cdf <- integrate(function(y) {
              inner_integral(y)
            }, lower = -Inf, upper = 1)$value
            est_cdf <- hermite_est %>% cum_prob(c(0.5, 1))
            expect_equal(est_cdf, target_cdf, tolerance = get_eps())
            
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
            hermite_est <-
              hermite_estimator(N = 10,
                                standardize = TRUE,
                                est_type = "bivariate", 
                                observations = test_observations_mat)
            inner_func <- function(x, y) {
              integrand_pdf_vec(hermite_est, x, y)
            }
            inner_integral <- Vectorize(function(y) {
              integrate(function(x) {
                inner_func(x, y)
              }, lower = -Inf, upper = 0.5)$value
            })
            target_cdf <- integrate(function(y) {
              inner_integral(y)
            }, lower = -Inf, upper = 1)$value
            est_cdf <- hermite_est %>% cum_prob(c(0.5, 1))
            expect_equal(est_cdf, target_cdf, tolerance = get_eps())
            
            hermite_est <-
              hermite_estimator(
                N = 10,
                standardize = TRUE,
                exp_weight_lambda = 0.1,
                est_type = "bivariate"
              )
            for (idx in seq_len(nrow(test_observations_mat))) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations_mat[idx, ])
            }
            inner_func <- function(x, y) {
              integrand_pdf_vec(hermite_est, x, y)
            }
            inner_integral <- Vectorize(function(y) {
              integrate(function(x) {
                inner_func(x, y)
              }, lower = -Inf, upper = 0.5)$value
            })
            target_cdf <- integrate(function(y) {
              inner_integral(y)
            }, lower = -Inf, upper = 1)$value
            est_cdf <- hermite_est %>% cum_prob(c(0.5, 1))
            expect_equal(est_cdf, target_cdf, tolerance = get_eps())
            
            hermite_est <-
              hermite_estimator(
                N = 10,
                standardize = FALSE,
                exp_weight_lambda = 0.1,
                est_type = "bivariate"
              )
            for (idx in seq_len(nrow(test_observations_mat))) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations_mat[idx, ])
            }
            inner_func <- function(x, y) {
              integrand_pdf_vec(hermite_est, x, y)
            }
            inner_integral <- Vectorize(function(y) {
              integrate(function(x) {
                inner_func(x, y)
              }, lower = -Inf, upper = 0.5)$value
            })
            target_cdf <- integrate(function(y) {
              inner_integral(y)
            }, lower = -Inf, upper = 1)$value
            est_cdf <- hermite_est %>% cum_prob(c(0.5, 1))
            expect_equal(est_cdf, target_cdf, tolerance = get_eps())
            
            x_vals <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
            x_mat <- matrix(x_vals,
                            nrow = 5,
                            ncol = 2,
                            byrow = FALSE)
            hermite_est <-
              hermite_estimator(N = 10,
                                standardize = FALSE,
                                est_type = "bivariate", 
                                observations = test_observations_mat)
            cdf_vals <- hermite_est %>% cum_prob(x_mat)
            target_cdf_vals_unstandardized <-
              c(
                0.870862110571656,
                0.867881158736368,
                1.0209208483946,
                1.02791959978022,
                1.02900881440701
              )
            expect_equal(cdf_vals, target_cdf_vals_unstandardized, 
                         tolerance = get_eps())
            hermite_est <-
              hermite_estimator(N = 10,
                                standardize = TRUE,
                                est_type = "bivariate", 
                                observations = test_observations_mat)
            cdf_vals <- hermite_est %>% cum_prob(x_mat)
            target_cdf_vals_standardized <-
              c(
                0.87449356687884,
                0.871314320667304,
                1.01860395621831,
                1.02209440814599,
                1.02819373498816
              )
            expect_equal(cdf_vals, target_cdf_vals_standardized, 
                         tolerance = get_eps())
            hermite_est <-
              hermite_estimator(N = 10,
                                standardize = TRUE,
                                est_type = "bivariate", 
                                observations = test_observations_mat)
            cdf_vals <- hermite_est %>% cum_prob(x_mat, clipped = TRUE)
            target_cdf_vals_standardized <-
              c(0.87449356687884, 0.871314320667304, 1, 1, 1)
            expect_equal(cdf_vals, target_cdf_vals_standardized, 
                         tolerance = get_eps())
            hermite_est <-
              hermite_estimator(
                N = 10,
                standardize = FALSE,
                exp_weight_lambda = 0.1,
                est_type = "bivariate"
              )
            for (idx in seq_len(nrow(test_observations_mat))) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations_mat[idx, ])
            }
            cdf_vals <- hermite_est %>% cum_prob(x_mat)
            target_cdf_vals_unstandardized <-
              c(
                0.855913097861038,
                0.856235166092503,
                0.969885108048889,
                0.979648555295636,
                0.986646514574038
              )
            expect_equal(cdf_vals, target_cdf_vals_unstandardized, 
                         tolerance = get_eps())
            hermite_est <-
              hermite_estimator(
                N = 10,
                standardize = TRUE,
                exp_weight_lambda = 0.1,
                est_type = "bivariate"
              )
            for (idx in seq_len(nrow(test_observations_mat))) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations_mat[idx, ])
            }
            cdf_vals <- hermite_est %>% cum_prob(x_mat)
            target_cdf_vals_standardized <-
              c(
                0.932106793550637,
                0.933634171610799,
                0.965211973901539,
                0.963558227610818,
                0.963606285800644
              )
            expect_equal(cdf_vals, target_cdf_vals_standardized, 
                         tolerance = get_eps())
            hermite_est <-
              hermite_estimator(N = 10, est_type = "bivariate")
            cdf_vals <- hermite_est %>% cum_prob(x_mat)
            expect_equal(length(cdf_vals), nrow(x_mat))
            expect_true(all(is.na(cdf_vals)))
          })

test_that("bivariate Spearman's correlation estimation works as expected",
          {
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
            hermite_est <-
              hermite_estimator(N = 10,
                                standardize = FALSE,
                                est_type = "bivariate", 
                                observations = test_observations_mat)
            est_spear <- hermite_est %>% spearmans()
            expect_equal(0.5486201, est_spear, tolerance = get_eps())
            hermite_est <-
              hermite_estimator(N = 10,
                                standardize = TRUE,
                                est_type = "bivariate", 
                                observations = test_observations_mat)
            est_spear <- hermite_est %>% spearmans()
            expect_equal(0.5639884, est_spear, tolerance = get_eps())
            hermite_est <-
              hermite_estimator(
                N = 10,
                standardize = FALSE,
                exp_weight_lambda = 0.1,
                est_type = "bivariate"
              )
            for (idx in seq_len(nrow(test_observations_mat))) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations_mat[idx, ])
            }
            est_spear <- hermite_est %>% spearmans()
            expect_equal(0.4455953, est_spear, tolerance = get_eps())
            hermite_est <-
              hermite_estimator(
                N = 10,
                standardize = TRUE,
                exp_weight_lambda = 0.1,
                est_type = "bivariate"
              )
            for (idx in seq_len(nrow(test_observations_mat))) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations_mat[idx, ])
            }
            est_spear <- hermite_est %>% spearmans()
            expect_equal(0.4494662, est_spear, tolerance = get_eps())
            hermite_est <-
              hermite_estimator(N = 10,
                                standardize = FALSE,
                                est_type = "bivariate", 
                                observations = test_observations_mat)
            spear_est <- hermite_est %>% spearmans()
            expect_equal(spear_est, 0.5486199, tolerance = get_eps())
            hermite_est <-
              hermite_estimator(N = 10,
                                standardize = TRUE,
                                est_type = "bivariate", 
                                observations = test_observations_mat)
            spear_est <- hermite_est %>% spearmans()
            expect_equal(spear_est, 0.5639886, tolerance = get_eps())
            
            hermite_est <-
              hermite_estimator(
                N = 10,
                standardize = FALSE,
                exp_weight_lambda = 0.1,
                est_type = "bivariate"
              )
            for (idx in seq_len(nrow(test_observations_mat))) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations_mat[idx, ])
            }
            spear_est <- hermite_est %>% spearmans()
            expect_equal(spear_est, 0.4455951, tolerance = get_eps())
            
            hermite_est <-
              hermite_estimator(
                N = 10,
                standardize = TRUE,
                exp_weight_lambda = 0.1,
                est_type = "bivariate"
              )
            for (idx in seq_len(nrow(test_observations_mat))) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations_mat[idx, ])
            }
            spear_est <- hermite_est %>% spearmans()
            expect_equal(spear_est, 0.4494662, tolerance = get_eps())
            
            hermite_est <-
              hermite_estimator(
                N = 10,
                standardize = TRUE,
                exp_weight_lambda = 0.1,
                est_type = "bivariate"
              )
            for (idx in seq_len(nrow(test_observations_mat))) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations_mat[idx, ])
            }
            spear_est <- hermite_est %>% spearmans(clipped = TRUE)
            expect_equal(spear_est, 0.4494662, tolerance = get_eps())
            
            hermite_est <-
              hermite_estimator(N = 10, est_type = "bivariate")
            spear_est <- hermite_est %>% spearmans()
            expect_equal(length(spear_est), 1)
            expect_true(all(is.na(spear_est)))
          })


test_that("bivariate Kendall correlation estimation works as expected",
          {
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
            hermite_est <-
              hermite_estimator(N = 10,
                                standardize = FALSE,
                                est_type = "bivariate", 
                                observations = test_observations_mat)
            est_kendall <- hermite_est %>% kendall()
            expect_equal(0.4500084, est_kendall, tolerance = get_eps())
            
            hermite_est <-
              hermite_estimator(N = 10,
                                standardize = TRUE,
                                est_type = "bivariate", 
                                observations = test_observations_mat)
            est_kendall <- hermite_est %>% kendall()
            expect_equal(0.4585442, est_kendall, tolerance = get_eps())
            hermite_est <-
              hermite_estimator(
                N = 10,
                standardize = FALSE,
                exp_weight_lambda = 0.1,
                est_type = "bivariate"
              )
            for (idx in seq_len(nrow(test_observations_mat))) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations_mat[idx, ])
            }
            est_kendall <- hermite_est %>% kendall()
            expect_equal(0.3050382, est_kendall, tolerance = get_eps())
            hermite_est <-
              hermite_estimator(
                N = 10,
                standardize = TRUE,
                exp_weight_lambda = 0.1,
                est_type = "bivariate"
              )
            for (idx in seq_len(nrow(test_observations_mat))) {
              hermite_est <-
                hermite_est %>% update_sequential(test_observations_mat[idx, ])
            }
            est_kendall <- hermite_est %>% kendall()
            expect_equal(0.2602773, est_kendall, tolerance = get_eps())
            est_kendall <- hermite_est %>% kendall(clipped = TRUE)
            expect_equal(0.2602773, est_kendall, tolerance = get_eps())
            hermite_est <-
              hermite_estimator(N = 10, est_type = "bivariate")
            kendall_est <- hermite_est %>% kendall()
            expect_equal(length(kendall_est), 1)
            expect_true(all(is.na(kendall_est)))
          })

test_that("Print and Summary work as expected", {
  hermite_est <- hermite_estimator(est_type = "bivariate")
  expect_equal(capture.output(print(hermite_est)), 
               c("Bivariate Hermite Estimator:",
                 "N = 30", "Standardize observations = TRUE",
                 "Exponential weighting for coefficents = FALSE",
                 "Number of observations = 0"))
  hermite_est <- hermite_estimator(est_type = "bivariate", 
                                   observations = matrix(c(1, 2, 3, 4,5, 6), 
                                                nrow=3, ncol=2, byrow = TRUE))
  expect_equal(capture.output(summary(hermite_est)), 
               c('Bivariate Hermite Estimator:','N = 30',
               'Standardize observations = TRUE',
               'Exponential weighting for coefficents = FALSE',
               'Number of observations = 3','','Mean x = 3','Mean y = 4',
               'Standard Deviation x = 2','Standard Deviation y = 2',
               'Spearman\'s Rho = 1.2154','Kendall Tau = 0.7662'))
})
