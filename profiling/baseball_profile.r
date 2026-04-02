library(profvis)
library(MMLN)

script_dir <- dirname(sys.frame(1)$ofile)
baseball_example <- clean_Lahman_data()

fmln_res_profile <- profvis({
  mlb_f <<- FMLN(
    Y = baseball_example$Y,
    X = cbind(1, baseball_example$X),
    n_iter = 100,
    burn_in = 30,
    proposal = "normbeta",
    mh_scale = .2,
    verbose = TRUE
  )
})

htmlwidgets::saveWidget(fmln_res_profile, file.path(script_dir, "fmln_baseball_profile.html"), selfcontained = TRUE)

mmln_res_profile <- profvis({
  mlb_m <<- MMLN(
    Y = baseball_example$Y,
    X = cbind(1, baseball_example$X),
    Z = baseball_example$Z,
    n_iter = 100,
    burn_in = 30,
    proposal = "normbeta",
    mh_scale = .2,
    verbose = TRUE
  )
})

htmlwidgets::saveWidget(mmln_res_profile, file.path(script_dir, "mmln_baseball_profile.html"), selfcontained = TRUE)


Y_pred_list_f <- lapply(seq_along(mlb_f$w_chain), function(i) {
  sample_posterior_predictive(X = cbind(1, baseball_example$X),
                              beta = mlb_f$beta_chain[[i]],
                              Sigma = mlb_f$sigma_chain[[i]],
                              n = baseball_example$PA,
                              mixed = FALSE,
                              verbose = FALSE
  )
})
mdres_fmln_baseball_profile <- profvis({
  resids_f <- MDres(baseball_example$Y, Y_pred_list_f)
})

htmlwidgets::saveWidget(mdres_fmln_baseball_profile, file.path(script_dir, "mdres_fmln_baseball_profile.html"), selfcontained = TRUE)


Y_pred_list_m <- lapply(seq_along(mlb_m$w_chain), function(i) {
  sample_posterior_predictive(X = cbind(1, baseball_example$X),
                              beta = mlb_m$beta_chain[[i]],
                              Sigma = mlb_m$sigma_chain[[i]],
                              n = baseball_example$PA,
                              Z = baseball_example$Z,
                              psi = mlb_m$psi_chain[[i]],
                              mixed = TRUE,
                              verbose = FALSE
  )
})
mdres_mmln_baseball_profile <- profvis({
resids_m <- MDres(baseball_example$Y, Y_pred_list_m)
})

htmlwidgets::saveWidget(mdres_mmln_baseball_profile, file.path(script_dir, "mdres_mmln_baseball_profile.html"), selfcontained = TRUE)


