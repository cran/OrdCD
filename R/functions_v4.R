admissible = function(i, j, gam_short_old) {
  if (gam_short_old[i, j]) {
    #delete an edge is always admissible
    return(TRUE)
  } else{
    gam_short_old[i, j] = 1
    return(gRbase::is.DAG(igraph::graph_from_adjacency_matrix(gam_short_old)))
  }
}

admissible_rev = function(i, j, gam_short_old) {
  if (gam_short_old[i, j] == gam_short_old[j, i]) {
    return(FALSE)
  } else{
    tmp = gam_short_old[i, j]
    gam_short_old[j, i] = tmp
    gam_short_old[i, j] = !tmp
    return(gRbase::is.DAG(igraph::graph_from_adjacency_matrix(gam_short_old)))
  }
}

mypolr = function(formula, data, ic, method, nq_y, nq_x) {
  boolFalse = FALSE
  if (nq_y > 2) {
    if (ic == "bic") {
      tryCatch({
        #sometimes MASS::polr fails to initialize
        IC = stats::BIC(MASS::polr(formula, data = data, method = method))
        boolFalse <- TRUE
      }, error = function(e) {

      }, finally = {

      })
      while (!boolFalse) {
        tryCatch({
          IC = stats::BIC(MASS::polr(
            formula,
            data = data,
            start = sort(stats::rnorm(nq_y - 1 + sum(nq_x - 1))),
            method = method
          ))
          boolFalse <- TRUE
        }, error = function(e) {

        }, finally = {

        })
      }
    } else if (ic == "aic") {
      tryCatch({
        IC = stats::AIC(MASS::polr(formula, data = data, method = method))
        boolFalse <- TRUE
      }, error = function(e) {

      }, finally = {

      })
      while (!boolFalse) {
        tryCatch({
          IC = stats::AIC(MASS::polr(
            formula,
            data = data,
            start = sort(stats::rnorm(nq_y - 1 + sum(nq_x - 1))),
            method = method
          ))
          boolFalse <- TRUE

        }, error = function(e) {

        }, finally = {

        })
      }
    }
  } else{
    if (ic == "bic") {
      tryCatch({
        IC = stats::BIC(stats::glm(
          formula,
          data = data ,
          family = stats::binomial(link = method)
        ))
        boolFalse <- TRUE
      }, error = function(e) {

      }, finally = {

      })
      while (!boolFalse) {
        tryCatch({
          IC = stats::BIC(stats::glm(
            formula,
            data = data,
            start = sort(stats::rnorm(nq_y - 1 + sum(nq_x - 1))),
            family = stats::binomial(link = method)
          ))
          boolFalse <- TRUE
        }, error = function(e) {

        }, finally = {

        })
      }
    } else if (ic == "aic") {
      tryCatch({
        IC = stats::AIC(stats::glm(
          formula,
          data = data ,
          family = stats::binomial(link = method)
        ))
        boolFalse <- TRUE
      }, error = function(e) {

      }, finally = {

      })
      while (!boolFalse) {
        tryCatch({
          IC = stats::AIC(stats::glm(
            formula,
            data = data,
            start = sort(stats::rnorm(nq_y - 1 + sum(nq_x - 1))),
            family = stats::binomial(link = method)
          ))
          boolFalse <- TRUE

        }, error = function(e) {

        }, finally = {

        })
      }
    }
  }
  return(IC)
}



oBN_greedy = function(y,
                      gam = NULL,
                      ic = "bic",
                      method = "probit",
                      verbose = verbose,
                      maxit = maxit) {
  #hill-climbing


  n = nrow(y)
  q = ncol(y)
  nq = rep(0, q)
  for (i in 1:q) {
    nq[i] = nlevels(y[, i])
  }
  if (is.null(gam)) {
    gam = matrix(FALSE, q, q)
  } else{
    gam = (gam != 0)
  }


  ind_q = matrix(0, q, q - 1)
  for (i in 1:q) {
    if (i == 1) {
      ind_noi = 2:q
    } else if (i == q) {
      ind_noi = 1:(q - 1)
    } else{
      ind_noi = c(1:(i - 1), (i + 1):q)
    }
    ind_q[i,] = ind_noi
  }

  iter = 0
  ic_improv = 1
  act_ind = c(NA, NA)
  state = "add" # or "del"
  if (ic == "bic") {
    ic_best = rep(0, q)
    for (i in 1:q) {
      if (sum(gam[i, ]) > 0) {
        ic_best[i] = mypolr(
          y[, i] ~ .,
          data = y[, gam[i, ]],
          ic = ic,
          method = method,
          nq_y = nq[i],
          nq_x = nq[gam[i, ]]
        )
      } else{
        if (nq[i] > 2) {
          ic_best[i] = stats::BIC(MASS::polr(y[, i] ~ 1, method = method))
        } else{
          ic_best[i] = stats::BIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
        }
      }
    }
    while (ic_improv > 0 && iter < maxit) {
      iter = iter + 1
      ic_improv = -Inf
      ic_improv_rev = rep(-Inf, 2)
      gam_new = gam
      ic_improv_new = -Inf
      ic_improv_rev_new = rep(-Inf, 2)
      ic_best_new = -Inf
      ic_rev_best_new = rep(-Inf, 2)
      for (i in 1:q) {
        for (j in 1:(q - 1)) {
          if (admissible(i, ind_q[i, j], gam)) {
            if (gam[i, ind_q[i, j]]) {
              #delete
              gam_new[i, ind_q[i, j]] = FALSE
              if (sum(gam_new[i, ]) > 0) {
                ic_best_new = mypolr(
                  y[, i] ~ .,
                  data = y[, gam_new[i, ]],
                  ic = ic,
                  method = method,
                  nq_y = nq[i],
                  nq_x = nq[gam_new[i, ]]
                )
              } else{
                if (nq[i] > 2) {
                  ic_best_new = stats::BIC(MASS::polr(y[, i] ~ 1, method = method))
                } else{
                  ic_best_new = stats::BIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
                }
              }
              ic_improv_new = ic_best[i] - ic_best_new
              if (ic_improv_new > ic_improv) {
                ic_improv = ic_improv_new
                act_ind = c(i, ind_q[i, j])
                state = "del"
              }
              gam_new[i, ind_q[i, j]] = TRUE
            } else{
              #add
              gam_new[i, ind_q[i, j]] = TRUE
              if (sum(gam_new[i, ]) > 0) {
                ic_best_new = mypolr(
                  y[, i] ~ .,
                  data = y[, gam_new[i, ]],
                  ic = ic,
                  method = method,
                  nq_y = nq[i],
                  nq_x = nq[gam_new[i, ]]
                )
              } else{
                if (nq[i] > 2) {
                  ic_best_new = stats::BIC(MASS::polr(y[, i] ~ 1, method = method))
                } else{
                  ic_best_new = stats::BIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
                }
              }
              ic_improv_new = ic_best[i] - ic_best_new
              if (ic_improv_new > ic_improv) {
                ic_improv = ic_improv_new
                act_ind = c(i, ind_q[i, j])
                state = "add"
              }
              gam_new[i, ind_q[i, j]] = FALSE
            }
          }
        }
      }
      #reverse edge
      for (i in 1:q) {
        for (j in 1:(q - 1)) {
          if (admissible_rev(i, ind_q[i, j], gam)) {
            tmp = gam_new[i, ind_q[i, j]]
            gam_new[ind_q[i, j], i] = tmp
            gam_new[i, ind_q[i, j]] = !tmp
            if (sum(gam_new[i, ]) > 0) {
              ic_rev_best_new[1] = mypolr(
                y[, i] ~ .,
                data = y[, gam_new[i, ]],
                ic = ic,
                method = method,
                nq_y = nq[i],
                nq_x = nq[gam_new[i, ]]
              )
            } else{
              if (nq[i] > 2) {
                ic_rev_best_new[1] = stats::BIC(MASS::polr(y[, i] ~ 1, method = method))
              } else{
                ic_rev_best_new[1] = stats::BIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
              }
            }
            if (gam_new[ind_q[i, j], i] > 0) {
              ic_rev_best_new[2] = mypolr(
                y[, ind_q[i, j]] ~ .,
                data = y[, gam_new[ind_q[i, j], ]],
                ic = ic,
                method = method,
                nq_y = nq[ind_q[i, j]],
                nq_x = nq[gam_new[ind_q[i, j], ]]
              )
            } else{
              if (nq[i] > 2) {
                ic_rev_best_new[2] = stats::BIC(MASS::polr(y[, ind_q[i, j]] ~ 1, method = method))
              } else{
                ic_rev_best_new[2] = stats::BIC(stats::glm(y[, ind_q[i, j]] ~ 1, family = stats::binomial(link = method)))
              }
            }
            ic_improv_rev_new[1] = ic_best[i] - ic_rev_best_new[1]
            ic_improv_rev_new[2] = ic_best[ind_q[i, j]] - ic_rev_best_new[2]
            ic_improv_new = ic_improv_rev_new[1] + ic_improv_rev_new[2]
            if (ic_improv_new > ic_improv) {
              ic_improv = ic_improv_new
              ic_improv_rev = ic_improv_rev_new
              act_ind = c(i, ind_q[i, j])
              state = "rev"
            }
            gam_new[i, ind_q[i, j]] = tmp
            gam_new[ind_q[i, j], i] = !tmp
          }
        }
      }

      if (ic_improv > 0) {
        if (state == "add") {
          gam[act_ind[1], act_ind[2]] = TRUE
          ic_best[act_ind[1]] = ic_best[act_ind[1]] - ic_improv
        } else if (state == "del") {
          gam[act_ind[1], act_ind[2]] = FALSE
          ic_best[act_ind[1]] = ic_best[act_ind[1]] - ic_improv
        } else if (state == "rev") {
          tmp = gam[act_ind[1], act_ind[2]]
          gam[act_ind[2], act_ind[1]] = tmp
          gam[act_ind[1], act_ind[2]] = !tmp
          ic_best[act_ind[1]] = ic_best[act_ind[1]] - ic_improv_rev[1]
          ic_best[act_ind[2]] = ic_best[act_ind[2]] - ic_improv_rev[2]
        }
      }
      if (verbose && iter %% 1 == 0) {
        print(paste(iter, " iterations have completed", sep = ""))
        print("The current DAG adjacency matrix is")
        print(gam)
        print(paste("with ",  ic, " = ", sum(ic_best), sep = ""))
      }
    }
  } else if (ic == "aic") {
    ic_best = rep(0, q)
    for (i in 1:q) {
      if (sum(gam[i, ]) > 0) {
        ic_best[i] = mypolr(
          y[, i] ~ .,
          data = y[, gam[i, ]],
          ic = ic,
          method = method,
          nq_y = nq[i],
          nq_x = nq[gam[i, ]]
        )
      } else{
        if (nq[i] > 2) {
          ic_best[i] = stats::AIC(MASS::polr(y[, i] ~ 1, method = method))
        } else{
          ic_best[i] = stats::AIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
        }
      }
    }
    while (ic_improv > 0 & iter < maxit) {
      iter = iter + 1
      ic_improv = -Inf
      ic_improv_rev = rep(-Inf, 2)
      gam_new = gam
      ic_improv_new = -Inf
      ic_improv_rev_new = rep(-Inf, 2)
      ic_best_new = -Inf
      ic_rev_best_new = rep(-Inf, 2)
      for (i in 1:q) {
        for (j in 1:(q - 1)) {
          if (admissible(i, ind_q[i, j], gam)) {
            if (gam[i, ind_q[i, j]]) {
              #delete
              gam_new[i, ind_q[i, j]] = FALSE
              if (sum(gam_new[i, ]) > 0) {
                ic_best_new = mypolr(
                  y[, i] ~ .,
                  data = y[, gam_new[i, ]],
                  ic = ic,
                  method = method,
                  nq_y = nq[i],
                  nq_x = nq[gam_new[i, ]]
                )
              } else{
                if (nq[i] > 2) {
                  ic_best_new = stats::AIC(MASS::polr(y[, i] ~ 1, method = method))
                } else{
                  ic_best_new = stats::AIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
                }
              }
              ic_improv_new = ic_best[i] - ic_best_new
              if (ic_improv_new > ic_improv) {
                ic_improv = ic_improv_new
                act_ind = c(i, ind_q[i, j])
                state = "del"
              }
              gam_new[i, ind_q[i, j]] = TRUE
            } else{
              #add
              gam_new[i, ind_q[i, j]] = TRUE
              if (sum(gam_new[i, ]) > 0) {
                ic_best_new = mypolr(
                  y[, i] ~ .,
                  data = y[, gam_new[i, ]],
                  ic = ic,
                  method = method,
                  nq_y = nq[i],
                  nq_x = nq[gam_new[i, ]]
                )
              } else{
                if (nq[i] > 2) {
                  ic_best_new = stats::AIC(MASS::polr(y[, i] ~ 1, method = method))
                } else{
                  ic_best_new = stats::AIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
                }
              }
              ic_improv_new = ic_best[i] - ic_best_new
              if (ic_improv_new > ic_improv) {
                ic_improv = ic_improv_new
                act_ind = c(i, ind_q[i, j])
                state = "add"
              }
              gam_new[i, ind_q[i, j]] = FALSE
            }
          }
        }
      }
      #reverse edge
      for (i in 1:q) {
        for (j in 1:(q - 1)) {
          if (admissible_rev(i, ind_q[i, j], gam)) {
            tmp = gam_new[i, ind_q[i, j]]
            gam_new[ind_q[i, j], i] = tmp
            gam_new[i, ind_q[i, j]] = !tmp
            if (sum(gam_new[i, ]) > 0) {
              ic_rev_best_new[1] = mypolr(
                y[, i] ~ .,
                data = y[, gam_new[i, ]],
                ic = ic,
                method = method,
                nq_y = nq[i],
                nq_x = nq[gam_new[i, ]]
              )
            } else{
              if (nq[i] > 2) {
                ic_rev_best_new[1] = stats::AIC(MASS::polr(y[, i] ~ 1, method = method))
              } else{
                ic_rev_best_new[1] = stats::AIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
              }
            }
            if (gam_new[ind_q[i, j], i] > 0) {
              ic_rev_best_new[2] = mypolr(
                y[, ind_q[i, j]] ~ .,
                data = y[, gam_new[ind_q[i, j], ]],
                ic = ic,
                method = method,
                nq_y = nq[ind_q[i, j]],
                nq_x = nq[gam_new[ind_q[i, j], ]]
              )
            } else{
              if (nq[i] > 2) {
                ic_rev_best_new[2] = stats::AIC(MASS::polr(y[, ind_q[i, j]] ~ 1, method = method))
              } else{
                ic_rev_best_new[2] = stats::AIC(stats::glm(y[, ind_q[i, j]] ~ 1, family = stats::binomial(link = method)))
              }
            }
            ic_improv_rev_new[1] = ic_best[i] - ic_rev_best_new[1]
            ic_improv_rev_new[2] = ic_best[ind_q[i, j]] - ic_rev_best_new[2]
            ic_improv_new = ic_improv_rev_new[1] + ic_improv_rev_new[2]
            if (ic_improv_new > ic_improv) {
              ic_improv = ic_improv_new
              ic_improv_rev = ic_improv_rev_new
              act_ind = c(i, ind_q[i, j])
              state = "rev"
            }
            gam_new[i, ind_q[i, j]] = tmp
            gam_new[ind_q[i, j], i] = !tmp
          }
        }
      }

      if (ic_improv > 0) {
        if (state == "add") {
          gam[act_ind[1], act_ind[2]] = TRUE
          ic_best[act_ind[1]] = ic_best[act_ind[1]] - ic_improv
        } else if (state == "del") {
          gam[act_ind[1], act_ind[2]] = FALSE
          ic_best[act_ind[1]] = ic_best[act_ind[1]] - ic_improv
        } else if (state == "rev") {
          tmp = gam[act_ind[1], act_ind[2]]
          gam[act_ind[2], act_ind[1]] = tmp
          gam[act_ind[1], act_ind[2]] = !tmp
          ic_best[act_ind[1]] = ic_best[act_ind[1]] - ic_improv_rev[1]
          ic_best[act_ind[2]] = ic_best[act_ind[2]] - ic_improv_rev[2]
        }
      }
      if (verbose && iter %% 1 == 0) {
        print(paste(iter, " iterations have completed", sep = ""))
        print("The current DAG adjacency matrix is")
        print(gam)
        print(paste("with ",  ic, " = ", sum(ic_best), sep = ""))
      }
    }
  }
  if (iter == maxit) {
    warning("The maximum number of iterations was reached. The algorithm has not converged.")
  }
  return(list(gam = gam + 0, ic_best = sum(ic_best)))
}

oBN_greedy_CPDAG = function(y,
                            gam = NULL,
                            ic = "bic",
                            edge_list = NULL,
                            method = "probit",
                            verbose = verbose,
                            maxit = maxit) {
  #hill-climbing


  n = nrow(y)
  q = ncol(y)
  nq = rep(0, q)
  for (i in 1:q) {
    nq[i] = nlevels(y[, i])
  }
  if (is.null(gam)) {
    gam = matrix(FALSE, q, q)
  } else{
    gam = (gam != 0)
  }

  ind_q = vector("list", q)
  nind_q = rep(0, q)
  for (e in 1:nrow(edge_list)) {
    i = edge_list[e, 1]
    ind_q[[i]] = c(ind_q[[i]], edge_list[e, 2])
    nind_q[i] = nind_q[i] + 1
  }


  iter = 0
  ic_improv = 1
  act_ind = c(NA, NA)
  if (ic == "bic") {
    ic_best = rep(0, q)
    for (i in 1:q) {
      if (sum(gam[i,]) > 0) {
        ic_best[i] = mypolr(
          y[, i] ~ .,
          data = y[, gam[i,]],
          ic = ic,
          method = method,
          nq_y = nq[i],
          nq_x = nq[gam[i,]]
        )
      } else{
        if (nq[i] > 2) {
          ic_best[i] = stats::BIC(MASS::polr(y[, i] ~ 1, method = method))
        } else{
          ic_best[i] = stats::BIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
        }
      }
    }
    while (ic_improv > 0 && iter < maxit) {
      iter = iter + 1
      ic_improv = -Inf
      ic_improv_rev = rep(-Inf, 2)
      gam_new = gam
      ic_improv_new = -Inf
      ic_improv_rev_new = rep(-Inf, 2)
      ic_best_new = -Inf
      ic_rev_best_new = rep(-Inf, 2)
      #reverse edge
      for (i in 1:q) {
        if (nind_q[i] > 0) {
          for (j in 1:(nind_q[i])) {
            if (admissible_rev(i, ind_q[[i]][j], gam)) {
              tmp = gam_new[i, ind_q[[i]][j]]
              gam_new[ind_q[[i]][j], i] = tmp
              gam_new[i, ind_q[[i]][j]] = !tmp
              if (sum(gam_new[i,]) > 0) {
                ic_rev_best_new[1] = mypolr(
                  y[, i] ~ .,
                  data = y[, gam_new[i,]],
                  ic = ic,
                  method = method,
                  nq_y = nq[i],
                  nq_x = nq[gam_new[i,]]
                )
              } else{
                if (nq[i] > 2) {
                  ic_rev_best_new[1] = stats::BIC(MASS::polr(y[, i] ~ 1, method = method))
                } else{
                  ic_rev_best_new[1] = stats::BIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
                }
              }
              if (gam_new[ind_q[[i]][j], i] > 0) {
                ic_rev_best_new[2] = mypolr(
                  y[, ind_q[[i]][j]] ~ .,
                  data = y[, gam_new[ind_q[[i]][j],]],
                  ic = ic,
                  method = method,
                  nq_y = nq[ind_q[[i]][j]],
                  nq_x = nq[gam_new[ind_q[[i]][j],]]
                )
              } else{
                if (nq[i] > 2) {
                  ic_rev_best_new[2] = stats::BIC(MASS::polr(y[, ind_q[[i]][j]] ~ 1, method = method))
                } else{
                  ic_rev_best_new[2] = stats::BIC(stats::glm(y[, ind_q[[i]][j]] ~ 1, family = stats::binomial(link = method)))
                }
              }
              ic_improv_rev_new[1] = ic_best[i] - ic_rev_best_new[1]
              ic_improv_rev_new[2] = ic_best[ind_q[[i]][j]] - ic_rev_best_new[2]
              ic_improv_new = ic_improv_rev_new[1] + ic_improv_rev_new[2]
              if (ic_improv_new > ic_improv) {
                ic_improv = ic_improv_new
                ic_improv_rev = ic_improv_rev_new
                act_ind = c(i, ind_q[[i]][j])
              }
              gam_new[i, ind_q[[i]][j]] = tmp
              gam_new[ind_q[[i]][j], i] = !tmp
            }
          }
        }
      }

      if (ic_improv > 0) {
        tmp = gam[act_ind[1], act_ind[2]]
        gam[act_ind[2], act_ind[1]] = tmp
        gam[act_ind[1], act_ind[2]] = !tmp
        ic_best[act_ind[1]] = ic_best[act_ind[1]] - ic_improv_rev[1]
        ic_best[act_ind[2]] = ic_best[act_ind[2]] - ic_improv_rev[2]
      }
      if (verbose && iter %% 1 == 0) {
        print(paste(iter, " iterations have completed", sep = ""))
        print("The current DAG adjacency matrix is")
        print(gam)
        print(paste("with ",  ic, " = ", sum(ic_best), sep = ""))
      }
    }
  } else if (ic == "aic") {
    ic_best = rep(0, q)
    for (i in 1:q) {
      if (sum(gam[i,]) > 0) {
        ic_best[i] = mypolr(
          y[, i] ~ .,
          data = y[, gam[i,]],
          ic = ic,
          method = method,
          nq_y = nq[i],
          nq_x = nq[gam[i,]]
        )
      } else{
        if (nq[i] > 2) {
          ic_best[i] = stats::AIC(MASS::polr(y[, i] ~ 1, method = method))
        } else{
          ic_best[i] = stats::AIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
        }
      }
    }
    while (ic_improv > 0 && iter < maxit) {
      iter = iter + 1
      ic_improv = -Inf
      ic_improv_rev = rep(-Inf, 2)
      gam_new = gam
      ic_improv_new = -Inf
      ic_improv_rev_new = rep(-Inf, 2)
      ic_best_new = -Inf
      ic_rev_best_new = rep(-Inf, 2)
      #reverse edge
      for (i in 1:q) {
        if (nind_q[i] > 0) {
          for (j in 1:(nind_q[i])) {
            if (admissible_rev(i, ind_q[[i]][j], gam)) {
              tmp = gam_new[i, ind_q[[i]][j]]
              gam_new[ind_q[[i]][j], i] = tmp
              gam_new[i, ind_q[[i]][j]] = !tmp
              if (sum(gam_new[i,]) > 0) {
                ic_rev_best_new[1] = mypolr(
                  y[, i] ~ .,
                  data = y[, gam_new[i,]],
                  ic = ic,
                  method = method,
                  nq_y = nq[i],
                  nq_x = nq[gam_new[i,]]
                )
              } else{
                if (nq[i] > 2) {
                  ic_rev_best_new[1] = stats::AIC(MASS::polr(y[, i] ~ 1, method = method))
                } else{
                  ic_rev_best_new[1] = stats::AIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
                }
              }
              if (gam_new[ind_q[[i]][j], i] > 0) {
                ic_rev_best_new[2] = mypolr(
                  y[, ind_q[[i]][j]] ~ .,
                  data = y[, gam_new[ind_q[[i]][j],]],
                  ic = ic,
                  method = method,
                  nq_y = nq[ind_q[[i]][j]],
                  nq_x = nq[gam_new[ind_q[[i]][j],]]
                )
              } else{
                if (nq[i] > 2) {
                  ic_rev_best_new[2] = stats::AIC(MASS::polr(y[, ind_q[[i]][j]] ~ 1, method = method))
                } else{
                  ic_rev_best_new[2] = stats::AIC(stats::glm(y[, ind_q[[i]][j]] ~ 1, family = stats::binomial(link = method)))
                }
              }
              ic_improv_rev_new[1] = ic_best[i] - ic_rev_best_new[1]
              ic_improv_rev_new[2] = ic_best[ind_q[[i]][j]] - ic_rev_best_new[2]
              ic_improv_new = ic_improv_rev_new[1] + ic_improv_rev_new[2]
              if (ic_improv_new > ic_improv) {
                ic_improv = ic_improv_new
                ic_improv_rev = ic_improv_rev_new
                act_ind = c(i, ind_q[[i]][j])
              }
              gam_new[i, ind_q[[i]][j]] = tmp
              gam_new[ind_q[[i]][j], i] = !tmp
            }
          }
        }
      }

      if (ic_improv > 0) {
        tmp = gam[act_ind[1], act_ind[2]]
        gam[act_ind[2], act_ind[1]] = tmp
        gam[act_ind[1], act_ind[2]] = !tmp
        ic_best[act_ind[1]] = ic_best[act_ind[1]] - ic_improv_rev[1]
        ic_best[act_ind[2]] = ic_best[act_ind[2]] - ic_improv_rev[2]
      }
      if (verbose && iter %% 1 == 0) {
        print(paste(iter, " iterations have completed", sep = ""))
        print("The current DAG adjacency matrix is")
        print(gam)
        print(paste("with ",  ic, " = ", sum(ic_best), sep = ""))
      }
    }
  }
  if (iter == maxit) {
    warning("The maximum number of iterations was reached. The algorithm has not converged.")
  }
  return(list(gam = gam + 0, ic_best = sum(ic_best)))
}


#allow multiple runs of the greedy search with different initial DAGs
oBN_greedy_wrap = function(y,
                           ic = "bic",
                           edge_list = NULL,
                           method = "probit",
                           nstart = 1,
                           verbose = verbose,
                           maxit = maxit) {
  q = ncol(y)
  gam_list  = vector("list", nstart)
  ic_best_list = rep(NA, nstart)
  if (nstart == 1) {
    if (is.null(edge_list)) {
      fit = oBN_greedy(
        y,
        gam = NULL,
        ic = ic,
        method = method,
        verbose = verbose,
        maxit = maxit
      )
    } else{
      gam = matrix(0, q, q)
      gam[edge_list] = 1
      und_edge = which(as.matrix(Matrix::tril(gam * t(gam))) == 1, arr.ind = TRUE)
      for (i in 1:nrow(und_edge)) {
        if (stats::rbinom(1, 1, .5) == 1) {
          gam[und_edge[i, 1], und_edge[i, 2]] = 0
        } else{
          gam[und_edge[i, 2], und_edge[i, 1]] = 0
        }
      }
      fit = oBN_greedy_CPDAG(
        y,
        gam = gam,
        ic = ic,
        edge_list = und_edge,
        method = method,
        verbose = verbose,
        maxit = maxit
      )
    }
    gam_list[[1]] = fit$gam
    ic_best_list[1] = fit$ic_best
  } else{
    if (is.null(edge_list)) {
      netlist = bnlearn::random.graph(
        nodes = as.character(1:q),
        method = "ordered",
        num = nstart - 1,
        prob = 1 / q
      )
      if (nstart == 2) {
        netlist = list(netlist)
      }
      for (i in 1:nstart) {
        gam = matrix(FALSE, q, q)
        if (i != 1) {
          gam[apply(netlist[[i - 1]]$arcs, 2, as.numeric)] = TRUE
        }
        fit = oBN_greedy(
          y,
          gam = gam,
          ic = ic,
          method = method,
          verbose = verbose,
          maxit = maxit
        )

        gam_list[[i]] = fit$gam
        ic_best_list[i] = fit$ic_best
      }
    } else{
      gam = matrix(0, q, q)
      gam[edge_list] = 1
      und_edge = which(as.matrix(Matrix::tril(gam * t(gam))) == 1, arr.ind = TRUE)
      for (i in 1:nstart) {
        gam_ini = gam

        for (ii in 1:nrow(und_edge)) {
          if (stats::rbinom(1, 1, .5) == 1) {
            gam_ini[und_edge[ii, 1], und_edge[ii, 2]] = 0
          } else{
            gam_ini[und_edge[ii, 2], und_edge[ii, 1]] = 0
          }
        }
        fit = oBN_greedy_CPDAG(
          y,
          gam = gam_ini,
          ic = ic,
          edge_list = und_edge,
          method = method,
          verbose = verbose,
          maxit = maxit
        )
        gam_list[[i]] = fit$gam
        ic_best_list[i] = fit$ic_best
      }
    }
  }
  i = which.min(ic_best_list)
  gam = gam_list[[i]]
  ic_best = ic_best_list[i]
  return(list(gam = gam, ic_best = ic_best))
  # }
}

oBN_exhaust = function(y,
                       gam_list = NULL,
                       ic = "bic",
                       method = "probit") {
  n = nrow(y)
  q = ncol(y)
  if (is.null(gam_list)) {
    if (q == 2) {
      gam_list = array(0, c(q, q, 3))
      gam_list[1, 2, 2] = 1
      gam_list[2, 1, 3] = 1
    } else if (q == 3) {
      gam_list = array(0, c(q, q, 25))
      gam_list[1, 2, 2] = 1
      gam_list[1, 3, 3] = 1
      gam_list[2, 3, 4] = 1
      gam_list[2, 1, 5] = 1
      gam_list[3, 1, 6] = 1
      gam_list[3, 2, 7] = 1
      gam_list[1, 2, 8] = gam_list[1, 3, 8] = 1
      gam_list[2, 3, 9] = gam_list[2, 1, 9] = 1
      gam_list[3, 1, 10] = gam_list[3, 2, 10] = 1

      gam_list[2, 1, 11] = gam_list[3, 1, 11] = 1
      gam_list[3, 2, 12] = gam_list[1, 2, 12] = 1
      gam_list[1, 3, 13] = gam_list[2, 3, 13] = 1

      gam_list[1, 2, 14] = gam_list[2, 3, 14] = 1
      gam_list[1, 3, 15] = gam_list[3, 2, 15] = 1
      gam_list[2, 3, 16] = gam_list[3, 1, 16] = 1
      gam_list[2, 1, 17] = gam_list[1, 3, 17] = 1
      gam_list[3, 2, 18] = gam_list[2, 1, 18] = 1
      gam_list[3, 1, 19] = gam_list[1, 2, 19] = 1

      gam_list[1, 2, 20] = gam_list[1, 3, 20] = gam_list[2, 3, 20] = 1
      gam_list[1, 2, 21] = gam_list[1, 3, 21] = gam_list[3, 2, 21] = 1
      gam_list[2, 1, 22] = gam_list[2, 3, 22] = gam_list[1, 3, 22] = 1
      gam_list[2, 1, 23] = gam_list[2, 3, 23] = gam_list[3, 1, 23] = 1
      gam_list[3, 1, 24] = gam_list[3, 2, 24] = gam_list[1, 2, 24] = 1
      gam_list[3, 1, 25] = gam_list[3, 2, 25] = gam_list[2, 1, 25] = 1
    } else{
      stop("The number of nodes must be 2 or 3")
    }
  }
  IC = rep(0, dim(gam_list)[3])

  nl = rep(0, q)
  for (i in 1:q) {
    nl[i] = nlevels(y[, i])
  }

  if (ic == "bic") {
    for (m in 1:length(IC)) {
      gam = gam_list[, , m]
      for (i in 1:q) {
        gam_tmp = gam[i,]
        if (sum(gam_tmp) == 0) {
          if (nl[i] == 2) {
            IC[m] = IC[m] + stats::BIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
          } else{
            IC[m] = IC[m] + stats::BIC(MASS::polr(y[, i] ~ 1, method = method))
          }
        } else{
          if (nl[i] == 2) {
            IC[m] = IC[m] + stats::BIC(stats::glm(
              y[, i] ~ .,
              data = y[, gam_tmp == 1],
              family = stats::binomial(link = method)
            ))
          } else{
            boolFalse = FALSE
            tryCatch({
              #sometimes MASS::polr fails to initialize
              IC[m] = IC[m] + stats::BIC(MASS::polr(y[, i] ~ ., data = y[, gam_tmp ==
                                                                           1], method = method))
              boolFalse <- TRUE
            }, error = function(e) {

            }, finally = {

            })

            if (!boolFalse) {
              IC[m] = IC[m] + stats::BIC(MASS::polr(
                y[, i] ~ .,
                data = y[, gam_tmp == 1],
                start = stats::rnorm(nl[i] - 1 + sum(nl[gam_tmp == 1] - 1)),
                method = method
              ))
            }
          }
        }
      }
    }
  } else if (ic == "aic") {
    for (m in 1:length(IC)) {
      gam = gam_list[, , m]
      for (i in 1:q) {
        gam_tmp = gam[i,]
        if (sum(gam_tmp) == 0) {
          if (nl[i] == 2) {
            IC[m] = IC[m] + stats::AIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
          } else{
            IC[m] = IC[m] + stats::AIC(MASS::polr(y[, i] ~ 1, method = method))
          }
        } else{
          if (nl[i] == 2) {
            IC[m] = IC[m] + stats::AIC(stats::glm(
              y[, i] ~ .,
              data = y[, gam_tmp == 1],
              family = stats::binomial(link = method)
            ))
          } else{
            boolFalse = FALSE
            tryCatch({
              #sometimes MASS::polr fails to initialize
              IC[m] = IC[m] + stats::AIC(MASS::polr(y[, i] ~ ., data = y[, gam_tmp ==
                                                                           1], method = method))
              boolFalse <- TRUE
            }, error = function(e) {

            }, finally = {

            })

            if (!boolFalse) {
              IC[m] = IC[m] + stats::AIC(MASS::polr(
                y[, i] ~ .,
                data = y[, gam_tmp == 1],
                start = stats::rnorm(nl[i] - 1 + sum(nl[gam_tmp == 1] - 1)),
                method = method
              ))
            }
          }
        }
      }
    }
  }

  mi = which.min(IC)
  ic_best = ic[mi]
  gam = gam_list[, , mi]
  return(list(gam = gam + 0, ic_best = ic_best))
}


#Almost the same as the main function OrdCD but without bootstrapping
OCD = function(y,
               search = "greedy",
               ic = "bic",
               edge_list = NULL,
               link = "probit",
               G = NULL,
               nstart = 1,
               verbose = FALSE,
               maxit = 50) {
  if (search == "exhaust") {
    G = oBN_exhaust(y, G, ic, link)
  } else{
    G = oBN_greedy_wrap(y, ic, edge_list, link, nstart, verbose, maxit)
  }
  return(G)
}

#' @title Causal Discovery for Ordinal Categorical Data
#'
#' @description Estimate a causal directed acyclic graph (DAG) for ordinal categorical data with greedy or exhaustive search.
#'
#' @param y a data frame with each column being an ordinal categorical variable, which must be a factor.
#' @param search the search method used to find the best-scored DAG. The default search method is "greedy". When the number of nodes is less than 4, "exhaust" search is available.
#' @param ic the information criterion (AIC or BIC) used to score DAGs. The default is "bic".
#' @param edge_list an edge list of a CPDAG, which may contain both directed and undirected edges. This option can significantly speed up the algorithm by restricting the search space to DAGs that are Markov equivalent to the input CPDAG. Such CPDAG may be obtained by e.g., the PC algorithm; see the example below.
#' @param link the link function for ordinal regression. The default is "probit". Other choices are "logistic", "loglog", "cloglog", and "cauchit".
#' @param G a list of DAG adjacency matrices that users want to restrict their search on for the "exhaust" search. The default is "NULL" meaning no restriction imposed on the search.
#' @param nstart number of random graph initializations for the "greedy" search.
#' @param verbose if TRUE, messages are printed during the run of the greedy search algorithm.
#' @param maxit the maximum number of iterations for the greedy search algorithm. The default is 50. When the maximum number of iteration is achieved, a warning message will be generated to caution the user that the algorithm has not converged.
#' @param boot the number of bootstrap samples. Default is no bootstrapping.
#' @return A list with "boot" elements. Each element is a list with two elements, gam and ic_best. gam is an estimated DAG adjacency matrix whose (i,j)th entry is 1 if j->i is present in the graph and 0 otherwise. ic_best is the corresponding information criterion value.
#'
#' @export
#'
#' @examples
#' set.seed(2020)
#' n=1000 #sample size
#' q=3 #number of nodes
#' y = u = matrix(0,n,q)
#' u[,1] = 4*rnorm(n)
#' y[,1] = (u[,1]>1) + (u[,1]>2)
#' for (j in 2:q){
#'   u[,j] = 2*y[,j-1] + rnorm(n)
#'   y[,j]=(u[,j]>1) + (u[,j]>2)
#' }
#' A=matrix(0,q,q) #true DAG adjacency matrix
#' A[2,1]=A[3,2]=1
#' y=as.data.frame(y)
#' for (j in 1:q){
#'   y[,j]=as.factor(y[,j])
#' }
#'
#' time=proc.time()
#' G=OrdCD(y) #estimated DAG adjacency matrix
#' time=proc.time() - time
#' print(A) #display the true adjacency
#' print(G) #display the estimated adjacency
#' print(time[3]) #elapsed time
#'
#' time2=proc.time()
#' colnames(y)=1:ncol(y)
#' PC=bnlearn::pc.stable(y,test="mi-sh",alpha=0.01)
#' edge_list=matrix(as.numeric(PC$arcs),ncol=2)
#' G2=OrdCD(y, edge_list=edge_list) #estimated DAG with an input CPDAG from PC algorithm
#' time2=proc.time()-time2
#' print(G2) #display the estimated adjacency
#' print(time2[3]) #elapsed time
#'
#'\dontrun{
#' time3=proc.time()
#' G3=OrdCD(y,boot=10) #estimated DAG adjacency matrix with bootstrapping
#' time3=proc.time() - time3
#' #average over bootstrap samples
#' G_avg = matrix(colMeans(matrix(unlist(G3),nrow=q^2+1,byrow=TRUE))[1:(q^2)],q,q)
#' print(G_avg) #display the estimated adjacency averaged over 10 bootstrap samples
#' print(time3[3]) #elapsed time
#'}

OrdCD = function(y,
               search = "greedy",
               ic = "bic",
               edge_list = NULL,
               link = "probit",
               G = NULL,
               nstart = 1,
               verbose = FALSE,
               maxit = 50,
               boot = NULL) {
  if (is.null(boot)){
    G = OCD(y, search, ic, edge_list, link, G, nstart, verbose, maxit)
  }else{
    G = vector("list",boot)
    for (b in 1:boot){
      G[[b]] = OCD(y[sample(nrow(y),replace = TRUE),], search, ic, edge_list, link, G, nstart, verbose, maxit)
    }
  }
  return(G)
}

