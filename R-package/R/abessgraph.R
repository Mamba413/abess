# #' Title
# #'
# #' @inheritParams abess.default
# #'
# #' @return
# #' @export
# #'
# #' @examples
# #' library(abess)
# #' p <- 16
# #' n <- 1e4
# #' train <- generate.bmn.freq.data(n, p, type = 10, graph.seed = 1, seed = 1, beta = 0.4)
# #' valid <- generate.bmn.freq.data(n, p, type = 10, graph.seed = 1, seed = 2, beta = 0.4)
# #' pool_data <- rbind(train[["data"]], valid[["data"]])
# #' res <- abessbmn(pool_data[, -1], weight = pool_data[, 1], support.size = 27:37, tune.type = "cv", nfolds = 2, foldid = rep(1:2, each = nrow(pool_data) / 2), graph.threshold = 0.2)
# #' all((res[["omega"]][, , which.min(res[["tune.value"]])] != 0) == (train[["theta"]] != 0))
# #' 
# #' res <- abessbmn(pool_data[, -1], weight = pool_data[, 1], tune.type = "gic", support.size = 32, c.max = round(2 * p / 3))
# #' all((res[["omega"]][, , 1] != 0) == (train[["theta"]] != 0))
# #' 
# abessbmn <- function(x,
#                      tune.type = c("gic", "bic", "cv", "ebic", "aic"),
#                      weight = NULL,
#                      c.max = 5,
#                      support.size = NULL,
#                      gs.range = NULL,
#                      always.include = NULL,
#                      splicing.type = 0,
#                      max.splicing.iter = 20,
#                      warm.start = TRUE,
#                      nfolds = 5,
#                      foldid = NULL,
#                      newton.thresh = 1e-6,
#                      max.newton.iter = 500,
#                      ic.scale = 1.0,
#                      num.threads = 0,
#                      seed = 1,
#                      graph.threshold = 0.0, 
#                      ...)
# {
#   early.stop <- FALSE
#   screening.num <- NULL
#   important.search <- NULL
#   lambda <- 0
#   group.index <- NULL
  
#   nobs <- nrow(x)
#   nvars <- ncol(x)
#   maximum_degree <- nvars * (nvars - 1) / 2
#   y <- matrix(0, nrow = nobs, ncol = 1)
  
#   if (is.null(weight)) {
#     weight <- rep(1, nobs)
#   }
  
#   if (is.null(support.size)) {
#     support_size <- seq.int(0, maximum_degree, length.out = maximum_degree)
#   } else {
#     support_size <- support.size
#   }
  
#   tune.path <- c("sequence", "gsection")
#   tune.path <- tune.path[1]
#   # tune_path <- match.arg(tune.path)
#   tune_path <- ifelse(tune.path == "sequence", 1, 2)
#   tune_path <- as.integer(tune_path)
  
#   tune.type <- match.arg(tune.type)
#   ic_type <- switch(tune.type,
#                     "aic" = 1,
#                     "bic" = 2,
#                     "gic" = 3,
#                     "ebic" = 4,
#                     "cv" = 1
#   )
#   is_cv <- ifelse(tune.type == "cv", TRUE, FALSE)
#   if (is_cv) {
#     stopifnot(is.numeric(nfolds) & nfolds >= 2)
#     check_integer_warning(
#       nfolds,
#       "nfolds should be an integer value. It is coerced to be as.integer(nfolds). "
#     )
#     nfolds <- as.integer(nfolds)
    
#     if (is.null(foldid)) {
#       cv_fold_id <- integer(0)
#     } else {
#       stopifnot(is.vector(foldid))
#       stopifnot(is.numeric(foldid))
#       stopifnot(length(foldid) == nobs)
#       check_integer_warning(
#         foldid,
#         "nfolds should be an integer value. It is coerced to be as.integer(foldid). "
#       )
#       foldid <- as.integer(foldid)
#       cv_fold_id <- foldid
#     }
#   } else {
#     cv_fold_id <- integer(0)
#   }
  
#   ## group variable:
#   group_select <- FALSE
#   if (is.null(group.index)) {
#     g_index <- 1:nvars - 1
#     ngroup <- 1
#     max_group_size <- 1
#     # g_df <- rep(1, nvars)
#   } else {
#     g_index <- group.index
#   }
  
#   # check always included variables:
#   if (is.null(always.include)) {
#     always_include <- integer(0)
#   } else {
#     always_include <- always.include - 1
#   }
  
#   # newton <- c("approx", "exact")
#   # newton <- match.arg(newton)
#   # newton_type <- switch(newton,
#   #                       "exact" = 0,
#   #                       "approx" = 1,
#   #                       "auto" = 2
#   # )
#   # approximate_newton <- ifelse(newton_type == 1, TRUE, FALSE)
#   approximate_newton <- FALSE
  
#   result <- abessCpp2(
#     x = x,
#     y = y,
#     n = nobs,
#     p = nvars,
#     data_type = as.integer(1),
#     weight = as.double(weight),
#     sigma = matrix(-1),
#     is_normal = FALSE,
#     algorithm_type = as.integer(6),
#     model_type = as.integer(8),
#     max_iter = as.integer(max.splicing.iter),
#     exchange_num = as.integer(c.max),
#     path_type = tune_path,
#     is_warm_start = warm.start,
#     ic_type = as.integer(ic_type),
#     ic_coef = as.double(ic.scale),
#     is_cv = is_cv,
#     Kfold = nfolds,
#     status = c(0),
#     sequence = support_size,
#     lambda_seq = lambda,
#     s_min = as.integer(0),
#     s_max = as.integer(0),
#     K_max = as.integer(0),
#     epsilon = as.double(0.0001),
#     lambda_min = as.double(0),
#     lambda_max = as.double(0),
#     nlambda = as.integer(0),
#     is_screening = FALSE,
#     screening_size = -1,
#     powell_path = as.integer(1),
#     g_index = g_index,
#     always_select = always_include,
#     tau = 0,
#     primary_model_fit_max_iter = as.integer(max.newton.iter),
#     primary_model_fit_epsilon = as.double(newton.thresh),
#     early_stop = FALSE,
#     approximate_Newton = approximate_newton,
#     thread = num.threads,
#     covariance_update = FALSE,
#     sparse_matrix = FALSE,
#     splicing_type = as.integer(splicing.type),
#     sub_search = as.integer(0),
#     cv_fold_id = cv_fold_id
#   )
  
#   omega <- lapply(result[["beta_all"]], recovery_adjacent_matrix, p = nvars)
#   omega <- simplify2array(omega)
  
#   if (is_cv) {
#     names(result)[which(names(result) == "test_loss_all")] <- "tune.value"
#     result[["ic_all"]] <- NULL
#   } else {
#     names(result)[which(names(result) == "ic_all")] <- "tune.value"
#     result[["test_loss_all"]] <- NULL
#   }
#   result[["tune.value"]] <- as.vector(result[["tune.value"]])
  
#   names(result)[which(names(result) == "train_loss_all")] <- "pseudo.loglik"
#   result[["pseudo.loglik"]] <- as.vector(result[["pseudo.loglik"]])
  
#   optimal_omega <- omega[, , which.min(result[["tune.value"]])]
#   if (graph.threshold > 0.0 && length(support_size) > 1) {
#     optimal_omega <- thres_bmn_est(optimal_omega, graph.threshold)
#   }
  
#   res_out <- list(
#     omega = omega,
#     support.size = support_size,
#     pseudo.loglik = result[["pseudo.loglik"]],
#     tune.value = result[["tune.value"]],
#     nobs = nobs,
#     nvars = nvars,
#     tune.type = tune.type, 
#     optimal.omega = optimal_omega
#   )
#   class(res_out) <- "abessbmn"
  
#   return(res_out)
# }

# recovery_adjacent_matrix <- function(x, p) {
#   zero_mat <- matrix(data = 0, nrow = p, ncol = p)
#   # zero_mat[lower.tri(zero_mat)] <- x
#   # zero_mat <- zero_mat + t(zero_mat)
#   # diag(zero_mat) <- diag(zero_mat) / 2
#   i <- 1
#   j <- 1
#   for (k in 1:as.integer(p * (p - 1) / 2)) {
#     if (i == j) {
#       i <- 1
#       j <- j + 1
#     }
#     zero_mat[j, i] <- zero_mat[i, j] <- x[k]
#     i <- i + 1
#   }
#   zero_mat
# }


#' 
#' Nodewise logistic regression for inverse Ising problem
#' 
#' @inheritParams abess.default
#'
#' @param max.support.size 
#'
#' @return
#' @export
#'
#' @examples
#' p <- 16
#' n <- 1e3
#' library(abess)
#' train <- generate.bmn.data(n, p, type = 10, graph.seed = 1, seed = 1, beta = 0.4)
#' res <- nodewise_L0(train[["data"]], train[["weight"]], tune.type = "gic", 
#'                    max.support.size = rep(4, p), support.size = rep(4, p))
#' all((res[[1]] != 0) == (train[["theta"]] != 0))
#' 
#' ## use cross validation to nodewisely select support.size
#' valid <- generate.bmn.data(n, p, type = 10, graph.seed = 1, seed = 10000, beta = 0.4)
#' all(train[["theta"]] == valid[["theta"]])
#' x <- rbind(train[["data"]], valid[["data"]])
#' sample_weight <- c(train[["weight"]], valid[["weight"]])
#' fold_id <- c(rep(1, length(train[["weight"]])), rep(2, length(valid[["weight"]])))
#' res <- nodewise_L0(x, sample_weight, tune.type = "cv", foldid = fold_id, graph.threshold = 0.2)
#' all((res[[1]] != 0) == (train[["theta"]] != 0))
#' 
#' ## use IC to nodewisely select support.size (without post-thresholding)
#' res <- nodewise_L0(x, sample_weight, tune.type = "gic", ic.scale = 2)
#' all((res[[1]] != 0) == (train[["theta"]] != 0))
#' 
#' res <- nodewise_L0(x, sample_weight, tune.type = "gic", ic.scale = 1, graph.threshold = 0.2)
#' all((res[[1]] != 0) == (train[["theta"]] != 0))
#' 
#' train <- generate.bmn.data(10000 * n, p, type = 8, alpha = 0.9, beta = 0.9, degree = 3, seed = 1)
#' system.time(res <- nodewise_L0(train[["data"]], train[["weight"]], tune.type = "gic", ic.scale = 1, graph.threshold = 0.2, magnetic = TRUE))
#' mean(res[[1]] - train[["theta"]])
#' mean(diag(res[[1]]) - diag(train[["theta"]]))
#' diag(res[[1]]) <- 0
#' all((res[[1]] != 0) == (train[["theta"]] != 0))
#' \dontrun{
#' ## simple test suggest that sparsity is not useful
#' system.time(res <- nodewise_L0(train[["data"]], train[["weight"]], tune.type = "gic", ic.scale = 1, graph.threshold = 0.2, sparse = TRUE))
#' mean(res[[1]] - train[["theta"]])
#' mean(diag(res[[1]]) - diag(train[["theta"]]))
#' diag(res[[1]]) <- 0
#' all((res[[1]] != 0) == (train[["theta"]] != 0))
#' }
#' 
#' train <- generate.bmn.data(9000, 48, type = 8, alpha = 0.9, beta = 0.9, degree = 3, method = "gibbs", seed = 1)
#' system.time(res <- nodewise_L0(train[["data"]], train[["weight"]], tune.type = "gic", ic.scale = 1, graph.threshold = 0.2))
#' mean(res[[1]] - train[["theta"]])
#' mean(diag(res[[1]]) - diag(train[["theta"]]))
#' diag(res[[1]]) <- 0
#' all((res[[1]] != 0) == (train[["theta"]] != 0))
#' 
nodewise_L0 <- function(x,
                        weight = NULL, 
                        max.support.size = NULL,
                        tune.type = "cv",
                        c.max = NULL, 
                        foldid = NULL, 
                        support.size = NULL,
                        ic.scale = 1, 
                        graph.threshold = 0.0, 
                        newton = c("approx", "exact"), 
                        max.newton.iter = 100, 
                        sparse = FALSE, 
                        magnetic = FALSE, 
                        ...) 
{
  p <- ncol(x)
  n <- nrow(x)
  if (is.null(max.support.size)) {
    max.support.size <- min(c(p - 2, 20))
    max.support.size <- rep(max.support.size, p)
  }
  if (is.null(foldid) && tune.type == "cv") {
    foldid <- c()
    nfolds <- 2
  } else if (tune.type == "cv") {
    nfolds <- length(unique(foldid))
  } else {
    nfolds <- 1
  }
  
  ## simple test suggest that sparsity is not useful:
  # sparse <- FALSE
  if (sparse) {
    non_zero_ind <- which(x == 1, arr.ind = TRUE)
    x <- sparseMatrix(i = non_zero_ind[, 1], j = non_zero_ind[, 2], x = 1, dims = c(n, p))
  }
  
  newton_method <- match.arg(newton)
  
  if (is.null(c.max)) {
    c_max_value <- round(max(max.support.size) / 2)
  } else {
    c_max_value <- c.max
  }
  
  theta <- matrix(0, p, p)
  tune_value_list <- list()
  for (node in 1:p) {
    no_fit_flag <- FALSE
    min_nobs <- min(table(x[, node]))
    if (min_nobs > 1) {
      init_active_set <- which(theta[-node, node] != 0)
      model_node <-
        abess::abess(
          x = x[, -node],
          y = x[, node],
          weight = weight,
          family = "binomial",
          tune.path = "sequence",
          support.size = 0:max.support.size[node],
          init.active.set = init_active_set, 
          tune.type = tune.type, 
          ic.scale = ic.scale, 
          nfolds = nfolds,
          foldid = foldid,
          c.max = c_max_value,
          max.splicing.iter = 100,
          newton = newton_method,
          newton.thresh = log(p) * log(log(n)) / n,
          max.newton.iter = as.integer(max.newton.iter),
          num.threads = nfolds, 
          seed = 1
        )
      tune_value_list[[node]] <- model_node[["tune.value"]]
      if (is.null(support.size)) {
        est_theta_node <- as.vector(coef(model_node, support.size = model_node[["best.size"]]))[-1]
      } else {
        est_theta_node <- as.vector(coef(model_node, support.size = support.size[node]))[-1]
      }
      theta[node, -node] <- (est_theta_node / 2)
      if (magnetic) {
        theta[node, node] <- as.vector(coef(model_node, support.size = support.size[node], sparse = FALSE))[1] / 2
      }
      if (anyNA(theta[node, ])) {
        theta[node, ] <- 0.0
        no_fit_flag <- TRUE
      }
    } else {
      no_fit_flag <- TRUE
    }
    # additional handling not fitted case:
    if (no_fit_flag) {
      if (magnetic & min_nobs > 1) {
        theta[node, node] <- -0.5 * log(1 / mean(x[, node] == 1, na.rm=TRUE) - 1)
      }
    }
  }
  theta <- (t(theta) + theta) / 2
  if (sparse) {
    # For {0, 1} dataset, the pseudo likelihood no need to divide 2. 
    # Here, 2 is multiply to compensate the divide 2 operation before. 
    theta <- 2 * theta 
    diag(theta) <- 0.25 * (colSums(theta) + diag(theta))
    theta <- theta / 4
    diag(theta) <- diag(theta) * 4
  }
  
  if (graph.threshold > 0.0 && is.null(support.size)) {
    theta <- thres_bmn_est(theta, graph.threshold)
  }
  
  names(tune_value_list) <- NULL
  res_out <- list(
    omega = theta, 
    tune.value = tune_value_list
  )
  class(res_out) <- "abessbmn"
  res_out
}


thres_bmn_est <- function(theta, thres) {
  theta_diag <- diag(theta)
  if (thres > 0) {
    theta[abs(theta) <= thres] <- 0
  } else if (thres < 0) {
    theta_vec <- as.vector(theta)
    ## TODO: use finite mixture model to cluster
  } 
  diag(theta) <- theta_diag
  theta
}
