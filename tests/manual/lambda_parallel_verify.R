# Stress harness for parallel_lambda branch.
# NOT in tests/testthat/ to avoid CRAN check time.
#
# Usage:
#   # Snapshot the sequential baseline (run once on Step 1 commit):
#   VERIFY_MODE=snapshot BASELINE_TAG=clean Rscript tests/manual/lambda_parallel_verify.R
#
#   # Verify after each subsequent step:
#   VERIFY_MODE=verify BASELINE_TAG=clean Rscript tests/manual/lambda_parallel_verify.R

suppressPackageStartupMessages({
    library(PAGFL)
})

configs <- list(
    list(name = "pls_small",    N = 20,  T = 50,  p = 2, K = 3, n_lambda = 5,  method = "PLS",  fuse = FALSE),
    list(name = "pls_typical",  N = 50,  T = 60,  p = 2, K = 3, n_lambda = 8,  method = "PLS",  fuse = FALSE),
    list(name = "pls_large",    N = 80,  T = 80,  p = 3, K = 4, n_lambda = 10, method = "PLS",  fuse = FALSE),
    list(name = "pls_edge1",    N = 30,  T = 50,  p = 2, K = 2, n_lambda = 1,  method = "PLS",  fuse = FALSE),  # inner-only path
    list(name = "pls_edge3",    N = 30,  T = 50,  p = 2, K = 2, n_lambda = 3,  method = "PLS",  fuse = FALSE),  # at threshold
    list(name = "pls_bias",     N = 40,  T = 60,  p = 2, K = 3, n_lambda = 6,  method = "PLS",  fuse = FALSE, bias_correc = TRUE),
    list(name = "pgmm_small",   N = 20,  T = 50,  p = 2, K = 3, n_lambda = 5,  method = "PGMM", fuse = FALSE),
    list(name = "pgmm_bias",    N = 50,  T = 40,  p = 2, K = 3, n_lambda = 6,  method = "PGMM", fuse = FALSE, bias_correc = TRUE),
    list(name = "fuse_small",   N = 30,  T = 60,  p = 2, K = 3, n_lambda = 5,  fuse = TRUE),
    list(name = "fuse_typical", N = 40,  T = 80,  p = 2, K = 3, n_lambda = 6,  fuse = TRUE)
)

run_one <- function(cfg, parallel) {
    set.seed(1)
    lam <- exp(seq(log(1e-3), log(1e-1), length.out = cfg$n_lambda))
    bias_correc <- isTRUE(cfg$bias_correc)
    if (isTRUE(cfg$fuse)) {
        dgp <- sim_tv_DGP(N = cfg$N, n_periods = cfg$T, p = cfg$p, n_groups = cfg$K)
        fuse_time(y ~ ., data = dgp$data, n_periods = cfg$T,
                  lambda = lam, parallel = parallel, verbose = FALSE)
    } else if (identical(cfg$method, "PGMM")) {
        dgp <- sim_DGP(N = cfg$N, n_periods = cfg$T, p = cfg$p,
                       n_groups = cfg$K, q = cfg$p + 1L)
        pagfl(y ~ ., data = dgp$data, n_periods = cfg$T,
              lambda = lam, method = "PGMM", Z = dgp$Z,
              bias_correc = bias_correc,
              parallel = parallel, verbose = FALSE)
    } else {
        dgp <- sim_DGP(N = cfg$N, n_periods = cfg$T, p = cfg$p, n_groups = cfg$K)
        pagfl(y ~ ., data = dgp$data, n_periods = cfg$T,
              lambda = lam, method = cfg$method,
              bias_correc = bias_correc,
              parallel = parallel, verbose = FALSE)
    }
}

# Strip fields that vary harmlessly between runs (call object captures the literal call)
canonicalize <- function(res) {
    res$call <- NULL
    res
}

compare <- function(a, b, tol = 1e-10) {
    isTRUE(all.equal(canonicalize(a), canonicalize(b), tolerance = tol,
                     check.attributes = FALSE))
}

mode    <- Sys.getenv("VERIFY_MODE", "snapshot")
tag     <- Sys.getenv("BASELINE_TAG", "clean")
outfile <- file.path("tests", "manual", paste0("baseline_", tag, ".rds"))

if (mode == "snapshot") {
    cat("Snapshotting baseline (parallel=FALSE) ...\n")
    results <- lapply(configs, function(cfg) {
        cat(sprintf("  %s ...\n", cfg$name))
        run_one(cfg, parallel = FALSE)
    })
    names(results) <- vapply(configs, function(x) x$name, character(1))
    saveRDS(results, outfile)
    cat("Saved baseline:", outfile, "\n")
} else if (mode == "verify") {
    cat("Verifying against baseline:", outfile, "\n")
    baseline <- readRDS(outfile)
    n_pass <- 0L
    n_fail <- 0L
    for (cfg in configs) {
        seq_res <- run_one(cfg, parallel = FALSE)
        par_res <- run_one(cfg, parallel = TRUE)
        seq_ok <- compare(seq_res, baseline[[cfg$name]])
        par_ok <- compare(par_res, baseline[[cfg$name]])
        sp_ok  <- compare(seq_res, par_res)
        all_ok <- seq_ok && par_ok && sp_ok
        if (all_ok) n_pass <- n_pass + 1L else n_fail <- n_fail + 1L
        cat(sprintf("[%-13s] seq~base=%-5s  par~base=%-5s  par~seq=%-5s  %s\n",
                    cfg$name, seq_ok, par_ok, sp_ok,
                    if (all_ok) "OK" else "FAIL"))
    }
    cat(sprintf("\n%d / %d configurations pass\n", n_pass, n_pass + n_fail))
    if (n_fail > 0L) quit(status = 1L)
} else {
    stop("VERIFY_MODE must be 'snapshot' or 'verify'")
}
