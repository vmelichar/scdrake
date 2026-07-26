cli::cli_h1("{.file tests/testthat/test-cc-scoring.R}")

## -- Regression test for the cell-cycle scoring nbin retry (sce_cc_fn). A permissive gene filter
## -- leaves a tie-pile of near-zero genes; Seurat::CellCycleScoring() (via AddModuleScore) then
## -- fails to bin all genes at the default nbin = 24. sce_cc_fn() must retry with smaller nbin and
## -- still produce phases instead of silently setting everything to NA.

test_that("sce_cc_fn() recovers when CellCycleScoring fails at the default nbin = 24", {
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("scuttle")
  suppressMessages(library(SingleCellExperiment))

  set.seed(42)
  ng <- 300L
  nc <- 200L
  ## -- 21/300 = 7% of genes all-zero -> tie-pile at the minimum mean expression, which trips the
  ## -- nbin = 24 binning; a lower nbin still bins fine.
  rates <- c(rep(0, 21), runif(ng - 21, 0.1, 5))
  counts <- matrix(rpois(ng * nc, lambda = rep(rates, times = nc)), nrow = ng)
  rownames(counts) <- paste0("ENSG", seq_len(ng))
  colnames(counts) <- paste0("cell", seq_len(nc))
  sce <- scuttle::logNormCounts(SingleCellExperiment(assays = list(counts = counts)))
  sce$batch <- "b1" # second colData column so sce_add_colData() keeps a DataFrame, not a vector
  cc_genes <- data.frame(
    phase = c(rep("S", 20), rep("G2M", 20)),
    ENSEMBL = rownames(sce)[22:61]
  )

  res <- suppressMessages(
    scdrake::sce_cc_fn(sce, cc_genes, data = NULL, cc_eval = TRUE, spatial = FALSE)
  )

  testthat::expect_false(all(is.na(res$phase)))
  testthat::expect_false(isTRUE(S4Vectors::metadata(res)$cc_phase_failed))
})
