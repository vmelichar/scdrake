cli::cli_h1("{.file tests/testthat/test-dummy-plots.R}")

## -- Unit tests for the dummy-plot fallbacks (create_dummy_plot) and the code paths that
## -- produce or save them. These exercise the edge cases handled in R/visualization.R and
## -- R/cell_annotation.R without running the full drake pipeline.
## -- Non-scdrake calls are fully qualified so the file also runs via
## -- testthat::test_file() after devtools::load_all().

test_that("create_dummy_plot() returns a ggplot", {
  p <- scdrake::create_dummy_plot("hello")
  testthat::expect_s3_class(p, "ggplot")
})

test_that("plot_clustree() returns a dummy plot when there are fewer than two resolutions", {
  p <- scdrake::plot_clustree(
    cluster_list = list(res_0.1 = c(1, 1, 1)),
    params = 0.1,
    prefix = "res_"
  )
  testthat::expect_s3_class(p, "ggplot")
})

test_that("plot_clustree() returns a dummy plot when every resolution has a single cluster", {
  ## -- Two resolutions, but each is a single cluster -> no tree to draw -> dummy plot.
  ## -- This also guards against the old unique() desync between cluster_list and params.
  p <- scdrake::plot_clustree(
    cluster_list = list(res_0.1 = c(1, 1, 1), res_0.2 = c(1, 1, 1)),
    params = c(0.1, 0.2),
    prefix = "res_"
  )
  testthat::expect_s3_class(p, "ggplot")
})

test_that("plot_clustree() builds a real tree when there are >=2 resolutions with >=2 clusters", {
  testthat::skip_if_not_installed("clustree")
  ## -- The guard must NOT over-trigger: with genuine splits we get a clustree plot, not a dummy.
  p <- scdrake::plot_clustree(
    cluster_list = list(res_0.1 = c(1, 1, 2, 2), res_0.2 = c(1, 2, 3, 4)),
    params = c(0.1, 0.2),
    prefix = "res_"
  )
  testthat::expect_s3_class(p, "ggplot")
})

test_that("cell_annotation_diagnostic_plots_files_fn() saves present plots and skips NA/missing ones", {
  tmp <- withr::local_tempdir()
  delta_file <- as.character(fs::path(tmp, "delta.pdf"))
  mh_file <- as.character(fs::path(tmp, "marker_heatmaps.pdf"))

  ## -- Mimics the tibble produced by cell_annotation_diagnostic_plots_fn(): a delta plot and
  ## -- a (single-cell-type) placeholder marker heatmap to save, plus an NA score heatmap to skip.
  df <- tibble::tibble(
    name = "test_ref",
    score_heatmaps = list(list(NULL)),
    score_heatmaps_out_file = NA_character_,
    delta_distribution_plot = list(scdrake::create_dummy_plot("delta")),
    delta_distribution_plot_out_file = delta_file,
    marker_heatmaps = list(magrittr::set_names(list(scdrake::create_dummy_plot("mh")), "A")),
    marker_heatmaps_out_file = mh_file
  )

  out <- scdrake::cell_annotation_diagnostic_plots_files_fn(df)

  testthat::expect_true(fs::file_exists(delta_file))
  testthat::expect_true(fs::file_exists(mh_file))
  ## -- The NA score-heatmap path is neither written nor returned.
  testthat::expect_false(fs::file_exists(fs::path(tmp, "score_heatmaps.pdf")))
  testthat::expect_setequal(out, c(delta_file, mh_file))
})

test_that("cell_annotation_diagnostic_plots_files_fn() does not error when out-file columns are absent", {
  tmp <- withr::local_tempdir()
  delta_file <- as.character(fs::path(tmp, "delta.pdf"))

  ## -- Only the delta plot columns are present (score_heatmaps_out_file / marker_heatmaps_out_file
  ## -- missing). The %in% names(row) guards added in commit 0fdbfe1 must handle this gracefully.
  df <- tibble::tibble(
    name = "test_ref",
    delta_distribution_plot = list(scdrake::create_dummy_plot("delta")),
    delta_distribution_plot_out_file = delta_file
  )

  ## -- Absent *_out_file columns make the direct `$` accesses in `to_return` warn (benign); the
  ## -- point of the test is that the %in% names(row) guards prevent an error.
  out <- suppressWarnings(scdrake::cell_annotation_diagnostic_plots_files_fn(df))

  testthat::expect_true(fs::file_exists(delta_file))
  testthat::expect_setequal(out, delta_file)
})

test_that("cell_annotation_diagnostic_plots_fn() uses a placeholder marker heatmap for a single cell type", {
  testthat::skip_if_not_installed("SingleR")
  testthat::skip_if_not_installed("scuttle")

  set.seed(1)
  ref <- scuttle::mockSCE(ncells = 40, ngenes = 200)
  ref <- scuttle::logNormCounts(ref)
  ref$label <- rep(c("A", "B"), length.out = ncol(ref))
  test_sce <- scuttle::mockSCE(ncells = 20, ngenes = 200)
  test_sce <- scuttle::logNormCounts(test_sce)

  pred <- SingleR::SingleR(test = test_sce, ref = ref, labels = ref$label)
  ## -- Force a single assigned cell type to hit the placeholder branch.
  pred$labels <- rep("A", nrow(pred))
  testthat::expect_lt(length(unique(pred$labels)), 2)

  ca <- tibble::tibble(
    name = "test_ref",
    cell_annotation = list(pred),
    train_params = list(list(genes = "de")),
    diagnostics_params = list(list(heatmap_n_top_markers = 5))
  )

  res <- scdrake::cell_annotation_diagnostic_plots_fn(
    cell_annotation = ca,
    cell_data = NULL,
    sce = NULL,
    base_out_dir = withr::local_tempdir(),
    do_heatmaps = FALSE
  )

  ## -- A placeholder marker heatmap is produced (a ggplot), the delta plot is built, and the row
  ## -- gets a non-NA output file. `lapply_rows()` unwraps the length-1 marker list to a bare
  ## -- ggplot, so accept either shape here.
  mh <- res$marker_heatmaps[[1]]
  testthat::expect_true(inherits(mh, "ggplot") || (is.list(mh) && inherits(mh[[1]], "ggplot")))
  testthat::expect_s3_class(res$delta_distribution_plot[[1]], "ggplot")
  testthat::expect_false(is.na(res$marker_heatmaps_out_file[[1]]))

  ## -- End-to-end: the placeholder must actually save without error (see the dedicated save_pdf
  ## -- regression test below for the content check).
  scdrake::cell_annotation_diagnostic_plots_files_fn(res)
  testthat::expect_true(fs::file_exists(res$marker_heatmaps_out_file[[1]]))
})

test_that("save_pdf() renders a bare ggplot identically to a one-element list", {
  ## -- Regression guard: lapply_rows() can hand save_pdf() a bare ggplot (the single-cell-type
  ## -- placeholder). Before the fix, save_pdf() iterated the ggplot's internal slots instead of
  ## -- printing it, yielding a broken/blank PDF that differs from the canonical list form.
  p <- scdrake::create_dummy_plot("only one cell type")
  f_bare <- withr::local_tempfile(fileext = ".pdf")
  f_list <- withr::local_tempfile(fileext = ".pdf")

  scdrake::save_pdf(p, output_file = f_bare)
  scdrake::save_pdf(list(p), output_file = f_list)

  testthat::expect_true(fs::file_exists(f_bare))
  ## -- Same single-page render -> identical content length (PDF dates are fixed-width, so size is stable).
  testthat::expect_equal(file.info(f_bare)$size, file.info(f_list)$size)
})
