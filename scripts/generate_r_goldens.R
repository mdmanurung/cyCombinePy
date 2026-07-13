#!/usr/bin/env Rscript

# Generate small R cyCombine golden fixtures for Python parity tests.
# This script must not download packages or data from the network.

required_packages <- c("cyCombine", "jsonlite", "digest")
missing_packages <- required_packages[!vapply(
  required_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required R package(s): ",
      paste(missing_packages, collapse = ", "),
      ". Install them locally before generating cyCombine goldens; ",
      "this script does not download packages from the network."
    ),
    call. = FALSE
  )
}

seed <- 13L
set.seed(seed)

script_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
  error = function(...) file.path("scripts", "generate_r_goldens.R")
)
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
cycombine_version <- as.character(utils::packageVersion("cyCombine"))
golden_dir <- file.path(
  repo_root,
  "tests",
  "golden",
  paste0("cycombine_r_", cycombine_version)
)

paths <- list(
  manifest = file.path(golden_dir, "manifest.json"),
  input = file.path(golden_dir, "input.csv"),
  normalize_scale = file.path(golden_dir, "normalize_scale.npz"),
  normalize_rank = file.path(golden_dir, "normalize_rank.npz"),
  emd = file.path(golden_dir, "emd.csv"),
  mad = file.path(golden_dir, "mad.csv"),
  corrected_fixed_labels = file.path(golden_dir, "corrected_fixed_labels.npz")
)

dir.create(golden_dir, recursive = TRUE, showWarnings = FALSE)

# TODO: Build a fixed small input data frame and write paths$input.
# TODO: Run cyCombine normalization variants and write paths$normalize_scale and
#       paths$normalize_rank in NumPy-compatible .npz form.
# TODO: Run EMD and MAD evaluation on the fixed data and write paths$emd and
#       paths$mad.
# TODO: Run fixed-label correction and write paths$corrected_fixed_labels.
# TODO: Compute SHA256 hashes for all generated fixture files and write
#       paths$manifest using jsonlite::write_json().

hash_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE)
}

manifest_template <- list(
  schema_version = 1L,
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  cycombine = list(
    version = cycombine_version,
    commit = Sys.getenv("CYCOMBINE_R_COMMIT", unset = "unknown")
  ),
  package_versions = stats::setNames(
    lapply(required_packages, function(package) {
      as.character(utils::packageVersion(package))
    }),
    required_packages
  ),
  random_seed = seed,
  provenance = list(
    generator = "scripts/generate_r_goldens.R",
    network = FALSE
  ),
  sha256 = list(
    input.csv = "<sha256>",
    normalize_scale.npz = "<sha256>",
    normalize_rank.npz = "<sha256>",
    emd.csv = "<sha256>",
    mad.csv = "<sha256>",
    corrected_fixed_labels.npz = "<sha256>"
  )
)

message("Golden fixture output directory: ", golden_dir)
message("Manifest template path: ", paths$manifest)
message("Generation TODOs remain; no fixture files were written.")
