#!/usr/bin/env Rscript

required_packages <- c(
  "DBI", "RPostgres", "data.table", "WeightIt", "cobalt",
  "ggplot2", "sandwich", "Hmisc", "glmnet"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    ". Install them before rerunning the pipeline."
  )
}

package_versions <- vapply(
  required_packages,
  function(package) as.character(utils::packageVersion(package)),
  character(1)
)

message("R version: ", R.version.string)
message("Required packages available:")
for (package in required_packages) {
  message("  ", package, " ", package_versions[[package]])
}
