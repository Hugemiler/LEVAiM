#####
# LEVAiM serialization
#####

#' Save a LEVAiM object
#'
#' Serializes a LEVAiM object in a versioned envelope using R's native RDS
#' format. The envelope records the serialization schema, package version,
#' object class, and creation time.
#'
#' @param x A LEVAiM object.
#' @param file Path to the output `.rds` file.
#' @param ... Additional arguments passed to [saveRDS()].
#'
#' @return Invisibly returns `file`.
#' @export
save_levaim_object <- function(x, file, ...) {
  if (!is_levaim_serializable(x)) {
    cli::cli_abort(
      "`x` must be a LEVAiM fit, results, model frame, cross-validation, dataset, specification, trajectory, or control object."
    )
  }

  envelope <- list(
    format = "LEVAiM_serialized_object",
    schema_version = levaim_serialization_schema(),
    package_version = levaim_package_version(),
    object_class = class(x),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    object = x
  )

  saveRDS(envelope, file = file, ...)

  invisible(file)
}


#' Load a LEVAiM object
#'
#' Reads a LEVAiM object saved with [save_levaim_object()], validates its schema
#' and class, and attaches the envelope metadata as `levaim_serialization`.
#' Legacy raw-object RDS files remain readable as schema version 0.
#'
#' @param file Path to an `.rds` file produced by [save_levaim_object()].
#' @param ... Additional arguments passed to [readRDS()].
#'
#' @return The deserialized LEVAiM object.
#' @export
load_levaim_object <- function(file, ...) {
  stored <- readRDS(file = file, ...)
  if (is_levaim_serialization_envelope(stored)) {
    if (stored$schema_version > levaim_serialization_schema()) {
      cli::cli_abort(c(
        "The file uses unsupported LEVAiM serialization schema {stored$schema_version}.",
        "i" = "This installation supports schema {levaim_serialization_schema()}."
      ))
    }
    metadata <- stored[setdiff(names(stored), "object")]
    x <- stored$object
  } else {
    metadata <- list(
      format = "LEVAiM_legacy_raw_object",
      schema_version = 0L,
      package_version = NA_character_,
      object_class = class(stored),
      created_utc = NA_character_
    )
    x <- stored
  }

  if (!is_levaim_serializable(x)) {
    cli::cli_abort(
      "The file does not contain a recognized LEVAiM object."
    )
  }

  attr(x, "levaim_serialization") <- metadata
  x
}


#' @keywords internal
levaim_serialization_schema <- function() 1L


#' @keywords internal
levaim_package_version <- function() {
  tryCatch(
    as.character(utils::packageVersion("LEVAiM")),
    error = function(e) "unknown"
  )
}


#' @keywords internal
is_levaim_serialization_envelope <- function(x) {
  is.list(x) && identical(x$format, "LEVAiM_serialized_object") &&
    is.numeric(x$schema_version) && length(x$schema_version) == 1L &&
    "object" %in% names(x)
}


#' Check whether an object uses a serializable LEVAiM class
#'
#' @param x An object.
#'
#' @return `TRUE` when `x` has a recognized LEVAiM class.
#' @keywords internal
is_levaim_serializable <- function(x) {
  inherits(
    x,
    c(
      "LongitudinalDataset",
      "levaim_model_spec",
      "levaim_trajectory",
      "levaim_trajectory_control",
      "levaim_model_frame",
      "levaim_fit",
      "levaim_results",
      "levaim_cv",
      "levaim_trajectory_moments",
      "levaim_moment_comparisons"
    )
  )
}
