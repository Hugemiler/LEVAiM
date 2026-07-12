#####
# Time-varying exposure derivation
#####

#' Derive analysis features from a time-varying exposure
#'
#' Adds subject-level metadata columns derived from a longitudinal exposure
#' column. The derived columns are repeated on each sample row for the subject
#' and recorded in an exposure-derivation audit table.
#'
#' @param ds A `LongitudinalDataset`.
#' @param exposure Metadata column containing the time-varying exposure.
#' @param summaries Named list of exposure summary specifications created with
#'   `exposure_window()`, `exposure_ever()`, `exposure_last_observed()`, or
#'   `exposure_transition()`.
#' @param subject Optional subject column. Defaults to `ds$subject_col`.
#' @param time Optional time column. Defaults to `ds$time_col`.
#'
#' @return A `LongitudinalDataset` with additional metadata columns.
#' @export
derive_time_varying_features <- function(
    ds,
    exposure,
    summaries,
    subject = NULL,
    time = NULL
) {
  if (!inherits(ds, "LongitudinalDataset")) {
    cli::cli_abort("`ds` must be a `LongitudinalDataset` object.")
  }

  subject <- subject %||% ds$subject_col
  time <- time %||% ds$time_col

  metadata <- ds$metadata

  required <- c(exposure, subject, time)
  missing <- setdiff(required, names(metadata))
  if (length(missing) > 0L) {
    cli::cli_abort("Missing metadata column(s): {paste(missing, collapse = ', ')}")
  }

  if (!is.list(summaries) || length(summaries) == 0L) {
    cli::cli_abort("`summaries` must be a non-empty named list of exposure summary specs.")
  }

  if (is.null(names(summaries)) || any(names(summaries) == "")) {
    cli::cli_abort("`summaries` must be a named list.")
  }

  bad <- !vapply(
    summaries,
    inherits,
    logical(1),
    what = "levaim_exposure_summary_spec"
  )
  if (any(bad)) {
    cli::cli_abort("All `summaries` entries must be created with exposure summary constructors.")
  }

  duplicated_names <- intersect(names(summaries), names(metadata))
  if (length(duplicated_names) > 0L) {
    cli::cli_abort(
      "Derived exposure column(s) already exist in metadata: {paste(duplicated_names, collapse = ', ')}"
    )
  }

  subject_values <- as.character(metadata[[subject]])
  subjects <- unique(subject_values)
  derived <- vector("list", length(summaries))
  names(derived) <- names(summaries)

  for (summary_name in names(summaries)) {
    spec <- summaries[[summary_name]]
    subject_summary <- vapply(
      subjects,
      function(subject_id) {
        rows <- metadata[subject_values == subject_id, , drop = FALSE]
        rows <- rows[order(rows[[time]]), , drop = FALSE]
        summarize_exposure(rows, exposure = exposure, time = time, spec = spec)
      },
      character(1)
    )
    names(subject_summary) <- subjects
    derived[[summary_name]] <- unname(subject_summary[subject_values])
  }

  for (summary_name in names(derived)) {
    metadata[[summary_name]] <- derived[[summary_name]]
  }

  audit <- exposure_derivation_audit(
    summaries = summaries,
    exposure = exposure,
    subject = subject,
    time = time
  )

  previous_audit <- attr(ds, "exposure_derivations", exact = TRUE)
  if (is.data.frame(previous_audit)) {
    audit <- rbind(previous_audit, audit)
  }

  ds$metadata <- metadata
  attr(ds, "exposure_derivations") <- audit
  ds
}


#' Summarize exposure values within a time window
#'
#' @param start Window start, inclusive.
#' @param end Window end, inclusive.
#' @param method Summary method. `"modal"` returns the most frequent observed
#'   value; `"first"` and `"last"` use chronological order; `"any"` returns
#'   whether any exposure value, or any matching `value`, occurs in the window.
#' @param value Optional value used by `method = "any"`.
#'
#' @return An exposure summary specification.
#' @export
exposure_window <- function(
    start,
    end,
    method = c("modal", "first", "last", "any"),
    value = NULL
) {
  method <- match.arg(method)

  if (!is.numeric(start) || length(start) != 1L || !is.finite(start)) {
    cli::cli_abort("`start` must be a finite numeric scalar.")
  }
  if (!is.numeric(end) || length(end) != 1L || !is.finite(end)) {
    cli::cli_abort("`end` must be a finite numeric scalar.")
  }
  if (end < start) {
    cli::cli_abort("`end` must be greater than or equal to `start`.")
  }

  exposure_summary_spec(
    type = "window",
    start = start,
    end = end,
    method = method,
    value = value
  )
}


#' Derive whether an exposure value ever occurs
#'
#' @param value Exposure value to query.
#'
#' @return An exposure summary specification.
#' @export
exposure_ever <- function(value) {
  exposure_summary_spec(type = "ever", value = value)
}


#' Derive the last observed exposure value
#'
#' @return An exposure summary specification.
#' @export
exposure_last_observed <- function() {
  exposure_summary_spec(type = "last_observed")
}


#' Derive the observed exposure transition path
#'
#' @return An exposure summary specification.
#' @export
exposure_transition <- function() {
  exposure_summary_spec(type = "transition")
}


#' Return the exposure derivation audit table
#'
#' @param ds A `LongitudinalDataset`.
#'
#' @return A data.frame.
#' @export
exposure_derivations <- function(ds) {
  audit <- attr(ds, "exposure_derivations", exact = TRUE)
  if (is.null(audit)) {
    return(data.frame())
  }
  audit
}


#' @keywords internal
exposure_summary_spec <- function(type, ...) {
  structure(
    c(list(type = type), list(...)),
    class = "levaim_exposure_summary_spec"
  )
}


#' @keywords internal
summarize_exposure <- function(rows, exposure, time, spec) {
  values <- as.character(rows[[exposure]])
  values[!nzchar(values)] <- NA_character_

  if (spec$type == "window") {
    keep <- rows[[time]] >= spec$start & rows[[time]] <= spec$end
    values <- values[keep]
    values <- values[!is.na(values)]

    if (length(values) == 0L) {
      return(NA_character_)
    }

    if (spec$method == "modal") {
      counts <- sort(table(values), decreasing = TRUE)
      return(names(counts)[[1L]])
    }
    if (spec$method == "first") {
      return(values[[1L]])
    }
    if (spec$method == "last") {
      return(values[[length(values)]])
    }
    if (spec$method == "any") {
      if (is.null(spec$value)) {
        return(as.character(any(!is.na(values))))
      }
      return(as.character(any(values == spec$value)))
    }
  }

  observed <- values[!is.na(values)]

  if (spec$type == "ever") {
    return(as.character(any(observed == spec$value)))
  }

  if (spec$type == "last_observed") {
    if (length(observed) == 0L) {
      return(NA_character_)
    }
    return(observed[[length(observed)]])
  }

  if (spec$type == "transition") {
    if (length(observed) == 0L) {
      return(NA_character_)
    }
    transitions <- observed[c(TRUE, observed[-1L] != observed[-length(observed)])]
    return(paste(transitions, collapse = "->"))
  }

  cli::cli_abort("Unknown exposure summary type: {.val {spec$type}}.")
}


#' @keywords internal
exposure_derivation_audit <- function(summaries, exposure, subject, time) {
  rows <- lapply(names(summaries), function(name) {
    spec <- summaries[[name]]
    data.frame(
      derived_column = name,
      exposure = exposure,
      subject = subject,
      time = time,
      type = spec$type,
      method = spec$method %||% NA_character_,
      value = spec$value %||% NA_character_,
      window_start = spec$start %||% NA_real_,
      window_end = spec$end %||% NA_real_,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
