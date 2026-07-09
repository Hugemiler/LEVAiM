#####
# 02. Data IO - HUMAnN functional profiles
#####

## PUBLIC functions
# usethis::use_package("data.table")

#' Read HUMAnN functional profiles
#'
#' Reads HUMAnN 3/4-style functional output from a single profile, a merged
#' table, or a vector of single-profile paths.
#'
#' HUMAnN outputs may represent different feature spaces, including UniRefs,
#' KOs, ECs, MetaCyc pathways, or other gene/pathway families. They may also
#' use different units, such as RPK, CPM, relative abundance, or raw coverage.
#'
#' @param path Character scalar or vector. Path to one HUMAnN file, one merged
#'   HUMAnN table, or multiple single-profile files.
#' @param input_type Input type. One of `"auto"`, `"single"`, or `"merged"`.
#' @param feature_space Functional feature space. One of `"auto"`, `"uniref"`,
#'   `"ko"`, `"ec"`, `"pathway"`, `"genefamily"`, or `"other"`.
#' @param units Feature units. One of `"auto"`, `"cpm"`, `"rpk"`,
#'   `"relative_abundance"`, `"coverage"`, or `"other"`.
#' @param stratification Which rows to retain. One of `"community"`,
#'   `"stratified"`, or `"all"`.
#' @param remove_unmapped Logical. Whether to remove `UNMAPPED` and
#'   `UNINTEGRATED` rows.
#' @param sample_names Optional sample names. For multiple paths, must have the
#'   same length as `path`.
#' @param clean_sample_names Logical. Whether to remove common HUMAnN suffixes
#'   from sample names.
#'
#' @return A data.frame with rows as samples and columns as functional features.
#' @export
read_humann <- function(
    path,
    input_type = c("auto", "single", "merged"),
    feature_space = c(
      "auto", "uniref", "ko", "ec", "pathway", "genefamily", "other"
    ),
    units = c(
      "auto", "cpm", "rpk", "relative_abundance", "coverage", "other"
    ),
    stratification = c("community", "stratified", "all"),
    remove_unmapped = TRUE,
    sample_names = NULL,
    clean_sample_names = TRUE
) {
  input_type <- match.arg(input_type)
  feature_space <- match.arg(feature_space)
  units <- match.arg(units)
  stratification <- match.arg(stratification)

  if (length(path) > 1L) {
    if (input_type == "auto") input_type <- "single"

    if (input_type != "single") {
      cli::cli_abort("Multiple paths currently imply `input_type = 'single'`.")
    }

    return(read_humann_many(
      paths = path,
      feature_space = feature_space,
      units = units,
      stratification = stratification,
      remove_unmapped = remove_unmapped,
      sample_names = sample_names,
      clean_sample_names = clean_sample_names
    ))
  }

  if (input_type == "auto") {
    input_type <- detect_humann_input_type(path)
  }

  switch(
    input_type,
    single = read_humann_single(
      path = path,
      feature_space = feature_space,
      units = units,
      stratification = stratification,
      remove_unmapped = remove_unmapped,
      sample_name = sample_names,
      clean_sample_names = clean_sample_names
    ),
    merged = read_humann_merged(
      path = path,
      feature_space = feature_space,
      units = units,
      stratification = stratification,
      remove_unmapped = remove_unmapped,
      clean_sample_names = clean_sample_names
    )
  )
}

## PRIVATE functions

#' Fast-read a HUMAnN-style table
#'
#' @keywords internal
read_humann_table_fast <- function(path) {
  data.table::fread(
    file = path,
    sep = "\t",
    header = TRUE,
    quote = "",
    check.names = FALSE,
    data.table = FALSE,
    showProgress = TRUE
    # showProgress = interactive()
  )
}

#' Detect whether HUMAnN input is single-profile or merged
#'
#' @keywords internal
detect_humann_input_type <- function(path) {
  lines <- readLines(path, n = 50L, warn = FALSE)
  nonempty <- lines[nzchar(lines)]

  if (length(nonempty) == 0L) {
    cli::cli_abort("Could not detect HUMAnN input type: file is empty.")
  }

  header <- strsplit(nonempty[[1]], "\t", fixed = TRUE)[[1]]

  # HUMAnN single-profile tables usually have two columns:
  # feature name + one abundance column.
  # Merged HUMAnN tables usually have feature name + many sample columns.
  if (length(header) > 2L) {
    return("merged")
  }

  "single"
}


#' Detect HUMAnN feature space
#'
#' @keywords internal
detect_humann_feature_space <- function(feature_ids, path = NULL) {
  ids <- feature_ids
  lower_path <- if (!is.null(path)) tolower(path) else ""

  if (any(grepl("^uniref", ids, ignore.case = TRUE)) ||
      grepl("uniref", lower_path)) {
    return("uniref")
  }

  if (any(grepl("^K[0-9]{5}$", ids)) ||
      grepl("ko|kegg", lower_path)) {
    return("ko")
  }

  if (any(grepl("^[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+$", ids)) ||
      grepl("ec", lower_path)) {
    return("ec")
  }

  if (any(grepl("^PWY-|^UNIPATHWAY|^\\S+:", ids)) ||
      grepl("pathabundance|pathcoverage|pathway", lower_path)) {
    return("pathway")
  }

  if (grepl("genefamilies|genefamily|gene_families", lower_path)) {
    return("genefamily")
  }

  "other"
}


#' Detect HUMAnN units
#'
#' @keywords internal
detect_humann_units <- function(path = NULL, header = character()) {
  lower_path <- if (!is.null(path)) tolower(path) else ""
  lower_header <- tolower(header)

  txt <- paste(c(lower_path, lower_header), collapse = " ")

  if (grepl("cpm|copies per million", txt)) {
    return("cpm")
  }

  if (grepl("rpk|reads per kilobase", txt)) {
    return("rpk")
  }

  if (grepl("relab|relative", txt)) {
    return("relative_abundance")
  }

  if (grepl("coverage|pathcoverage", txt)) {
    return("coverage")
  }

  "other"
}


#' Filter HUMAnN rows by stratification and special rows
#'
#' @keywords internal
filter_humann_features <- function(
    feature_ids,
    stratification = c("community", "stratified", "all"),
    remove_unmapped = TRUE
) {
  stratification <- match.arg(stratification)

  keep <- rep(TRUE, length(feature_ids))

  if (remove_unmapped) {
    keep <- keep & !grepl(
      "^UNMAPPED$|^UNINTEGRATED$",
      feature_ids,
      ignore.case = TRUE
    )
  }

  if (stratification == "community") {
    keep <- keep & !grepl("\\|", feature_ids)
  }

  if (stratification == "stratified") {
    keep <- keep & grepl("\\|", feature_ids)
  }

  keep
}


#' Read a single HUMAnN profile
#'
#' @keywords internal
read_humann_single <- function(
    path,
    feature_space = c(
      "auto", "uniref", "ko", "ec", "pathway", "genefamily", "other"
    ),
    units = c(
      "auto", "cpm", "rpk", "relative_abundance", "coverage", "other"
    ),
    stratification = c("community", "stratified", "all"),
    remove_unmapped = TRUE,
    sample_name = NULL,
    clean_sample_names = TRUE
) {
  feature_space <- match.arg(feature_space)
  units <- match.arg(units)
  stratification <- match.arg(stratification)

  raw <- read_humann_table_fast(path)
  
  if (ncol(raw) < 2L) {
    cli::cli_abort("Single HUMAnN profile must have at least two columns.")
  }
  
  feature_col <- colnames(raw)[[1]]
  value_col <- colnames(raw)[[2]]
  
  feature_ids <- raw[[feature_col]]

  if (feature_space == "auto") {
    feature_space <- detect_humann_feature_space(feature_ids, path)
  }

  if (units == "auto") {
    units <- detect_humann_units(path, colnames(raw))
  }

  keep <- filter_humann_features(
    feature_ids = feature_ids,
    stratification = stratification,
    remove_unmapped = remove_unmapped
  )

  values <- raw[[value_col]][keep]
  names(values) <- feature_ids[keep]

  if (is.null(sample_name)) {
    sample_name <- tools::file_path_sans_ext(basename(path))
  }

  if (clean_sample_names) {
    sample_name <- clean_humann_sample_names(sample_name)
  }

  out <- as.data.frame(t(values), check.names = FALSE)
  rownames(out) <- sample_name

  attr(out, "input_type") <- "single"
  attr(out, "feature_space") <- feature_space
  attr(out, "units") <- units
  attr(out, "stratification") <- stratification

  out
}


#' Read a merged HUMAnN table
#'
#' @keywords internal
read_humann_merged <- function(
    path,
    feature_space = c(
      "auto", "uniref", "ko", "ec", "pathway", "genefamily", "other"
    ),
    units = c(
      "auto", "cpm", "rpk", "relative_abundance", "coverage", "other"
    ),
    stratification = c("community", "stratified", "all"),
    remove_unmapped = TRUE,
    clean_sample_names = TRUE
) {
  feature_space <- match.arg(feature_space)
  units <- match.arg(units)
  stratification <- match.arg(stratification)

  df <- read_humann_table_fast(path)
  
  if (ncol(df) < 2L) {
    cli::cli_abort("Merged HUMAnN table must have at least two columns.")
  }
  
  feature_col <- colnames(df)[[1]]
  feature_ids <- df[[feature_col]]
  
  abund <- df[, -1, drop = FALSE]
  rownames(abund) <- feature_ids

  if (feature_space == "auto") {
    feature_space <- detect_humann_feature_space(feature_ids, path)
  }

  if (units == "auto") {
    units <- detect_humann_units(path, colnames(df))
  }

  keep <- filter_humann_features(
    feature_ids = feature_ids,
    stratification = stratification,
    remove_unmapped = remove_unmapped
  )

  out <- as.data.frame(t(abund[keep, , drop = FALSE]), check.names = FALSE)
  
  if (clean_sample_names) {
    rownames(out) <- clean_humann_sample_names(rownames(out))
  }

  attr(out, "input_type") <- "merged"
  attr(out, "feature_space") <- feature_space
  attr(out, "units") <- units
  attr(out, "stratification") <- stratification

  out
}


#' Read many single HUMAnN profiles
#'
#' @keywords internal
read_humann_many <- function(
    paths,
    feature_space = c(
      "auto", "uniref", "ko", "ec", "pathway", "genefamily", "other"
    ),
    units = c(
      "auto", "cpm", "rpk", "relative_abundance", "coverage", "other"
    ),
    stratification = c("community", "stratified", "all"),
    remove_unmapped = TRUE,
    sample_names = NULL,
    clean_sample_names = TRUE
) {
  feature_space <- match.arg(feature_space)
  units <- match.arg(units)
  stratification <- match.arg(stratification)

  if (!is.null(sample_names) && length(sample_names) != length(paths)) {
    cli::cli_abort("`sample_names` must have the same length as `paths`.")
  }

  pieces <- lapply(seq_along(paths), function(i) {
    read_humann_single(
      path = paths[[i]],
      feature_space = feature_space,
      units = units,
      stratification = stratification,
      remove_unmapped = remove_unmapped,
      sample_name = if (!is.null(sample_names)) sample_names[[i]] else NULL,
      clean_sample_names = clean_sample_names
    )
  })

  all_features <- Reduce(union, lapply(pieces, colnames))

  pieces <- lapply(pieces, function(x) {
    missing <- setdiff(all_features, colnames(x))

    if (length(missing) > 0L) {
      x[missing] <- 0
    }

    x[, all_features, drop = FALSE]
  })

  out <- do.call(rbind, pieces)

  if (clean_sample_names) {
    rownames(out) <- clean_humann_sample_names(rownames(out))
  }

  attr(out, "input_type") <- "many_single"
  attr(out, "feature_space") <- feature_space
  attr(out, "units") <- units
  attr(out, "stratification") <- stratification

  out
}


#' Clean HUMAnN sample names
#'
#' @keywords internal
clean_humann_sample_names <- function(x) {
  x <- sub("_genefamilies$", "", x)
  x <- sub("_pathabundance$", "", x)
  x <- sub("_pathcoverage$", "", x)
  x <- sub("_profile$", "", x)
  x <- sub("\\.tsv$", "", x)
  x <- sub("\\.txt$", "", x)
  x
}
