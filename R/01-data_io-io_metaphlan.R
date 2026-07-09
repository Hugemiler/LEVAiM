#####
# 01. Data IO - MetaPhlAn 3/4
#####

## PUBLIC functions

#' Read MetaPhlAn profiles
#'
#' Reads MetaPhlAn 3/4 output from a single profile, a merged table, or a
#' vector of single-profile paths.
#'
#' @param path Character scalar or vector. Path to one MetaPhlAn file, one
#'   merged MetaPhlAn table, or multiple single-profile files.
#' @param tax_level Taxonomic level to retain. One of `"kingdom"`, `"phylum"`,
#'   `"class"`, `"order"`, `"family"`, `"genus"`, `"species"`, or `"strain"`.
#' @param input_type Input type. One of `"auto"`, `"single"`, or `"merged"`.
#' @param metaphlan_version MetaPhlAn version. One of `"auto"`, `"3"`, or `"4"`.
#' @param include_unclassified Logical. Whether to retain `UNCLASSIFIED` /
#'   `UNIDENTIFIED` rows when present.
#' @param sample_names Optional sample names. For multiple paths, must have the
#'   same length as `path`.
#' @param clean_sample_names Logical. Whether to remove common MetaPhlAn suffixes
#'   from sample names.
#'
#' @return A data.frame with rows as samples and columns as taxa.
#' @export
read_metaphlan <- function(
    path,
    tax_level = "species",
    input_type = c("auto", "single", "merged"),
    metaphlan_version = c("auto", "3", "4"),
    include_unclassified = TRUE,
    sample_names = NULL,
    clean_sample_names = TRUE
) {
  input_type <- match.arg(input_type)
  metaphlan_version <- match.arg(metaphlan_version)

  if (length(path) > 1L) {
    if (input_type == "auto") input_type <- "single"

    if (input_type != "single") {
      cli::cli_abort("Multiple paths currently imply `input_type = 'single'`.")
    }

    return(read_metaphlan_many(
      paths = path,
      tax_level = tax_level,
      metaphlan_version = metaphlan_version,
      include_unclassified = include_unclassified,
      sample_names = sample_names,
      clean_sample_names = clean_sample_names
    ))
  }

  if (input_type == "auto") {
    input_type <- detect_metaphlan_input_type(path)
  }

  switch(
    input_type,
    single = read_metaphlan_single(
      path = path,
      tax_level = tax_level,
      metaphlan_version = metaphlan_version,
      include_unclassified = include_unclassified,
      sample_name = sample_names,
      clean_sample_names = clean_sample_names
    ),
    merged = read_metaphlan_merged(
      path = path,
      tax_level = tax_level,
      metaphlan_version = metaphlan_version,
      include_unclassified = include_unclassified,
      clean_sample_names = clean_sample_names
    )
  )
}

## PRIVATE functions

detect_metaphlan_input_type <- function(path) {
  lines <- readLines(path, n = 100L, warn = FALSE)

  nonempty <- lines[nzchar(lines)]

  if (length(nonempty) == 0L) {
    cli::cli_abort("Could not detect MetaPhlAn input type: file is empty.")
  }

  header_candidates <- nonempty[
    grepl("clade_name|relative_abundance|NCBI_tax_id", nonempty, ignore.case = TRUE)
  ]

  if (length(header_candidates) > 0L) {
    header <- strsplit(header_candidates[[1]], "\t", fixed = TRUE)[[1]]

    if (any(grepl("relative_abundance", header, ignore.case = TRUE))) {
      return("single")
    }
  }

  data_lines <- nonempty[!startsWith(nonempty, "#")]

  if (length(data_lines) == 0L) {
    cli::cli_abort("Could not detect MetaPhlAn input type: no data lines found.")
  }

  first_fields <- strsplit(data_lines[[1]], "\t", fixed = TRUE)[[1]]

  if (length(first_fields) > 2L) {
    return("merged")
  }

  "single"
}

detect_metaphlan_version <- function(feature_ids, header = character()) {
  feature_ids_upper <- toupper(feature_ids)
  header_lower <- tolower(header)

  if (any(feature_ids_upper == "UNCLASSIFIED")) {
    return("4")
  }

  if (any(grepl("sgb", feature_ids, ignore.case = TRUE))) {
    return("4")
  }

  if (any(grepl("estimated_number_of_reads", header_lower))) {
    return("4")
  }

  "3"
}

filter_metaphlan_taxa <- function(
    feature_ids,
    tax_level = "species",
    include_unclassified = TRUE
) {
  level_prefix <- c(
    kingdom = "k__",
    phylum  = "p__",
    class   = "c__",
    order   = "o__",
    family  = "f__",
    genus   = "g__",
    species = "s__",
    strain  = "t__"
  )

  tax_level <- match.arg(tax_level, names(level_prefix))
  pfx <- level_prefix[[tax_level]]
  rank_i <- match(tax_level, names(level_prefix))

  mask <- grepl(paste0("\\|", pfx), feature_ids) |
    startsWith(feature_ids, pfx)

  if (rank_i < length(level_prefix)) {
    next_levels <- level_prefix[(rank_i + 1L):length(level_prefix)]
    deeper_pattern <- paste0("\\|(", paste(next_levels, collapse = "|"), ")")
    mask <- mask & !grepl(deeper_pattern, feature_ids)
  }

  if (include_unclassified) {
    mask <- mask | toupper(feature_ids) %in% c("UNCLASSIFIED", "UNIDENTIFIED")
  }

  mask
}

read_metaphlan_single <- function(
    path,
    tax_level = "species",
    metaphlan_version = c("auto", "3", "4"),
    include_unclassified = TRUE,
    sample_name = NULL,
    clean_sample_names = TRUE
) {
  metaphlan_version <- match.arg(metaphlan_version)

  raw <- utils::read.table(
    path,
    sep = "\t",
    header = TRUE,
    comment.char = "#",
    check.names = FALSE,
    quote = "",
    stringsAsFactors = FALSE
  )

  clade_col <- intersect(
    c("clade_name", "clade", "taxon", "#clade_name"),
    colnames(raw)
  )[1]

  abundance_col <- grep(
    "relative_abundance|abundance",
    colnames(raw),
    ignore.case = TRUE,
    value = TRUE
  )[1]

  if (is.na(clade_col) || is.na(abundance_col)) {
    cli::cli_abort("Could not identify clade and abundance columns in single MetaPhlAn profile.")
  }

  feature_ids <- raw[[clade_col]]

  if (metaphlan_version == "auto") {
    metaphlan_version <- detect_metaphlan_version(feature_ids, colnames(raw))
  }

  keep <- filter_metaphlan_taxa(
    feature_ids = feature_ids,
    tax_level = tax_level,
    include_unclassified = include_unclassified
  )

  values <- raw[[abundance_col]][keep]
  names(values) <- feature_ids[keep]

  if (is.null(sample_name)) {
    sample_name <- tools::file_path_sans_ext(basename(path))
    sample_name <- sub("_profile$", "", sample_name)
  }

  if (clean_sample_names) {
    sample_name <- clean_metaphlan_sample_names(sample_name)
  }

  out <- as.data.frame(t(values), check.names = FALSE)
  rownames(out) <- sample_name

  attr(out, "metaphlan_version") <- metaphlan_version
  attr(out, "input_type") <- "single"
  attr(out, "tax_level") <- tax_level

  out
}

read_metaphlan_merged <- function(
    path,
    tax_level = "species",
    metaphlan_version = c("auto", "3", "4"),
    include_unclassified = TRUE,
    clean_sample_names = TRUE
) {
  metaphlan_version <- match.arg(metaphlan_version)

  df <- utils::read.table(
    path,
    sep = "\t",
    header = TRUE,
    row.names = 1,
    comment.char = "#",
    check.names = FALSE,
    quote = "",
    stringsAsFactors = FALSE
  )

  feature_ids <- rownames(df)

  if (metaphlan_version == "auto") {
    metaphlan_version <- detect_metaphlan_version(feature_ids, colnames(df))
  }

  keep <- filter_metaphlan_taxa(
    feature_ids = feature_ids,
    tax_level = tax_level,
    include_unclassified = include_unclassified
  )

  out <- as.data.frame(t(df[keep, , drop = FALSE]), check.names = FALSE)

  if (clean_sample_names) {
    rownames(out) <- clean_metaphlan_sample_names(rownames(out))
  }

  attr(out, "metaphlan_version") <- metaphlan_version
  attr(out, "input_type") <- "merged"
  attr(out, "tax_level") <- tax_level

  out
}

clean_metaphlan_sample_names <- function(x) {
  x <- sub("_profile$", "", x)
  x <- sub("\\.txt$", "", x)
  x <- sub("\\.tsv$", "", x)
  x
}

read_metaphlan_many <- function(
    paths,
    tax_level = "species",
    metaphlan_version = c("auto", "3", "4"),
    include_unclassified = TRUE,
    sample_names = NULL,
    clean_sample_names = TRUE
) {
  metaphlan_version <- match.arg(metaphlan_version)

  if (!is.null(sample_names) && length(sample_names) != length(paths)) {
    cli::cli_abort("`sample_names` must have the same length as `paths`.")
  }

  pieces <- lapply(seq_along(paths), function(i) {
    read_metaphlan_single(
      path = paths[[i]],
      tax_level = tax_level,
      metaphlan_version = metaphlan_version,
      include_unclassified = include_unclassified,
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
    rownames(out) <- clean_metaphlan_sample_names(rownames(out))
  }

  attr(out, "input_type") <- "many_single"
  attr(out, "tax_level") <- tax_level

  out
}
