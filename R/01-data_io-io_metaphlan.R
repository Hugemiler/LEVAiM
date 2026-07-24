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
  x <- sub("_mpa.*$", "", x)
  x <- sub("_mpa$", "", x)
  x <- sub("\\.txt$", "", x)
  x <- sub("\\.tsv$", "", x)
  x
}


#' Aggregate a taxonomic assay to a broader rank
#'
#' Sums MetaPhlAn lineage columns at a declared taxonomic rank. Composition is
#' handled explicitly: values can be retained as supplied or renormalized to a
#' constant row total of 100 after aggregation.
#'
#' @param assay A matrix or data.frame with samples in rows and MetaPhlAn
#'   lineage identifiers in columns.
#' @param tax_level Target rank. One of `"kingdom"`, `"phylum"`, `"class"`,
#'   `"order"`, `"family"`, or `"genus"`.
#' @param composition_policy One of `"as_is"` or `"renormalize"`.
#'
#' @return A data.frame with samples in rows and aggregated taxa in columns.
#' @export
aggregate_taxa <- function(
    assay,
    tax_level = c("genus", "family", "order", "class", "phylum", "kingdom"),
    composition_policy = c("as_is", "renormalize")
) {
  tax_level <- match.arg(tax_level)
  composition_policy <- match.arg(composition_policy)
  assay <- as.data.frame(assay, check.names = FALSE)

  if (ncol(assay) == 0L || is.null(colnames(assay))) {
    cli::cli_abort("`assay` must have at least one named taxonomic feature.")
  }
  numeric_columns <- vapply(assay, is.numeric, logical(1))
  if (!all(numeric_columns)) {
    cli::cli_abort("Every column in `assay` must be numeric.")
  }
  if (any(as.matrix(assay) < 0, na.rm = TRUE)) {
    cli::cli_abort("Taxonomic abundances must be non-negative.")
  }

  prefix <- c(
    kingdom = "k__", phylum = "p__", class = "c__", order = "o__",
    family = "f__", genus = "g__"
  )[[tax_level]]
  ranks <- vapply(
    strsplit(colnames(assay), "|", fixed = TRUE),
    function(parts) {
      hit <- parts[startsWith(parts, prefix)]
      if (length(hit)) hit[[1]] else NA_character_
    },
    character(1)
  )
  keep <- !is.na(ranks) & nzchar(ranks)
  if (!any(keep)) {
    cli::cli_abort(
      "No {.val {tax_level}} identifiers were found in the assay column names."
    )
  }

  values <- as.matrix(assay[, keep, drop = FALSE])
  groups <- factor(ranks[keep], levels = unique(ranks[keep]))
  aggregated <- vapply(
    levels(groups),
    function(group) rowSums(values[, groups == group, drop = FALSE], na.rm = TRUE),
    numeric(nrow(values))
  )
  if (is.null(dim(aggregated))) {
    aggregated <- matrix(aggregated, ncol = 1L)
  }
  colnames(aggregated) <- levels(groups)
  rownames(aggregated) <- rownames(assay)

  if (composition_policy == "renormalize") {
    totals <- rowSums(aggregated, na.rm = TRUE)
    usable <- is.finite(totals) & totals > 0
    aggregated[usable, ] <- aggregated[usable, , drop = FALSE] /
      totals[usable] * 100
  }

  out <- as.data.frame(aggregated, check.names = FALSE)
  attr(out, "tax_level") <- tax_level
  attr(out, "composition_policy") <- composition_policy
  attr(out, "source_row_totals") <- rowSums(as.matrix(assay), na.rm = TRUE)
  out
}


#' Derive community-level features from a taxonomic assay
#'
#' Computes prespecified sample-level summaries that can be modeled through the
#' same longitudinal workflow as individual taxa. Diversity summaries use
#' within-sample proportions over the supplied feature universe; profiled mass
#' retains the original row total.
#'
#' @param assay A non-negative numeric matrix or data.frame with samples in rows.
#' @param detection_threshold Abundance above which a feature contributes to
#'   observed richness.
#'
#' @return A data.frame containing `richness`, `shannon`, `simpson`,
#'   `dominance`, and `profiled_mass`.
#' @export
community_features <- function(assay, detection_threshold = 0) {
  assay <- as.data.frame(assay, check.names = FALSE)
  values <- as.matrix(assay)
  if (ncol(values) == 0L || !is.numeric(values) || any(values < 0, na.rm = TRUE)) {
    cli::cli_abort("`assay` must contain non-negative numeric abundances.")
  }
  if (!is.numeric(detection_threshold) || length(detection_threshold) != 1L ||
      !is.finite(detection_threshold) || detection_threshold < 0) {
    cli::cli_abort("`detection_threshold` must be one finite non-negative number.")
  }

  profiled_mass <- rowSums(values, na.rm = TRUE)
  proportions <- matrix(0, nrow = nrow(values), ncol = ncol(values))
  usable <- is.finite(profiled_mass) & profiled_mass > 0
  proportions[usable, ] <- values[usable, , drop = FALSE] / profiled_mass[usable]
  log_proportions <- matrix(0, nrow = nrow(values), ncol = ncol(values))
  positive <- proportions > 0
  log_proportions[positive] <- log(proportions[positive])

  out <- data.frame(
    richness = rowSums(values > detection_threshold, na.rm = TRUE),
    shannon = -rowSums(proportions * log_proportions, na.rm = TRUE),
    simpson = 1 - rowSums(proportions^2, na.rm = TRUE),
    dominance = apply(proportions, 1L, max, na.rm = TRUE),
    profiled_mass = profiled_mass,
    row.names = rownames(assay),
    check.names = FALSE
  )
  out$shannon[!usable] <- NA_real_
  out$simpson[!usable] <- NA_real_
  out$dominance[!usable] <- NA_real_
  attr(out, "detection_threshold") <- detection_threshold
  out
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
