#' Balanced Allocation of Subjects to Treatments
#'
#' Convenience wrapper function for \code{\link{anticlustering}} using different defaults and additional functionality, provided by Dimitry Wintermantel as used in Wintermantel et al (2026; <doi:10.48550/arXiv.2607.07543>) used in . This function assigns subjects to treatment groups while balancing specified covariates and ensuring replicates per group. If a table of groups within treatments is supplied via `group_data`, the subjects (already assigned to treatments) are instead allocated to those groups within each treatment, again balancing the covariates.
#'
#' @param data A data frame containing the subjects to be allocated, including columns for the covariates. If `group_data` is supplied, it must also contain the treatment column named in `treatment_var`.
#' @param covariates A vector of column names in `data` representing covariates to balance.
#' @param treatments A vector specifying the treatment group labels. Required unless `group_data` is supplied, in which case the target groups are taken from `group_data`.
#' @param treatment_var A string specifying the name of the treatment column. If `group_data` is not supplied, this column is created in the output (defaults to `"Treatment"`). If `group_data` is supplied, it names the existing treatment column, which must be present in both `data` and `group_data`.
#' @param group_data An optional data frame describing the available groups within each treatment. Must include columns named in `group_var` and `treatment_var` indicating which groups belong to which treatment. Any further columns (e.g. group-level covariates) are merged into the output but are not used for balancing. Defaults to `NULL`.
#' @param group_var A string specifying the name of the group column in `group_data` to which subjects will be allocated within each treatment (e.g., site, cage, plot). Required when `group_data` is supplied.
#' @param n_replicates Integer. Specifies the number of replicates per group. If `NULL`, either all subjects are allocated or the highest possible equal number per group is allocated depending on `ensure_equal_n_replicates`. Defaults to `NULL`.
#' @param ensure_equal_n_replicates Logical. If `TRUE`, ensures equal numbers of replicates across groups. Defaults to `TRUE`. If `FALSE`, all subjects are allocated.
#' @param keep_excluded_data Logical. If `TRUE`, retains subjects that were not allocated within the specified constraints. Defaults to `FALSE`.
#' @param objective A string specifying the objective function for anticlustering. Options include `"variance"` (default), `"diversity"`, `"average-diversity"`, `"kplus"`, and `"dispersion"`.
#' @param method A string specifying the optimization method for anticlustering. Options include `"local-maximum"` (default), `"exchange"`, `"brusco"`, `"ilp"`, and `"2PML"`.
#' @param repetitions Integer. Specifies the number of times the optimization is repeated when using heuristic methods (`"exchange"`, `"local-maximum"`, `"brusco"`, or `"2PML"`). The best solution is selected. Defaults to `10`.
#' @param match_within A column name in `data` (optional). Specifies a variable within which matching should occur, ensuring that subjects are grouped within subsets defined by this variable. Defaults to `NULL`.
#' @param standardize Logical. If `TRUE`, covariates are standardized via `scale()` before optimization starts. Defaults to `TRUE`.
#'
#' @return
#' A data frame with subjects assigned to treatments, in the column named by `treatment_var`. If `group_data` is supplied, the output instead contains the subjects assigned to groups within their treatment (in the column named by `group_var`), merged with the corresponding group information from `group_data`. If `keep_excluded_data = TRUE`, excluded subjects are included in the output with `NA` in the allocation column.
#'
#' @details
#' The function balances covariates among groups by creating sets of similar individuals based on the specified covariates, and then assigning the target labels to these sets so that similar individuals are spread across different groups.
#'
#' If `group_data` is supplied, the allocation is carried out separately for each treatment: for every treatment level present in `data[[treatment_var]]`, the subjects of that treatment are allocated to the groups that belong to it in `group_data`, and the result is merged with the corresponding rows of `group_data`.
#'
#' - If `ensure_equal_n_replicates = TRUE`, the number of replicates per group is enforced.
#' - If `keep_excluded_data = TRUE`, subjects that cannot be allocated under the constraints are retained in the output with `NA` in the allocation column.
#' - Covariates can be scaled to standardize their ranges if `standardize = TRUE`.
#' - If `match_within` is specified, subjects are matched within levels of the specified variable, ensuring allocations respect the structure defined by this variable.
#'
#' @section Validations:
#' The function validates the following conditions:
#' - The required arguments are provided (`data`, `covariates`, and either `treatments` or, when `group_data` is supplied, `group_var`).
#' - Covariates exist in `data`, are numeric, and do not take the same value across all rows.
#' - When `group_data` is supplied, `treatment_var` exists in both `data` and `group_data`, and `group_var` exists in `group_data`.
#' - `match_within`, if specified, exists in `data`.
#'
#' @section How to cite:
#' If you use this function in academic work, please cite:
#' Wintermantel, D., Osterman, J., Mair, M. M., & Hartig, F. (2026).
#' *Equivalence testing in pesticide risk assessment – Evaluation and practical
#' https://doi.org/10.48550/arXiv.2607.07543
#'
#' The treatment allocation implemented here relies on anticlustering algorithms.
#' Users are encouraged to also cite the `anticlust` package where appropriate
#' (see `citation("anticlust")`).
#'
#' @examples
#'
#' # Example dataset: Bee subjects
#' example_bee_data <- data.frame(
#'   ID = as.factor(as.character(seq(1, 100, 1))),
#'   Bee_count = rnorm(100, mean = 1000, sd = 200),
#'   Colony_weight = rnorm(100, mean = 500, sd = 100)
#' )
#'
#' # Example dataset: Sites (groups) with a site-level covariate and a blocking variable
#' example_site_data <- data.frame(
#'   Site = as.factor(as.character(seq(1, 16, 1))),
#'   Irrigated = as.factor(as.character(c(rep("yes", 4), rep("no", 12)))),
#'   Field_quality = rnorm(16, mean = 3, sd = 0.5)
#' )
#'
#' treatments <- c("Control", "Pesticide")
#' bee_covariates <- c("Bee_count", "Colony_weight")
#'
#' # 1) Allocate subjects to treatments (balancing colony-level covariates)
#' allocated_bee_data <- experimental_allocation(
#'   data = example_bee_data,
#'   treatments = treatments,
#'   covariates = bee_covariates
#' )
#'
#' # 2) Allocate sites to treatments (balancing site quality, matching within irrigation status)
#' allocated_site_data <- experimental_allocation(
#'   data = example_site_data,
#'   treatments = treatments,
#'   covariates = "Field_quality",
#'   match_within = "Irrigated"
#' )
#'
#' # 3) Allocate subjects to specific sites within each treatment
#' allocated_data <- experimental_allocation(
#'   data = allocated_bee_data,
#'   covariates = bee_covariates,
#'   treatment_var = "Treatment",
#'   group_data = allocated_site_data,
#'   group_var = "Site"
#' )
#' head(allocated_data)
#'
#' @export
experimental_allocation <- function(data,
                                covariates,
                                treatments = NULL,
                                treatment_var = "Treatment",
                                group_data = NULL,
                                group_var = NULL,
                                n_replicates = NULL,
                                ensure_equal_n_replicates = TRUE,
                                keep_excluded_data = FALSE,
                                objective = "variance",
                                method = "local-maximum",
                                repetitions = 10,
                                match_within = NULL,
                                standardize = TRUE) {

  # ---- argument checks -----------------------------------------------------
  require_arg <- function(arg, name) {
    if (missing(arg) || is.null(arg)) {
      stop("The required argument '", name, "' is missing or NULL.", call. = FALSE)
    }
  }
  require_arg(data, "data")
  require_arg(covariates, "covariates")

  nested <- !is.null(group_data)

  if (!nested) {
    require_arg(treatments, "treatments")
  } else {
    require_arg(group_var, "group_var")
    if (!treatment_var %in% names(data)) {
      stop("When 'group_data' is supplied, 'data' must contain the treatment column '",
           treatment_var, "'.", call. = FALSE)
    }
    if (!treatment_var %in% names(group_data)) {
      stop("'group_data' must contain the treatment column '",
           treatment_var, "'.", call. = FALSE)
    }
    if (!group_var %in% names(group_data)) {
      stop("'group_data' must contain the group column '",
           group_var, "'.", call. = FALSE)
    }
  }

  missing_cov <- setdiff(covariates, colnames(data))
  if (length(missing_cov) > 0) {
    stop("Covariates missing in 'data': ", paste(missing_cov, collapse = ", "),
         call. = FALSE)
  }
  non_numeric <- covariates[!vapply(data[, covariates, drop = FALSE], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop("Covariates are not numeric: ", paste(non_numeric, collapse = ", "),
         call. = FALSE)
  }
  constant <- covariates[vapply(data[, covariates, drop = FALSE],
                                function(x) length(unique(x[!is.na(x)])) == 1, logical(1))]
  if (length(constant) > 0) {
    stop("Covariates take the same value for all observations: ",
         paste(constant, collapse = ", "), call. = FALSE)
  }
  if (!is.null(match_within) && !match_within %in% colnames(data)) {
    stop("'match_within' column '", match_within, "' is not present in 'data'.",
         call. = FALSE)
  }
  if (!is.null(n_replicates) && !ensure_equal_n_replicates) {
    ensure_equal_n_replicates <- TRUE
    warning("'ensure_equal_n_replicates' was changed to TRUE as 'n_replicates' was specified.",
            call. = FALSE)
  }

  # ---- allocate one set of subjects to one set of target labels ------------
  allocate_one <- function(df, targets, out_var, label = NULL) {
    ctx <- if (is.null(label)) "" else paste0(" (", treatment_var, " = ", label, ")")
    n_targets <- length(targets)
    if (n_targets < 2) {
      stop("Fewer than two target groups available", ctx, ".", call. = FALSE)
    }
    if (nrow(df) < n_targets) {
      stop("Fewer subjects than target groups", ctx, ".", call. = FALSE)
    }
    const_here <- covariates[vapply(df[, covariates, drop = FALSE],
                                    function(x) length(unique(x[!is.na(x)])) == 1, logical(1))]
    if (length(const_here) > 0) {
      stop("Covariate(s) take the same value for all observations", ctx, ": ",
           paste(const_here, collapse = ", "), call. = FALSE)
    }

    targets_shuffled <- sample(targets)
    n_rep <- if (is.null(n_replicates)) floor(nrow(df) / n_targets) else n_replicates

    covariate_data <- if (standardize) {
      scale(df[, covariates, drop = FALSE])
    } else {
      as.matrix(df[, covariates, drop = FALSE])
    }
    match_within_vector <- if (!is.null(match_within)) df[[match_within]] else NULL

    set_col <- ".set_similar"

    if (ensure_equal_n_replicates || n_targets * n_rep == nrow(df)) {
      # Create sets of similar individuals based on covariates
      df[[set_col]] <- matching(covariate_data, p = n_targets,
                                match_within = match_within_vector)
      keep <- !is.na(df[[set_col]]) & df[[set_col]] <= n_rep
      included <- df[keep, , drop = FALSE]

      # Assign individuals of the same set to different groups
      included[[out_var]] <- anticlustering(
        included[, covariates, drop = FALSE],
        K = n_targets,
        objective = objective,
        method = method,
        categories = included[[set_col]],
        repetitions = repetitions
      )
      included[[set_col]] <- NULL

      if (keep_excluded_data) {
        excluded <- df[!keep, , drop = FALSE]
        excluded[[set_col]] <- NULL
        excluded[[out_var]] <- NA
        out <- rbind(included, excluded)
      } else {
        out <- included
      }
    } else {
      out <- df
      out[[out_var]] <- anticlustering(
        covariate_data,
        K = n_targets,
        objective = objective,
        method = method,
        categories = match_within_vector,
        repetitions = repetitions
      )
    }

    # Map the anticlustering output to the (shuffled) target labels
    target_map <- stats::setNames(targets_shuffled, as.character(seq_len(n_targets)))
    out[[out_var]] <- unname(target_map[as.character(out[[out_var]])])
    out[[out_var]] <- factor(out[[out_var]], levels = targets)
    out
  }

  # ---- allocate subjects to treatments -------------------------------------
  if (!nested) {
    out <- allocate_one(data, treatments, treatment_var)
    out <- out[, c(treatment_var, setdiff(names(out), treatment_var)), drop = FALSE]
    return(out)
  }

  # ---- allocate subjects to groups within each treatment -------------------
  data[[treatment_var]] <- as.character(data[[treatment_var]])
  group_data[[treatment_var]] <- as.character(group_data[[treatment_var]])
  group_data[[group_var]] <- as.character(group_data[[group_var]])
  treatment_levels <- unique(data[[treatment_var]])

  allocate_in_treatment <- function(trt) {
    grp_rows <- group_data[group_data[[treatment_var]] == trt, , drop = FALSE]
    targets <- unique(grp_rows[[group_var]])
    subjects <- data[data[[treatment_var]] == trt, , drop = FALSE]

    if (group_var %in% names(subjects)) subjects[[group_var]] <- NULL
    subjects <- allocate_one(subjects, targets, out_var = group_var, label = trt)
    subjects[[group_var]] <- as.character(subjects[[group_var]])

    # Merge subject and group data for this treatment
    merge(grp_rows, subjects, by = intersect(names(grp_rows), names(subjects)),
          all = TRUE)
  }

  allocated_list <- lapply(treatment_levels, allocate_in_treatment)
  allocated <- do.call(rbind, allocated_list)
  rownames(allocated) <- NULL

  front <- intersect(c(treatment_var, group_var), names(allocated))
  allocated[, c(front, setdiff(names(allocated), front)), drop = FALSE]
}