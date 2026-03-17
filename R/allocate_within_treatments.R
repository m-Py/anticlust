#' Allocate Subjects to Groups *Within* Treatments
#'
#' This function assigns subjects to experimental groups *within each treatment level* while balancing specified covariates and (optionally) ensuring equal replicates per group. Internally, it applies [`allocate_treatments()`] to each treatment subset, using the available groups as the target labels.
#'
#' @param data_subjects A data frame containing the subjects to be allocated. Must include a column named in `treatment_var` as well as the specified `covariates`.
#' @param data_groups A data frame describing the available groups within each treatment. Must include columns named in `group_var` and `treatment_var` indicating which groups belong to which treatment.
#' @param treatment_var A string specifying the name of the treatment column present in both `data_subjects` and `data_groups`.
#' @param group_var A string specifying the name of the group column in `data_groups` to which subjects will be allocated within each treatment (e.g., site, cage, plot).
#' @param covariates A vector of column names in `data_subjects` representing covariates to balance.
#' @param ensure_equal_n_replicates Logical. If `TRUE`, ensures equal numbers of replicates across groups within each treatment (handled inside `allocate_treatments()`). Defaults to `TRUE`. If `FALSE`, all subjects are allocated without enforcing equal replicates.
#' @param keep_excluded_data Logical. If `TRUE`, retains subjects that were not allocated within the specified constraints (as determined by `allocate_treatments()`). Defaults to `FALSE`.
#' @param objective A string specifying the objective function for anticlustering. Options include `"variance"` (default), `"diversity"`, `"average-diversity"`, `"kplus"`, and `"dispersion"`.
#' @param method A string specifying the optimization method for anticlustering. Options include `"local-maximum"` (default), `"exchange"`, `"brusco"`, `"ilp"`, and `"2PML"`.
#' @param repetitions Integer. Specifies the number of times the optimization is repeated when using heuristic methods (`"exchange"`, `"local-maximum"`, `"brusco"`, or `"2PML"`). The best solution is selected. Defaults to `10`.
#' @param match_within A column name in `data_subjects` (optional). Specifies a variable within which matching should occur, ensuring that subjects are grouped within subsets defined by this variable before allocation. Defaults to `NULL`.
#' @param standardize Logical. If `TRUE`, covariates are standardized via `scale()` before optimization starts. Defaults to `TRUE`.
#'
#' @return
#' A data frame containing the merged group information from `data_groups` and the subject allocations from `data_subjects`. The output includes:
#' - the original treatment column `treatment_var`,
#' - the allocated group column `group_var` (created by assigning subjects to groups within each treatment),
#' - all remaining columns from both inputs (joined per treatment).
#' If `keep_excluded_data = TRUE`, excluded subjects are included with `NA` in the group column.
#'
#' @details
#' The function operates treatment-by-treatment:
#' 1. For each treatment level present in `data_subjects[[treatment_var]]`, it identifies the corresponding set of groups in `data_groups` that belong to that treatment.
#' 2. It subsets the subjects to the current treatment and calls [`allocate_treatments()`], passing the available groups as the `treatments` argument, to balance the specified `covariates`.
#' 3. It restores the original treatment column and merges the per-treatment allocations with the corresponding rows in `data_groups`.
#'
#' - If `ensure_equal_n_replicates = TRUE`, the number of replicates per group (within a treatment) is enforced by the internal call to `allocate_treatments()`.
#' - If `keep_excluded_data = TRUE`, subjects that cannot be allocated under the constraints are retained in the output with `NA` in `group_var`.
#' - Covariates can be scaled to standardize their ranges if `standardize = TRUE`.
#' - If `match_within` is specified, subjects are matched within levels of the specified variable prior to allocation, respecting that structure.
#'
#' @section Validations:
#' The function validates the following conditions:
#' - All required arguments (`data_subjects`, `data_groups`, `group_var`, `covariates`) are provided.
#' - `treatment_var` exists in `data_subjects`; `group_var` exists in `data_groups`.
#' - Types are harmonized internally (e.g., factors/characters) to ensure consistent merging and labeling.
#' Additional checks on covariates (presence, numeric type, non-constant) are performed inside [`allocate_treatments()`].
#'
#' @section How to cite:
#' If you use this function in academic work, please cite:
#' Wintermantel, D., Osterman, J., Mair, M. M., & Hartig, F. (in preparation).
#' *Equivalence testing in pesticide risk assessment – Evaluation and practical guidance
#' for design, analysis and interpretation*.
#'
#' The treatment allocation implemented here relies on anticlustering functions
#' Users are encouraged to also cite the `anticlust` package where appropriate
#' (see `citation("anticlust")`).
#'
#' @examples
#' set.seed(123)
#'
#' # Example dataset: Bee subjects (to be allocated to treatments, then to sites within treatments)
#' example_bee_data <- data.frame(
#'   ID = as.factor(as.character(seq(1, 100, 1))),
#'   Bee_count = rnorm(100, mean = 1000, sd = 200),
#'   Colony_weight = rnorm(100, mean = 500, sd = 100)
#' )
#'
#' # Example dataset: Sites (groups) – with site-level covariates and a blocking variable
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
#' allocated_bee_data <- allocate_treatments(
#'   data = example_bee_data,
#'   treatments = treatments,
#'   covariates = bee_covariates,
#'   ensure_equal_n_replicates = TRUE,
#'   keep_excluded_data = FALSE
#' )
#' nrow(allocated_bee_data)
#'
#' # 2) Allocate sites to treatments (balancing site quality and matching within irrigation status)
#' allocated_site_data <- allocate_treatments(
#'   data = example_site_data,
#'   treatments = treatments,
#'   covariates = "Field_quality",
#'   ensure_equal_n_replicates = TRUE,
#'   keep_excluded_data = FALSE,
#'   match_within = "Irrigated"
#' )
#' nrow(allocated_site_data)
#'
#' # 3) Allocate subjects to specific sites *within each treatment*
#' allocated_data <- allocate_within_treatments(
#'   data_subjects = allocated_bee_data,
#'   data_groups   = allocated_site_data,
#'   treatment_var = "Treatment",
#'   group_var     = "Site",
#'   covariates    = bee_covariates,
#'   ensure_equal_n_replicates = TRUE,
#'   keep_excluded_data = FALSE,
#'   standardize = TRUE
#' )
#' head(allocated_data)
#'
#' @seealso
#' [`allocate_treatments()`] for the underlying allocation within each treatment using anticlustering.
#'
#' @section Required Packages:
#' This function relies on:
#' - `anticlust`: for anticlustering optimization (used within `allocate_treatments()`).
#'
#' @export

allocate_within_treatments <- function(data_subjects, data_groups, treatment_var, group_var, covariates, 
                                       ensure_equal_n_replicates = TRUE, keep_excluded_data = FALSE,  
                                       objective = "variance", method = "local-maximum",
                                       repetitions = 10, match_within = NULL, standardize = TRUE) {
  
  # Helper function to validate arguments
  validate_args <- function(arg, arg_name) {
    if (missing(arg) || is.null(arg)) {
      stop(paste0("Error: The required argument '", arg_name, "' is missing or NULL. Please provide a valid input."))
    }
  }
  
  # Validate necessary arguments
  validate_args(data_subjects, "data_subjects")
  validate_args(data_groups, "data_groups")
  validate_args(group_var, "group_var")
  validate_args(covariates, "covariates")
  
  if (!treatment_var %in% names(data_subjects)) {
    stop("The treatment_var column is missing in data_subjects.")
  }
  if (!group_var %in% names(data_groups)) {
    stop("The group_var column is missing in data_groups.")
  }
  
  treatments <- unique(as.character(data_subjects[[treatment_var]]))
  
  # Ensure correct data types
  data_subjects[[treatment_var]] <- as.character(data_subjects[[treatment_var]])
  data_groups[[group_var]] <- as.factor(as.character(data_groups[[group_var]]))
  data_groups[[treatment_var]] <- as.character(data_groups[[treatment_var]])
  
  # Function to allocate subjects to groups within a single treatment
  allocate_to_groups_within_treatment <- function(treatment) {
    sub_group_data <- data_groups[data_groups[[treatment_var]] == treatment, , drop = FALSE]
    groups <- unique(as.character(sub_group_data[[group_var]]))
    
    sub_subject_data <- data_subjects[data_subjects[[treatment_var]] == treatment, , drop = FALSE]
    
    allocated_sub_subject_data <- allocate_treatments(
      data = sub_subject_data,
      treatments = groups,
      covariates = covariates,
      ensure_equal_n_replicates = ensure_equal_n_replicates,
      keep_excluded_data = keep_excluded_data,
      objective = objective, 
      method = method,
      repetitions = repetitions,
      match_within = match_within,
      standardize = standardize
    ) 
    
    # Rename the column that was overwritten by group assignment
    colnames(allocated_sub_subject_data)[colnames(allocated_sub_subject_data) == treatment_var] <- group_var
    
    # Add the original treatment column back
    allocated_sub_subject_data[[treatment_var]] <- treatment
    
    # Merge subject and group data for this treatment
    merged_data <- merge(sub_group_data, allocated_sub_subject_data, all = TRUE)
    
    return(merged_data)
  }
  
  # Apply the allocation for each treatment
  allocated_data_list <- lapply(treatments, allocate_to_groups_within_treatment)
  allocated_data <- do.call(rbind, allocated_data_list)
  
  return(allocated_data)
}