library(anticlust)

#' Allocate Treatments to Subjects
#'
#' This function assigns subjects to treatment groups while balancing specified covariates and ensuring replicates per group.
#'
#' @param data A data frame containing the subjects to be allocated, including columns for the covariates.
#' @param treatments A vector specifying the treatment group labels.
#' @param treatment_var A string specifying the name of the treatment column to be created. Defaults to `"Treatment"`. 
#' @param covariates A vector of column names in `data` representing covariates to balance.
#' @param n_replicates Integer. Specifies the number of replicates per treatment group. If `NULL`, either all subjects are allocated or the highest possible equal number per treatment is allocated depending on `ensure_equal_n_replicates`. Defaults to `NULL`.
#' @param ensure_equal_n_replicates Logical. If `TRUE`, ensures equal numbers of replicates across treatment groups. Defaults to `TRUE`. If `FALSE` all subjects are allocated. 
#' @param keep_excluded_data Logical. If `TRUE`, retains subjects that were not allocated within the specified constraints. Defaults to `FALSE`.
#' @param objective A string specifying the objective function for anticlustering. Options include `"variance"` (default), `"diversity"`, `"average-diversity"`, `"kplus"`, and `"dispersion"`. 
#' @param method. A string specifying the optimization method for anticlustering. Options include `"local-maximum"` (default), `"exchange"`, `"brusco"`, `"ilp"`, and `"2PML"`.
#' @param repetitions Integer. Specifies the number of times the optimization is repeated when using heuristic methods (`"exchange"`, `"local-maximum"`, `"brusco"`, or `"2PML"`). The best solution is selected. Defaults to `10`.
#' @param match_within A column name in `data` (optional). Specifies a variable within which matching should occur, ensuring that subjects are grouped within subsets defined by this variable. Defaults to `NULL`.
#' @param standardize Logical. If `TRUE`, covariates are standardized via `scale()` before optimization starts. Defaults to `TRUE`.
#' 
#'
#' @return A data frame with subjects assigned to treatments. The column name for treatment allocation is specified by `treatment_var`. If `keep_excluded_data = TRUE`, excluded subjects are included in the output with `NA` in the treatment column.
#'
#' @details
#' The function balances covariates among treatment groups by creating groups of similar individuals based on the specified covariates. The treatments are then assigned to these groups in a way that minimizes variance within groups.
#'
#' - If `ensure_equal_n_replicates = TRUE`, the number of replicates per treatment is enforced.
#' - If `keep_excluded_data = TRUE`, subjects that cannot be allocated to treatments are retained in the output.
#' - Covariates can be scaled to standardize their ranges if `standardize = TRUE`.
#' - If `match_within` is specified, subjects are grouped within levels of the specified variable, ensuring allocations respect the structure defined by this variable.
#'
#' @section Validations:
#' The function validates the following conditions:
#' - All required arguments (`data`, `treatments`, `covariates`) are provided.
#' - Covariates exist in the `data` frame and are numeric.
#' - Covariates do not have constant values across all rows.
#'
#' @examples
#' # Example dataset: Bee data
#' set.seed(123)
#' example_bee_data <- data.frame(
#'   ID = as.factor(as.character(seq(1, 100, 1))),
#'   Bee_count = rnorm(100, mean = 1000, sd = 200),
#'   Colony_weight = rnorm(100, mean = 500, sd = 100))
#' 
#' # Example dataset: Site data
#' example_site_data <- data.frame(
#'   Site = as.factor(as.character(seq(1, 16, 1))),
#'   Organic = as.factor(as.character(c(rep("yes", 4), rep("no", 12)))), 
#'   Field_quality = rnorm(16, mean = 3, sd = 0.5))
#' 
#' treatments <- c("Control", "Pesticide")
#' bee_covariates <- c("Bee_count", "Colony_weight")
#' 
#' # Allocate bee data to treatments
#' allocated_bee_data <- allocate_treatments(
#'   example_bee_data, 
#'   treatments = treatments, 
#'   covariates = bee_covariates,
#'   ensure_equal_n_replicates = TRUE, 
#'   keep_excluded_data = FALSE)
#' 
#' # Allocate site data to treatments while matching within "Organic"
#' allocated_site_data <- allocate_treatments(
#'   example_site_data, 
#'   treatments = treatments, 
#'   covariates = "Field_quality", 
#'   ensure_equal_n_replicates = TRUE, 
#'   keep_excluded_data = FALSE,
#'   match_within = "Organic")
#'
#' @section Required Packages:
#' This function requires the following packages:
#' - `anticlust`: for anticlustering optimization.
#'
#' @export
#' 

allocate_treatments <- function(data, treatments, treatment_var = "Treatment",
                                covariates, n_replicates = NULL, 
                                ensure_equal_n_replicates = TRUE, 
                                keep_excluded_data = FALSE,
                                objective = "variance", 
                                method = "local-maximum",
                                repetitions = 10, 
                                match_within = NULL,
                                standardize = TRUE) {
  
  # Use validate_args for argument validation
  validate_args <- function(arg, arg_name) {
    if (missing(arg) || is.null(arg)) {
      stop(paste0("Error: The required argument '", arg_name, "' is missing or NULL. Please provide a valid input."))
    }
  }
  
  # Validate required arguments
  validate_args(data, "data")
  validate_args(treatments, "treatments")
  validate_args(covariates, "covariates")
  
  # Check that specified covariates are present in the data
  missing_covariates <- setdiff(covariates, colnames(data))
  if (length(missing_covariates) > 0) {
    stop("The following covariates are missing in the data: ", paste(missing_covariates, collapse = ", "), ". Ensure these covariates are included in the dataset.")
  }
  
  # Check if any covariate is not numeric 
  non_numeric_covariates <- covariates[!sapply(data[, covariates, drop = FALSE], is.numeric)]
  if (length(non_numeric_covariates) > 0) {
    stop("The following covariates are not numeric: ", paste(non_numeric_covariates, collapse = ", "), ". Convert these covariates to numeric before running the function.")
  }
  
  # Check if any covariate takes the same value for all rows
  constant_covariates <- covariates[sapply(data[, covariates, drop = FALSE], function(x) length(unique(x[!is.na(x)])) == 1)]
  if (length(constant_covariates) > 0) {
    stop("The following covariates take the same value for all observations: ", paste(constant_covariates, collapse = ", "), ". Remove these covariates or provide covariates with variability.")
  }
  
  # Validate match_within
  if (!is.null(match_within)) {
    if (!match_within %in% colnames(data)) {
      stop("The column specified in match_within ('", match_within, "') is not present in the data frame. Please ensure it is a valid column name.")
    }
    match_within_vector <- data[[match_within]]
  } else {
    match_within_vector <- NULL
  }
  
  # Shuffle treatments for random assignment
  treatments_shuffled <- sample(treatments)
  
  # Calculate the number of replicates per treatment if not provided
  n_treatments <- length(treatments)
  
  if (is.null(n_replicates)) {
    n_replicates <- floor(nrow(data) / n_treatments)
  } else {
    if (!ensure_equal_n_replicates) {
      ensure_equal_n_replicates <- TRUE
      warning("'ensure_equal_n_replicates' was changed to TRUE as 'n_replicates' was manually specified.")
    }
  }
  
  # Scale covariates if specified
  covariate_data <- if (standardize) {
    scale(data[, covariates, drop = FALSE])
  } else {
    data[, covariates, drop = FALSE]
  }
  
  if (ensure_equal_n_replicates | n_treatments * n_replicates == nrow(data)) {
    # Create groups with similar individuals based on covariates
    data$set_similar <- matching(
      covariate_data, 
      p = n_treatments,
      match_within = match_within_vector
    )
    
    # Include individuals based on the number of replicates
    included <- data[!is.na(data$set_similar) & data$set_similar <= n_replicates, ]
    
    # Assign individuals of the same set to different treatment groups
    included[[treatment_var]] <- anticlustering(
      included[, covariates, drop = FALSE], 
      K = n_treatments,
      objective = objective, 
      method = method,
      categories = included$set_similar,
      repetitions = repetitions
    ) 
    
    included$set_similar <- NULL
    
    if (keep_excluded_data) {
      # Extract data that was excluded from the experiment
      excluded <- data[data$set_similar > n_replicates | is.na(data$set_similar), ]
      excluded$set_similar <- NULL
      excluded[[treatment_var]] <- NA
      
      # Combine included and excluded data
      output <- rbind(included, excluded)
    } else {
      output <- included
    }
  } else {
    # Directly assign treatments without ensuring equal replicates
    data[[treatment_var]] <- anticlustering(
      covariate_data, 
      K = n_treatments,
      objective = objective, 
      method = method,
      categories = match_within_vector,
      repetitions = repetitions
    )
    output <- data
  }
  
  # Map treatment values to shuffled treatments
  treatment_map <- setNames(treatments_shuffled, as.character(1:n_treatments))
  output[[treatment_var]] <- unname(treatment_map[as.character(output[[treatment_var]])])
  
  # Convert treatment variable to factor
  output[[treatment_var]] <- factor(output[[treatment_var]], levels = treatments)
  
  # Reorder columns to have treatment_var first
  output <- output[, c(treatment_var, setdiff(names(output), treatment_var))]
  
  return(output)
}
