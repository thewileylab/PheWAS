#' db_restrictPhecodesBySex
#' 
#' @details 
#' Add PheWAS code descriptions to existing data.
#' 
#' @description \code{restrictPhecodesBySex} alters a table for PheWAS with phecodes, as from \code{\link[PheWAS:createPhewasTable]{createPhewasTable}}, 
#' to exclude individuals with non-applicable genders from certain phenotypes.
#' 
#' @keywords utilities
#'
#' @param phenotypes The PheWAS table to have restrictions applied. The first column should be the id.
#' @param id.sex A `tbl()` source connection with the first column being the id and the second the gender, "M" or "F", of the 
#' individual. Individuals with any other specification will have all gender specific phenotypes set to NA.
#'
#' @importFrom dplyr case_when collect distinct filter left_join mutate rename select
#' @importFrom rlang .data
#' 
#' @return The \code{phenotypes} `tbl()` source connection with NA values for individuals that do not match the gender for 
#' gender-specific codes.
#' @export
#'

db_restrictPhecodesBySex <- function(phenotypes,
                                     id.sex) 
  {
  ## Add gender information to data ---- 
  data <- phenotypes %>%
    left_join(id.sex, by = "id")
  
  ## Identify gender restrictions found in the phenotype data ----
  current_phens <- phenotypes %>%
    distinct(.data$code) %>%
    rename(phecode = .data$code) %>%
    collect()
  current_gender_restriction <- PheWAS::gender_restriction %>%
    left_join(current_phens, by = "phecode") %>%
    filter(.data$phecode %in% current_phens$phecode)
  
  ## Get male and female-only phenotypes ----
  male_only <- current_gender_restriction %>%
    filter(.data$male_only == TRUE) %>%
    select(.data$phecode)
  female_only <- current_gender_restriction %>%
    filter(.data$female_only == TRUE) %>%
    select(.data$phecode)
  
  ## Set row column matches to NA where ids of a gender meet restricted phenotypes ----
  data <- data %>%
    mutate(case_status = case_when(.data$code %in% !!female_only$phecode & sex != 'F' ~ NA,
                                   .data$code %in% !!male_only$phecode & sex != 'M' ~ NA,
                                   TRUE ~ .data$count)
    )
  
  ## Return everything, sans gender column ----
  data <- data %>% select(-.data$sex, -.data$count)
  data
}
