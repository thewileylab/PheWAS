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
#' @param id.sex A DBI tbl connection with the first column being the id and the second the gender, "M" or "F", of the 
#' individual. Individuals with any other specification will have all gender specific phenotypes set to NA.
#'
#' @importFrom dplyr case_when collect distinct filter left_join mutate rename select
#' @importFrom rlang .data
#' 
#' @return The \code{phenotypes} tbl connection with NA values for individuals that do not match the gender for 
#' gender-specific codes.
#' @export
#'
#' @examples

db_restrictPhecodesBySex <- function(phenotypes,
                                      id.sex) 
  {
  ## Add gender information to data ---- 
  data <- phenotypes %>%
    left_join(id.sex)
  
  ## Identify column containing gender ----
  # g=dim(data)[2]
  # g=length( colnames(test) )
  
  ## Identify gender restrictions found in the phenotype data ----
  # current_gender_restriction=PheWAS::gender_restriction[PheWAS::gender_restriction$phecode %in% colnames(phenotypes)[-1],]
  # current_phens <- tibble(new_phe = colnames(phenotypes)[-1]) %>%
  #   mutate(phecode = str_replace(new_phe, 'code_', ''),
  #          phecode = str_replace(phecode, '_', '\\.')
  #          )
  current_phens <- phenotypes %>%
    distinct(.data$code) %>%
    rename(phecode = .data$code) %>%
    collect()
  current_gender_restriction <- PheWAS::gender_restriction %>%
    left_join(current_phens) %>%
    filter(.data$phecode %in% current_phens$phecode)
  
  ## Get male and female-only phenotypes ----
  # male_only=current_gender_restriction[current_gender_restriction$male_only,"phecode"]
  male_only <- current_gender_restriction %>%
    filter(.data$male_only == TRUE) %>%
    select(.data$phecode)
  # female_only=current_gender_restriction[current_gender_restriction$female_only,"phecode"]
  female_only <- current_gender_restriction %>%
    filter(.data$female_only == TRUE) %>%
    select(.data$phecode)
  
  ## Set row column matches to NA where ids of a gender meet restricted phenotypes ----
  # data[!is.na(data[,g])&data[,g]!="F",female_only]=NA
  # data[!is.na(data[,g])&data[,g]!="M",male_only]=NA
  
  # data <- data %>%
  #   mutate(across(any_of(!!female_only$new_phe), ~case_when(sex != 'F' ~ NA,
  #                                                             TRUE ~ .x)),
  #          across(any_of(!!male_only$new_phe), ~case_when(sex != 'M' ~ NA,
  #                                                           TRUE ~ .x))
  #          )
  data <- data %>%
    mutate(case_status = case_when(.data$code %in% !!female_only$phecode & sex != 'F' ~ NA,
                                   .data$code %in% !!male_only$phecode & sex != 'M' ~ NA,
                                   TRUE ~ .data$count)
    )
  
  ## Return everything, sans gender column ----
  # data[,-g]
  data <- data %>% select(-.data$sex, -.data$count)
  data
}
