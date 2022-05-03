#' db_createPhenotypes
#' 
#' @details 
#' A database optimized version of (\code{\link[PheWAS:createPhenotypes]{PheWAS::createPhenotypes}}). Creates a phenotype table 
#' from id, ICD9CM, ICD10CM (or phecode, etc), data.
#' 
#' By default, this function returns a long format data frame with boolean phenotypes suitable for PheWAS analysis. Specifying 
#' a \code{min.code.count=NA} will permit continuous code count phenotypes.
#' 
#' The default exclusions can be skipped with \code{add.exclusions=F}. In conjuntion with \code{translate=F} 
#' (and optionally adjusting \code{min.code.count} and \code{aggregate.fun}), one can use this function as a simple reshaping wrapper.
#' 
#' @description 
#' This function takes a data frame with four columns: id, vocabulary_id, code, and index. It returns a wide table with phecodes as 
#' TRUE/FALSE/NA. It can optionally use the PheWAS exclusion criteria.
#' 
#' @keywords utilities
#'
#' @param id.vocab.code.index Data frame with four columns of information: id, vocabulary_id, code, and index. The id and index 
#' columns can have other names, but id must be consistent among input files. The vocabulary_id and code must match up with the 
#' vocabulary.map file. The default supports the vocabularies "ICD9CM" and "ICD10CM". Code contains the raw code value.
#' @param min.code.count The minimum code count to be considered a case. NA results in a continuous output.
#' @param add.phecode.exclusions Apply PheWAS exclusions to phecodes.
#' @param translate Should the input be translated to phecodes? Defaults to TRUE. Generally recommended, though can be skipped 
#' if phecodes are provided.
#' @param id.sex If supplied, restrict the phecodes by sex. This should be a data frame with the first column being the id and 
#' the second the sex, "M" or "F", of the individual. Individuals with any other specification will have all sex specific phenotypes 
#' set to NA.
#' @param full.population.ids List of IDs in the "complete" population. This allows for individuals with no observed codes to have 
#' appropriate "control" status, eg 0s or FALSE in every field.
#' @param vocabulary.map Map between supplied vocabularies and phecodes. Allows for custom phecode maps. By default uses 
#' \code{\link[PheWAS:phecode_map]{PheWAS::phecode_map}}, which supports ICD9CM (v1.2) and ICD10CM (beta-2018). The package also 
#' includes the ICD10 beta map (\code{\link[PheWAS:phecode_map_icd10]{PheWAS::phecode_map_icd10}}), which can be used in this parameter.
#' @param rollup.map Map between phecodes and all codes that they expand to, eg parent codes. By default uses the PheWAS::phecode_rollup_map.
#' @param exclusion.map Map between phecodes and their exclusions. By default uses the PheWAS::phecode_exclude.
#'
#' @importFrom dplyr collect compute count distinct group_by filter inner_join left_join pull rename select summarise transmute type_sum union_all
#' @importFrom magrittr %>% extract2
#' @importFrom rlang .data abort format_error_bullets inform warn
#' @importFrom utils head
#' 
#' @return A DBI tbl connection. The first column contains the supplied id for each individual (preserving the name of the original column). 
#' The following columns contain the phecode and case_status of all present phewas codes. They contain T/F/NA for case/control/exclude or 
#' continuous/NA if min.code.count was NA.
#' @export
#'
#' @examples
#' \dontrun{
#' library(PheWAS)
#' library(DBI)
#' library(dplyr)
#' library(RSQLite)
#' library(tidyr)
#' ## Create DBI connection object to target database
#' con <- dbConnect(RSQLite::SQLite(), ":memory:")
#' ## Generate Example Data
#' set.seed(2020)
#' example_data <- PheWAS::generateExample(n = 50)
#' ## Create Example Database
#' dbWriteTable(con, 'id.vocab.code.count', example_data$id.vocab.code.count)
#' dbWriteTable(con, 'id.sex', example_data$id.sex)
#' ## Upload PheWAS package data
#' dbWriteTable(con, 'phecode_map', PheWAS::phecode_map)
#' dbWriteTable(con, 'phecode_rollup_map', PheWAS::phecode_rollup_map)
#' dbWriteTable(con, 'phecode_exclude', PheWAS::phecode_exclude)
#' ## Define tbl connections to example data
#' db_icd_summary <- tbl(con, 'id.vocab.code.count')
#' db_id_sex <- tbl(con, 'id.sex')
#' ## Define tbl connections to PheWAS package Datasets
#' db_phecode_map <- tbl(con, 'phecode_map')
#' db_phecode_rollup_map <- tbl(con, 'phecode_rollup_map')
#' db_phecode_exclude <- tbl(con, 'phecode_exclude')
#' ## Create Phenotypes in Database
#' phenotype_table <- db_createPhenotypes(id.vocab.code.index = db_icd_summary, 
#'                                        min.code.count = 2, 
#'                                        add.phecode.exclusions = T, 
#'                                        id.sex = db_id_sex,
#'                                        vocabulary.map = db_phecode_map, 
#'                                        rollup.map = db_phecode_rollup_map, 
#'                                        exclusion.map = db_phecode_exclude
#' )
#' local_phenotype_table <- phenotype_table %>% 
#'   collect()
#' local_phenotype_table %>% 
#'   pivot_wider(names_from = code, values_from = case_status)
#' }

db_createPhenotypes <- function(id.vocab.code.index, 
                                min.code.count=2, 
                                add.phecode.exclusions=T, 
                                translate=T, 
                                id.sex,
                                full.population.ids=id.vocab.code.index %>% distinct(!!as.name(colnames(id.vocab.code.index)[1])),
                                vocabulary.map=PheWAS::phecode_map,
                                rollup.map=PheWAS::phecode_rollup_map,
                                exclusion.map=PheWAS::phecode_exclude
)
{
  ## Setup ----
  ### Obtain ID column name ----
  id.name=colnames(id.vocab.code.index)[1]
  
  ### Warn if id.sex is not provided ----
  if( missing(id.sex) ) { warn(format_error_bullets(c(' '='','!'="It is recommended to provide id.sex information to help address spurious sex-specific associations.")) ) }
  
  ## Translate to Phecode ----  
  if( !translate ) {
    ### Warn about exclusions if input is not translated and not phecodes. 
    if(add.phecode.exclusions & id.vocab.code.index %>% count(.data$code) %>% pull(.data$n) %>% sum() !=  id.vocab.code.index %>% count() %>% pull(n)) {
      abort("Codes are not translated and vocab is not 'phecode' for every row, but exclusions are to be applied. Ensure that the code column has only phecodes or disable add.phecode.exclusions for accurate results.")
    }
    ### Warn about exclusions if input is not translated and not phecodes.
    if( !missing(id.sex)  & id.vocab.code.index %>% count(.data$code) %>% pull(.data$n) %>% sum() !=  id.vocab.code.index %>% count() %>% pull(.data$n)) {
      abort("Codes are not translated and vocab is not 'phecode' for every row, but id.sex is supplied for sex-based exclusions. Ensure that the code column has only phecodes or omit id.sex for accurate results.")
    }
    phemapped=id.vocab.code.index
  } else {
    ### check to make sure numeric codes were not passed in
    if( !id.vocab.code.index %>% head() %>% collect() %>% lapply(type_sum) %>% magrittr::extract2(3) %in% c("chr","fct") ) {
      abort("Please ensure character or factor code representation. Some vocabularies, eg ICD9CM, require strings to be represented accurately: E.G.: 250, 250.0, and 250.00 are different codes and necessitate string representation")
    }
    id.vocab.code.index_cols <- colnames(id.vocab.code.index)
    id.vocab.code.index <- id.vocab.code.index %>%
      rename('id' = id.vocab.code.index_cols[[1]], 
             'vocabulary_id' = id.vocab.code.index_cols[[2]], 
             'code' = id.vocab.code.index_cols[[3]], 
             'index' = id.vocab.code.index_cols[[4]]
      )
    inform(format_error_bullets(c('i' = "Mapping codes to phecodes...")) )
    phemapped <- db_mapCodesToPhecodes(input = id.vocab.code.index, vocabulary.map=vocabulary.map, rollup.map=rollup.map) %>%
      transmute(.data$id, code=.data$phecode, .data$index)
  }
  
  ## Aggregate Codes ----  
  inform(format_error_bullets(c('i' = "Aggregating codes...")) )
  if( is.numeric(phemapped %>% head() %>% collect() %>% pull(.data$index)) ) {
    phecode <- phemapped %>%
      group_by(.data$id,.data$code) %>%
      summarise(count = sum(.data$index, na.rm = T), .groups = 'drop')
  } else {
    phecode <- phemapped %>%
      group_by(.data$id,.data$code) %>%
      summarise(count = length(distinct(.data$index)), .groups = 'drop' )
  }
  phecode <- phecode %>%
    filter(.data$count > 0)
  
  ## Map Exclusions ----
  ### Check exclusions, and add them to the list
  if( add.phecode.exclusions ) {
    inform(format_error_bullets(c('i'= "Mapping exclusions...")) )
    exclusions <- phecode %>%
      rename(exclusion_criteria=.data$code) %>%
      inner_join(exclusion.map, by = "exclusion_criteria")
    exclusions <- exclusions %>%
      transmute(id, code, count = -1) %>%
      distinct()
    phecode <- union_all(phecode, exclusions)%>% 
      dplyr::compute()
  }
  
  ### If there is request for a min code count, adjust counts to -1 if needed
  if(!is.na(min.code.count)) {
    phecode <- phecode %>%
      mutate(count = case_when(!is.na(.data$count) & .data$count < min.code.count ~ -1,
                               TRUE ~ .data$count)
      )%>% 
      dplyr::compute()
  }
  
  if( !is.na(min.code.count) | add.phecode.exclusions ) {
    inform(format_error_bullets(c('i' = "Coalescing exclusions and min.code.count as applicable...")) )
    phecode <- phecode %>%
      group_by(.data$id, .data$code) %>%
      summarise(count = max(.data$count, na.rm = T), .groups = 'drop') %>% ##???????????
      dplyr::compute()
  }
  
  ## Add Population Controls ----  
  # rlang::inform(rlang::format_error_bullets(c('i' = glue::glue('Long data is {crayon::red("l")}{crayon::green("o")}{crayon::yellow("o")}{crayon::blue("o")}{crayon::magenta("o")}{crayon::cyan("o")}{crayon::white("o")}{crayon::silver("o")}{crayon::red("o")}{crayon::green("o")}{crayon::yellow("o")}{crayon::blue("o")}{crayon::magenta("o")}{crayon::cyan("o")}{crayon::white("n")}{crayon::silver("g")}...!'))))
  inform(format_error_bullets(c('i' = 'Add population controls...')) )
  distinct_phe_population_ids <- phecode %>%
    distinct(.data$id)
  distinct_phe_population <- phecode %>%
    distinct(.data$code)
  control_setup <- distinct_phe_population_ids %>%
    full_join(distinct_phe_population, by = character()) %>%
    mutate(count_2 = NA)
  phecode <- phecode %>%
    full_join(control_setup, by = c("id", "code"))
  
  ### Set exclusions to NA, controls to 0, preserving IDs just in case one is -1 ----
  phens <- phecode %>%
    mutate(count = case_when(is.na(.data$count) ~ 0,
                             .data$count == -1 ~ NA_real_,
                             TRUE ~ .data$count
    )
    ) %>%
    select(-.data$count_2)
  
  ### Add back in ids present in input or the full population list, but without mapped phecodes ----
  missing_ids <- full.population.ids %>%
    anti_join(phens %>% select("id"), by = "id" )
  
  if( count(missing_ids) %>% pull(.data$n) > 0 ) {
    empty_records <- missing_ids %>%
      full_join(distinct_phe_population, by = character()) %>%
      mutate(count = 0)
    phens <- union_all(phens, empty_records)
  }
  
  ### Change to logical if there is a min code count ----
  if( !is.na(min.code.count) ) {
    phens <- phens %>%
      mutate(count = case_when(is.na(.data$count) ~ NA,
                               .data$count > 0 ~ TRUE,
                               TRUE ~ FALSE
      )
      )
  }
  
  ## Gender Restrictions ----
  ### If there are sex restrictions, set them to NA
  inform(format_error_bullets(c('i' = 'Mapping gender restrictions...')) )
  if( !missing(id.sex) ) {
    phens <- db_restrictPhecodesBySex(phens,id.sex)
  }
  
  ## Limit to full population ids ----
  phens <- full.population.ids %>%
    left_join(phens, by = "id")
  
  ## Rename the ID column ----
  ### Use the original input ID column name 
  phens <- phens %>%
    rename(!!id.name := 1)
  
  ## Return the output ----  
  inform(format_error_bullets(c('v' = 'Complete!')) )
  phens
}
