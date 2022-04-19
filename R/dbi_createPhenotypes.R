#' dbi_createPhenotypes
#' 
#' @details 
#' Creates a phenotype table from id, ICD9CM, ICD10CM (or phecode, etc), data.
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
#' @param aggregate.fun Aggregate function for duplicated phenotypes (phecodes, etc) in an individual. The default supports will 
#' use \code{sum} for numeric values, otherwise it will count the distinct values, eg, for dates.
#' @param vocabulary.map Map between supplied vocabularies and phecodes. Allows for custom phecode maps. By default uses 
#' \code{\link[PheWAS:phecode_map]{PheWAS::phecode_map}}, which supports ICD9CM (v1.2) and ICD10CM (beta-2018). The package also 
#' includes the ICD10 beta map (\code{\link[PheWAS:phecode_map_icd10]{PheWAS::phecode_map_icd10}}), which can be used in this parameter.
#' @param rollup.map Map between phecodes and all codes that they expand to, eg parent codes. By default uses the PheWAS::phecode_rollup_map.
#' @param exclusion.map Map between phecodes and their exclusions. By default uses the PheWAS::phecode_exclude.
#'
#' @importFrom DBI collect
#' @importFrom dplyr count distinct group_by inner_join left_join pull rename summarise transmute union_all
#' @importFrom magrittr %>% extract2
#' 
#' @return A DBI tbl connection. The first column contains the supplied id for each individual (preserving the name of the original column). 
#' The following columns contain the phecode and case_status of all present phewas codes. They contain T/F/NA for case/control/exclude or 
#' continuous/NA if min.code.count was NA.
#' @export
#'
#' @examples

dbi_createPhenotypes <- function(id.vocab.code.index, 
                                 min.code.count=2, 
                                 add.phecode.exclusions=T, 
                                 translate=T, 
                                 id.sex,
                                 full.population.ids=id.vocab.code.index %>% distinct(!!as.name(colnames(id.vocab.code.index)[1])),
                                 aggregate.fun=PheWAS:::default_code_agg,
                                 vocabulary.map=PheWAS::phecode_map,
                                 rollup.map=PheWAS::phecode_rollup_map,
                                 exclusion.map=PheWAS::phecode_exclude
                                 )
  {
  ## Obtain ID column name
  id.name=colnames(id.vocab.code.index)[1]
  ## Warn if id.sex information is not provided.
    if( missing(id.sex) ) { warning("It is recommended to provide id.sex information to help address spurious sex-specific associations.") }
    
    if( !translate ) {
      # Warn about exclusions if input is not translated and not phecodes. Same with id.sex
      # if(add.phecode.exclusions & sum(tolower(id.vocab.code.index[[2]])=='phecode')!=nrow(id.vocab.code.index)){stop("Codes are not translated and vocab is not 'phecode' for every row, but exclusions are to be applied. Ensure that the code column has only phecodes or disable add.phecode.exclusions for accurate results.")}
      if(add.phecode.exclusions & id.vocab.code.index %>% count(code) %>% pull(n) %>% sum() !=  id.vocab.code.index %>% count() %>% pull(n)) {
        stop("Codes are not translated and vocab is not 'phecode' for every row, but exclusions are to be applied. Ensure that the code column has only phecodes or disable add.phecode.exclusions for accurate results.")
      }
      # if(!missing(id.sex) & sum(tolower(id.vocab.code.index[[2]])=='phecode')!=nrow(id.vocab.code.index)){stop("Codes are not translated and vocab is not 'phecode' for every row, but id.sex is supplied for sex-based exclusions. Ensure that the code column has only phecodes or omit id.sex for accurate results.")}
      if( !missing(id.sex )  & id.vocab.code.index %>% count(code) %>% pull(n) %>% sum() !=  id.vocab.code.index %>% count() %>% pull(n)) {
        stop("Codes are not translated and vocab is not 'phecode' for every row, but id.sex is supplied for sex-based exclusions. Ensure that the code column has only phecodes or omit id.sex for accurate results.")
      }
      # phemapped=tbl_df(data.frame(id=id.vocab.code.index[[1]],code=id.vocab.code.index[[3]],index=id.vocab.code.index[[4]],stringsAsFactors = F))
      phemapped=id.vocab.code.index
    } else {
      #check to make sure numeric codes were not passed in
      # if(!class(id.vocab.code.index[[3]]) %in% c("character","factor")) {stop("Please ensure character or factor code representation. Some vocabularies, eg ICD9CM, require strings to be represented accurately: E.G.: 250, 250.0, and 250.00 are different codes and necessitate string representation")}
      if( !id.vocab.code.index %>% head() %>% collect() %>% lapply(type_sum) %>% magrittr::extract2(3) %in% c("chr","fct") ) {
        stop("Please ensure character or factor code representation. Some vocabularies, eg ICD9CM, require strings to be represented accurately: E.G.: 250, 250.0, and 250.00 are different codes and necessitate string representation")
      }
      id.vocab.code.index_cols <- colnames(id.vocab.code.index)
      id.vocab.code.index <- id.vocab.code.index %>%
        rename('id' = id.vocab.code.index_cols[[1]], 'vocabulary_id' = id.vocab.code.index_cols[[2]], 'code' = id.vocab.code.index_cols[[3]], 'index' = id.vocab.code.index_cols[[4]])
      message("Mapping codes to phecodes...")
      phemapped <- dbi_mapCodesToPhecodes(input = id.vocab.code.index, vocabulary.map=vocabulary.map, rollup.map=rollup.map) %>%
        transmute(id, code=phecode, index)
    }
    
    message("Aggregating codes...")
    # phecode=ungroup(summarize(group_by(phemapped,id,code),count=aggregate.fun(index)))
    if( is.numeric(phemapped %>% head() %>% collect() %>% pull(index)) ) {
      phecode <- phemapped %>%
        group_by(id,code) %>%
        summarise(count = sum(index, na.rm = T), .groups = 'drop')
    } else {
      phecode <- phemapped %>%
        group_by(id,code) %>%
        summarise(count = length(distinct(index)), .groups = 'drop' )
    }
    
    # phecode=phecode[phecode$count>0,]
    phecode <- phecode %>%
      filter(count > 0)
    
    #Check exclusions, and add them to the list
    if( add.phecode.exclusions ) {
      message("Mapping exclusions...")
      exclusions <- phecode %>%
        rename(exclusion_criteria=code) %>%
        inner_join(exclusion.map, by = "exclusion_criteria")
      exclusions <- exclusions %>%
        transmute(id, code, count= -1) %>%
        distinct()
      phecode <- union_all(phecode, exclusions)
    }
    
    #If there is request for a min code count, adjust counts to -1 if needed
    # if(!is.na(min.code.count) & ( max(!is.na(phecode$count) & phecode$count<min.code.count) ) ) {
    #   phecode[!is.na(phecode$count)&phecode$count<min.code.count,]$count=-1
    # }
    if(!is.na(min.code.count)) {
      phecode <- phecode %>%
        mutate(count = case_when(!is.na(count) & count < min.code.count ~ -1,
                                 TRUE ~ count)
        )
    }
    
    if( !is.na(min.code.count) | add.phecode.exclusions ) {
      message("Coalescing exclusions and min.code.count as applicable...")
      # phecode=ungroup(summarize(group_by(phecode,id,code),count=max(count)))
      phecode <- phecode %>%
        group_by(id, code) %>%
        summarise(count = max(count, na.rm = T), .groups = 'drop') ##???????????
    }
    
    # message("Reshaping data...")
    # # phens=spread(phecode,code,count,fill=0)
    # phens <- phecode %>%
    #   mutate(code = str_replace(code, '\\.','_') ) %>% ## BQ Col Restrictions
    #   pivot_wider(names_from = code, names_prefix = 'code_', names_sort = T, values_from = count) %>% ## values_fill doesn't work right
    #   mutate(across(starts_with('code_'), ~case_when(is.na(.x) ~ 0,
    #                                                  TRUE ~ .x)
    #                 )
    #          ) %>%
    #   compute()
    
    rlang::inform(rlang::format_error_bullets(c('i' = glue::glue('Long data is {crayon::red("l")}{crayon::green("o")}{crayon::yellow("o")}{crayon::blue("o")}{crayon::magenta("o")}{crayon::cyan("o")}{crayon::white("o")}{crayon::silver("o")}{crayon::red("o")}{crayon::green("o")}{crayon::yellow("o")}{crayon::blue("o")}{crayon::magenta("o")}{crayon::cyan("o")}{crayon::white("n")}{crayon::silver("g")}...!'))))
    ## Add Controls
    distinct_phe_population_ids <- phecode %>%
      distinct(id)
    distinct_phe_population <- phecode %>%
      distinct(code)
    control_setup <- distinct_phe_population_ids %>%
      full_join(distinct_phe_population, by = character()) %>%
      mutate(count_2 = NA)
    phecode <- phecode %>%
      full_join(control_setup)
    
    # #Set exclusions to NA, preserving IDs just in case one is -1
    # # tmp_id=phens[,1]
    # # phens[phens==-1]=NA
    # # phens[,1]=tmp_id
    # phens <- phens %>%
    #   mutate(across(starts_with('code_'), ~case_when(.x == -1 ~ NA,
    #                                                  TRUE ~ .x)
    #                 )
    #          )
    ## Exclusions to NA, controls to 0
    phens <- phecode %>%
      mutate(count = case_when(is.na(count) ~ 0,
                               count == -1 ~ NA,
                               TRUE ~ count)
      ) %>%
      select(-count_2)
    
    #Add in inds present in input or the full population list, but without mapped phecodes
    # missing_ids=setdiff(full.population.ids, phens %>% select("id") )
    missing_ids <- full.population.ids %>%
      anti_join(phens %>% select("id") )
    
    if( count(missing_ids) %>% pull(n) > 0 ) {
      # empty_record=phens[1,-1]
      # empty_record[]=0
      # empty_record <- phens %>%
      #   head(1) %>%
      #   select(-1) %>%
      #   mutate(across(everything(), ~0))
      # missing_ids <- union_all(missing_ids, empty_record) %>%
      #   filter(!is.na(id)) %>%
      #   mutate(across(starts_with('code'), ~ 0))
      empty_records <- missing_ids %>%
        full_join(distinct_phe_population, by = character()) %>%
        mutate(count = 0)
      # phens=rbind(phens,data.frame(id=missing_ids,empty_record,check.names=F))
      phens <- union_all(phens, empty_records)
    }
    
    #Change to logical if there is a min code count
    # if(!is.na(min.code.count)) {phens[,-1]=phens[,-1]>0}
    if( !is.na(min.code.count) ) {
      # phens <- phens %>%
      #   mutate(across(starts_with('code_'), ~case_when(is.na(.x) ~ NA,
      #                                                  .x > 0 ~ TRUE,
      #                                                  TRUE ~ FALSE
      #                                                  )
      #                 )
      #          )
      phens <- phens %>%
        mutate(count = case_when(is.na(count) ~ NA,
                                 count > 0 ~ TRUE,
                                 TRUE ~ FALSE)
        )
    }
    
    #If there are sex restrictions, set them to NA
    if( !missing(id.sex) ) {
      phens <- dbi_restrictPhecodesBySex(phens,id.sex)
    }
    
    #Limit to full population ids
    # phens = filter(phens, id %in% full.population.ids)
    phens <- full.population.ids %>%
      left_join(phens)
    
    #Rename the ID column to the input ID column name
    # names(phens)[1]=id.name
    phens <- phens %>%
      rename(!!id.name := 1)
    
    phens
  }
