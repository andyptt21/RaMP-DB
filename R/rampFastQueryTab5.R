#' Function that query database to find ontology information based on 
#' the given list of analytes
#' @param analytes a vector of analytes or a analytes delimited by new line character
#' @param dbname database name for this database
#' @param host host for this database
#' @param username username for the database
#' @param conpass password for the database
#' @param sourceOrName specify the type of given data
rampFastOntoFromMeta <- function(analytes,conpass = NULL,
                         dbname = 'ramp',
                         host = 'localhost',
                         username = 'root',
                         sourceOrName = 'source'){
  if(!(sourceOrName %in% c('source','name')))
    stop("Specifiy the type of given data to 'source' or 'name'")
  
  if(is.null(conpass)) {
    stop("Please define the password for the mysql connection")
  }
  now <- proc.time()
  if(is.character(analytes)){
    if(grepl("\n",analytes)[1]){
      list_metabolite <- strsplit(analytes,"\n")
      list_metabolite <- unlist(list_metabolite)
    } else if(grepl(",",analytes)[1]){
      list_metabolite <- strsplit(analytes,"\n")
      list_metabolite <- unlist(list_metabolite)
    } else {
      list_metabolite <- analytes
    }
  } else if(is.data.frame(analytes)){
    list_metabolite <- unlist(analytes)
  }
  list_metabolite <- unique(list_metabolite)
  list_metabolite <- sapply(list_metabolite,shQuote)
  list_metabolite <- paste(list_metabolite,collapse = ",")
  
  con <- DBI::dbConnect(RMySQL::MySQL(),
                        dbname = dbname,
                        username = username,
                        host = host,
                        password = conpass)
  if(sourceOrName == 'source'){
    sql <- paste0('select * from source where sourceId in (',list_metabolite,');')
  } else if (sourceOrName == 'name'){
    sql <- paste0('select * from source where rampId in (select rampId from analytesynonym where Synonym in (',
                  list_metabolite,'));')
  }
  df <- DBI::dbGetQuery(con,sql)
  print(colnames(df))
  DBI::dbDisconnect(con)
  if(nrow(df) == 0) {
    message('No searching result since this source id 
            is not existed in source table')
    return(NULL)
  }
  rampid <- unique(df$rampId)
  rampid <- sapply(rampid,shQuote)
  rampid <- paste(rampid,collapse = ',')
  con <- DBI::dbConnect(RMySQL::MySQL(),
                        dbname = dbname,
                        username = username,
                        host = host,
                        password = conpass)
  
  sql <- paste0('select * from analytehasontology where rampCompoundId in (',
                rampid,');')
  df2 <- DBI::dbGetQuery(con,sql)
  if(nrow(df2) == 0) {
    message('No searching result because these metabolites are not linked to ontology')
    return(NULL)
  }
  DBI::dbDisconnect(con)
  rampontoid <- unique(df2$rampOntologyIdLocation)
  rampontoid <- sapply(rampontoid,shQuote)
  rampontoid <- paste(rampontoid,collapse = ',')
  sql <- paste0('select * from ontology where rampOntologyIdLocation in (',
                rampontoid,');')
  con <- DBI::dbConnect(RMySQL::MySQL(),
                        dbname = dbname,
                        username = username,
                        host = host,
                        password = conpass)
  df3 <- DBI::dbGetQuery(con,sql)
  print(colnames(df3))
  DBI::dbDisconnect(con)
  mdf <- unique(merge(df3,df2,all.x=T))
  mdf <- unique(merge(mdf,df,all.x = T,by.x = 'rampCompoundId',
                      by.y= 'rampId'))
  
  mdf <- mdf[c('sourceId','commonName.x','biofluidORcellular')]
  
  colnames(mdf) <- c('Metabolites','Ontology','Ontology_Type')
  return(mdf)
}
#' function that query database to find analytes in given ontologies 
#' @param ontology a vector of ontology or ontologies delimited by new line character
#' @param dbname a database name for the database connected
#' @param conpass a database password
#' @param host host name for the database
#' @param username user name for the database
rampFastMetaFromOnto <- function(ontology,conpass = NULL,
                                 dbname = 'ramp',
                                 host = 'localhost',
                                 username = 'root'){
  if(is.null(conpass)) {
    stop("Please define the password for the mysql connection")
  }
  now <- proc.time()
  if(is.character(ontology)){
    if(grepl("\n",ontology)[1]){
      list_ontology <- strsplit(ontology,"\n")
      list_ontology <- unlist(list_ontology)
    } else if(grepl(",",ontology)[1]){
      list_ontology <- strsplit(ontology,"\n")
      list_ontology <- unlist(list_ontology)
    } else {
      list_ontology <- ontology
    }
  } else if(is.data.frame(ontology)){
    list_ontology <- unlist(ontology)
  }
  list_ontology <- unique(list_ontology)
  list_ontology <- sapply(list_ontology,shQuote)
  list_ontology <- paste(list_ontology,collapse = ",")
  
  con <- DBI::dbConnect(RMySQL::MySQL(),
                        dbname = dbname,
                        username = username,
                        host = host,
                        password = conpass)
  sql <- paste0('select * from ontology where commonName in (',
                list_ontology,');')
  
  df <- DBI::dbGetQuery(con,sql)
  DBI::dbDisconnect(con)
  print(colnames(df))
  rampontoid <- paste(sapply(unique(df$rampOntologyIdLocation),shQuote),
                      collapse = ',')
  con <- DBI::dbConnect(RMySQL::MySQL(),
                        dbname = dbname,
                        username = username,
                        host = host,
                        password = conpass)
  sql <- paste0('select * from analytehasontology where rampOntologyIdLocation in (',
                rampontoid,');')
  
  df2 <- DBI::dbGetQuery(con,sql)
  DBI::dbDisconnect(con)
  
  print(colnames(df2))
  rampid <- paste(sapply(unique(df2$rampCompoundId),shQuote),collapse = ',')
  con <- DBI::dbConnect(RMySQL::MySQL(),
                        dbname = dbname,
                        username = username,
                        host = host,
                        password = conpass)
  sql <- paste0('select * from source where rampId in (',rampid,');')
  df3 <- DBI::dbGetQuery(con,sql)
  df3 <- unique(df3)
  print(colnames(df3))
  DBI::dbDisconnect(con)
  if(nrow(df3) == 0) {
    message('No searching result since this ramp id 
            is not existed in source table')
    return(NULL)
  }
  print('Merging 1...')
  mdf <- unique(merge(df3,df2,by.x = 'rampId',by.y = 'rampCompoundId',all.x = T))
  print('Merging 2...')
  mdf <- unique(merge(mdf,df,by.x = 'rampOntologyIdLocation',
                      by.y = 'rampOntologyIdLocation',
                      all.x = T))
  mdf <- mdf[c('sourceId','commonName.y','biofluidORcellular')]
  colnames(mdf) <- c('Metabolites','Ontology','Ontology_Type')
  return(mdf)
  
}