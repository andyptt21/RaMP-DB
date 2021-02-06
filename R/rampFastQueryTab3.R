#' Do fisher test for only one pathway from search result
#' clicked on highchart
#' @param pathwaydf a data frame resulting from getPathwayFromAnalyte
#' @param backgrounddf optional dataframe resulting from getPathwayFromAnalyte run on all analytes detected in study. This will be used as the background for the Fisher's contingency table. If left NULL, all analytes in RaMP are used as background
#' @param total_metabolites number of metabolites analyzed in the experiment (e.g. background) (default is 1000; set to 'NULL' to retrieve total number of metabolites that map to any pathway in RaMP). Assumption that analyte_type is "metabolite")
#' @param total_genes number of genes analyzed in the experiment (e.g. background) (default is 20000, with assumption that analyte_type is "genes")
#' @param analyte_type "metabolites" or "genes" (default is "metabolites")
#' @param holmadjust boolean indicating if holm-adjusted p values should be included in output
#' @param conpass password for database access (string)
#' @param dbname name of the mysql database (default is "ramp")
#' @param username username for database access (default is "root")
#' @param host host name for database access (default is "localhost")
#' @param socket (optional) location of mysql.sock file
#' @return a dataframe with columns containing pathway ID, fisher's p value, user analytes in pathway, and total analytes in pathway
runFisherTest <- function(pathwaydf,backgrounddf=NULL,
                          total_metabolites=NULL,total_genes=20000,
                          analyte_type="metabolites",conpass=NULL,
                          holmadjust=TRUE,
                          dbname="ramp",username="root",
                          host = "localhost",
                          socket = NULL){
    now <- proc.time()
    print("Fisher Testing ......")
    if(is.null(conpass)) {
        stop("Please define the password for the mysql connection")
    }
    if(analyte_type=="metabolites") {total_analytes=total_metabolites
    } else if (analyte_type=="genes") {
	total_analytes=total_genes
    } else {
	stop("Please define the analyte_type variable as 'metabolites' or 'genes'")
    }
    
    contingencyTb <- matrix(0,nrow = 2,ncol = 2)
    colnames(contingencyTb) <- c("In Pathway","Not In Pathway")
    rownames(contingencyTb) <- c("All Metabolites","User's Metabolites")
    
    ## Get pathway ids that contain the user analytes
    pid <- unique(pathwaydf$pathwayRampId);
    list_pid <- sapply(pid,shQuote)
    list_pid <- paste(list_pid,collapse = ",")
    
  ## Get the total number of metabolites that are mapped to pathways in RaMP (that's the default background)
    query <- "select * from analytehaspathway"
    con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
                          password = conpass,
                          dbname = dbname,
                          host = host,
                          unix.socket = socket)
    allids <- DBI::dbGetQuery(con,query)
    DBI::dbDisconnect(con)
    allids <- allids[!duplicated(allids),]
    
    if((analyte_type == "metabolites") && (is.null(total_metabolites))) {
	wiki_totanalytes <- length(unique(allids$rampId[grep("RAMP_C",allids[which(allids$pathwaySource=="wiki"),"rampId"])]))
	react_totanalytes <- length(unique(allids$rampId[grep("RAMP_C",allids[which(allids$pathwaySource=="reactome"),"rampId"])]))
	kegg_totanalytes <- length(unique(allids$rampId[grep("RAMP_C",allids[which(allids$pathwaySource=="kegg"),"rampId"])]))
    }
    if(analyte_type=="genes") {
	wiki_totanalytes <- react_totanalytes <- kegg_totanalytes <- total_genes
    }
    
    print("Calculating p-values for pathways in input")
    ## Retrieve the Ramp compound ids associated with the ramp pathway id and count them:

    query1 <- paste0("select rampId,pathwayRampId from analytehaspathway where pathwayRampId in (",
                     list_pid,")")
    
    con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
                          password = conpass,
                          dbname = dbname,
                          host = host,
                          unix.socket = socket)
    cids <- DBI::dbGetQuery(con,query1)#[[1]]
    
    DBI::dbDisconnect(con)
    
    if(!is.null(backgrounddf)){
        ## Get background pathway ids that contain the user analytes
        pid_bg <- unique(backgrounddf$pathwayRampId);
        list_pid_bg <- sapply(pid_bg,shQuote)
        list_pid_bg <- paste(list_pid_bg,collapse = ",")
        
        ## Retrieve compound ids associated with background pathways and count 
        query1 <- paste0("select rampId,pathwayRampId from analytehaspathway where pathwayRampId in (",
                         list_pid_bg,")")
        
        con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
                              password = conpass,
                              dbname = dbname,
                              host = host,
                              unix.socket = socket)
        cids_bg <- DBI::dbGetQuery(con,query1)#[[1]]
        
        DBI::dbDisconnect(con)
    }
    
   # Loop through each user pathway, build the contingency table, and calculate Fisher's Exact
   # test p-value
   pval=totinpath=userinpath=pidused=c()
   for (i in pid) {
        curpathcids <- unique(cids[which(cids[,"pathwayRampId"]==i),"rampId"])
        if(analyte_type=="metabolites") {
                tot_in_pathway <- length(grep("RAMP_C",curpathcids))
        }else {
                tot_in_pathway <- length(grep("RAMP_G",curpathcids))
        }
	# Check that the pathway being considered has your analyte type, if not, move on
	if(tot_in_pathway==0) {
		next;
	} else if(is.null(backgrounddf)){
		if(allids$pathwaySource[which(allids$pathwayRampId==i)[1]] == "wiki") {
			total_analytes <- wiki_totanalytes
		} else if (allids$pathwaySource[which(allids$pathwayRampId==i)[1]] == "reactome") {
	                total_analytes <- react_totanalytes
	        } else if (allids$pathwaySource[which(allids$pathwayRampId==i)[1]] == "kegg") {
	                total_analytes <- kegg_totanalytes
	        } else {stop("Couldn't find pathway type for current pathway!")}

	 	tot_out_pathway <- total_analytes - tot_in_pathway
		  # fill the rest of the table out
		  user_in_pathway <- length(unique(pathwaydf[which(pathwaydf$pathwayRampId==i),"rampId"]))
		  user_out_pathway <- length(unique(pathwaydf$rampId)) - user_in_pathway
		  contingencyTb[1,1] <- tot_in_pathway - user_in_pathway
		  contingencyTb[1,2] <- tot_out_pathway - user_out_pathway
		  contingencyTb[2,1] <- user_in_pathway
		  contingencyTb[2,2] <- user_out_pathway
		  result <- stats::fisher.test(contingencyTb)
		  pval <- c(pval,result$p.value )
		  userinpath<-c(userinpath,user_in_pathway)
		  totinpath<-c(totinpath,tot_in_pathway)
		  pidused <- c(pidused,i)
	}else{
            user_in_pathway <- length(unique(pathwaydf[which(pathwaydf$pathwayRampId==i),"rampId"]))
            user_out_pathway <- length(unique(pathwaydf$rampId)) - user_in_pathway
            
            bg_in_pathway <- length(unique(backgrounddf[which(backgrounddf$pathwayRampId==i),"rampId"]))
            bg_out_pathway <- length(unique(backgrounddf$rampId)) - user_in_pathway
            
            contingencyTb[1,1] <- user_in_pathway
            contingencyTb[1,2] <- ifelse(bg_in_pathway - user_in_pathway < 0,
                                         0, bg_in_pathway - user_in_pathway)
            contingencyTb[2,1] <- user_out_pathway
            contingencyTb[2,2] <- ifelse(bg_out_pathway - user_out_pathway < 0,
                                         0, bg_out_pathway - user_out_pathway)
            result <- stats::fisher.test(contingencyTb)
            pval <- c(pval,result$p.value )
            userinpath<-c(userinpath,user_in_pathway)
            totinpath<-c(totinpath,tot_in_pathway)
            pidused <- c(pidused,i)
        }
   } # end for loop
    if(holmadjust){
        if(is.null(backgrounddf)){
            ## Now run fisher's tests for all other pids
            query <- "select distinct(pathwayRampId) from analytehaspathway where pathwaySource != 'hmdb';"
            con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
                                  password = conpass,
                                  dbname = dbname,
                                  host = host,
                                  unix.socket = socket)
            allpids <- DBI::dbGetQuery(con,query)
            DBI::dbDisconnect(con)
            pidstorun <- setdiff(allpids[,1],pid)
            pidstorunlist <- sapply(pidstorun,shQuote)
            pidstorunlist <- paste(pidstorunlist,collapse = ",")
            
            query2 <- paste0("select rampId,pathwayRampId from analytehaspathway where pathwayRampId in (",
                             pidstorunlist,")")
            
            con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
                                  password = conpass,
                                  dbname = dbname,
                                  host = host,
                                  unix.socket = socket)
            restcids <- DBI::dbGetQuery(con,query2)#[[1]]
            DBI::dbDisconnect(con)
            
            query1 <- paste0("select rampId,pathwayRampId from analytehaspathway;")
            
            con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
                                  password = conpass,
                                  dbname = dbname,
                                  host = host,
                                  unix.socket = socket)
            allcids <- DBI::dbGetQuery(con,query1)#[[1]]
            DBI::dbDisconnect(con)
            
            print("Calculating p-values for all other pathways")
            count=1;
            pval2=userinpath2=totinpath2=c()
            for (i in pidstorun) {
                if(( count %% 100) ==0) {print(paste0("Processed ",count))}
                count=count+1
                user_in_pathway=0
                curpathcids <- unique(allcids[which(allcids[,"pathwayRampId"]==i),"rampId"])
                if(analyte_type=="metabolites") {
                    tot_in_pathway <- length(grep("RAMP_C",curpathcids))
                }else {
                    tot_in_pathway <- length(grep("RAMP_G",curpathcids))
                }
                ## Check that the pathway being considered has your analyte type, if not, move on
                if(tot_in_pathway==0) {next;
                } else {
                    if(allids$pathwaySource[which(allids$pathwayRampId==i)[1]] == "wiki") {
                        total_analytes <- wiki_totanalytes
                    } else if (allids$pathwaySource[which(allids$pathwayRampId==i)[1]] == "reactome") {
                        total_analytes <- react_totanalytes
                    } else if (allids$pathwaySource[which(allids$pathwayRampId==i)[1]] == "kegg") {
                        total_analytes <- kegg_totanalytes
                    } else if (allids$pathwaySource[which(allids$pathwayRampId==i)[1]] == "hmdb") {
                        total_analytes <- NULL
                    } else {stop("Couldn't find pathway type for current pathway!")}
                    
                    if(is.null(total_analytes)) {next;}
                    tot_out_pathway <- total_analytes - tot_in_pathway
                    ## fill the rest of the table out
                    user_out_pathway <- length(unique(pathwaydf$rampId))
                    contingencyTb[1,1] <- tot_in_pathway - user_in_pathway
                    contingencyTb[1,2] <- tot_out_pathway - user_out_pathway
                    contingencyTb[2,1] <- user_in_pathway
                    contingencyTb[2,2] <- user_out_pathway
                    result <- stats::fisher.test(contingencyTb)
                    pval2 <- c(pval2,result$p.value )
                    userinpath2<-c(userinpath2,user_in_pathway)
                    totinpath2<-c(totinpath2,tot_in_pathway)
                }
            } # end for loop
            
        }else{ ## Holm adjustment pvals with a custom background
            pidstorun <- setdiff(pid_bg,pid)
        
            pval2=totinpath2=userinpath2=pidused2=c()
            for (i in pidstorun) {
                curpathcids <- unique(cids_bg[which(cids_bg[,"pathwayRampId"]==i),"rampId"])
                if(analyte_type=="metabolites") {
                    tot_in_pathway <- length(grep("RAMP_C",curpathcids))
                }else {
                    tot_in_pathway <- length(grep("RAMP_G",curpathcids))
                }
                ## Check that the pathway being considered has your analyte type, if not, move on
                if(tot_in_pathway==0) {
                    next;
                } else{
                    user_in_pathway <- length(unique(pathwaydf[which(pathwaydf$pathwayRampId==i),"rampId"]))
                    user_out_pathway <- length(unique(pathwaydf$rampId)) - user_in_pathway
                    
                    bg_in_pathway <- length(unique(backgrounddf[which(backgrounddf$pathwayRampId==i),"rampId"]))
                    bg_out_pathway <- length(unique(backgrounddf$rampId)) - user_in_pathway
                    
                    contingencyTb[1,1] <- user_in_pathway
                    contingencyTb[1,2] <- ifelse(bg_in_pathway - user_in_pathway < 0,
                                                 0, bg_in_pathway - user_in_pathway)
                    contingencyTb[2,1] <- user_out_pathway
                    contingencyTb[2,2] <- ifelse(bg_out_pathway - user_out_pathway < 0,
                                                 0, bg_out_pathway - user_out_pathway)
                    result <- stats::fisher.test(contingencyTb)
                    pval2 <- c(pval2,result$p.value )
                    userinpath2<-c(userinpath2,user_in_pathway)
                    totinpath2<-c(totinpath2,tot_in_pathway)
                    pidused2 <- c(pidused2,i)
                } 
            }
        }
        print(paste0("Calculated p-values for ",length(c(pval,pval2))," pathways"))
        ## only keep pathways that have > 8 or < 100 compounds
        keepers <- intersect(which(c(totinpath,totinpath2)>=8),
                             which(c(totinpath,totinpath2)<100))
        print(paste0("Keeping ",length(keepers)," pathways"))
        
        ## format output (retrieve pathway name for each unique source id first
        out <- data.frame(pathwayRampId=c(pidused,pidstorun)[keepers],
                          Pval=c(pval,pval2)[keepers],
                          Num_In_Path=c(userinpath,userinpath2)[keepers],
                          Total_In_Path=c(totinpath,totinpath2)[keepers])
        print(dim(out))
        out = out[!duplicated(out),]
        print(colnames(out))
        return(out)
    }else{
        keepers <- intersect(which(c(totinpath)>=8),
                             which(c(totinpath)<100))
        out <- data.frame(pathwayRampId=c(pidused)[keepers],
                          Pval=c(pval)[keepers],
                          Num_In_Path=c(userinpath)[keepers],
                          Total_In_Path=c(totinpath)[keepers])
        
        out = out[!duplicated(out),]
        return(out)
    }
}



#' Do fisher test for only one pathway from search result
#' clicked on highchart
#' @param pathwaydf a data frame resulting from getPathwayFromAnalyte
#' @param holmadjust boolean indicating if holm-adjusted p values should be included in output
#' @param analyte_type "metabolites" or "genes" (default is "metabolites")
#' @param backgrounddf optional dataframe resulting from getPathwayFromAnalyte run on all analytes detected in study. This will be used as the background for the Fisher's contingency table. If left NULL, all analytes in RaMP are used as background
#' @param total_metabolites number of metabolites analyzed in the experiment (e.g. background) (default is 1000; set to 'NULL' to retrieve total number of metabolites that map to any pathway in RaMP). Assumption that analyte_type is "metabolite")
#' @param total_genes number of genes analyzed in the experiment (e.g. background) (default is 20000, with assumption that analyte_type is "genes")
#' @param min_analyte if the number of analytes (gene or metabolite) in a pathway is
#' < min_analyte, do not report
#' @param conpass password for database access (string)
#' @param dbname name of the mysql database (default is "ramp")
#' @param username username for database access (default is "root")
#' @param host host name for database access (default is "localhost")
#' @param socket (optional) location of mysql.sock file
#' @return a list containing two entries: [[1]] fishresults, a dataframe containing pathways with Fisher's p values (raw and with FDR and Holm adjustment), number of user analytes in pathway, total number of analytes in pathway, and pathway source ID/database. [[2]] analyte_type, a string specifying the type of analyte input into the function ("genes", "metabolites", or "both")
#'@examples
#'\dontrun{
#' pathwaydf<-getPathwayFromAnalyte(c("MDM2","TP53","glutamate","creatinine"),
#'                 NameOrIds="names", conpass=conpass)
#' fisher.results <- runCombinedFisherTest(pathwaydf=pathwaydf,conpass=conpass)
#'}
#' @export
runCombinedFisherTest <- function(pathwaydf,backgrounddf=NULL,
                          total_metabolites=NULL,total_genes=20000,
                          analyte_type="metabolites",conpass=NULL,
                          min_analyte=2,
                          holmadjust=TRUE,
                          dbname="ramp",username="root",
                          host = "localhost",
                          socket = NULL){

    if(is.null(conpass)) {
    stop("Please define the password for the mysql connection")
  }

    G <- M <- 0

    # Grab pathways that contain metabolites to run Fisher on metabolites
    # This will return all pathways that have at 8-120 metabolites/genes in them
    fishmetab <- pathwaydf[grep("RAMP_C_",pathwaydf$rampId),]
    if(nrow(fishmetab) == 0) {
        outmetab=NULL
    } else{
        M=1
        print("Running Fisher's tests on metabolites")
        outmetab <- runFisherTest(pathwaydf=fishmetab,backgrounddf=backgrounddf,
                                  analyte_type="metabolites",
                                  total_metabolites=total_metabolites,total_genes=total_genes,
                                  holmadjust=holmadjust,
                                  conpass=conpass,dbname=dbname,
                                  username=username,host=host,
                                  socket=socket)
    }

    # Grab pathways that contain genes to run Fisher on genes
    fishgene <- pathwaydf[grep("RAMP_G_",pathwaydf$rampId),]
    if(nrow(fishgene) == 0) {
        outgene=NULL
    } else{
	G=1
	print("Running Fisher's tests on genes")
        outgene <- runFisherTest(pathwaydf=fishgene, backgrounddf=backgrounddf,
                                 analyte_type="genes",
                                 total_metabolites=total_metabolites,total_genes=total_genes,
                                 holmadjust=holmadjust,
                                 conpass=conpass,dbname=dbname,
                                 username=username,host=host,
                                 socket=socket)
    }
    
    if(is.null(outgene) & !is.null(outmetab)) {
	out <- outmetab
	fdr <- stats::p.adjust(out$Pval,method="fdr")
	out<-cbind(out,fdr);colnames(out)[ncol(out)]="Pval_FDR"
        if(holmadjust){
            holm <- stats::p.adjust(out$Pval,method="holm")
            out<-cbind(out,holm);colnames(out)[ncol(out)]="Pval_Holm"
        }
	keepers <- which(out$Num_In_Path>=min_analyte)
	out2 <- merge(out[keepers,],
        	pathwaydf[,c("pathwayName","pathwayRampId","pathwaysourceId",
        	"pathwaysource")],by="pathwayRampId")
     } else if (!is.null(outgene) & is.null(outmetab)) {
        out <- outgene
        fdr <- stats::p.adjust(out$Pval,method="fdr")
        out<-cbind(out,fdr);colnames(out)[ncol(out)]="Pval_FDR"
        if(holmadjust){
            holm <- stats::p.adjust(out$Pval,method="holm")
            out<-cbind(out,holm);colnames(out)[ncol(out)]="Pval_Holm"
        }
        keepers <- which(out$Num_In_Path>=min_analyte)
        out2 <- merge(out[keepers,],
                      pathwaydf[,c("pathwayName","pathwayRampId","pathwaysourceId",
                                   "pathwaysource")],by="pathwayRampId")
    } else {
        ## merge the results if both genes and metabolites were run
        G = M = 1
        allfish <- merge(outmetab,outgene,
                             by="pathwayRampId",all.x=T,all.y=T)
        colnames(allfish)[which(colnames(allfish)=="Pval.x")]="Pval.Metab"
        colnames(allfish)[which(colnames(allfish)=="Pval.y")]="Pval.Gene"
        colnames(allfish)[which(colnames(allfish)=="Total_In_Path.x")]="Total_In_Path.Metab"
        colnames(allfish)[which(colnames(allfish)=="Total_In_Path.y")]="Total_In_Path.Gene"
        colnames(allfish)[which(colnames(allfish)=="Num_In_Path.x")]="Num_In_Path.Metab"
        colnames(allfish)[which(colnames(allfish)=="Num_In_Path.y")]="Num_In_Path.Gene"
        
        ## Calculate combined p-values for pathways that have both genes and metabolites
        gm <- intersect(which(!is.na(allfish$Pval.Metab)),which(!is.na(allfish$Pval.Gene)))
        combpval <- stats::pchisq(-2 * (log(allfish$Pval.Metab[gm])+log(allfish$Pval.Gene[gm])),
                                  df=2,lower.tail=FALSE)
        
        g <- which(is.na(allfish$Pval.Metab))
        gpval <- allfish$Pval.Gene[g]
        m <- which(is.na(allfish$Pval.Gene))
        mpval <- allfish$Pval.Metab[m]
        
        out <- rbind(allfish[gm,],allfish[g,],allfish[m,])
        out <- cbind(out,c(combpval,gpval,mpval))
        colnames(out)[ncol(out)]="Pval_combined"
        fdr <- stats::p.adjust(out$Pval_combined,method="fdr")
        out <- cbind(out,fdr)
        colnames(out)[ncol(out)]="Pval_combined_FDR"
        if(holmadjust){
            holm <- stats::p.adjust(out$Pval_combined,method="holm")
            out <- cbind(out,holm)
            colnames(out)[ncol(out)]="Pval_combined_Holm"
        }
        keepers <- intersect(c(which(out$Num_In_Path.Metab>=min_analyte),
                               which(is.na(out$Num_In_Path.Metab))),
                             c(which(out$Num_In_Path.Gene>=min_analyte),
                               which(is.na(out$Num_In_Path.Gene)))
                             )


	    # Now that p-values are calculated, only return pathways that are in the list
	    # of pathways that contain user genes and metabolites
	    out2 <- merge(out[keepers,],
		pathwaydf[,c("pathwayName","pathwayRampId","pathwaysourceId",
		"pathwaysource")],by="pathwayRampId")
	} # end merging when genes and metabolites were run
    out2 <- out2[!duplicated(out2),]

    analyte_type=c()
    if(G==1 && M==1) {
	analyte_type="both"
    } else if (G==1 && M==0) {
		analyte_type="genes"
    } else if (G==0 && M==1) {
	analyte_type="metabolites"
    }

    return(list(fishresults=out2,analyte_type=analyte_type))
}

#' Function that search analytes (gene or compounds)  or a list of analytes and
#' returns associated pathways
#'
#' @param analytes a vector of analytes (genes or metabolites) that need to be searched
#' @param find_synonym find all synonyms or just return same synonym (T/F)
#' @param conpass password for database access (string)
#' @param NameOrIds whether input is "names" or "ids" (default is "ids")
#' @param host host name for database access (default is "localhost")
#' @param dbname name of the mysql database (default is "ramp")
#' @param socket (optional) location of mysql.sock file
#' @param username username for database access (default is "root")
#' @return a list contains all metabolites as name and pathway inside.
#' @examples
#' \dontrun{
#' mypath <- getPathwayFromAnalyte(analytes=c("2-hydroxyglutarate","glutamate"), conpass="mypassword")
#' }
#' @export
getPathwayFromAnalyte<- function (
                                 analytes = NULL,
                                 find_synonym = FALSE,
                                 conpass = NULL, 
                                 host = "localhost",
                                 dbname = "ramp",
                                 username = "root",
                                 NameOrIds = "ids",
                                 socket=NULL) {
    if (is.null(conpass)) {
        stop("Please define the password for the mysql connection")
    }
    now <- proc.time()
    if (is.null(analytes)) {
        return(NULL)
    }
    if (NameOrIds == "names") {
        synonym <- rampFindSynonymFromSynonym(synonym = analytes, 
            find_synonym = find_synonym, conpass = conpass, host = host, 
            dbname = dbname, username = username, socket = socket)
        colnames(synonym)[1] = "commonName"
        synonym$commonName <- tolower(synonym$commonName)
        if (nrow(synonym) == 0) {
            stop("Could not find any matches to the analytes entered.  If pasting, please make sure the names are delimited by end of line (not analyte per line)\nand that you are selecting 'names', not 'ids'")
        }
        list_metabolite <- unique(synonym$rampId)
        list_metabolite <- sapply(list_metabolite, shQuote)
        list_metabolite <- paste(list_metabolite, collapse = ",")
    }
    else if (NameOrIds == "ids") {
        sourceramp <- rampFindSourceRampId(sourceId = analytes, 
            conpass = conpass, host = host, dbname = dbname, 
            username = username, socket=socket)
        if (nrow(sourceramp) == 0) {
            stop("Make sure you are actually inputting ids and not names (you have NameOrIds set to 'ids'. If you are, then no ids were matched in the RaMP database.")
        }
        list_metabolite <- unique(sourceramp$rampId)
        list_metabolite <- sapply(list_metabolite, shQuote)
        list_metabolite <- paste(list_metabolite, collapse = ",")
    }
    else {
        stop("Make sure NameOrIds is set to 'names' or 'ids'")
    }
    if (list_metabolite == "") {
        warning("Unable to retrieve metabolites")
        return(NULL)
    }
    query2 <- paste0("select pathwayRampId,rampId from analytehaspathway where\n                      rampId in (", 
        list_metabolite, ");")
    con <- connectToRaMP(dbname = dbname, username = username, 
        conpass = conpass, host = host, socket = socket)
    df2 <- DBI::dbGetQuery(con, query2)
    DBI::dbDisconnect(con)
    pathid_list <- df2$pathwayRampId
    pathid_list <- sapply(pathid_list, shQuote)
    pathid_list <- paste(pathid_list, collapse = ",")
    if (pathid_list == "") {
        warning("The input list of analytes do not map to any pathways")
        return(NULL)
    }
    query3 <- paste0("select pathwayName,sourceId as pathwaysourceId,type as pathwaysource,pathwayRampId from pathway where pathwayRampId in (", 
        pathid_list, ");")
    con <- connectToRaMP(dbname = dbname, username = username, 
        conpass = conpass, host = host, socket = socket)
    df3 <- DBI::dbGetQuery(con, query3)
    DBI::dbDisconnect(con)
    mdf <- merge(df3, df2, all.x = T)
    if (NameOrIds == "ids") {
        list_analytes <- sapply(analytes, shQuote)
        list_analytes <- paste(list_analytes, collapse = ",")
        query4 <- paste0("select sourceId,commonName,rampId from source where sourceId in (", 
            list_analytes, ");")
        con <- connectToRaMP(dbname = dbname, username = username, 
            conpass = conpass, host = host, socket = socket)
        df4 <- DBI::dbGetQuery(con, query4)
        DBI::dbDisconnect(con)
        df4$commonName <- sapply(as.character(df4$commonName), 
            function(x) if (stringi::stri_enc_mark(x) == "native") {
                x <- iconv(x, "latin1", "UTF-8")
            }
            else {
                x
            })
        mdf <- merge(mdf, df4, all.x = T, by.y = "rampId")
        mdf$commonName = tolower(mdf$commonName)
    }
    else {
        mdf <- merge(mdf, synonym, all.x = T, by.y = "rampId")
    }
    out <- mdf[!duplicated(mdf), ]
    return(out[which(out$pathwaysource != "hmdb"), c("rampId", 
        "pathwayRampId", "pathwayName", "pathwaysourceId", "pathwaysource", 
        "commonName")])
}
#' Perform fuzzy multiple linkage partitioning clustering on pathways identified by
#' Fisher's test
#'
#' @param fishers_df The data frame generated by runFisherTest
#' @param perc_analyte_overlap Minimum overlap for pathways to be considered similar
#' (Default = 0.5)
#' @param min_pathway_tocluster Minimum number of 'similar' pathways required to start
#' a cluster (medoid) (Default = 3)
#' @param perc_pathway_overlap Minimum overlap for clusters to merge (Default = 0.5)
#'
#' @return list:[[1]] Pathway enrichment result dataframe with cluster assignment column added
#' [[2]] analyte type
#' [[3]] cluster assignment in the list form
#'@examples
#'\dontrun{
#' pathwaydf<-getPathwayFromAnalyte(c("MDM2","TP53","glutamate","creatinine"),
#'                 NameOrIds="names", conpass=conpass)
#' fisher.results <- runCombinedFisherTest(pathwaydf=pathwaydf,conpass=conpass)
#' filtered.fisher.results <- FilterFishersResults(fisher.results,p_holmadj_cutoff=0.05)
#' filteredclust.fisher.results <- findCluster(filtered.fisher.results)
#'}
#' @export
findCluster <- function(fishers_df,perc_analyte_overlap = 0.5,
                        min_pathway_tocluster = 2,perc_pathway_overlap = 0.5){
  if(perc_analyte_overlap <= 0 || perc_analyte_overlap >= 1 ||
     perc_pathway_overlap <= 0 || perc_pathway_overlap >= 1){
    return(NULL)
  }

  analyte_type=fishers_df$analyte_type
  fishers_df=fishers_df$fishresults
  if(nrow(fishers_df)==0){
    return(NULL)
  }else if(nrow(fishers_df)==1){
    fishers_df$cluster_assignment="Did not cluster"
    fishers_df$rampids<-fishers_df$pathwayRampId
    fishers_df$pathwayRampId<-NULL
    output<-list(fishresults=fishers_df,analyte_type=analyte_type,cluster_list="Did not cluster")
    return(output)
  } else {
    #similarity_matrix_list<-loadOverlapMatrices()
    similarity_matrix_gene <- genes_result
    similarity_matrix_analyte <- analyte_result
    similarity_matrix_metab <- metabolites_result
    if(analyte_type=="both"){
      #similarity_matrix = similarity_matrix_list[["analyte"]]
      similarity_matrix = similarity_matrix_analyte
    }else if(analyte_type=="metabolites"){
      #similarity_matrix = similarity_matrix_list[["metab"]]
      similarity_matrix = similarity_matrix_metab
    } else if(analyte_type=="genes"){
      #similarity_matrix = similarity_matrix_list[["gene"]]
      similarity_matrix = similarity_matrix_gene
    } else {
      stop("analyte_type should be 'genes' or metabolites'")
    }
    pathway_list<-fishers_df[,"pathwayRampId"]

    pathway_indices<-match(pathway_list,rownames(similarity_matrix))

    if(length(which(is.na(pathway_indices)))>0){
      pathway_indices<-pathway_indices[-which(is.na(pathway_indices))]
    }

    pathway_matrix<-similarity_matrix[pathway_indices,pathway_indices]
    unmerged_clusters<-apply(pathway_matrix, 1, function(x){
      # if(length(which(x>=perc_analyte_overlap))>(min_pathway_tocluster+1)){
      if(length(which(x>=perc_analyte_overlap))>(min_pathway_tocluster-1)){
        return(colnames(pathway_matrix)[which(x>=perc_analyte_overlap)])
      } else {
        return(NA)
      }
    })
    # Remove the unmerged clusters
    if(length(which(is.na(unmerged_clusters)))>0){
      unmerged_clusters<-unmerged_clusters[-which(is.na(unmerged_clusters))]
    }

    if(length(unmerged_clusters)==0){
      #stop("No medoids found, make perc_analyte_overlap or min_pathway_tocluster smaller")
      cluster_list<-rep("Did not cluster",times = nrow(fishers_df))
    }else{
      # Evaluate similarity between clusters
      cluster_similarity<-matrix(0,ncol = length(unmerged_clusters),nrow = length(unmerged_clusters))
      for(i in 1:length(unmerged_clusters)){
        for(j in 1:length(unmerged_clusters)){
          cluster_similarity[i,j]<-length(intersect(unmerged_clusters[[i]],unmerged_clusters[[j]]))/
            length(unique(c(unmerged_clusters[[i]],unmerged_clusters[[j]])))
        }
      }
      colnames(cluster_similarity)<-rownames(cluster_similarity)<-names(unmerged_clusters)
      unmerged_cluster_similarity<-cluster_similarity

      cluster_list<-unmerged_clusters

      # Merge Clusters
      count = 1
      while(length(which(cluster_similarity >= perc_pathway_overlap)) > nrow(cluster_similarity)){
        cluster_similarity_mod<-cluster_similarity
        for(i in 1:nrow(cluster_similarity_mod)){
          cluster_similarity_mod[i,i]<-0
        }

        clusters_to_merge<-which(cluster_similarity_mod == max(cluster_similarity_mod), arr.ind = TRUE)
        clusters_to_merge<-unique(t(apply(clusters_to_merge, 1, sort)))

        for(i in 1:nrow(clusters_to_merge)){
          if(!is.na(cluster_list[[clusters_to_merge[i,1]]])&&!is.na(cluster_list[[clusters_to_merge[i,2]]])){
            cluster_list[[clusters_to_merge[i,1]]]<-unique(unlist(cluster_list[c(clusters_to_merge[i,1],clusters_to_merge[i,2])]))
            cluster_list[[clusters_to_merge[i,2]]]<-NA
          }
        }

        if(length(which(is.na(cluster_list)))>0){
          cluster_list<-cluster_list[-which(is.na(cluster_list))]
        }

        cluster_similarity<-matrix(0,ncol = length(cluster_list),nrow = length(cluster_list))
        for(i in 1:length(cluster_list)){
          for(j in 1:length(cluster_list)){
            cluster_similarity[i,j]<-length(intersect(cluster_list[[i]],cluster_list[[j]]))/
              length(unique(c(cluster_list[[i]],cluster_list[[j]])))
          }
        }

        if(nrow(cluster_similarity)==1){
          #stop("Clusters converged, use larger perc_pathway_overlap")
          #return(rep(1,times = nrow(fishers_df)))
          cluster_list<-rep("Did not cluster",times = nrow(fishers_df))
        }
        count = count + 1
        if(count == length(unmerged_clusters)+1){
          stop("ERROR: while loop failed to terminate")
          #return(rep(1,times = nrow(fishers_df)))
          #cluster_list<-rep("Did not cluster",times = nrow(fishers_df))
        }
      }
      if(length(unique(cluster_list))!=1){
        colnames(cluster_similarity) = rownames(cluster_similarity) = paste0("cluster_",c(1:length(cluster_list)))
      }
    }
    #return(cluster_list)

    # Reformat cluster list to embed into results file
    rampids<-as.vector(fishers_df$pathwayRampId)
    fishers_df$pathwayRampId<-NULL

    if(length(cluster_list)>1){
      cluster_assignment<-sapply(rampids,function(x){
        pathway<-x
        clusters<-""
        for(i in 1:length(cluster_list)){
          if(pathway %in% cluster_list[[i]]){
            clusters<-paste0(clusters,i,sep = ", ",collapse = ", ")
          }
        }
        if(clusters!=""){
          clusters=substr(clusters,1,nchar(clusters)-2)
        }else{
          clusters = "Did not cluster"
        }
        return(clusters)
      })
      fishers_df<-cbind(fishers_df,cluster_assignment)
    }else{
      fishers_df<-cbind(fishers_df,rep("Did not cluster",times=nrow(fishers_df)))
    }
    fishers_df$rampids<-rampids
    output<-list(fishresults=fishers_df,analyte_type=analyte_type,cluster_list=cluster_list)
    return(output)
  }
}

#' Filter pathways by p-value cutoff for display and clustering
#' @param fishers_df The data frame generated by runFisherTest
#' @param p_holmadj_cutoff return pathways where Holm adjusted pvalues are < p_holmadj_cutoff
#' @param p_fdradj_cutoff return pathways where FDR adjusted pvalues are < p_fdradj_cutoff
#' @return list:[[1]]Dataframe with pathway enrichment results, only significant pathways
#' [[2]]analyte type
#'@examples
#'\dontrun{
#' pathwaydf<-getPathwayFromAnalyte(c("MDM2","TP53","glutamate","creatinine"),
#'                 NameOrIds="names", conpass=conpass)
#' fisher.results <- runCombinedFisherTest(pathwaydf=pathwaydf,conpass=conpass)
#' filtered.fisher.results <- FilterFishersResults(fisher.results,p_holmadj_cutoff=0.05)
#'}
#' @export
FilterFishersResults<-function(fishers_df,p_holmadj_cutoff=NULL,
                               p_fdradj_cutoff=NULL){

  # Check to see whether the output is from ORA performed on genes and metabolites
  # or genes or metabolites
  analyte_type=fishers_df$analyte_type
  fishers_df=fishers_df$fishresults

  if(length(grep("Pval_combined",colnames(fishers_df)))==0) {
    if(!is.null(p_holmadj_cutoff)) {
      return(list(fishresults=fishers_df[which(fishers_df[,"Pval_Holm"] <=
                                                 p_holmadj_cutoff),],analyte_type=analyte_type))
    } else if (!is.null(p_fdradj_cutoff)) {
      return(list(fishresults=fishers_df[which(fishers_df[,"Pval_FDR"] <=
                                                 p_fdradj_cutoff),],analyte_type=analyte_type))
    } else {
      stop("Please set a cutoff for Holm Adjusted pvalues
			(p_holmadj_cutoff paramter) or FDR Adjusted pvalues
			(p_fdradj_cutoff)")
    }
  }  else { # ORA was performed on both genes and metabolites:
    if(!is.null(p_holmadj_cutoff)) {
      return(list(fishresults=fishers_df[which(fishers_df[,"Pval_combined_Holm"] <=
                                                 p_holmadj_cutoff),],analyte_type=analyte_type))
    } else if (!is.null(p_fdradj_cutoff)) {
      return(list(fishresults=fishers_df[which(fishers_df[,"Pval_combined_FDR"] <=
                                                 p_fdradj_cutoff),],analyte_type=analyte_type))
    } else {
      stop("Please set a cutoff for Holm Adjusted pvalues
                        (p_holmadj_cutoff paramter) or FDR Adjusted pvalues
                        (p_fdradj_cutoff)")
    }
  }
}

#' Fisher's test with the p value permutation approach developed by Westfall and Young to account for overlapping pathway annotations
#' @param pathwaydf a data frame resulting from getPathwayFromAnalyte
#' @param backgrounddf optional dataframe resulting from getPathwayFromAnalyte run on all analytes detected in study. This will be used as the background for the Fisher's contingency table. If left NULL, all analytes in RaMP are used as background
#' @param total_metabolites number of metabolites analyzed in the experiment (e.g. background) (default is 1000; set to 'NULL' to retrieve total number of metabolites that map to any pathway in RaMP). Assumption that analyte_type is "metabolite")
#' @param total_genes number of genes analyzed in the experiment (e.g. background) (default is 20000, with assumption that analyte_type is "genes")
#' @param analyte_type "metabolites" or "genes" (default is "metabolites")
#' @param conpass password for database access (string)
#' @param dbname name of the mysql database (default is "ramp")
#' @param username username for database access (default is "root")
#' @param host host name for database access (default is "localhost")
#' @param socket (optional) location of mysql.sock file
#' @param repetitions number of times p values are permuted
#' @param metaboliteClusters Optional output of cluster_metabolites
#' @param ncpus optional parameter indicating number of cpus available for parallel computation
#' @return a dataframe with columns containing pathway ID, fisher's p value, user analytes in pathway, total analytes in pathway, and Westfall-Young permuted p-values
#' @export
WY_pathway_adjustment<-
    function(pathwaydf,backgrounddf,
             total_metabolites=NULL,total_genes=20000,
             analyte_type="metabolites",conpass=NULL,
             dbname="ramp",username="root",
             host = "localhost",
             socket = NULL,
             repetitions=100,metaboliteClusters=NULL,ncpus=1){
        now <- proc.time()
        print("Fisher Testing ......")
        if(is.null(conpass)) {
            stop("Please define the password for the mysql connection")
        }
        if(analyte_type=="metabolites") {
            total_analytes=total_metabolites
        } else if (analyte_type=="genes") {
            total_analytes=total_genes
        } else {
            stop("Please define the analyte_type variable as 'metabolites' or 'genes'")
        }
        
        contingencyTb <- matrix(0,nrow = 2,ncol = 2)
        colnames(contingencyTb) <- c("In Pathway","Not In Pathway")
        rownames(contingencyTb) <- c("All Metabolites","User's Metabolites")
        
        ## Get pathway ids that contain the user analytes
        pid <- unique(pathwaydf$pathwayRampId);
        list_pid <- sapply(pid,shQuote)
        list_pid <- paste(list_pid,collapse = ",")
        
        ## Get the total number of metabolites that are mapped to pathways in RaMP (that's the default background)
        query <- "select * from analytehaspathway"
        con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
                              password = conpass,
                              dbname = dbname,
                              host = host,
                              unix.socket = socket)
        allids <- DBI::dbGetQuery(con,query)
        DBI::dbDisconnect(con)
        allids <- allids[!duplicated(allids),]
        
        ## Get background pathway ids that contain the user analytes
        pid_bg <- unique(backgrounddf$pathwayRampId);
        list_pid_bg <- sapply(pid_bg,shQuote)
        list_pid_bg <- paste(list_pid_bg,collapse = ",")
        
        ## Retrieve compound ids associated with background pathways and count 
        query1 <- paste0("select rampId,pathwayRampId from analytehaspathway where pathwayRampId in (",
                         list_pid_bg,")")
        
        con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
                              password = conpass,
                              dbname = dbname,
                              host = host,
                              unix.socket = socket)
        cids_bg <- DBI::dbGetQuery(con,query1)#[[1]]
        
        DBI::dbDisconnect(con)
        
        fast_fishers<-function(cids_bg,pathwaydf,backgrounddf){
            ## Loop through each user pathway, build the contingency table, and calculate Fisher's Exact
        ## test p-value
            pval=totinpath=userinpath=pidused=c()
            for (i in pid) {
                ##curpathcids <- unique(cids[which(cids[,"pathwayRampId"]==i),"rampId"])
                curpathcids<-unique(intersect(cids_bg[which(cids_bg[,"pathwayRampId"]==i),"rampId"],
                                              pathwaydf$rampId))
                if(analyte_type=="metabolites") {
                    tot_in_pathway <- length(grep("RAMP_C",curpathcids))
                }else {
                    tot_in_pathway <- length(grep("RAMP_G",curpathcids))
                }
                ## Check that the pathway being considered has your analyte type, if not, move on
                if(tot_in_pathway==0) {
                    next;
                }else{
                    user_in_pathway <- length(unique(pathwaydf[which(pathwaydf$pathwayRampId==i),"rampId"]))
                    user_out_pathway <- length(unique(pathwaydf$rampId)) - user_in_pathway
                    
                    bg_in_pathway <- length(unique(backgrounddf[which(backgrounddf$pathwayRampId==i),"rampId"]))
                    bg_out_pathway <- length(unique(backgrounddf$rampId)) - user_in_pathway
                    
                    contingencyTb[1,1] <- user_in_pathway
                    contingencyTb[1,2] <- ifelse(bg_in_pathway - user_in_pathway < 0,
                                                 0, bg_in_pathway - user_in_pathway)
                    contingencyTb[2,1] <- user_out_pathway
                    contingencyTb[2,2] <- ifelse(bg_out_pathway - user_out_pathway < 0,
                                                 0, bg_out_pathway - user_out_pathway)
                    result <- stats::fisher.test(contingencyTb)
                    pval <- c(pval,result$p.value )
                    userinpath<-c(userinpath,user_in_pathway)
                    totinpath<-c(totinpath,tot_in_pathway)
                    pidused <- c(pidused,i)
                }
            }
            
            keepers <- intersect(which(c(totinpath)>=8),
                                 which(c(totinpath)<100))
            out <- data.frame(pathwayRampId=c(pidused)[keepers],
                              Pval=c(pval)[keepers],
                              Num_In_Path=c(userinpath)[keepers],
                              Total_In_Path=c(totinpath)[keepers])    
            out = out[!duplicated(out),]
            return(out)        
        }
        if(is.null(metaboliteClusters)){
            true_results<-fast_fishers(cids_bg,pathwaydf,backgrounddf) %>% dplyr::arrange("Pval")
            true_pathway_count<-nrow(true_results)
        }else{
            ## So cluster numbering is consistent between network and pathway visual
            metaboliteClusters <-metaboliteClusters %>%
                dplyr::arrange("cluster")
            
            true_results<-lapply(unique(metaboliteClusters$cluster),function(x){
                metabolites<-metaboliteClusters %>%
                    dplyr::filter("cluster"==x) %>%
                    dplyr::select("metabolite") %>%
                    as.matrix() %>% as.vector() %>%
                    rampFindSourceRampId(username = username,
                              conpass = conpass,
                              dbname = dbname,
                              host = host,
                              socket = socket) %>%
                    dplyr::select("rampId")
                pathwaydf<-pathwaydf %>%
                    dplyr::filter("rampId" %in% as.matrix(as.vector(metabolites)))
                true_results<-fast_fishers(cids_bg,pathwaydf,backgrounddf) %>% dplyr::arrange("Pval")
                return(true_results)
            })
            true_pathway_count<-max(sapply(true_results,nrow))
        }    
        metabolites_of_interest<-unique(pathwaydf$rampId)
        background_metabolites<-unique(backgrounddf$rampId)
        
        ## Make a matrix of pathway p vals using set of interest the same size as the true set
        empirical_pval_matrix<-parallel::mclapply(1:repetitions,function(x){
            fake_metabolites_of_interest<-sample(background_metabolites,
                                                 length(metabolites_of_interest))
            pathwaydf_fake<-backgrounddf %>%
                dplyr::filter(backgrounddf$rampId %in% fake_metabolites_of_interest)
            fake_result<-fast_fishers(cids_bg,pathwaydf_fake,backgrounddf)
            if(x %% 10 == 0){
                print(paste0(x, " out of ", repetitions, " pathway sets simulated"))
            }
            return(sort(fake_result$Pval))
        }, mc.cores=ncpus
        )
    
        ## empirical_pval_matrix<-do.call(cbind,
        ##                       empirical_pval_matrix)
        empirical_pval_matrix<-t(plyr::ldply(empirical_pval_matrix,rbind)) %>% as.matrix()
        empirical_pval_matrix[is.na(empirical_pval_matrix)]<-1
        if(nrow(empirical_pval_matrix)<true_pathway_count){
            rowdiff = true_pathway_count - nrow(empirical_pval_matrix)
            empirical_pval_matrix<-rbind(empirical_pval_matrix,
                                         matrix(1,nrow=rowdiff,ncol=ncol(empirical_pval_matrix)))
        }
        ## Find percentile of true value in distribution of false values
        find_percentile <- function(true_results){
            true_value_perc=c()
            for(i in 1:nrow(true_results)){
                percentile <- stats::ecdf(empirical_pval_matrix[i,])
                true_value_perc=c(true_value_perc,
                                  percentile(as.matrix(true_results$Pval)[i]))
            }
            true_results$WY_Pval=true_value_perc
            return(true_results)
        }
        
        if(is.null(metaboliteClusters)){
            final_results<-find_percentile(true_results)
        }else{
            output<-list()
            for(i in 1:length(true_results)){
                if(nrow(true_results[[i]])!=0){
                    output=rlist::list.append(output,find_percentile(true_results[[i]]))
                }else{
                    output=rlist::list.append(output,NULL)
                }
            }
            empty_clusters<-which(sapply(true_results,nrow)==0)
            final_results<-dplyr::bind_rows(output, .id="index") %>% dplyr::rename("cluster"="index")
            for(i in empty_clusters){
                final_results<-final_results %>%
                    dplyr::mutate("cluster"=sapply("cluster",
                                          function(x) ifelse(as.numeric(x)>=i,return(as.numeric(x)+1),return(x))))
            }
        }    
        final_results$pathwaySourceId=pathway_ramp_to_source(
            final_results$pathwayRampId,
            conpass=conpass,
            dbname=dbname,
            username=username,
            host=host,
            socket=socket)
        final_results$pathwayCommonName=pathway_source_to_common(
            final_results$pathwaySourceId,
            conpass=conpass,
            dbname=dbname,
            username=username,
            host=host,
            socket=socket)
        return(final_results)
    }
