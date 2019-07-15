#' Bacterial Data
#'
#' Generate and download a csv file for either Zone or Minimum Inhibitory Concentration (MIC) Data.
#'
#' @param mic logical; FALSE (default) Zone Data are generated
#'
#'
#' @return Dowload a csv file with either MIC or Zone Data.
#'
#'
#' @details
#' The function scrapeBacteria webscrapes the antimicrobial wild type distributions of microorganisms from the EUCAST-webpage.
#' It provides the collected data for all combinations of antibiotics and bacteria species in a CSV file for either MIC Data or Zone Data.
#'
#' MIC Data:
#' \itemize{
#'  \item{Antimicrobial (the name of the antibiotic)}
#'  \item{Bacterium (the name of the bacterium)}
#'  \item{M0.002, ... ,M64 (19 variables, each of which contains the number of observations where the value specified in the name was observed)}
#'  \item{Distributions(the number of sources for the observations)}
#'  \item{Observations(the total number of observations, so the sum of M0.002, ... ,M64 )}
#' }
#'
#' Zone Data:
#' \itemize{
#'  \item{Antimicrobial (the name of the antibiotic)}
#'  \item{Bacterium (the name of the bacterium)}
#'  \item{DiskContent (numeric)}
#'  \item{Z6, ... ,Z50 (45 variables, each of which contains the number of observations where the value (in millimeters) specified in the name was observed)}
#'  \item{ECOFF (the offcial Epidemiological Cutoff Value (short: ECOFF), which was determined by EUCAST)}
#'  \item{Distributions(the number of sources for the observations)}
#'  \item{Observations(the total number of observations, so the sum of M0.002, ... ,M64 )}
#' }
#'
#'
#' @source \url{https://mic.eucast.org/Eucast2/SearchController/search.jsp?action=init}
#'
#' @seealso \code{\link{MIC}} \code{\link{ZD}}
#'
#' @export
#'
#' @importFrom utils write.csv2



scrapeBacteria <- function(mic=T){
  stopifnot(is.logical(mic))
  #### setup webscraping ####
  url <- 'https://mic.eucast.org/Eucast2/SearchController/search.jsp?action=init'
  session1 <- rvest::html_session(url)
  f <- rvest::html_form(session1)[[1]]
  antibiotics <- f$fields[[6]]$options #get options for antibiotics
  antibiotics <- antibiotics[antibiotics!="-1"]

  #There is no submit-button in the html, therefore we create one, so submit_form() works
  #source: https://stackoverflow.com/questions/33885629/submit-form-with-no-submit-button-in-rvest
  fake_submit <- list(name = NULL,type = "submit",value = NULL,checked = NULL,
                      disabled = NULL,readonly = NULL,required = FALSE)
  attr(fake_submit, "class") <- "input"
  #### setup output df####
  ifelse(mic,
         DF <- data.frame("Antimicrobial"=character(),"Bacterium"=character(),"M0.002"=numeric(),"M0.004"=numeric(),	"M0.008"=numeric(),
                          "M0.016"=numeric(),	"M0.032"=numeric(),	"M0.064"=numeric(),	"M0.125"=numeric(), "M0.25"=numeric()	,"M0.5"=numeric(),
                          "M1"=numeric(),	"M2"=numeric(),	"M4"=numeric(),	"M8"=numeric(),	"M16"=numeric(),	"M32"=numeric(),	"M64"=numeric(),
                          "M128"=numeric(),"M256"=numeric(),	"M512"=numeric(),
                          "ECOFF"=numeric(),"Distributions"=numeric(),"Observations"=numeric()),
         DF <- data.frame( "Antimicrobial"=character(),"Bacterium"=character(),
                           "DiskContent"=numeric(),
                           "Z6"=numeric(),  "Z7" =numeric(), "Z8" =numeric(), "Z9" =numeric(), "Z10" =numeric(),"Z11"=numeric(),
                           "Z12"=numeric(), "Z13"=numeric(), "Z14" =numeric(),"Z15"=numeric(), "Z16"=numeric(), "Z17" =numeric(),
                           "Z18"=numeric(), "Z19" =numeric(),"Z20" =numeric(),"Z21"=numeric(), "Z22"=numeric(), "Z23"=numeric(),
                           "Z24" =numeric(),"Z25"=numeric(), "Z26"=numeric(), "Z27"=numeric(), "Z28" =numeric(),"Z29" =numeric(),
                           "Z30"=numeric(),"Z31" =numeric(),"Z32" =numeric(),"Z33"=numeric(), "Z34"=numeric(), "Z35" =numeric(),
                           "Z36"=numeric(), "Z37" =numeric(),"Z38" =numeric(),"Z39" =numeric(),"Z40"=numeric(), "Z41" =numeric(),
                           "Z42"=numeric(), "Z43"=numeric(), "Z44"=numeric(), "Z45"=numeric(), "Z46" =numeric(),"Z47" =numeric(),
                           "Z48"=numeric(), "Z49"=numeric(), "Z50"=numeric(),
                           "ECOFF"=numeric(), 	"Distributions"=numeric(),"Observations"=numeric())
  )


  #### loop through options ####
  suppressMessages(suppressWarnings( # Message= "Submitting with NumberIndex", Warning when empty cell in table on website
    for(i in 1:length(antibiotics)){
      f$fields[[3]]$checked <- NULL
      f$fields[[4]]$checked <- "checked"
      f$fields[[5]]$value <- "1000"
      f$fields[[6]]$value <- antibiotics[i]
      f$fields[["submit"]] <- fake_submit
      session1 <- rvest::submit_form(session1, f)

      if(!mic){session1 <- rvest::html_session(stringr::str_replace(session1$url, regex("=mic"), "=dif"))}
      ifelse(mic,nc <- 23,nc<- 50)
      d <- tryCatch( # sometimes no data available => table doesn't exist => tryCatch necessary
        d <- rvest::html_table(rvest::html_node(xml2::read_html(session1), xpath="/html/body/table[3]"), header =T),
        error = function(e){NA})

      D <-data.frame(as.list(rep(NA,nc)))
      if( !is.null(dim(d)) && ncol(d)>1){
        colnames(d)[1]<-"Bacterium"
        if (any(d$Bacterium =="")){d <- d[-which(d$Bacterium==""),]}
        d <- d[,-which(colnames(d)=="")]
        ifelse(nrow(d)>1,
               D <- data.frame(d[,1],as.data.frame(apply(d[,-1],2, as.numeric), stringsAsFactors = F), stringsAsFactors = F),
               D <- data.frame(Antimicrobial=d[,1], as.list(as.numeric(d[,-1])), stringsAsFactors = F)
        )}
      #add data to DF
      df <- data.frame(names(antibiotics)[i],as.list(D), stringsAsFactors = F)
      colnames(df) <- colnames(DF)
      DF <- rbind(DF, df)
      #prepare for next loop
      d <- NA
      f <- rvest::html_form(session1)[[1]]
    }))
  #### write csv ####
  filename <- ifelse(mic, "MicData.csv", "ZoneData.csv")
  write.csv2(x=DF, file = filename,row.names=FALSE, na="")
}




