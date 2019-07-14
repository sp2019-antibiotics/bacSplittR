#' Visualization of the results of estimateCutoff
#'
#' Plotting the nls function with the estimated ECOFF.
#' Visualiziation through a pdf.
#'
#'
#' @param m matrix of all possible ECOFF estimates, the corresponding estimates for parameters mean, sd and n  and fit values
#' @param startP numeric vector of start parameters
#' @param DF data frame, containing columns "diam" and "obs" which specify the Zone Data
#'
#' @seealso \code{\link{estimateCutoff}}
#' @import ggplot2
#' @importFrom stats dnorm density

plotDistr <- function(m,startP,DF){
  # check input
  stopifnot(
    is.matrix(m) & is.numeric(m) & ncol(m)==5 & nrow(m)>=1
    & is.numeric(startP) & length(startP)==3
    & is.numeric(DF$diam) & is.numeric(DF$obs))

  best <- which.min(m[,5]) #= estimate for the cutoff
  diam <- obs <- NULL #for package, otherwise Note, because variables are undefined
  if(best==1){
    ggplot2::ggplot()+
      geom_col(data=DF,mapping=aes(x=diam,y=obs))+ #use for i = 1
      geom_line(mapping=aes(x=seq(min(DF$diam),max(DF$diam),by=0.2),y=dnorm(seq(min(DF$diam),max(DF$diam),by=0.2),startP[1],startP[2])*startP[3],colour="Starting Distribution"), linetype="dashed", size=1)+
      geom_line(mapping=aes(x=seq(min(DF$diam),max(DF$diam),by=0.2),y=dnorm(seq(min(DF$diam),max(DF$diam),by=0.2),m[,2][best],m[,3][best])*m[,4][best],colour="Fitted Distribution"), size=1)+
      geom_vline(xintercept=m[,1][best], colour="indianred1", size=1)+
      ggtitle("ECOFFinder")+
      guides(colour=guide_legend(title=NULL))+
      theme(legend.position = c(0.8,0.8))+
      ylab("Observed Bacteria") +
      xlab("Zone Diameter")
  }else{
    ggplot2::ggplot()+
      geom_col(data=DF[-(1:(best-1)),],mapping=aes(x=diam,y=obs))+
      geom_col(data=DF[1:(best-1),],mapping=aes(x=diam,y=obs), fill="gray")+
      geom_line(mapping=aes(x=seq(min(DF$diam),max(DF$diam),by=0.2),y=dnorm(seq(min(DF$diam),max(DF$diam),by=0.2),startP[1],startP[2])*startP[3],colour="Starting Distribution"), linetype="dashed", size=1)+
      geom_line(mapping=aes(x=seq(min(DF$diam),max(DF$diam),by=0.2),y=dnorm(seq(min(DF$diam),max(DF$diam),by=0.2),m[,2][best],m[,3][best])*m[,4][best],colour="Fitted Distribution"), size=1)+
      geom_vline(xintercept=m[,1][best], colour="indianred1", size=1)+
      ggtitle("ECOFFinder")+
      guides(colour=guide_legend(title=NULL))+
      theme(legend.position = c(0.8,0.8))+
      ylab("Observed Bacteria") +
      xlab("Zone Diameter")
  }
}
