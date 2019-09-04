colorvis <- function(COL){
  plot(NULL, xlim=c(0,length(COL)), ylim=c(0,1), 
       xlab="", ylab="", xaxt="n", yaxt="n")
  rect(0:(length(COL)-1), 0, 1:length(COL), 1, col=COL)
}