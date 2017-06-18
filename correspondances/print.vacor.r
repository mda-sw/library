print.vacor <- function(vacorres) {

    centc <- vacorres$vals
    labsc <- vacorres$labsc
    labsr <- vacorres$labsr

    n <- nrow(centc)
    m <- ncol(centc)

    cat("Projections of clusters of obs. vs. clusters of vbes. \n")
    cat("Projections on scale of 0 to 9.\n")
    #asimple <- round( 9.499*(centc-min(centc))/(max(centc)-min(centc)) )
    #asimple <- round ( 10*centc )
    asimple <- centc  

    cat("Var. labs.: ")
    for (j in 1:m) cat(labsc[[j]], " ")
    cat("\n") 
    for (i in 1:n) {
        # cat(as.matrix(asimple[i,]), " Cluster ", i, ": ", labsr[[i]], "\n")
        # cat(i, "&", labsr[[i]], " & ", asimple[i,1], " & ", asimple[i,2], 
        #    " \\\\ \n")
        cat(i, " & ")
        for (j in 1:m) {
            cat(asimple[i,j], " & ")
        }
        cat("\\\\ \n")
    }   

}

