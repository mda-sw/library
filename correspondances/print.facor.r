print.facor <- function(facorres) {

    centr <- facorres$centersr
    labsr <- facorres$cluslabrow
    centc <- facorres$centersc
    labsc <- facorres$cluslabcol

    nclr <- nrow(centr)
    nfac <- ncol(centr)
    nclc <- nrow(centc)
    # nfac <- ncol(centc)

    cat("Projections of observations x variables (rows x columns) on \n")
    cat(nfac, " factors.  Projections on scale of 0 to 9.\n")

    # asimple <- round( 9.499*(centr-min(centr))/(max(centr)-min(centr)) )
    asimple <- centr
    cat("Obs. (rows) x number of factors retained (on scale of 0 to 9).\n")
    for (i in 1:nclr) {
        cat(asimple[i,], " Cluster ", i, ":", labsr[[i]], "\n")
    }   

    # asimple <- round( 9.499*(centc-min(centc))/(max(centc)-min(centc)) )
    asimple <- round( centr )
    cat("Vbes. (columns) x number of factors retained (on scale of 0 to 9).\n")
    for (j in 1:nclc) { 
        cat(asimple[j,], " Cluster ", j, ":", labsc[[j]], "\n")
    }   
    
}

