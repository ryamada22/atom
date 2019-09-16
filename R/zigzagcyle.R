my.zigzag.cycle <- function(f){
	WcWd <- my.W(f)
	series <- my.zigzag(WcWd$Wc,WcWd$Wd)
	unique.series <- my.odd.unique(series)
	cycles <- apply(unique.series,1,my.odd.return)
	if(is.matrix(cycles)){
		ret <- list()
		for(i in 1:length(cycles[1,])){
			ret[[i]] <- cycles[,i]
		}
		cycles <- ret
	}
	return(cycles)
}

cy.tetra <- my.zigzag.cycle(tetra)
cy.antiprism6 <- my.zigzag.cycle(antiprism6)
tetantiprism6 <- my.paste.tri(tetra,antiprism6,1,3)
cy.tetantiprism6 <- my.zigzag.cycle(tetantiprism6)

