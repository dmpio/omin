try:
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr

    RPY2_INSTALL=True

except Exception as err:
    RPY2_INSTALL=False
    msg = """rpi2 does not appear to be instaled on your system.
    It can be installed by following the directions found here:

    https://rpy2.readthedocs.io/en/version_2.8.x/overview.html#install-from-source
    """
    print(msg)


# Adding a function to the R library limma.
reduce_fit_script = """
reduceFit <- function(fit, results=NULL, file, digits=NULL, adjust="none", method="separate", F.adjust="none", quote=FALSE, sep="\t", row.names=TRUE, ...)

{
	if(!is(fit, "MArrayLM")) stop("fit should be an MArrayLM object")
	if(!is.null(results) && !is(results,"TestResults")) stop("results should be a TestResults object")
	if(is.null(fit$t) || is.null(fit$p.value)) fit <- eBayes(fit)
	method <- match.arg(method, c("separate","global"))

	p.value <- as.matrix(fit$p.value)
	if(adjust=="none") {
		p.value.adj <- NULL
	} else {
		p.value.adj <- p.value
		if(method=="separate") for (j in 1:ncol(p.value)) p.value.adj[,j] <- p.adjust(p.value[,j],method=adjust)
		if(method=="global") p.value.adj[] <- p.adjust(p.value,method=adjust)
	}
	if(F.adjust=="none" || is.null(fit$F.p.value))
		F.p.value.adj <- NULL
	else
		F.p.value.adj <- p.adjust(fit$F.p.value,method=F.adjust)

#	Optionally, round results for easy reading
	if(is.null(digits)) {
		rn <- function(x,digits=digits) x
	} else {
		rn <- function(x,digits=digits)
				if(is.null(x))
					NULL
				else {
					if(is.matrix(x) && ncol(x)==1) x <- x[,1]
					round(x,digits=digits)
				}
	}

#	Prepare output data.frame
	tab <- list()
	tab$A <- rn(fit$Amean,digits=digits-1)
	tab$Coef <- rn(fit$coef,digits=digits)
	tab$t <- rn(fit$t,digits=digits-1)
	tab$p.value <- rn(p.value,digits=digits+2)
	tab$p.value.adj <- rn(p.value.adj,digits=digits+3)
	tab$F <- rn(fit$F,digits=digits-1)
	tab$F.p.value <- rn(fit$F.p.value,digits=digits+2)
	tab$F.p.value.adj <- rn(F.p.value.adj,digits=digits+3)
	tab$Res <- unclass(results)
	tab$Genes <- fit$genes
	tab <- data.frame(tab,check.names=FALSE)

	if(is.null(row.names(fit))) {
		row.names <- FALSE
	} else {
		row.names(tab) <- row.names(fit)
	}

	#write.table(tab,quote=quote,row.names=row.names,sep=sep,...)
    return(tab)
}

"""

def install_limma():
    base = importr('base')
    base.source("http://www.bioconductor.org/biocLite.R")
    biocinstaller = importr("BiocInstaller")
    biocinstaller.biocLite("limma")


if RPY2_INSTALL:

    class RTools:
        """
        """
        # Set robjects as ro.
        ro = ro
        # Set the robjects.r as R
        R = ro.r
        # Set the rprint function.
        rprint = ro.globalenv.get("print")
        # Initialize the base package.
        base = importr('base')

        try:
            # Initialize the limma module.
            limma = importr('limma')
            # Initialize the reducefit function.
            limma.reduce_fit = ro.r(reduce_fit_script)

        except Exception as err:
            print(err)
            print("Attempting to install limma")

            try:
                install_limma()
                # Initialize the limma module.
                limma = importr('limma')
                # Initialize the reducefit function.
                limma.reduce_fit = ro.r(reduce_fit_script)

            except Exception as err:
                print(err)


else:# If rpy2 is not installed an empty class is given.
    class RTools(object):
        """
        """
        pass
