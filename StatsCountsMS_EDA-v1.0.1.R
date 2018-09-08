###################################
###  MS/MS EDA PIPE-LINE v 1.55 ###
###################################
library(gtools)
library(gplots)
library(MASS)
library(RColorBrewer)
library(msmsEDA)
cls <- brewer.pal(8,"Dark2")

###  Widget to show status and output
output <- NULL

###  LC-MS/MS Exploratory Data Analysis
#########################################
msms.EDA <- function()
{
  etit <- tit <- get.expnm()
  pdf.flnm <- tit
  txt.flnm <- paste(tit,"_EDA.txt",sep="")

  insert(output,"\n\n   Exploratory Data Analysis")
  insert(output,"-------------------------------")

  ###  Read samples description
  samples.flnm <- svalue(g.smpls.flnm)
  if( !file.exists(samples.flnm) )
  { galert("Samples file not found!",cont=window)
    reset.smpls()
    return()
  }
  samples <- read.table(file=samples.flnm,sep="\t",header=TRUE,dec=",")
  cnms <- colnames(samples)
  if(cnms[1]!="Sample")
  { galert("Wrong samples file: No 'Sample' column!",cont=window)
    reset.smpls()
    return()
  }
  samples$Sample <- as.character(samples$Sample)
  if(!is.null(samples$offsets))
    samples$offsets <- as.numeric(as.character(samples$offsets))
    
  ###  Check if normalizing offsets given
  jdx <- 1
  cell.j <- which(cnms=="offsets")
  if(length(cell.j)>0) jdx <- c(jdx,cell.j)

  ###  Build factors matrix
  facs <- data.frame(samples[,-jdx,drop=FALSE])
  rownames(facs) <- as.character(samples$Sample)

  ###  Read dataset
  data.flnm <- svalue(g.data.flnm)
  msms <- read.table(file=data.flnm,header=TRUE,sep="\t",
                   stringsAsFactors=FALSE)

  ###  Check for data consistency
  if( length(intersect(colnames(msms),samples$Sample)) < nrow(samples) )
    stop("\nFATAL ERROR: Samples description file does not match data file.\n")

  if( ! "Accession" %in% colnames(msms) )
    stop("\nFATAL ERROR: no 'Accession' column in data file.\n")

  ### Main factor and treatment and control levels
  fnm <- colnames(facs)[1]
  condA <- levels(facs[,fnm])[1]
  condB <- levels(facs[,fnm])[2]
  fl <- (facs[,fnm]==condA | facs[,fnm]==condB)

  ### Subset to needed samples
  msms.counts <- data.matrix(msms[ ,samples$Sample[fl]])
  rownames(msms.counts) <- msms$Accession
  facs <- facs[fl,,drop=FALSE]
  nsmpl <- sum(fl)
  smpls <- samples[fl,,drop=FALSE]
  mnf <- facs[,1]
  
  ###  Remove rows of all equal values
  all.eq <- which(apply(msms.counts,1,function(x) all(x==x[1])))
  if(length(all.eq))
    msms.counts <- msms.counts[-all.eq,,drop=FALSE]

  ###  Short sample names to be used in graphics
  short.nms <- paste("X",1:nsmpl,sep="")

  ###  Create MSnSet object, and preprocess it by setting NAs to 0,
  ###    and removing all 0 rows and -R proteins.
  e <- SpCm2MSnSet(msms.counts,facs)
  if(!validObject(e))
   { galert("\nInvalid MSnSet object.\nPlease check.",
           cont=window)
    return()
  }
  msms.counts <- exprs(e)
  
  ###  Check for a minimum dimensionality (just warn)
  if( nrow(exprs(e)) < 100 )
    if ( ! gconfirm("\nExpression matrix with less than 100 rows.\nProceed ?",
              title="WARNING",icon=c("question"),cont=window) ) 
                     return()
    
  ###  Start EDA report
  #######################
  sink(file=txt.flnm)
  cat("MS/MS EXPLORATORI DATA ANALYSIS\n\n")
  cat("Dataset: ",etit,"\n\n")
  cat("Samples description:\n\n")
  print(samples)

  ###  Start HTML index file
  html.all <- html.tail

  ###  RAW SpC STATISTICS & PLOTS
  insert(output,"\nDistribution of SpC by sample.")
  pdf(file=paste(pdf.flnm,"-RAW-SpCDistrib.pdf",sep=""),width=6,height=11,
      paper="a4")
  layout(mat=matrix(1:2,ncol=1),widths=1,heights=c(0.55,0.45))
  spc.distribs(e,paste(tit,"- raw SpC"))
  dev.off()
  insert(output,paste("  See: ",pdf.flnm,"-RAW-SpCDistrib.pdf",sep=""))
  tfvnm <- count.stats(e)
  cat("\nSample statistics after removing NAs and -R:\n\n")
  cat("SpC matrix dimension:",dim(msms.counts),"\n\n")
  print(tfvnm)
  pca.lst <- msms.eda.plots(e,fnm=NULL,tit=tit,sbtit="RAW")
  
  ###  EXPLORE SAMPLE SIZES
  tcnts <- apply(msms.counts,2,sum)
  medt <- median(tcnts)
  div <- tcnts/medt
  names(div) <- colnames(exprs(e))
  cat("\nSize normalizing divisors scaled to the median sample SpC:\n")
  print(round(div[mnf==condA],2))
  print(round(div[mnf==condB],2))
  insert(output,"\nComputing total SpC normalizing divisors.")
  pdf(file=paste(pdf.flnm,"SizeDivsBarplots.pdf",sep="-"),width=7.5,height=5,
      paper="a4")
  par(cex.main=1,cex.lab=0.8,cex.axis=1)
  omar <- par(mar=c(7,3,2,2))
  barplot(div,las=2,col=cls[facs[,1]])
  abline(h=1,lty=2,col="blue")
  title(main="Size normalizing divisors scaled to the median sample SpC",cex=1,
        line=1)
  dev.off()
  insert(output,paste("  See: ",pdf.flnm,"SizeDivsBarplots.pdf",sep="-"))
  szdiv <- div
  
  ###  Add raw SpC files to file index
  html.all <- add.html(tit,html.all,html.raw)
  
  ###  NORMALIZE BY SIZE, TO TOTAL MEDIAN COUNTS
  if(is.null(samples$offsets))
  { e.n <- norm.counts(e,div)
    ###  NORMALIZED SpC STATISTICS & PLOTS
    tfvnm <- count.stats(e.n)
    cat("\nSample statistics after size normalization:\n\n")
    print(round(tfvnm,1))
    pca.lst <- msms.eda.plots(e.n,fnm=NULL,tit=tit,sbtit="SizeNorm")
    ###  Save normalized SpC matrix
    insert(output,"\nSaving SpC normalized expression data.")
    tblflnm <- paste(tit,"SizeNorm.csv",sep="-")
    write.table(exprs(e.n),file=tblflnm,sep="\t",row.names=FALSE,dec=",")	
    insert(output,paste("  See: ",tblflnm,sep=""))
    ###  Add size normalized SpC files to file index
    html.all <- add.html(tit,html.all,html.sizediv)
  
  ###  NORMALIZE BY GIVEN DIVISORS AND SIZE, SCALED TO THE MEAN
  } else {
    ###  Scale given divisors to the mean
    div <- samples$offsets/mean(samples$offsets)
    names(div) <- colnames(exprs(e))
    cat("\nBio normalizing divisors, scaled to the mean:\n")
    print(round(div[mnf==condA],2))
    print(round(div[mnf==condB],2))
    insert(output,"\nBio normalizing divisors.")
    pdf(file=paste(pdf.flnm,"BioDivsBarplots.pdf",sep="-"),width=7.5,height=5,
        paper="a4")
    par(cex.main=1,cex.lab=0.8,cex.axis=1)
    omar <- par(mar=c(7,3,2,2))
    barplot(div,las=2,col=cls[facs[,1]])
    abline(h=1,lty=2,col="blue")
    title(main="Bio normalizing divisors scaled to the mean",cex=1,line=1)
    dev.off()
    insert(output,paste("  See: ",pdf.flnm,"BioDivsBarplots.pdf",sep="-"))
    ###  Scale Bio + Size divisors 
    div <- samples$offsets*tcnts
    div <- div/mean(div)
    names(div) <- colnames(exprs(e))
    cat("\nBio + size normalizing divisors, scaled to the mean:\n")
    print(round(div[mnf==condA],2))
    print(round(div[mnf==condB],2))
    insert(output,"\nBio + size normalizing divisors.")
    pdf(file=paste(pdf.flnm,"BioSzDivsBarplots.pdf",sep="-"),width=7.5,height=5,
        paper="a4")
    par(cex.main=1,cex.lab=0.8,cex.axis=1)
    omar <- par(mar=c(7,3,2,2))
    barplot(div,las=2,col=cls[facs[,1]])
    abline(h=1,lty=2,col="blue")
    title(main="Bio + size normalizing divisors scaled to the mean",cex=1,line=1)
    dev.off()
    insert(output,paste("  See: ",pdf.flnm,"BioSzDivsBarplots.pdf",sep="-"))

    ###  NORMALIZE SpC MATRIX
    e.n <- norm.counts(e,div)
     
    ###  NORMALIZED SpC STATISTICS 
    tfvnm <- count.stats(e.n)
    cat("\nSample statistics after normalization:\n\n")
    print(round(tfvnm,1))
    pca.lst <- msms.eda.plots(e.n,fnm=NULL,tit=tit,sbtit="BioSzNorm")

    ###  Save normalized SpC matrix
    insert(output,"\nSaving normalized expression data.")
    tblflnm <- paste(tit,"BioSzNorm.csv",sep="-")
    write.table(exprs(e.n),file=tblflnm,sep="\t",row.names=FALSE,dec=",")	
    insert(output,paste("  See: ",tblflnm,sep=""))

    ###  Add bio + size normalized SpC files to file index
    html.all <- add.html(tit,html.all,html.biodiv)
  }
 
  ###  CORRECT BATCH EFFECTS
  if(ncol(facs)>1)
  { batch <- as.factor(facs[,2])
    bmc.counts <- batch.neutralize(exprs(e.n),batch,half=TRUE,
                      sqrt.trans=TRUE)
    e.nbc <- SpCm2MSnSet(bmc.counts,facs[,fnm,drop=FALSE])
    ###  SpC STATISTICS 
    tfvnm <- count.stats(e.nbc)
    cat("\nSample statistics after batch correction:\n\n")
    print(round(tfvnm,1))
    pca.lst <- msms.eda.plots(e.nbc,fnm=NULL,tit=tit,sbtit="BatchCorrected")
    ###  Save batch corrected SpC matrix
    insert(output,"\nSaving batch corrected expression data.")
    tblflnm <- paste(tit,"BatchCorrected.csv",sep="-")
    write.table(exprs(e.nbc),file=tblflnm,sep="\t",row.names=FALSE,dec=",")	
    insert(output,paste("  See: ",tblflnm,sep=""))
    e.n <- e.nbc

    ###  Add batch corrected SpC files to file index
    html.all <- add.html(tit,html.all,html.batch)
  }
		 
  ###   DISPERSION ESTIMATES
  insert(output,"\nComputing dispersion estimates.")
  dsp <- disp.estimates(e.n,facs=facs[,fnm,drop=FALSE],etit=tit,to.pdf=TRUE)
  cat("\nObserved residual dispersion quantiles:\n")
  print(round(dsp,3))
  insert(output,paste("  See: ",tit,"-DispPlots.pdf",sep=""))
  insert(output,"\n   = = =  D O N E  = = =\n")
  sink()

  ###  Add last files and header to HTML file and save it
  html.all <- add.html(tit,html.all,html.last)
  html.all <- c(html.head,html.all)
  writeLines(html.all,paste(tit,"html",sep="."))
  
}   ###  END OF msms.EDA()  

###  GUI - msmsEDA
require(gWidgets)
options(guiToolkit="RGtk2")

##  Data global to controls
tit="MSMS"

##  Functions used by handlers
show.file <- function(flnm,n=-1L)
{ insert(output,readLines(flnm,n),
     font.attr=c(family="monospace",style="tty"))
}

##  Get current contents of exp.nm gedit control
get.expnm <- function()
{ svalue(exp.nm) }

##  Handlers
h.exec <- function(h,...)
{ msms.EDA() }

h.about <- function(h,...)
{ cms <- c("Lable-free LC-MS/MS SpC-based Exploratory Data Analysis     ",
           "",
           "                     StatsCountsMS_EDA v 1.44",
           "",
           "See:",
		   "   Gregori J. et al. (2012) J. Proteomics 75, 3938-3951",
		   "   Gregori J. et al. (2013) J. Proteomics 95, 55-65",
		   "   Gregori J. et al. (2014) J. Proteome Res. Jun 4",
		   "",
		   "(C) J. Gregori, J. Villanueva, A. Sanchez, UB - VHIO, 2014-2018")
  amsg <- paste(cms,collapse="\r\n")		   
  gmessage(amsg,title="About",icon="info") 
}

h.smpls <- function(h,...)
{ samples.flnm <- svalue(h$obj)
  if( file.exists(samples.flnm) )
  { insert(output,"\n\nSamples description:\n")
    show.file(samples.flnm)
  }
  data.flnm <- svalue(g.data.flnm)
  if(!is.null(samples.flnm) && file.exists(samples.flnm) &&
     !is.null(data.flnm) && file.exists(data.flnm) )
    enabled(exec) <- TRUE
}

h.msms <-function(h,...)
{ data.flnm <- svalue(h$obj)
  if( file.exists(data.flnm) )
  { insert(output,"\n\nExpression matrix column names:\n")
    show.file(data.flnm,1)
  }
  samples.flnm <- svalue(g.smpls.flnm)
  if(!is.null(samples.flnm) && file.exists(samples.flnm) &&
     !is.null(data.flnm) && file.exists(data.flnm) )
    enabled(exec) <- TRUE
}    

h.clear <- function(h,...)
{ svalue(output) <- ""
  svalue(g.smpls.flnm) <- ""
  svalue(g.data.flnm) <- ""
  enabled(exec) <- FALSE
}


##  Main window
window <- gwindow("StatsCountsMS - Exploratory Data Analysis", 
                   visible=FALSE)

maing <- ggroup(cont = window, horizontal = FALSE)

## top group organized horizontally
group <- ggroup(cont = maing, horizontal = TRUE)

## Group for file selection organized vertically
grfl <- ggroup(cont = group, horizontal = FALSE, spacing=5)


glabel("Samples file name:",cont = grfl, anchor = c(-1,0))
g.smpls.flnm <- gfilebrowse(text = "Select a file ...",
                        quote = FALSE, handler = h.smpls,
                        type = "open", cont = grfl)
addSpace(grfl,12)
glabel("Counts file name:",cont = grfl, anchor = c(-1,0))
g.data.flnm <- gfilebrowse(text = "Select a file ...",
                        quote = FALSE, handler = h.msms,
                        type = "open", cont = grfl)
addSpace(grfl,18)

##  Separator
addSpace(group,20)
gseparator(horizontal=FALSE,cont=group)
addSpace(group,20)
#addSpring(group)

## Group for buttons and experiment name
grrg <- ggroup(cont = group, horizontal = FALSE)

## Group for the experiment name
addSpace(grrg,25)
grenm <- ggroup(cont = grrg, horizontal = TRUE)
glabel("Experiment name:",cont = grenm, anchor = c(-1,0))
exp.nm <- gedit(text=tit,width=15,cont=grenm)
##  el valor es recuperarà dins l'execució com svalue(exp.nm)

## Group for clear, exec and about buttons
addSpring(grrg)
grbt <- ggroup(cont = grrg, horizontal = TRUE)
addSpace(grbt,10)
## A button to clear dialog
clr_button <- gbutton("Clear",handler=h.clear, cont = grbt)

## A button to initiate the execution
addSpace(grbt,20)
exec <- gaction("Execute EDA",handler=h.exec)
exec <- gbutton(action=exec, cont = grbt)
enabled(exec) <- FALSE

## A button to initiate the execution
addSpace(grbt,20)
about <- gaction("About",handler=h.about)
about <- gbutton(action=about, cont = grbt)
addSpace(grbt,10)


## Area for output
frame <- gframe(" Output & Status ", cont = maing,  
                horizontal = FALSE)
output <- gtext("", cont = frame, expand = TRUE,
                 font.attr=c(family="monospace"))
size(output) <- c(500, 500)

###  Copyright message
glabel("(C) Josep Gregori, UB - VHIO, 2014-2018",cont = maing, anchor = c(-1,0))

visible(window) <- TRUE


###  Construct a MSnSet object from SpC matrix and factors
SpCm2MSnSet <- function(counts,facts)
{ e <- new("MSnSet", exprs=counts)
  pData(e) <- facts
  e <- pp.msms.data(e)
  return(e)
}  

###  Plots of SpC distribution by sample
spc.distribs <- function(e,tit)
{ msms.counts <- exprs(e)
  facs <- pData(e) 
  fact <- as.factor(facs[,1])
  spc.boxplots(msms.counts,fact,minSpC=2)
  spc.densityplots(msms.counts,fact,minSpC=2,main=tit)
}


###  EXPANDED HEAT MAP ON THE NORMALIZED MATRIX
exp.heatmap <- function(e,fnm,h,tit)
{
  cls <- c("red","blue")[as.integer(as.factor(pData(e)[,fnm]))]
  ###  Filter dataset 
  data <- exprs(e)
  data <- data[filter.flags(exprs(e)),]
  ###  Expression to write data on the heatmap
  add.expr <- expression(
  { par(mar = c(margins[1], 0, 0, margins[2]))
    for(j in 1:ncol(SpC.expr))
     text(j,1:nrow(SpC.expr),labels=round(SpC.expr[rowInd,colInd[j]],1),
          cex=0.6,font=2,col=tracecol)
  })
  gtit <- paste("Heat map - ",tit,sep="")
  hm <- heatmap.2jg(t(scale(t(data))),col=greenred(255),trace="none",
            key=FALSE,cexRow=0.6,cexCol=0.7,margins=c(5,6),dendrogram="both",
            add.expr=add.expr,lhei=c(2,h-3),SpC.expr=data,
            ColSideColors=cls)
}

###  EDA PLOTS GIVEN A SpC MATRIX
msms.eda.plots <- function(e,fnm=NULL,tit="MSMS",sbtit="RAW")
{ pdf.flnm <- paste(tit,sbtit,sep="-")
  tt <- paste(tit,sbtit,sep=" - ")
  ###  MEAN SpC SCATTERPLOTS
  if(is.null(fnm))
    fnm <- colnames(pData(e))[1]
  insert(output,"\nPerforming scatterplots on raw SpC.")
  pdf(file=paste(pdf.flnm,fnm,"Scatterplot.pdf",sep="-"),width=5.5,height=11,
    paper="a4")
  par(mfrow=c(2,1),mar=c(5,4.5,4,2)+0.1)
  spc.scatterplot(exprs(e),factor(pData(e)[,fnm]),trans="sqrt",
          minSpC=2,minLFC=1,main=tt)
  spc.scatterplot(exprs(e),factor(pData(e)[,fnm]),trans="log2",
          minSpC=2,minLFC=1,main=tt)
  dev.off()
  insert(output,paste("   See: ",pdf.flnm,"-",fnm,"-Scatterplot.pdf",sep=""))
  ###  PRINCIPAL COMPONENTS ANALYSIS
  insert(output,"\nPerforming PCA on the expression data.")
  pdf(file=paste(pdf.flnm,"PCA.pdf",sep="-"),width=6,height=5,
    paper="a4")
  par(cex.main=1,cex.lab=0.8,cex.axis=0.8)
  fm <- pData(e)
  pca.lst <- counts.pca(e,facs=fm[,1,drop=FALSE],snms=NULL)
  if(ncol(fm)>1) counts.pca(e,facs=fm[,2,drop=FALSE],snms=NULL)
  dev.off()
  ###  Report on variances
  cat("\nPrincipal components analisis on the SpC matrix\n")
  cat("Relative variance of the first four principal components:\n\n")
  print(pca.lst$pc.vars[2,1:4])
  insert(output,paste("   See: ",pdf.flnm,"-PCA.pdf",sep=""))
  ###  SAMPLE HC DENDROGRAMS
  insert(output,"\nPerforming HC on expression data.")
  pdf(file=paste(pdf.flnm,"-HCD.pdf",sep=""),width=5,height=7.5,
      paper="a4")
  par(cex.main=1,cex.lab=0.8,cex.axis=0.8)
  counts.hc(e,facs=fm[,1,drop=FALSE])
  if(ncol(fm)>1) counts.hc(e,facs=fm[,2,drop=FALSE])
  dev.off()
  insert(output,paste("   See: ",pdf.flnm,"-HCD.pdf",sep=""))
  ###  HEAT MAP with filter dataset 
  data <- exprs(e)
  data <- data[filter.flags(exprs(e)),]
  insert(output,"\nHeatmap on the expression data.")
  pdf(file=paste(pdf.flnm,"-HeatMap.pdf",sep=""),paper="a4",
    width=7.5,height=11)
  gtit <- paste("Heat map - ",tit,sbtit,sep=" ")
  cls <- c("red","blue")[as.integer(as.factor(pData(e)[,fnm]))]
  hm <- heatmap(t(scale(t(data))),col=greenred(255),labRow=NA,
                cexCol=0.7,ColSideColors=cls)
  title(main=gtit,line=1,cex=1)
  dev.off()
  insert(output,paste("  See: ",pdf.flnm,"-HeatMap.pdf",
         sep=""))
  ###  EXPANDED HEATMAP
  h <- nrow(data)/(2.54/0.35)
  pdf(file=paste(pdf.flnm,"-ExpandedHeatmap.pdf",sep=""),
      width=7,height=h)
  exp.heatmap(e,fnm,h,tit=pdf.flnm)
  dev.off()
  insert(output,paste("  See: ",pdf.flnm,"-ExpandedHeatmap.pdf",
         sep=""))
  return(pca.lst)
}

### S'ha hagut de resoldre una mancança en heatmap.2
###   per executar 'add.expr'
heatmap.2jg <-
  function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
      distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), symm = FALSE, scale = c("none", 
        "row", "column"), na.rm = TRUE, revC = identical(Colv, 
        "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || 
        scale != "none", col = "heat.colors", colsep, rowsep, 
    sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, 
    notecol = "cyan", na.color = par("bg"), trace = c("column", 
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
    vline = median(breaks), linecol = tracecol, margins = c(5, 
        5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
    key = TRUE, keysize = 1.5, density.info = c("histogram", 
        "density", "none"), denscol = tracecol, symkey = min(x < 
        0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL, 
    xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, 
	SpC.expr = NULL, ...) 
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv)) 
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv)) 
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv)) 
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
        if (missing(col) || is.function(col)) 
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks) 
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) != 
                nc) 
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr")) 
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(add.expr)      ###  Hi havia eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, 
            length(csep)), xright = csep + 0.5 + sepwidth[1], 
            ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, 
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row") 
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column") 
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, "Value", line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}

###  Add HTML module
add.html <- function(tit,html.all,html.lns)
{ html.lns <- sub("#MSMS#",tit,html.lns)
  html.lns <- sub("#MSMS#",tit,html.lns)
  c(html.lns,html.all)
}

# HTML modules
 
# Els mòduls HTML
html.head <- c(
  "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">",                                                    
  "",
  "<html>",
  "<head>" ,
  "",
  "<STYLE TYPE=\"text/css\">",
  "<!--",
  "    Body",
  "     { font-family:Arial,Verdana,helvetica;",
  "       font-size:10pt;",
  "       color:#000000;",
  "       margin:1em;",
  "     }",
  "    P { margin-left:1em; text-align: justify; margin-right:1em;}",
  "    H1 {  text-decoration: underline  }",                                                                              
  "    A { color:maroon; }",
  "    A:hover { color:#ff4422; }",
  "    A.blue { color:blue; }",
  "    A.blue:hover { color:blue; text-decoration:underline; }",
  "    A.navy { color:navy; text-decoration:none; }",
  "    A.navy:hover { text-decoration:underline; }",
  "    A.tit { color:maroon; font-size:18pt; text-decoration:underline;",
  "            font-weight:bold; }",
  "    A.tit:hover { color:#FF6633; }",
  "-->",
  "</STYLE>",
  "",
  "<title>StatsCountsMS_EDA - Exploratory data analysis results</title>",
  "</head>",
  "<body>",
  "",
  "<table width=\"850\" cellpadding=\"15\">",
  "<tr><td align=\"center\"><H1>StatsCountsMS_EDA Graphical User Interface</H1></td></tr>",
  "<tr><td><center>",
  "",
  "<table border=\"1\" bordercolor=\"#00008B\" cellpadding=\"5\" cellspacing=\"1\" summary=\"\">",
  "<caption><b>EXPLORATORY DATA ANALYSIS RESULTS - Output files index</b></caption>",
  "<tr bgcolor=\"#D8CCB1\"><th>FILE / LINK</th><th>Contents</th></tr>",
  "")

html.last <- c(
  "<tr bgcolor=\"#e5ddcc\">",
  "<td><a target=\"_blank\" href=\"#MSMS#_EDA.txt\">#MSMS#_EDA.txt</a></td>",
  "<td><b>Text file with statistic summaries at each step</b></td></tr>",
  "<tr bgcolor=\"#e5ddcc\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-DispPlots.pdf\">#MSMS#-DispPlots.pdf</a></td>",
  "<td>PDF file with plots<br>Final SpC matrix<br><b>Residual dispersion plots</b></td></tr>",
  "")                                                                                                                     
  
html.batch <- c(
  "<tr bgcolor=\"#e5ddcc\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-BatchCorrected-PCA.pdf\">#MSMS#-BatchCorrected-PCA.pdf</a></td>",
  "<td>PDF file with plots<br>Normalized + batch corrected SpC matrix<br><b>Principal components plot</b></td></tr>",
  "<tr bgcolor=\"#e5ddcc\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-BatchCorrected-HCD.pdf\">#MSMS#-BatchCorrected-HCD.pdf</a></td>",
  "<td>PDF file with plots<br>Normalized + batch corrected SpC matrix<br><b>Hierarchical clustering plot</b></td></tr>",
  "<tr bgcolor=\"#e5ddcc\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-BatchCorrected-HeatMap.pdf\">#MSMS#-BatchCorrected-HeatMap.pdf</a></td>",
  "<td>PDF file with plots<br>Normalized + batch corrected SpC matrix<br><b>Heatmap</b></td></tr>",
  "<tr bgcolor=\"#e5ddcc\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-BatchCorrected-ExpandedHeatmap.pdf\">#MSMS#-BatchCorrected-ExpandedHeatmap.pdf</a></td>",
  "<td>PDF file with plots<br>Normalized + batch corrected SpC matrix<br><b>Expanded heatmap with SpC data</b></td></tr>",
  "<tr bgcolor=\"#e5ddcc\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-BatchCorrected-Treat-Scatterplot.pdf\">#MSMS#-BatchCorrected-Treat-Scatterplot.pdf</a></td>",
  "<td>PDF file with plots<br>Normalized + batch corrected SpC matrix<br><b>Treatment level scatterplots</b></td></tr>",
  "<tr bgcolor=\"#e5ddcc\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-BatchCorrected.csv\">#MSMS#-BatchCorrected.csv</a></td>",
  "<td>CSV tab separated file<br>Normalized + batch corrected SpC matrix<br><b>Expression matrix</b></td></tr>",
  "")                                                                                                                     

html.biodiv <- c(
  "<tr bgcolor=\"#dceaea\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-BioDivsBarplots.pdf\">#MSMS#-BioDivsBarplots.pdf</a></td>",
  "<td>PDF file with plots<br>RAW SpC matrix<br><b>Bio normalizing divisors barplot</b></td></tr>",
  "<tr bgcolor=\"#dceaea\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-BioSzDivsBarplots.pdf\">#MSMS#-BioSzDivsBarplots.pdf</a></td>",
  "<td>PDF file with plots<br>RAW SpC matrix<br><b>Bio + size normalizing divisors barplot</b></td></tr>",
  "<tr bgcolor=\"#dceaea\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-BioSzNorm-PCA.pdf\">#MSMS#-BioSzNorm-PCA.pdf</a></td>",
  "<td>PDF file with plots<br>Bio + size normalized SpC matrix<br><b>Principal components plot</b></td></tr>",
  "<tr bgcolor=\"#dceaea\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-BioSzNorm-HCD.pdf\">#MSMS#-BioSzNorm-HCD.pdf</a></td>",
  "<td>PDF file with plots<br>Bio + size normalized SpC matrix<br><b>Hierarchical clustering plot</b></td></tr>",
  "<tr bgcolor=\"#dceaea\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-BioSzNorm-HeatMap.pdf\">#MSMS#-BioSzNorm-HeatMap.pdf</a></td>",
  "<td>PDF file with plots<br>Bio + size normalized SpC matrix<br><b>Heatmap</b></td></tr>",
  "<tr bgcolor=\"#dceaea\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-BioSzNorm-ExpandedHeatmap.pdf\">#MSMS#-BioSzNorm-ExpandedHeatmap.pdf</a></td>",
  "<td>PDF file with plots<br>Bio + size normalized SpC matrix<br><b>Expanded heatmap with SpC data</b></td></tr>",
  "<tr bgcolor=\"#dceaea\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-BioSzNorm-Treat-Scatterplot.pdf\">#MSMS#-BioSzNorm-Treat-Scatterplot.pdf</a></td>",
  "<td>PDF file with plots<br>Bio + size normalized SpC matrix<br><b>Treatment level scatterplots</b></td></tr>",
  "<tr bgcolor=\"#dceaea\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-BioSzNorm.csv\">#MSMS#-BioSzNorm.csv</a></td>",
  "<td>CSV tab separated file<br>Bio + size normalized SpC matrix<br><b>Expression matrix</b></td></tr>",
  "")                                                                                                                     

html.sizediv <- c(
  "<tr bgcolor=\"#eae2d2\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-SizeNorm-PCA.pdf\">#MSMS#-SizeNorm-PCA.pdf</a></td>",
  "<td>PDF file with plots<br>Size normalized SpC matrix<br><b>Principal components plot<b></td></tr>",
  "<tr bgcolor=\"#eae2d2\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-SizeNorm-HCD.pdf\">#MSMS#-SizeNorm-HCD.pdf</a></td>",
  "<td>PDF file with plots<br>Size normalized SpC matrix<br><b>Hierarchical clustering plot<b></td></tr>",
  "<tr bgcolor=\"#eae2d2\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-SizeNorm-HeatMap.pdf\">#MSMS#-SizeNorm-HeatMap.pdf</a></td>",
  "<td>PDF file with plots<br>Size normalized SpC matrix<br><b>Heatmap</b></td></tr>",
  "<tr bgcolor=\"#eae2d2\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-SizeNorm-ExpandedHeatmap.pdf\">#MSMS#-SizeNorm-ExpandedHeatmap.pdf</a></td>",
  "<td>PDF file with plots<br>Size normalized SpC matrix<br><b>Expanded heatmap with SpC data</b></td></tr>",
  "<tr bgcolor=\"#eae2d2\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-SizeNorm-Treat-Scatterplot.pdf\">#MSMS#-SizeNorm-Treat-Scatterplot.pdf</a></td>",
  "<td>PDF file with plots<br>Size normalized SpC matrix<br><b>Treatment level scatterplots</b></td></tr>",
  "<tr bgcolor=\"#eae2d2\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-SizeNorm.csv\">#MSMS#-SizeNorm.csv</a></td>",
  "<td>CSV tab separated file<br>Size normalized SpC matrix<br><b>Expression matrix</b></td></tr>",
  "")                                                                                                                     

html.raw <- c(
  "<tr bgcolor=\"#dedbf4\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-SizeDivsBarplots.pdf\">#MSMS#-SizeDivsBarplots.pdf</a></td>",
  "<td>PDF file with plots<br>RAW SpC matrix<br><b>Size normalizing divisors barplot</b></td></tr>",
  "<tr bgcolor=\"#dedbf4\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-RAW-PCA.pdf\">#MSMS#-RAW-PCA.pdf</a></td>",
  "<td>PDF file with plots<br>RAW SpC matrix<br><b>Principal components plot</b></td></tr>",
  "<tr bgcolor=\"#dedbf4\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-RAW-HCD.pdf\">#MSMS#-RAW-HCD.pdf</a></td>",
  "<td>PDF file with plots<br>RAW SpC matrix<br><b>Hierarchical clustering plot</b></td></tr>",
  "<tr bgcolor=\"#dedbf4\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-RAW-HeatMap.pdf\">#MSMS#-RAW-HeatMap.pdf</a></td>",
  "<td>PDF file with plots<br>RAW SpC matrix<br><b>Heatmap</b></td></tr>",
  "<tr bgcolor=\"#dedbf4\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-RAW-ExpandedHeatmap.pdf\">#MSMS#-RAW-ExpandedHeatmap.pdf</a></td>",
  "<td>PDF file with plots<br>RAW SpC matrix<br><b>Expanded heatmap with SpC data</b></td></tr>",
  "<tr bgcolor=\"#dedbf4\">",
  "<td><a target=\"_blank\" href=\"#MSMS#-RAW-Treat-Scatterplot.pdf\">#MSMS#-RAW-Treat-Scatterplot.pdf</a></td>",
  "<td>PDF file with plots<br>Size normalized SpC matrix<br><b>Treatment level scatterplots</b></td></tr>",
  "")                                                                                                                     

html.tail <- c(
  "</table>",
  "",
  "</center></td></tr>",
  "<tr><td align=\"center\"><font size=\"1\">Josep Gregori - UB+VHIO Bcn - 2014</font></td></tr>",
  "</table>",
  "",
  "</body>",
  "</html>")
