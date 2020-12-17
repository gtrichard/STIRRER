#### Dependecies ####

library(ggplot2)
library(ggrastr)
library(cowplot)
library(devtools)
library(vroom)
library(limma) 
library(edgeR)
library(splitstackshape)
library(RColorBrewer)
library(kableExtra)
library(tidyr)
library(cowplot)
source_gist("524eade46135f6348140")

#### Compute Parental Mean ####

voom_lcpm = function (counts, design = NULL, lib.size = NULL, normalize.method = "none", 
                      block = NULL, correlation = NULL, weights = NULL, span = 0.5, 
                      plot = FALSE, save.plot = FALSE)
  {
  out <- list()
  if (is(counts, "DGEList")) {
    out$genes <- counts$genes
    out$targets <- counts$samples
    if (is.null(design) && diff(range(as.numeric(counts$sample$group))) > 
        0) 
      design <- model.matrix(~group, data = counts$samples)
    if (is.null(lib.size)) 
      lib.size <- counts$samples$lib.size * counts$samples$norm.factors
    counts <- counts$counts
  }
  else {
    isExpressionSet <- suppressPackageStartupMessages(is(counts, 
                                                         "ExpressionSet"))
    if (isExpressionSet) {
      if (length(Biobase::fData(counts))) 
        out$genes <- Biobase::fData(counts)
      if (length(Biobase::pData(counts))) 
        out$targets <- Biobase::pData(counts)
      counts <- Biobase::exprs(counts)
    }
    else {
      counts <- as.matrix(counts)
    }
  }
  n <- nrow(counts)
  if (n < 2L) 
    stop("Need at least two genes to fit a mean-variance trend")
  if (is.null(design)) {
    design <- matrix(1, ncol(counts), 1)
    rownames(design) <- colnames(counts)
    colnames(design) <- "GrandMean"
  }
  if (is.null(lib.size)) 
    lib.size <- colSums(counts)
  #    y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  y <- counts
  y <- normalizeBetweenArrays(y, method = normalize.method)
  fit <- lmFit(y, design, block = block, correlation = correlation, 
               weights = weights)
  if (is.null(fit$Amean)) 
    fit$Amean <- rowMeans(y, na.rm = TRUE)
  NWithReps <- sum(fit$df.residual > 0L)
  if (NWithReps < 2L) {
    if (NWithReps == 0L) 
      warning("The experimental design has no replication. Setting weights to 1.")
    if (NWithReps == 1L) 
      warning("Only one gene with any replication. Setting weights to 1.")
    out$E <- y
    out$weights <- y
    out$weights[] <- 1
    out$design <- design
    if (is.null(out$targets)) 
      out$targets <- data.frame(lib.size = lib.size)
    else out$targets$lib.size <- lib.size
    return(new("EList", out))
  }
  sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
  sy <- sqrt(fit$sigma)
  allzero <- rowSums(counts) == 0
  if (any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx, sy, f = span)
  if (plot) {
    plot(sx, sy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )", 
         pch = 16, cex = 0.25)
    title("voom: Mean-variance trend")
    lines(l, col = "red")
  }
  f <- approxfun(l, rule = 2, ties = list("ordered", mean))
  if (fit$rank < ncol(design)) {
    j <- fit$pivot[1:fit$rank]
    fitted.values <- fit$coef[, j, drop = FALSE] %*% t(fit$design[, 
                                                                  j, drop = FALSE])
  }
  else {
    fitted.values <- fit$coef %*% t(fit$design)
  }
  fitted.cpm <- 2^fitted.values
  fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
  fitted.logcount <- log2(fitted.count)
  w <- 1/f(fitted.logcount)^4
  dim(w) <- dim(fitted.logcount)
  out$E <- y
  out$weights <- w
  out$design <- design
  if (is.null(out$targets)) 
    out$targets <- data.frame(lib.size = lib.size)
  else out$targets$lib.size <- lib.size
  if (save.plot) {
    out$voom.xy <- list(x = sx, y = sy, xlab = "log2( count size + 0.5 )", 
                        ylab = "Sqrt( standard deviation )")
    out$voom.line <- l
  }
  new("EList", out)
}

#### Compute DEGs ####

getDEG <- function(path, files, group, columns, contr.matrix, lfc = 0.584962501, genome = ".*", normalisation = "RLE", col = 
                     c( brewer.pal(8, "Blues")[6:8], brewer.pal(8, "Greens")[6:8], brewer.pal(8, "Reds")[6:8], 
                        brewer.pal(8, "Oranges")[6:8], brewer.pal(8, "Greys")[6:8] ), keep_homeo = F, homeo_path = NULL)
  {
  # LOAD DATA

  out <- list()
  data_list = list()
  count = 0
  for (f in files){
    count = count + 1 
    file_full =  paste0(path, files[[count]])
    data_list[[f]] = suppressMessages( vroom(file_full, skip = 1, col_names = T) )
  }
  
  genes_expression =  do.call(cbind,data_list)
  to_keep = seq(from = 7, to = ncol(genes_expression), by = 7)
  genes_expression = genes_expression[,c(1:6,to_keep)]
  
  colnames(genes_expression) = columns
  genes_expression_basic = genes_expression
  rownames(genes_expression_basic) = genes_expression_basic$gene
  
  
  # FILTER GENES BY HOMEOLOGUES
  
  if (keep_homeo == T){
    homeo = vroom(homeo_path, col_names = F)
    homeoA = as.data.frame(homeo$X1)
    homeoC = as.data.frame(homeo$X2)
    colnames(homeoA) = "gene"
    colnames(homeoC) = "gene"
    homeo = as.data.frame(rbind(homeoA, homeoC))
    homeo = as.data.frame(homeo[order(homeo),])
    colnames(homeo) = "gene"
    homeo = unique(homeo)
    rownames(homeo) = homeo$gene
    
    genes_expression_basic = merge(genes_expression_basic, homeo, by = "row.names" )
    genes_expression_basic = genes_expression_basic[,3:ncol(genes_expression_basic)-1]
    colnames(genes_expression_basic)[1] = "gene"
    rownames(genes_expression_basic) = genes_expression_basic$gene
  }
  
  # GENES FILTRATION DEPENDING ON CHOSEN GENOME
  
  genes_expression_basic = genes_expression_basic[grep(genes_expression_basic$chr, pattern = genome),]
  
  # LOAD GENES ANNOTATION FROM FEATURECOUNTS OUTPUT
  
  genes             = genes_expression_basic[,c(1:6)]
  genes[,2]         = sub(x = genes[,2],pattern = "^.*;",replacement = "")
  genes[,3]         = sub(x = genes[,3],pattern = ";.*$",replacement = "")
  genes[,4]         = sub(x = genes[,4],pattern = "^.*;",replacement = "")
  genes[,5]         = sub(x = genes[,5],pattern = "^.*;",replacement = "")
  
  
  # CREATE DGE OBJECT
  
  counts            = genes_expression_basic[,c(7:ncol(genes_expression_basic))]
  row.names(counts) <- genes_expression_basic[,1]
  lib.size          = apply( genes_expression_basic[,7:ncol(genes_expression_basic)], 2, sum )
  norm.factors      = rep(1, length(lib.size))
  samplenames      = colnames(counts)
  samples           = as.data.frame( cbind(files, group) )
  rownames(samples) = samplenames
  
  x = DGEList(counts = counts,
              lib.size = lib.size,
              genes = genes,
              samples = samples, 
              remove.zeros = F)
  
  a <- as.data.frame(data.frame(cbind(samplenames, lib.size) ))
  a
  rownames(x$counts) = x$genes$gene
  
  
  # KEEP ONLY EXPRESSED GENES AND PLOT
  
  keep.exprs <- filterByExpr(x, group=group)
  x_expr <- x[keep.exprs,, keep.lib.sizes=FALSE]
  
  nsamples <- ncol(x)
  
  L <- mean(x$samples$lib.size) * 1e-6
  M <- median(x$samples$lib.size) * 1e-6
  lcpm.cutoff <- log2(10/M + 2/L)
  
  pdf(NULL)
  par(mfrow=c(1,2))
  lcpm <- cpm(x, log=TRUE)
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.40), las=2, main="", xlab="")
  title(main=paste("All Genes n =", nrow(x)), xlab="Log-cpm")
  abline(v=lcpm.cutoff, lty=3)
  for (i in 1:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", samplenames, text.col=col, bty="n")
  lcpm <- cpm(x_expr, log=TRUE)
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.40), las=2, main="", xlab="")
  title(main=paste("Expressed Genes n =", nrow(x_expr)), xlab="Log-cpm")
  abline(v=lcpm.cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", samplenames, text.col=col, bty="n")
  out$plot_expr_filt <- recordPlot()
  invisible(dev.off())
  
  # NORM FACTOR COMPUTATION AND PLOT MDS
  
  x_norm <- calcNormFactors(x_expr, method = normalisation)
  
  lcpm <- cpm(x_norm, log=TRUE)
  
  col.group <- as.factor(group)
  levels(col.group) <- col
  col.group <- as.character(col.group)
  
  pdf(NULL)
  par(mfrow=c(1,2))
  
  if (length(grep("A", genome)) == 1){
    lcpmA = lcpm[grep(row.names(lcpm), pattern = "^A"), ]  
    boxplot(lcpmA, las=2, col=col, main="A subgenome")
  }
  
  if (length(grep("C", genome)) == 1){
    lcpmC = lcpm[grep(row.names(lcpm), pattern = "^C"), ]
    boxplot(lcpmC, las=2, col=col, main="C subgenome")
  }
  
  par(mfrow=c(1,1))
  plotMDS(lcpm, labels=group, col=col)
  out$plot_mds <- recordPlot()
  invisible(dev.off())
  
  # DEG ANALYSIS
  
  pdf(NULL)
  par(mfrow=c(1,2))
  
  design <- model.matrix(~0+group)
  colnames(design) <- gsub("group", "", colnames(design))
  design
  contr.matrix
  
  v <- voom(x_norm, design, plot=TRUE)
  
  vfit <- lmFit(v, design)
  vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
  efit <- eBayes(vfit)
  plotSA(efit)
  out$plot_fit <- recordPlot()
  invisible(dev.off())
  
  tfit <- treat(vfit, lfc = lfc)
  dt <- decideTests(tfit)
  
  dtA = dt[grep(row.names(dt), pattern = "^A"), ]
  dtC = dt[grep(row.names(dt), pattern = "^C"), ]
  
  
#  par(mfrow=c(1,ncol(dt)))
#  for (i in 1:ncol(dt)){
#    plotMD(tfit, column=i, status=dt[,i], main=colnames(tfit)[i], xlim=c(-8,13))
#  }
  
  
  # OUTPUT RESULTS
  

  out$tfit <- tfit
  out$dt  <- dt
  out$dtA <- dtA
  out$dtC <- dtC
  out$x_norm <- x_norm
  out$cpm <- cpm(x_norm)
  out$lcpm <- cpm(x_norm, log = T)
  out$keep.exprs <- keep.exprs
  return(out)
}

#### Compute DEGs for parental mean ####

getDEG_pmix <- function(path, files, group, columns, contr.matrix, lfc, genome, normalisation, col, keep_homeo)
  {
  # LOAD DATA

  data_list = list()
  count = 0
  for (f in files){
    count = count + 1 
    file_full =  paste0(path, files[[count]])
    data_list[[f]] = suppressMessages( vroom(file_full, skip = 1, col_names = T) )
  }
  
  genes_expression =  do.call(cbind,data_list)
  to_keep = seq(from = 7, to = ncol(genes_expression), by = 7)
  genes_expression = genes_expression[,c(1:6,to_keep)]
  
  colnames(genes_expression) = columns
  genes_expression_basic = genes_expression
  rownames(genes_expression_basic) = genes_expression_basic$gene
  
  
  # FILTER GENES BY HOMEOLOGUES
  
  if (keep_homeo == T){
    homeo = vroom(homeo_path, col_names = F)
    homeoA = as.data.frame(homeo$X1)
    homeoC = as.data.frame(homeo$X2)
    colnames(homeoA) = "gene"
    colnames(homeoC) = "gene"
    homeo = as.data.frame(rbind(homeoA, homeoC))
    homeo = as.data.frame(homeo[order(homeo),])
    colnames(homeo) = "gene"
    homeo = unique(homeo)
    rownames(homeo) = homeo$gene
    
    genes_expression_basic = merge(genes_expression_basic, homeo, by = "row.names" )
    genes_expression_basic = genes_expression_basic[,3:ncol(genes_expression_basic)-1]
    colnames(genes_expression_basic)[1] = "gene"
    rownames(genes_expression_basic) = genes_expression_basic$gene
    
  }
  
  
  # GENES FILTRATION DEPENDING ON CHOSEN GENOME
  
  genes_expression_basic = genes_expression_basic[grep(genes_expression_basic$chr, pattern = genome),]
  
  # LOAD GENES ANNOTATION FROM FEATURECOUNTS OUTPUT
  
  genes             = genes_expression_basic[,c(1:6)]
  genes[,2]         = sub(x = genes[,2],pattern = "^.*;",replacement = "")
  genes[,3]         = sub(x = genes[,3],pattern = ";.*$",replacement = "")
  genes[,4]         = sub(x = genes[,4],pattern = "^.*;",replacement = "")
  genes[,5]         = sub(x = genes[,5],pattern = "^.*;",replacement = "")
  
  
  # CREATE DGE OBJECT
  
  counts            = genes_expression_basic[,c(7:ncol(genes_expression_basic))]
  row.names(counts) <- genes_expression_basic[,1]
  lib.size          = apply( genes_expression_basic[,7:ncol(genes_expression_basic)], 2, sum )
  norm.factors      = rep(1, length(lib.size))
  samplenames      = colnames(counts)
  samples           = as.data.frame( cbind(files, group) )
  rownames(samples) = samplenames
  
  x = DGEList(counts = counts,
              lib.size = lib.size,
              genes = genes,
              samples = samples, 
              remove.zeros = F)
  
  cat(cbind(samplenames,lib.size))
  
  rownames(x$counts) = x$genes$gene
  
  
  # KEEP ONLY EXPRESSED GENES AND PLOT
  
  keep.exprs <- filterByExpr(x, group=group)
  x_expr <- x[keep.exprs,, keep.lib.sizes=FALSE]
  
  nsamples <- ncol(x)
  
  L <- mean(x$samples$lib.size) * 1e-6
  M <- median(x$samples$lib.size) * 1e-6
  lcpm.cutoff <- log2(10/M + 2/L)
  
  par(mfrow=c(1,2))
  lcpm <- cpm(x, log=TRUE)
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.40), las=2, main="", xlab="")
  title(main=paste("All Genes n =", nrow(x)), xlab="Log-cpm")
  abline(v=lcpm.cutoff, lty=3)
  for (i in 1:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", samplenames, text.col=col, bty="n")
  lcpm <- cpm(x_expr, log=TRUE)
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.40), las=2, main="", xlab="")
  title(main=paste("Expressed Genes n =", nrow(x_expr)), xlab="Log-cpm")
  abline(v=lcpm.cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", samplenames, text.col=col, bty="n")
  
  
  # NORM FACTOR COMPUTATION AND PLOT MDS
  
  x_norm <- calcNormFactors(x_expr, method = normalisation)
  
  lcpm <- cpm(x_norm, log=TRUE)
  
  col.group <- as.factor(group)
  levels(col.group) <- col
  col.group <- as.character(col.group)
  
  par(mfrow=c(1,2))
  
  if (length(grep("A", genome)) == 1){
    lcpmA = lcpm[grep(row.names(lcpm), pattern = "^A"), ]  
    boxplot(lcpmA, las=2, col=col, main="A subgenome")
  }
  
  if (length(grep("C", genome)) == 1){
    lcpmC = lcpm[grep(row.names(lcpm), pattern = "^C"), ]
    boxplot(lcpmC, las=2, col=col, main="C subgenome")
  }
  
  
  par(mfrow=c(1,1))
  plotMDS(lcpm, labels=group, col=col)
  
  # PARENTAL MIX CALCULATION
  
  lcpm <- cpm(x_norm, log=TRUE)
  
  x_norm$counts = lcpm
  for (r in 1:nrow(x_norm$counts)){
    x_norm$counts[r,7] = mean( x_norm$counts[r,1],x_norm$counts[r,7])
    x_norm$counts[r,8] = mean( x_norm$counts[r,2],x_norm$counts[r,8])
    x_norm$counts[r,9] = mean( x_norm$counts[r,3],x_norm$counts[r,9])
  }
  
  colnames(x_norm$counts)[7] = "Parental_mix R1"
  colnames(x_norm$counts)[8] = "Parental_mix R2"
  colnames(x_norm$counts)[9] = "Parental_mix R3"
  
  x_norm$counts = x_norm$counts[, -c(1,2,3) ]
  
  x_norm$samples[7,2] = mean(c(x_norm$samples[1,"lib.size"],x_norm$samples[7,"lib.size"]))
  x_norm$samples[8,2] = mean(c(x_norm$samples[2,"lib.size"],x_norm$samples[8,"lib.size"]))
  x_norm$samples[9,2] = mean(c(x_norm$samples[3,"lib.size"],x_norm$samples[9,"lib.size"]))
  
  x_norm$samples[7,3] = mean(c(x_norm$samples[1,"norm.factors"],x_norm$samples[7,"norm.factors"]))
  x_norm$samples[8,3] = mean(c(x_norm$samples[2,"norm.factors"],x_norm$samples[8,"norm.factors"]))
  x_norm$samples[9,3] = mean(c(x_norm$samples[3,"norm.factors"],x_norm$samples[9,"norm.factors"]))
  
  x_norm$samples[,1] = as.character(x_norm$samples[,1])
  x_norm$samples[7,1] = "Parental_mix"
  x_norm$samples[8,1] = "Parental_mix"
  x_norm$samples[9,1] = "Parental_mix"
  x_norm$samples[,1] = as.factor(x_norm$samples[,1])
  
  row.names(x_norm$samples)[7] = "Parental_mix R1"
  row.names(x_norm$samples)[8] = "Parental_mix R2"
  row.names(x_norm$samples)[9] = "Parental_mix R3"
  
  x_norm$samples = x_norm$samples[ -c(1, 2, 3),]
  
  group2 = substr(colnames(x_norm$counts), 1, nchar(colnames(x_norm$counts))-3)
  
  # DEG ANALYSIS
  
  par(mfrow=c(1,2))
  
  design <- model.matrix(~0+group2)
  colnames(design) <- gsub("group2", "", colnames(design))
  
  contr.matrix <- makeContrasts(
    Hybrid_vs_Parental_mix = Triploid_Ch_AAC - Parental_mix,
    levels = colnames(design)
  )
  
  
  v <- voom_lcpm(x_norm, design, plot=TRUE)
  
  vfit <- lmFit(v, design)
  vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
  efit <- eBayes(vfit)
  
  plotSA(efit)
  
  tfit <- treat(vfit, lfc = lfc)
  dt <- decideTests(tfit)
  
  dtA = dt[grep(row.names(dt), pattern = "^A"), ]
  dtC = dt[grep(row.names(dt), pattern = "^C"), ]
  
  par(mfrow=c(1,ncol(dt)))
  for (i in 1:ncol(dt)){
    plotMD(tfit, column=i, status=dt[,i], main=colnames(tfit)[i], xlim=c(-8,13))
  }
  
  
  # OUTPUT RESULTS
  
  out <- list()
  out$tfit <- tfit
  out$dt  <- dt
  out$dtA <- dtA
  out$dtC <- dtC
  out$x_norm <- x_norm
  out$cpm <- 2^x_norm$counts
  out$lcpm <- x_norm$counts
  out$keep.exprs <- keep.exprs
  return(out)
}

#### Plot DEGs & correlation ####

plotDEG <- function(data, cols_1, cols_2, name_1, name_2, comp_name, genome)
  {
  
  
  comparison = data$lcpm[,c(1:2)]
  
  if(genome == "A"){
    sum = summary(data$dtA) 
    data$lcpm = data$lcpm[grep( x = row.names(data$lcpm), pattern = "^A"),]
    comparison = comparison[grep( x = row.names(comparison), pattern = "^A"),]
  }
  
  0.
  if(genome == "C"){
    sum = summary(data$dtC) 
    data$lcpm = data$lcpm[grep( x = row.names(data$lcpm), pattern = "^C"),]
    comparison = comparison[grep( x = row.names(comparison), pattern = "^C"),]
  }
  
  comparison = data$lcpm[,c(1:2)]
  
  for (r in 1:nrow(data$lcpm)){
    comparison[r,1] = mean(data$lcpm[r,cols_1])
    comparison[r,2] = mean(data$lcpm[r,cols_2])
  }
  
  comparison = as.data.frame(comparison)
  colnames(comparison) = c(name_1, name_2)
  
  if (genome == "A"){
    comparison[grep( x = row.names(data$lcpm), pattern = "^A"),"genome"] =  paste("A -", sum[3,comp_name], "Up &", sum[1,comp_name],"Down DEGs")
  }
  
  if (genome == "C"){
    comparison[grep( x = row.names(data$lcpm), pattern = "^C"),"genome"] =  paste("C -", sum[3,comp_name], "Up &", sum[1,comp_name],"Down DEGs")
  }
  
  comparison["differential"] = "No diff. expression"
  comparison[which(data$dt[,comp_name] ==  1),"differential"] = "DEG Up"
  comparison[which(data$dt[,comp_name] == -1),"differential"] = "DEG Down"
  
  comparison = na.omit(comparison)
  
  ggplot(comparison, aes(x = comparison[,2], y = comparison[,1]))+
    facet_wrap(genome ~ .) +
    geom_point_rast(aes( color = differential, alpha = differential)) +
    theme_cowplot() +
    scale_color_manual(paste0("Genes, n=",nrow(comparison)), values = c("blue","red","black")) +
    scale_alpha_manual(paste0("Genes, n=",nrow(comparison)), values = c(1,1,0.1)) +
    xlab(paste0(name_2, " log2(cpm)")) +
    ylab(paste0(name_1, " log2(cpm)"))+
    ggtitle(paste0(name_1, " vs ", name_2)) +
    stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE,ypos = max(comparison[,1:2])+0.5, xpos = min(comparison[,1:2])) +
    geom_smooth(method="lm",se=FALSE, color = "grey") +
    theme(plot.title = element_text(size = 12))
}