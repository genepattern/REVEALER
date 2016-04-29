  # Check if all neccessary packages are installed and if not install those missing

# list.of.packages <- c("parmigene", "verification", "sn", "e1071", "misc3d", "MASS", "PerformanceAnalytics", "smacof", "RColorBrewer", "ppcor", "optparse")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)

  # Check if all neccessary packages are installed and if not install those missing

   suppressPackageStartupMessages(library(parmigene))
   suppressPackageStartupMessages(library(verification))
   suppressPackageStartupMessages(library(sn))
   suppressPackageStartupMessages(library(e1071))
   suppressPackageStartupMessages(library(misc3d))
   suppressPackageStartupMessages(library(MASS))
   suppressPackageStartupMessages(library(PerformanceAnalytics))
   suppressPackageStartupMessages(library(smacof))
   suppressPackageStartupMessages(library(NMF))
   suppressPackageStartupMessages(library(RColorBrewer))
   suppressPackageStartupMessages(library(ppcor))
   suppressPackageStartupMessages(library(optparse))

### Revealer 2.0  February 6, 2013

  REVEALER.v1 <- function(
   # REVEALER (Repeated Evaluation of VariablEs conditionAL Entropy and Redundancy) is an analysis method specifically suited
   # to find groups of genomic alterations that match in a complementary way, a predefined functional activation, dependency of
   # drug response “target” profile. The method starts by considering, if available, already known genomic alterations (“seed”)
   # that are the known “causes” or are known or assumed “associated” with the target. REVEALER starts by finding the genomic
   # alteration that best matches the target profile “conditional” to the known seed profile using the conditional mutual information.
   # The newly discovered alteration is then merged with the seed to form a new “summary” feature, and then the process repeats itself
   # finding additional complementary alterations that explain more and more of the target profile.

   ds1,                                          # Dataset that contains the "target"
   target.name,                                  # Target feature (row in ds1)
   target.match = "positive",                    # Use "positive" to match the higher values of the target, "negative" to match the lower values
   ds2,                                          # Features dataset.
   seed.names = NULL,                            # Seed(s) name(s)
   exclude.features = NULL,                      # Features to exclude for search iterations
   max.n.iter = 5,                               # Maximun number of iterations
   pdf.output.file,                              # PDF output file
   count.thres.low = NULL,                       # Filter out features with less than count.thres.low 1's
   count.thres.high = NULL,                      # Filter out features with more than count.thres.low 1's

   n.markers = 30,                               # Number of top hits to show in heatmap for every iteration
   locs.table.file = NULL)                       # Table with chromosomal location for each gene symbol (optional)

{
   identifier = "REVEALER"                      # Documentation suffix to be added to output file    
   n.perm = 10                                  # Number of permutations (x number of genes) for computing p-vals and FRDs
   save_preprocessed_features_dataset = NULL    # Save preprocessed features dataset    
   seed.combination.op = "max"                  # Operation to consolidate and summarize seeds to one vector of values
   assoc.metric = "IC"                          # Assoc. Metric: "IC" information coeff.; "COR" correlation.
   normalize.features = F                       # Feature row normalization: F or "standardize" or "0.1.rescaling"
   top.n = 1                                    # Number of top hits in each iteration to diplay in Landscape plot
   max.n = 2                                    # Maximum number of iterations to diplay in Landscape plot
#   save.cmi = F                                # Save feature search cond. assoc. metric for each feature at every iteration
   phen.table = NULL                            # Table with phenotypes for each sample (optional)
   phen.column = NULL                           # Column in phen.table containing the relevant phenotype info
   phen.selected = NULL                         # Use only samples of these phenotypes in REVEALER analysis
   produce.lanscape.plot = F                    # Produce multi-dimensional scaling projection plot
   character.scaling = 1.25                      # Character scaling for heatmap
   r.seed = 34578                               # Random number generation seed
   consolidate.identical.features = F           # Consolidate identical features: F or "identical" or "similar" 
   cons.features.hamming.thres = NULL           # If consolidate.identical.features = "similar" then consolidate features within this Hamming dist. thres.

# ------------------------------------------------------------------------------------------------------------------------------
#   print(paste("ds1:", ds1))
#   print(paste("target.name:", target.name))
#   print(paste("target.match:", target.match))
#   print(paste("ds2:", ds2))                        
#   print(paste("seed.names:", seed.names))
#   print(paste("exclude.features:", exclude.features))
#   print(paste("max.n.iter:", max.n.iter))
#   print(paste("pdf.output.file:", pdf.output.file))
#   print(paste("n.perm:", n.perm))                  
#   print(paste("count.thres.low:", count.thres.low))
#   print(paste("count.thres.high:", count.thres.high))
#   print(paste("identifier:", identifier))                  
#   print(paste("n.markers:", n.markers))   
#   print(paste("locs.table.file:", locs.table.file))

# ------------------------------------------------------------------------------------------------------------------------------
   # Load libraries

   if (is.null(seed.names)) seed.names <- "NULLSEED"
  
#   pdf.file <- paste(output.dir, "REVEALER_", target.name, "_", paste(seed.names, collapse="."), "_", identifier, ".pdf", sep="")
#   gct.file <- paste(output.dir, "REVEALER_", target.name, "_", paste(seed.names, collapse="."), "_", identifier, ".gct", sep="")
   
   pdf(file=pdf.output.file, height=14, width=8.5)
   set.seed(r.seed)

   # Read table with HUGO gene symbol vs. chr location
   
   if (!is.null(locs.table.file)) {
      locs.table <- read.table(locs.table.file, header=T, sep="\t", skip=0, colClasses = "character")
    }
   
   # Define color map

   mycol <- vector(length=512, mode = "numeric")
   for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
   for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
   mycol <- rev(mycol)
   ncolors <- length(mycol)

   # Read datasets

   dataset.1 <- MSIG.Gct2Frame(filename = ds1)
   m.1 <- data.matrix(dataset.1$ds)
   row.names(m.1) <- dataset.1$row.names
   Ns.1 <- ncol(m.1)  
   sample.names.1 <- colnames(m.1) <- dataset.1$names

   dataset.2 <- MSIG.Gct2Frame(filename = ds2)
   m.2 <- data.matrix(dataset.2$ds)
   row.names(m.2) <- dataset.2$row.names
   Ns.2 <- ncol(m.2)  
   sample.names.2 <- colnames(m.2) <- dataset.2$names

    # exclude samples with target == NA

   target <- m.1[target.name,]
#   print(paste("initial target length:", length(target)))      
   locs <- seq(1, length(target))[!is.na(target)]
   m.1 <- m.1[,locs]
   sample.names.1 <- sample.names.1[locs]
#   print(paste("target length after excluding NAs:", ncol(m.1)))     

   overlap <- intersect(sample.names.1, sample.names.2)
   length(overlap)
   locs1 <- match(overlap, sample.names.1)
   locs2 <- match(overlap, sample.names.2)
   m.1 <- m.1[, locs1]
   m.2 <- m.2[, locs2]

#   print(dim(m.1))
#   print(dim(m.2))

#   print(colnames(m.1)[1:10])
#   print(colnames(m.2)[1:10])

   # Filter samples with only the selected phenotypes 

   if (!is.null(phen.selected)) {
      samples.table <- read.table(phen.table, header=T, row.names=1, sep="\t", skip=0)
#      print(colnames(samples.table)[1:10])
      table.sample.names <- row.names(samples.table)
      locs1 <- match(colnames(m.2), table.sample.names)
      phenotype <- as.character(samples.table[locs1, phen.column])
#      print("phenotype:")
#      print(phenotype)
      
      locs2 <- NULL
      for (k in 1:ncol(m.2)) {   
         if (!is.na(match(phenotype[k], phen.selected))) {
            locs2 <- c(locs2, k)
         }
      }
      length(locs2)
      m.1 <- m.1[, locs2]
      m.2 <- m.2[, locs2]
      phenotype <- phenotype[locs2]
      table(phenotype)
#      print(dim(m.1))
#      print(dim(m.2))
    }

   # Define target

   target <- m.1[target.name,]
   if (target.match == "negative") {
      ind <- order(target, decreasing=F)
   } else {
      ind <- order(target, decreasing=T)
   }
   target <- target[ind]
   m.2 <- m.2[, ind]

   if (!is.na(match(target.name, row.names(m.2)))) {
     loc <- match(target.name, row.names(m.2))
     m.2 <- m.2[-loc,]
   }

   MUT.count <- AMP.count <- DEL.count <- 0
   for (i in 1:nrow(m.2)) {
      temp <- strsplit(row.names(m.2)[i], split="_")
      temp <- strsplit(temp[[1]][length(temp[[1]])], split=" ")
      suffix <- temp[[1]][1]
      if (!is.na(suffix)) {
         if (suffix == "MUT") MUT.count <- MUT.count + 1
         if (suffix == "AMP") AMP.count <- AMP.count + 1
         if (suffix == "DEL") DEL.count <- DEL.count + 1
     }
   }
#   print(paste("Initial number of features ", nrow(m.2), " MUT:", MUT.count, " AMP:", AMP.count, " DEL:", DEL.count))
   
   # Eliminate flat, sparse or features that are too dense

   if (!is.null(count.thres.low) && !is.null(count.thres.high)) {
      sum.rows <- rowSums(m.2)
#      print(paste("rowSums:"))
#      print(sum.rows[1:100])
#      print(paste(sum.rows))

#      print(paste("m.2:"))
#      print(paste(m.2[40000,]))      
      
      seed.flag <- rep(0, nrow(m.2))
      if (seed.names != "NULLSEED") {
         locs <- match(seed.names, row.names(m.2))
         locs <- locs[!is.na(locs)]
         seed.flag[locs] <- 1
      }
      retain <- rep(0, nrow(m.2))
      for (i in 1:nrow(m.2)) {
         if ((sum.rows[i] >= count.thres.low) && (sum.rows[i] <= count.thres.high)) retain[i] <- 1
#         if ((sum.rows[i] >= 0) && (sum.rows[i] <= 100)) retain[i] <- 1          
#         if (i < 10) print(paste(sum.rows[i], count.thres.low, count.thres.high, retain[i]))
         if (seed.flag[i] == 1) retain[i] <- 1
      }

#      print(paste("retain[1:100]"))
#      print(retain[1:100])
      
      m.2 <- m.2[retain == 1,]
#      print(paste("Number of features kept:", sum(retain), "(", signif(100*sum(retain)/length(retain), 3), " percent)"))
  }
   
   # Normalize features and define seeds

  if (normalize.features == "standardized") {
      for (i in 1:nrow(m.2)) {
         mean.row <- mean(m.2[i,])
         sd.row <- ifelse(sd(m.2[i,]) == 0, 0.1*mean.row, sd(m.2[i,]))
         m.2[i,] <- (m.2[i,] - mean.row)/sd.row
       }
   } else if (normalize.features == "0.1.rescaling") {
      for (i in 1:nrow(m.2)) {
         max.row <- max(m.2[i,])
         min.row <- min(m.2[i,])
         range.row <- ifelse(max.row == min.row, 1, max.row - min.row)
         m.2[i,] <- (m.2[i,] - min.row)/range.row
       }
    }
   
  if (seed.names == "NULLSEED") {
#     seed <- as.vector(c(rep(0, ncol(m.2) -1), 0.000001))
     seed <- as.vector(rep(0, ncol(m.2)))      
     seed.vectors <- as.matrix(t(seed))
  } else {
#      print("Location(s) of seed(s):")
#      print(match(seed.names, row.names(m.2)))
      if (length(seed.names) > 1) {
         seed <- apply(m.2[seed.names,], MARGIN=2, FUN=seed.combination.op)
         seed.vectors <- as.matrix(m.2[seed.names,])
      } else {
         seed <- m.2[seed.names,]
         seed.vectors <- as.matrix(t(m.2[seed.names,]))
      }
      locs <- match(seed.names, row.names(m.2))
      locs
     m.2 <- m.2[-locs,]
     dim(m.2)
   }

  if (length(table(m.2[1,])) > ncol(m.2)*0.5) { # continuous target
     feature.type <- "continuous"
  } else {
     feature.type <- "discrete"
  }
    
  # Exclude user-specified features 
   
   if (!is.null(exclude.features)) {
      locs <- match(exclude.features, row.names(m.2))
      locs <- locs[!is.na(locs)]
#      print(paste("Excluding feature:", row.names(m.2)[locs]))
      m.2 <- m.2[-locs,]
    }

# Eliminate  features with .1 suffix (duplicates)
   
#   locs <- NULL
#   for (i in 1:nrow(m.2)) {
#      temp <- strsplit(row.names(m.2)[i], split="\\.")
#      suffix <- temp[[1]][length(temp[[1]])]
#      if (suffix == "1") locs <- c(locs, i)
#      if (suffix == "2") locs <- c(locs, i)      
#   }
#   if (!is.null(locs)) m.2 <- m.2[-locs,]
   
  #  Consolidate identical features

   if (consolidate.identical.features == "identical") {  # This is a very fast way to eliminate perfectly identical features compared with what we do below in "similar"
      dim(m.2)
      summary.vectors <- apply(m.2, MARGIN=1, FUN=paste, collapse="")
      ind <- order(summary.vectors)
      summary.vectors <- summary.vectors[ind]
      m.2 <- m.2[ind,]
      taken <- i.count <- rep(0, length(summary.vectors))
      i <- 1
      while (i <= length(summary.vectors)) {
        j <- i + 1
        while ((summary.vectors[i] == summary.vectors[j]) & (j <= length(summary.vectors))) {
            j <- j + 1
         }
        i.count[i] <- j - i
        if (i.count[i] > 1) taken[seq(i + 1, j - 1)] <- 1
        i <- j
      }
   
      if (sum(i.count) != length(summary.vectors)) stop("ERROR")     # Add counts in parenthesis
      row.names(m.2) <- paste(row.names(m.2), " (", i.count, ")", sep="")
      m.2 <- m.2[taken == 0,]
      dim(m.2)

   } else if (consolidate.identical.features == "similar") { # this uses the hamming distance to consolidate similar features up to the Hamming dist. threshold 
      hamming.matrix <- hamming.distance(m.2)
      taken <- rep(0, nrow(m.2))
      for (i in 1:nrow(m.2)) {
         if (taken[i] == 0) { 
            similar.features <- row.names(m.2)[hamming.matrix[i,] <= cons.features.hamming.thres]
            if (length(similar.features) > 1) {
               row.names(m.2)[i]  <- paste(row.names(m.2)[i], " [", length(similar.features), "]", sep="") # Add counts in brackets
               locs <- match(similar.features, row.names(m.2))
               taken[locs] <- 1
               taken[i] <- 0
            }
        }
      }
      m.2 <- m.2[taken == 0,]
     dim(m.2)
   }

   MUT.count <- AMP.count <- DEL.count <- 0
   for (i in 1:nrow(m.2)) {
      temp <- strsplit(row.names(m.2)[i], split="_")
      temp <- strsplit(temp[[1]][length(temp[[1]])], split=" ")      
      suffix <- temp[[1]][1]
      if (!is.na(suffix)) {
         if (suffix == "MUT") MUT.count <- MUT.count + 1
         if (suffix == "AMP") AMP.count <- AMP.count + 1
         if (suffix == "DEL") DEL.count <- DEL.count + 1
     }
   }
#   print(paste("Number of features (after filtering and consolidation)",
#   nrow(m.2), " MUT:", MUT.count, " AMP:", AMP.count, " DEL:", DEL.count))
   
   # Add location info

   if (!is.null(locs.table.file)) {
      gene.symbol <- row.names(m.2)
      chr <- rep(" ", length(gene.symbol))
      for (i in 1:length(gene.symbol)) {
        temp1 <- strsplit(gene.symbol[i], split="_")
        temp2 <- strsplit(temp1[[1]][1], split="\\.")
        gene.symbol[i] <- ifelse(temp2[[1]][1] == "", temp1[[1]][1], temp2[[1]][1])
        loc <- match(gene.symbol[i], locs.table[,"Approved.Symbol"])
        chr[i] <- ifelse(!is.na(loc), locs.table[loc, "Chromosome"], " ")
       }
      row.names(m.2)  <- paste(row.names(m.2), " ", chr, " ", sep="")
#      print(paste("Total unmatched to chromosomal locations:", sum(chr == " "), "out of ", nrow(m.2), "features"))
    }

   # save filtered and consolidated file

    if (!is.null(save_preprocessed_features_dataset)) {
       write.gct.2(gct.data.frame = data.frame(m.2), descs = row.names(m.2), filename = save_preprocessed_features_dataset)
   }
   
   # Compute MI and % explained with original seed(s)
   
   median_target <- median(target)
    if (target.match == "negative") {
      target.locs <- seq(1, length(target))[target <= median_target]
    } else {
      target.locs <- seq(1, length(target))[target > median_target]
    }

   cmi.orig.seed <- cmi.orig.seed.cum <- pct_explained.orig.seed <- pct_explained.orig.seed.cum <- vector(length=length(seed.names), mode="numeric")
   if (length(seed.names) > 1) {
      seed.cum <- NULL
      for (i in 1:nrow(seed.vectors)) {
         y <- seed.vectors[i,]
         cmi.orig.seed[i] <- assoc(target, y, assoc.metric)
         pct_explained.orig.seed[i] <- sum(y[target.locs])/length(target.locs)
         seed.cum <- apply(rbind(seed.vectors[i,], seed.cum), MARGIN=2, FUN=seed.combination.op)
         cmi.orig.seed.cum[i] <- assoc(target, seed.cum, assoc.metric)
         pct_explained.orig.seed.cum[i] <- sum(seed.cum[target.locs])/length(target.locs)
      }
   } else {
       y <- as.vector(seed.vectors)
       seed.cum <- y
       cmi.orig.seed <- cmi.orig.seed.cum <- assoc(target, y, assoc.metric)
       pct_explained.orig.seed <- sum(y[target.locs])/length(target.locs)
   }
    cmi.seed.iter0 <- assoc(target, seed, assoc.metric)
    pct_explained.seed.iter0 <- sum(seed[target.locs])/length(target.locs) 

   # CMI iterations

   cmi <- pct_explained <- cmi.names <- matrix(0, nrow=nrow(m.2), ncol=max.n.iter)
   cmi.seed <- pct_explained.seed <- vector(length=max.n.iter, mode="numeric")
   seed.names.iter <- vector(length=max.n.iter, mode="character")
   seed.initial <- seed
   seed.iter <- matrix(0, nrow=max.n.iter, ncol=ncol(m.2))

   target.rand <- matrix(target, nrow=n.perm, ncol=ncol(m.2), byrow=TRUE)
   for (i in 1:n.perm) target.rand[i,] <- sample(target.rand[i,])

   for (iter in 1:max.n.iter) {

      cmi.rand <- matrix(0, nrow=nrow(m.2), ncol=n.perm)     
      for (k in 1:nrow(m.2)) {
#         if (k %% 100 == 0) print(paste("Iter:", iter, " feature #", k, " out of ", nrow(m.2)))
         y <- m.2[k,]
         cmi[k, iter] <- cond.assoc(target, y, seed, assoc.metric)
         for (j in 1:n.perm) {
            cmi.rand[k, j] <- cond.assoc(target.rand[j,], y, seed, assoc.metric)
          }
       }

      if (target.match == "negative") {
         ind <- order(cmi[, iter], decreasing=F)
      } else {
         ind <- order(cmi[, iter], decreasing=T)
      }
      cmi[, iter] <- cmi[ind, iter]
      cmi.names[, iter] <- row.names(m.2)[ind]
      pct_explained[iter] <- sum(m.2[cmi.names[1, iter], target.locs])/length(target.locs)

#      if (save.cmi == T) {
#          tab <- cbind(cmi.names[, iter], cmi[, iter])
#          colnames(tab) <- c("Feature", paste("Conditional ", assoc.metric))
#          file <-  paste(gct.file, ".ITER.CMI.", iter, ".txt", sep="")
#          col.names <- c("Rank", colnames(tab))
#          col.names <- paste(colnames(tab), collapse = "\t")
#          write(noquote(col.names), file = file, append = F, ncolumns = length(col.names))
#          write.table(tab, file=file, quote=F, col.names = F, row.names = T, append = T, sep="\t")
#      }
      
      # Estimate p-vals and FDRs

      p.val <- FDR <- FDR.lower <- rep(0, nrow(m.2))
      for (i in 1:nrow(m.2)) p.val[i] <- sum(cmi.rand >  cmi[i, iter])/(nrow(m.2)*n.perm)
      FDR <- p.adjust(p.val, method = "fdr", n = length(p.val))
      FDR.lower <- p.adjust(1 - p.val, method = "fdr", n = length(p.val))
    
      for (i in 1:nrow(m.2)) {
         if (cmi[i, iter] < 0) {
            p.val[i] <- 1 - p.val[i]
            FDR[i] <- FDR.lower[i]
         }
         p.val[i] <- signif(p.val[i], 2)
         FDR[i] <- signif(FDR[i], 2)
      }
      p.zero.val <- paste("<", signif(1/(nrow(m.2)*n.perm), 2), sep="")
      p.val <- ifelse(p.val == 0, rep(p.zero.val, length(p.val)), p.val)

      # Make a heatmap of the n.marker top hits in this iteration

      size.mid.panel <- length(seed.names) + iter
      pad.space <- 15
      
      nf <- layout(matrix(c(1, 2, 3, 4), 4, 1, byrow=T), 1, c(2, size.mid.panel, ceiling(n.markers/2) + 4, pad.space), FALSE)
      cutoff <- 2.5
      x <- as.numeric(target)         
      x <- (x - mean(x))/sd(x)         
      ind1 <- which(x > cutoff)
      ind2 <- which(x < -cutoff)
      x[ind1] <- cutoff
      x[ind2] <- -cutoff
      V1 <- ceiling(ncolors * (x + cutoff)/(cutoff*2))
      par(mar = c(1, 22, 2, 12))
      image(1:length(target), 1:1, as.matrix(V1), col=mycol, zlim=c(0, ncolors), axes=FALSE, main=paste("REVEALER - Iteration:", iter), sub = "", xlab= "", ylab="",
            font=2, family="")
      axis(2, at=1:1, labels=paste("TARGET: ", target.name), adj= 0.5, tick=FALSE,las = 1, cex=1, cex.axis=character.scaling, font.axis=1,
           line=0, font=2, family="")
      axis(4, at=1:1, labels="  IC", adj= 0.5, tick=FALSE, las = 1, cex.axis=character.scaling, font.axis=1, line=0, font=2, family="")
 
      if (iter == 1) {
            V0 <- rbind(seed.vectors, seed + 2)
            cmi.vals <- c(cmi.orig.seed, cmi.seed.iter0)
            cmi.vals <- signif(cmi.vals, 2)
            row.names(V0) <- c(paste("SEED: ", seed.names), "SUMMARY SEED:")
            V0 <- apply(V0, MARGIN=2, FUN=rev)
       } else {
         V0 <- rbind(seed.vectors, m.2[seed.names.iter[1:(iter-1)],], seed + 2)
         row.names(V0) <- c(paste("SEED:   ", seed.names), paste("ITERATION ", seq(1, iter-1), ":  ",
                                                                 seed.names.iter[1:(iter-1)], sep=""), "SUMMARY SEED:")
         cmi.vals <- c(cmi.orig.seed, cmi[1, 1:iter-1], cmi.seed[iter-1])
         cmi.vals <- signif(cmi.vals, 2)
         pct.vals <- c(signif(pct_explained.orig.seed, 2), signif(pct_explained[seq(1, iter - 1)], 2), signif(pct_explained.seed[iter - 1], 2))         
         V0 <- apply(V0, MARGIN=2, FUN=rev)
       }

      all.vals <- cmi.vals
      par(mar = c(1, 22, 0, 12))
      if (feature.type == "discrete") {
         image(1:ncol(V0), 1:nrow(V0), t(V0), zlim = c(0, 3), col=c(brewer.pal(9, "Blues")[3], brewer.pal(9, "Blues")[9],
                                                                    brewer.pal(9, "Greys")[2], brewer.pal(9, "Greys")[5]),
               axes=FALSE, main="", sub = "", xlab= "", ylab="")
      } else { # continuous
         for (i in 1:length(V0[,1])) {
            x <- as.numeric(V0[i,])
            V0[i,] <- (x - mean(x))/sd(x)
            max.v <- max(max(V0[i,]), -min(V0[i,]))
            V0[i,] <- ceiling(ncolors * (V0[i,] - (- max.v))/(1.001*(max.v - (- max.v))))
         }
         image(1:ncol(V0), 1:nrow(V0), t(V0), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
      }
      axis(2, at=1:nrow(V0), labels=row.names(V0), adj= 0.5, tick=FALSE, las = 1, cex.axis=character.scaling, line=0, font=2, family="")
      axis(4, at=1:nrow(V0), labels=rev(all.vals), adj= 0.5, tick=FALSE, las = 1, cex.axis=character.scaling, line=0, font=2, family="")

      V0 <- m.2[cmi.names[1:n.markers, iter],]
      V0 <- apply(V0, MARGIN=2, FUN=rev)
      par(mar = c(6, 22, 3, 12))
      if (feature.type == "discrete") {
         image(1:ncol(V0), 1:nrow(V0), t(V0), zlim = c(0, 1), col=c(brewer.pal(9, "Blues")[3], brewer.pal(9, "Blues")[9]),
               axes=FALSE, main=paste("Top", n.markers, "Matches"), sub = "", xlab= "", ylab="")
      } else { # continuous
         for (i in 1:length(V0[,1])) {
            cutoff <- 2.5
            x <- as.numeric(V0[i,])            
            x <- (x - mean(x))/sd(x)         
            ind1 <- which(x > cutoff)
            ind2 <- which(x < -cutoff)
            x[ind1] <- cutoff
            x[ind2] <- -cutoff
            V0[i,] <- ceiling(ncolors * (x + cutoff)/(cutoff*2))
         }
         image(1:ncol(V0), 1:nrow(V0), t(V0), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
      }
      axis(2, at=1:nrow(V0), labels=row.names(V0), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.9*character.scaling, line=0, font=2, family="")

      all.vals <- paste(signif(cmi[1:n.markers, iter], 2), p.val[1:n.markers], FDR[1:n.markers], sep="   ")

      axis(4, at=nrow(V0)+0.4, labels=" CIC    p-val   FDR", adj= 0.5, tick=FALSE, las = 1, cex.axis=0.8*character.scaling, line=0, font=2, family="")
      axis(4, at=c(seq(1, nrow(V0) - 1), nrow(V0) - 0.2), labels=rev(all.vals), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.9*character.scaling, line=0, font=2, family="")
      axis(1, at=1:ncol(V0), labels=colnames(V0), adj= 0.5, tick=FALSE,las = 3, cex=1, cex.axis=0.45*character.scaling,  line=0, font=2, family="")

     # second page shows the same markers clustered in groups with similar profiles

     tab <- m.2[cmi.names[1:n.markers, iter],]
     all.vals <- paste(signif(cmi[1:n.markers, iter], 2), p.val[1:n.markers], FDR[1:n.markers], sep="   ")

     # Cluster and make heatmap of n.markers top hits in groups

     tab2 <- tab + 0.001

     k.min <- 2
     k.max <- 10
     NMF.models <- nmf(tab2, seq(k.min, k.max), nrun = 25, method="brunet", seed=9876)
     plot(NMF.models)
     NMF.sum <- summary(NMF.models)

     k.vec <- seq(k.min, k.max, 1)
     cophen <- NMF.sum[, "cophenetic"]

#     plot(k.vec, cophen)
#     points(k.vec, cophen, type="l")

     peak <- c(0, rep(0, k.max-2), 0)
     for (h in 2:(length(cophen) - 1)) if (cophen[h - 1] < cophen[h] & cophen[h] > cophen[h + 1]) peak[h] <- 1

     if (sum(peak) == 0) {
        if (cophen[1] > cophen[length(cophen)]) {
           k <- k.min
         } else {
           k <- k.max
         }
     } else {
        k.peaks <- k.vec[peak == 1]
        k <- rev(k.peaks)[1]
     }
#     print(paste("Number of groups:", k))
     NMF.model <- nmf(tab2, k, method="brunet", seed=9876)
     classes <- predict(NMF.model, "rows")
     table(classes)
     lens <- table(classes)

     lens2 <- ifelse(lens <= 5, 5, lens)
     lens2[length(lens2)] <- lens2[length(lens2)] + 5


     def.par <- par(no.readonly = TRUE)       
     nf <- layout(matrix(seq(1, k+3), k+3, 1, byrow=T), 1, c(3.5, size.mid.panel, lens2, pad.space), FALSE)      

      cutoff <- 2.5
      x <- as.numeric(target)         
      x <- (x - mean(x))/sd(x)         
      ind1 <- which(x > cutoff)
      ind2 <- which(x < -cutoff)
      x[ind1] <- cutoff
      x[ind2] <- -cutoff
      V1 <- ceiling(ncolors * (x + cutoff)/(cutoff*2))
      par(mar = c(1, 22, 1, 12))
      image(1:length(target), 1:1, as.matrix(V1), col=mycol, zlim=c(0, ncolors), axes=FALSE, main=paste("REVEALER - Iteration:", iter), sub = "", xlab= "", ylab="",
            font=2, family="")
     axis(2, at=1:1, labels=paste("TARGET: ", target.name), adj= 0.5, tick=FALSE,las = 1, cex=1, cex.axis=character.scaling, font.axis=1,
           line=0, font=2, family="")
      axis(4, at=1:1, labels="  IC", adj= 0.5, tick=FALSE, las = 1, cex.axis=0.8*character.scaling, font.axis=1, line=0, font=2, family="")
 
      if (iter == 1) {
            V0 <- rbind(seed.vectors, seed + 2)
            cmi.vals <- c(cmi.orig.seed, cmi.seed.iter0)
            cmi.vals <- signif(cmi.vals, 2)
            row.names(V0) <- c(paste("SEED: ", seed.names), "SUMMARY SEED:")
            V0 <- apply(V0, MARGIN=2, FUN=rev)
       } else {
         V0 <- rbind(seed.vectors, m.2[seed.names.iter[1:(iter-1)],], seed + 2)
         row.names(V0) <- c(paste("SEED:   ", seed.names), paste("ITERATION ", seq(1, iter-1), ":  ",
                                                                 seed.names.iter[1:(iter-1)], sep=""), "SUMMARY SEED:")
         cmi.vals <- c(cmi.orig.seed, cmi[1, 1:iter-1], cmi.seed[iter-1])
         cmi.vals <- signif(cmi.vals, 2)
         pct.vals <- c(signif(pct_explained.orig.seed, 2), signif(pct_explained[seq(1, iter - 1)], 2), signif(pct_explained.seed[iter - 1], 2))         
         V0 <- apply(V0, MARGIN=2, FUN=rev)
       }

      all.vals <- cmi.vals
      par(mar = c(1, 22, 0, 12))
      if (feature.type == "discrete") {
         image(1:ncol(V0), 1:nrow(V0), t(V0), zlim = c(0, 3), col=c(brewer.pal(9, "Blues")[3], brewer.pal(9, "Blues")[9],
                                                                    brewer.pal(9, "Greys")[2], brewer.pal(9, "Greys")[5]),
               axes=FALSE, main="", sub = "", xlab= "", ylab="")
      } else { # continuous
         for (i in 1:length(V0[,1])) {
            x <- as.numeric(V0[i,])
            V0[i,] <- (x - mean(x))/sd(x)
            max.v <- max(max(V0[i,]), -min(V0[i,]))
            V0[i,] <- ceiling(ncolors * (V0[i,] - (- max.v))/(1.001*(max.v - (- max.v))))
         }
         image(1:ncol(V0), 1:nrow(V0), t(V0), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
      }
      axis(2, at=1:nrow(V0), labels=row.names(V0), adj= 0.5, tick=FALSE, las = 1, cex.axis=character.scaling, line=0, font=2, family="")
      axis(4, at=1:nrow(V0), labels=rev(all.vals), adj= 0.5, tick=FALSE, las = 1, cex.axis=character.scaling, line=0, font=2, family="")

    # groups of abnormalities
     all.vals <- paste(signif(cmi[1:n.markers, iter], 2), p.val[1:n.markers], FDR[1:n.markers], sep="   ")
      
#     for (h in 1:k) {
     for (h in sort(unique(classes))) {      
         if (lens[h] == 1) {
            V0 <- t(as.matrix(tab[classes == h,]))
          } else {
            V0 <- tab[classes == h,]
            V0 <- apply(V0, MARGIN=2, FUN=rev)
          }
         r.names <- row.names(tab)[classes == h]
         all.vals0 <- all.vals[classes == h]         
         if (h < k) {
            par(mar = c(0.5, 22, 1, 12))
          } else {
            par(mar = c(3, 22, 1, 12))
          }
         if (feature.type == "discrete") {
            if (lens[h] == 1) {           
               image(1:ncol(V0), 1, t(V0), zlim = c(0, 1), col=c(brewer.pal(9, "Blues")[3], brewer.pal(9, "Blues")[9]),
                     axes=FALSE, main=paste("Top Matches. Group:", h, "(iter ", iter, ")"), sub = "", xlab= "", ylab="", cex.main=0.8)
             } else {
               image(1:ncol(V0), 1:nrow(V0), t(V0), zlim = c(0, 1), col=c(brewer.pal(9, "Blues")[3], brewer.pal(9, "Blues")[9]),
                     axes=FALSE, main=paste("Top Matches. Group:", h, "(iter ", iter, ")"), sub = "", xlab= "", ylab="", cex.main=0.8)
             }
         } else { # continuous
            for (i in 1:length(V0[,1])) {
               cutoff <- 2.5
               x <- as.numeric(V0[i,])            
               x <- (x - mean(x))/sd(x)         
               ind1 <- which(x > cutoff)
               ind2 <- which(x < -cutoff)
               x[ind1] <- cutoff
               x[ind2] <- -cutoff
               V0[i,] <- ceiling(ncolors * (x + cutoff)/(cutoff*2))
            }
            image(1:ncol(V0), 1:nrow(V0), t(V0), zlim = c(0, ncolors), col=mycol, axes=FALSE,
               main=paste("Top Matches. Group:", h), sub = "", xlab= "", ylab="")
         }

         if (lens[h] == 1) {
           axis(2, at=1, labels=rev(r.names), adj= 0.5, tick=FALSE, las = 1, cex.axis=character.scaling, line=-0.7, font=2, family="")
           axis(4, at=1+0.4, labels=" CIC     p-val     FDR", adj= 0.5, tick=FALSE, las = 1, cex.axis=0.8*character.scaling, line=-0.7, font=2, family="")
           axis(4, at=1, labels=all.vals0, adj= 0.5, tick=FALSE, las = 1, cex.axis=character.scaling, line=-0.7, font=2, family="")
         } else {
            axis(2, at=1:nrow(V0), labels=rev(r.names), adj= 0.5, tick=FALSE, las = 1, cex.axis=character.scaling, line=-0.7, font=2, family="")
            axis(4, at=nrow(V0)+0.4, labels=" CIC     p-val     FDR", adj= 0.5, tick=FALSE, las = 1, cex.axis=0.8*character.scaling, line=-0.7, font=2, family="")            
            axis(4, at=c(seq(1, nrow(V0) - 1), nrow(V0) - 0.2), labels=rev(all.vals0), adj= 0.5, tick=FALSE, las = 1, cex.axis=character.scaling, line=-0.7, font=2, family="")
         }
#         axis(1, at=1:ncol(V0), labels=colnames(V0), adj= 0.5, tick=FALSE,las = 3, cex=1, cex.axis=0.35, line=-0.7, font=2, family="")
      }
      par(def.par)
      
#      V0 <- m.2[cmi.names[1:n.markers, iter],]
#      file <-  paste(gct.file, ".ITER.", iter, ".txt", sep="")
#      col.names <- paste(colnames(V0), collapse = "\t")
#      col.names <- paste("Feature", col.names, sep= "\t")
#      write(noquote(col.names), file = file, append = F, ncolumns = length(col.names))
#      write.table(V0, file=file, quote=F, col.names = F, row.names = T, append = T, sep="\t")
      
      # Update seed

      seed.names.iter[iter] <- cmi.names[1, iter] # top hit from this iteration
      seed <- apply(rbind(seed, m.2[seed.names.iter[iter],]), MARGIN=2, FUN=seed.combination.op)
      seed.iter[iter,] <- seed
      cmi.seed[iter] <- assoc(target, seed, assoc.metric)
      pct_explained.seed[iter] <- sum(seed[target.locs])/length(target.locs)
      
    } # end of iterations loop

   # Final summary figures -----------------------------------------------------------------------------

   summ.panel <- length(seed.names) + 2 * max.n.iter + 2
  
   legend.size <- 4
   pad.space <- 30 - summ.panel - legend.size

   nf <- layout(matrix(c(1, 2, 3, 0), 4, 1, byrow=T), 1, c(2, summ.panel, legend.size, pad.space), FALSE)

   cutoff <- 2.5
   x <- as.numeric(target)         
   x <- (x - mean(x))/sd(x)         
   ind1 <- which(x > cutoff)
   ind2 <- which(x < -cutoff)
   x[ind1] <- cutoff
   x[ind2] <- -cutoff
   V1 <- ceiling(ncolors * (x + cutoff)/(cutoff*2))

   par(mar = c(1, 22, 2, 12))
   image(1:length(target), 1:1, as.matrix(V1), zlim=c(0, ncolors), col=mycol, axes=FALSE, main=paste("REVEALER - Results"),
         sub = "", xlab= "", ylab="", font=2, family="")
  
   axis(2, at=1:1, labels=paste("TARGET:  ", target.name), adj= 0.5, tick=FALSE,las = 1, cex=1, cex.axis=character.scaling,  line=0, font=2, family="")
   axis(4, at=1:1, labels="  IC   ", adj= 0.5, tick=FALSE, las = 1, cex.axis=character.scaling, line=0, font=2, family="")

#   print(paste("max.n.iter:", max.n.iter))
#   print(paste("dim m.2:", dim(m.2)))
#   print(paste("dim seed.iter:", dim(seed.iter)))  

   V0 <- rbind(seed.vectors + 2, seed.cum) 
   for (i in 1:max.n.iter) {
      V0 <- rbind(V0,
                  m.2[seed.names.iter[i],] + 2,
                  seed.iter[i,])
   }

   row.names.V0 <- c(paste("SEED:   ", seed.names), "SUMMARY SEED:")
   for (i in 1:max.n.iter) {
      row.names.V0 <- c(row.names.V0, paste("ITERATION ", i, ":  ", seed.names.iter[i], sep=""), paste("SUMMARY FEATURE ", i, ":  ", sep=""))
   }
   row.names(V0) <- row.names.V0

   cmi.vals <- c(cmi.orig.seed, cmi.orig.seed.cum[length(seed.names)])                 
   for (i in 1:max.n.iter) {
      cmi.vals <- c(cmi.vals, as.vector(cmi[1, i]), cmi.seed[i])
    }
   cmi.vals <- signif(cmi.vals, 2)
   all.vals <-cmi.vals
      
   V0 <- apply(V0, MARGIN=2, FUN=rev)
#   par(mar = c(7, 22, 2, 12))
   par(mar = c(7, 22, 0, 12))
   if (feature.type == "discrete") {  
       image(1:ncol(V0), 1:nrow(V0), t(V0), zlim = c(0, 3),
       col=c(brewer.pal(9, "Greys")[2], brewer.pal(9, "Greys")[5],                          
             brewer.pal(9, "Blues")[3], brewer.pal(9, "Blues")[9]), axes=FALSE, main="", sub = "", xlab= "", ylab="")
   } else { # continuous
      for (i in 1:nrow(V0)) {
         cutoff <- 2.5
         x <- as.numeric(V0[i,])
         x <- (x - mean(x))/sd(x)         
         ind1 <- which(x > cutoff)
         ind2 <- which(x < -cutoff)
         x[ind1] <- cutoff
         x[ind2] <- -cutoff
         V0[i,] <- ceiling(ncolors * (x + cutoff)/(cutoff*2))
      }
      image(1:ncol(V0), 1:nrow(V0), t(V0), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
   }
   axis(2, at=1:nrow(V0), labels=row.names(V0), adj= 0.5, tick=FALSE, las = 1, cex.axis=character.scaling, line=0, font=2, family="")
   axis(4, at=1:nrow(V0), labels=rev(all.vals), adj= 0.5, tick=FALSE, las = 1, cex.axis=character.scaling, line=0, font=2, family="")
   axis(1, at=1:ncol(V0), labels=colnames(V0), adj= 0.5, tick=FALSE,las = 3, cex=1, cex.axis=0.4*character.scaling, line=0, font=2, family="")

        # Legend
      
#      nf <- layout(matrix(c(1, 2, 3), 3, 1, byrow=T), 1, c(1,2,2),  FALSE)
      par.mar <- par("mar")
      par(mar = c(3, 35, 8, 10))
      leg.set <- seq(-cutoff, cutoff, 2*cutoff/100)
      image(1:101, 1:1, as.matrix(leg.set), zlim=c(-cutoff, cutoff), col=mycol, axes=FALSE, main="",
          sub = "", xlab= "", ylab="",font=2, family="")
      ticks <- c(-2, -1, 0, 1, 2)
      tick.cols <- rep("black", 5)
      tick.lwd <- c(1,1,2,1,1)
      locs <- NULL
      for (k in 1:length(ticks)) locs <- c(locs, which.min(abs(ticks[k] - leg.set)))
      axis(1, at=locs, labels=ticks, adj= 0.5, tick=T, cex=0.8, cex.axis=1, line=0, font=2, family="")
      mtext("Standardized Target Profile", cex=0.8, side = 1, line = 3.5, outer=F)
      par(mar = par.mar)

#   print("end of heatmap")

   V0 <- rbind(target, seed.vectors, seed.cum) 
   for (i in 1:max.n.iter) {
      V0 <- rbind(V0,
                  m.2[seed.names.iter[i],],
                  seed.iter[i,])
   }
   V0.colnames <- colnames(V0)
   V0 <- cbind(V0, c(1, all.vals))
   colnames(V0) <- c(V0.colnames, "IC")

   row.names.V0 <- c(target.name, seed.names, "SUMMARY SEED:")
   for (i in 1:max.n.iter) {
      row.names.V0 <- c(row.names.V0, seed.names.iter[i], paste("SUMMARY FEATURE ", i, ":  ", sep=""))
   }
   row.names(V0) <- row.names.V0
  
#  write.gct.2(gct.data.frame = V0, descs = row.names(V0), filename = gct.file)
#  print("gct results file has been produced")


#  file <-  paste(gct.file, ".txt", sep="")
#  col.names <- paste(colnames(V0), collapse = "\t")
#  col.names <- paste("Feature", col.names, sep= "\t")
#  write(noquote(col.names), file = file, append = F, ncolumns = length(col.names))
#  write.table(V0, file=file, quote=F, col.names = F, row.names = T, append = T, sep="\t")

  # Version without summaries ----------------------------------------------------

     # Final summary figures -----------------------------------------------------------------------------

   summ.panel <- length(seed.names) + max.n.iter + 2
  
   legend.size <- 4
   pad.space <- 30 - summ.panel - legend.size
   
   nf <- layout(matrix(c(1, 2, 3, 0), 4, 1, byrow=T), 1, c(2, summ.panel, legend.size, pad.space), FALSE)

   cutoff <- 2.5
   x <- as.numeric(target)         
   x <- (x - mean(x))/sd(x)         
   ind1 <- which(x > cutoff)
   ind2 <- which(x < -cutoff)
   x[ind1] <- cutoff
   x[ind2] <- -cutoff
   V1 <- ceiling(ncolors * (x + cutoff)/(cutoff*2))

   par(mar = c(1, 22, 2, 12))
   image(1:length(target), 1:1, as.matrix(V1), zlim=c(0, ncolors), col=mycol, axes=FALSE, main=paste("REVEALER - Results"),
         sub = "", xlab= "", ylab="", font=2, family="")
  
   axis(2, at=1:1, labels=paste("TARGET:  ", target.name), adj= 0.5, tick=FALSE,las = 1, cex=1, cex.axis=character.scaling,  line=0, font=2, family="")
   axis(4, at=1:1, labels="  IC   ", adj= 0.5, tick=FALSE, las = 1, cex.axis=character.scaling, line=0, font=2, family="")

#   print(paste("max.n.iter:", max.n.iter))
#   print(paste("dim m.2:", dim(m.2)))
#   print(paste("dim seed.iter:", dim(seed.iter)))  

   V0 <- seed.vectors + 2
   for (i in 1:max.n.iter) {
      V0 <- rbind(V0,
                  m.2[seed.names.iter[i],] + 2)
   }
   V0 <- rbind(V0, seed.iter[max.n.iter,])

   row.names.V0 <- c(paste("SEED:   ", seed.names))
   for (i in 1:max.n.iter) {
      row.names.V0 <- c(row.names.V0, paste("ITERATION ", i, ":  ", seed.names.iter[i], sep=""))
   }
   row.names(V0) <- c(row.names.V0, "FINAL SUMMARY")

   cmi.vals <- cmi.orig.seed
   for (i in 1:max.n.iter) {
      cmi.vals <- c(cmi.vals, as.vector(cmi[1, i]))
    }
   cmi.vals <- c(cmi.vals, cmi.seed[max.n.iter])
   cmi.vals <- signif(cmi.vals, 2)
   all.vals <-cmi.vals
      
   V0 <- apply(V0, MARGIN=2, FUN=rev)
#    par(mar = c(7, 22, 2, 12))
   par(mar = c(7, 22, 0, 12))   
   
   if (feature.type == "discrete") {  
       image(1:ncol(V0), 1:nrow(V0), t(V0), zlim = c(0, 3),
       col=c(brewer.pal(9, "Greys")[2], brewer.pal(9, "Greys")[5],                          
              brewer.pal(9, "Blues")[3], brewer.pal(9, "Blues")[9]), axes=FALSE, main="", sub = "", xlab= "", ylab="")
   } else { # continuous
      for (i in 1:nrow(V0)) {
         cutoff <- 2.5
         x <- as.numeric(V0[i,])
         x <- (x - mean(x))/sd(x)         
         ind1 <- which(x > cutoff)
         ind2 <- which(x < -cutoff)
         x[ind1] <- cutoff
         x[ind2] <- -cutoff
         V0[i,] <- ceiling(ncolors * (x + cutoff)/(cutoff*2))
      }
      image(1:ncol(V0), 1:nrow(V0), t(V0), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
   }
   axis(2, at=1:nrow(V0), labels=row.names(V0), adj= 0.5, tick=FALSE, las = 1, cex.axis=character.scaling, line=0, font=2, family="")
   axis(4, at=1:nrow(V0), labels=rev(all.vals), adj= 0.5, tick=FALSE, las = 1, cex.axis=character.scaling, line=0, font=2, family="")
   axis(1, at=1:ncol(V0), labels=colnames(V0), adj= 0.5, tick=FALSE,las = 3, cex=1, cex.axis=0.4*character.scaling, line=0, font=2, family="")

        # Legend
      
#      nf <- layout(matrix(c(1, 2, 3), 3, 1, byrow=T), 1, c(1,2,2),  FALSE)
      par.mar <- par("mar")
      par(mar = c(3, 35, 8, 10))
      leg.set <- seq(-cutoff, cutoff, 2*cutoff/100)
      image(1:101, 1:1, as.matrix(leg.set), zlim=c(-cutoff, cutoff), col=mycol, axes=FALSE, main="",
          sub = "", xlab= "", ylab="",font=2, family="")
      ticks <- c(-2, -1, 0, 1, 2)
      tick.cols <- rep("black", 5)
      tick.lwd <- c(1,1,2,1,1)
      locs <- NULL
      for (k in 1:length(ticks)) locs <- c(locs, which.min(abs(ticks[k] - leg.set)))
      axis(1, at=locs, labels=ticks, adj= 0.5, tick=T, cex=0.8, cex.axis=0.8, line=0, font=2, family="")
      mtext("Standardized Target Profile", cex=0.8, side = 1, line = 3.5, outer=F)
      par(mar = par.mar)

#   print("end of heatmap")

  # Landscape plot ---------------------------------------------------------------

  if (produce.lanscape.plot == T) {
  
   nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(2, 1), FALSE)
 
#    print("cmi.names")
#    print(cmi.names)
#    print(is.matrix(cmi.names))
#    print(dim(cmi.names))
#    print(c(top.n, max.n))
#    print(cmi.names[1:top.n, 1:max.n])
#    print(dim(seed.vectors))
#    print(as.vector(cmi.names[1:top.n, 1:max.n]))
#    print(dim(as.matrix(m.2[as.vector(cmi.names[1:top.n, 1:max.n]),])))

    if (length(as.vector(cmi.names[1:top.n, 1:max.n])) > 1) {
       V0 <- rbind(seed.vectors, as.matrix(m.2[as.vector(cmi.names[1:top.n, 1:max.n]),]))
    } else {
       V0 <- rbind(seed.vectors, t(as.matrix(m.2[as.vector(cmi.names[1:top.n, 1:max.n]),])))
    }
   
   number.seq <- NULL
   for (i in 1:max.n) number.seq <- c(number.seq, rep(i, top.n))

   row.names(V0) <-  c(paste("SEED: ", seed.names, "(", signif(cmi.orig.seed, 2), ")"),
                       paste("ITER ", number.seq, ":", as.vector(cmi.names[1:top.n, 1:max.n]),
                             "(", signif(as.vector(cmi[1:top.n, 1:max.n]), 2), ")"))
 
    cmi.vals <- c(cmi.orig.seed, as.vector(cmi[1:top.n, 1:max.n]))

   total.points <- row(V0)
   V2 <- V0
   metric.matrix <- matrix(0, nrow=nrow(V2), ncol=nrow(V2))
   row.names(metric.matrix)  <- row.names(V2)
   colnames(metric.matrix) <- row.names(V2)
   MI.ref <- cmi.vals
   for (i in 1:nrow(V2)) {
      for (j in 1:i) {
           metric.matrix[i, j] <- assoc (V2[j,], V2[i,], assoc.metric)
      }
   }
   metric.matrix
   metric.matrix <- metric.matrix + t(metric.matrix)
   metric.matrix
   alpha <- 5
   metric.matrix2 <- 1 - ((1/(1+exp(-alpha*metric.matrix))))
   for (i in 1:nrow(metric.matrix2)) metric.matrix2[i, i] <- 0
   metric.matrix2
 
   smacof.map <- smacofSphere(metric.matrix2, ndim = 2, weightmat = NULL, init = NULL,
                                     ties = "primary", verbose = FALSE, modulus = 1, itmax = 1000, eps = 1e-6)
   x0 <- smacof.map$conf[,1]
   y0 <- smacof.map$conf[,2]
   r <- sqrt(x0*x0 + y0*y0)
   radius <-  1 - ((1/(1+exp(-alpha*MI.ref))))
   x <- x0*radius/r
   y <- y0*radius/r
   angles <- atan2(y0, x0)
   
   par(mar = c(4, 7, 4, 7))
 
   plot(x, y, pch=20, bty="n", xaxt='n', axes = FALSE, type="n", xlab="", ylab="", main=paste("REVEALER - Landscape for ", target.name),
        xlim=1.2*c(-max(radius), max(radius)), ylim=1.2*c(-max(radius), max(radius)))
   line.angle <- seq(0, 2*pi-0.001, 0.001)
   for (i in 1:length(x)) {
      line.max.x <- radius[i] * cos(line.angle)
      line.max.y <- radius[i] * sin(line.angle)
      points(line.max.x, line.max.y, type="l", col="gray80", lwd=1)
      points(c(0, x[i]), c(0, y[i]), type="l", col="gray80", lwd=1)
   }
   line.max.x <- 1.2*max(radius) * cos(line.angle)
   line.max.y <- 1.2*max(radius) * sin(line.angle)
   points(line.max.x, line.max.y, type="l", col="purple", lwd=2)
   points(0, 0, pch=21, bg="red", col="black", cex=2.5)   
   points(x, y, pch=21, bg="steelblue", col="black", cex=2.5)
 
#   print.names <- c(paste("TARGET:", target.name), row.names(V0))
   x <- c(0, x)
   y <- c(0, y)
 
   text(x[1], y[1], labels=print.names[1], pos=2, cex=0.85, col="red", offset=1, font=2, family="")   
   for (i in 2:length(x)) {
      pos <- ifelse(x[i] <= 0.25, 4, 2)
     text(x[i], y[i], labels=print.names[i], pos=pos, cex=0.50, col="darkblue", offset=1, font=2, family="")   
    }

  }
   dev.off()
}


assoc <- function(x, y, metric) { # Pairwise association of x and y
    if (length(unique(x)) == 1 || length(unique(y)) == 1) return(0)
    if (metric == "IC") {
       return(mutual.inf.v2(x = x, y = y, n.grid=25)$IC)
    } else if (metric == "COR") {
        return(cor(x, y))
    }
}

cond.assoc <-  function(x, y, z, metric) { # Association of a and y given z

    if (length(unique(x)) == 1 || length(unique(y)) == 1) return(0)

    if (length(unique(z)) == 1) {  # e.g. for NULLSEED
       if (metric == "IC") {
          return(mutual.inf.v2(x = x, y = y, n.grid = 25)$IC)
       } else if (metric == "COR") {
          return(cor(x, y))
       }
   } else {
       if (metric == "IC") {
          return(cond.mutual.inf(x = x, y = y, z = z, n.grid = 25)$CIC)
       } else if (metric == "COR") {
          return(pcor.test(x, y, z)$estimate)
       }
   }
}

cond.mutual.inf <- function(x, y, z, n.grid=25, delta = 0.25*c(bcv(x), bcv(y), bcv(z))) {
 # Computes the Conditional mutual imnformation: 
 # I(X, Y | X) = H(X, Z) + H(Y, Z) - H(X, Y, Z) - H(Z)
 # The 0.25 in front of the bandwidth is because different conventions between bcv and kde3d

   rho <- cor(x, y)
   rho2 <- ifelse(rho < 0, 0, rho)
   delta <- delta*(1 + (-0.75)*rho2)
  
   kde3d.xyz <- kde3d(x=x, y=y, z=z, h=delta, n = n.grid)
   X <- kde3d.xyz$x
   Y <- kde3d.xyz$y
   Z <- kde3d.xyz$z
   PXYZ <- kde3d.xyz$d + .Machine$double.eps

   # get grid spacing
   dx <- X[2] - X[1]
   dy <- Y[2] - Y[1]
   dz <- Z[2] - Z[1]

   # normalize density and calculate marginal densities and entropies
   PXYZ <- PXYZ/(sum(PXYZ)*dx*dy*dz)
   PXZ <- colSums(aperm(PXYZ, c(2,1,3)))*dy
   PYZ <- colSums(PXYZ)*dx
   PZ <- rowSums(aperm(PXYZ, c(3,1,2)))*dx*dy
   PXY <- colSums(aperm(PXYZ, c(3,1,2)))*dz
   PX <- rowSums(PXYZ)*dy*dz
   PY <- rowSums(aperm(PXYZ, c(2,1,3)))*dx*dz
   
   HXYZ <- - sum(PXYZ * log(PXYZ))*dx*dy*dz
   HXZ <- - sum(PXZ * log(PXZ))*dx*dz
   HYZ <- - sum(PYZ * log(PYZ))*dy*dz
   HZ <-  - sum(PZ * log(PZ))*dz
   HXY <- - sum(PXY * log(PXY))*dx*dy   
   HX <-  - sum(PX * log(PX))*dx
   HY <-  - sum(PY * log(PY))*dy

   MI <- HX + HY - HXY   
   CMI <- HXZ + HYZ - HXYZ - HZ

   SMI <- sign(rho) * MI
   SCMI <- sign(rho) * CMI

   IC <- sign(rho) * sqrt(1 - exp(- 2 * MI))
   CIC <- sign(rho) * sqrt(1 - exp(- 2 * CMI))
   
   return(list(CMI=CMI, MI=MI, SCMI=SCMI, SMI=SMI, HXY=HXY, HXYZ=HXYZ, IC=IC, CIC=CIC))
 }

mutual.inf.v2 <- function(x, y, n.grid=25, delta = c(bcv(x), bcv(y))) {

    # For definitions of mutual information and the universal metric (NMI) see the 
    # definition of "Mutual Information" in wikipedia and Thomas and Cover's book

   rho <- cor(x, y)
   rho2 <- abs(rho)
   delta <- delta*(1 + (-0.75)*rho2)

   kde2d.xy <- kde2d(x, y, n = n.grid, h = delta)
   FXY <- kde2d.xy$z + .Machine$double.eps
   dx <- kde2d.xy$x[2] - kde2d.xy$x[1]
   dy <- kde2d.xy$y[2] - kde2d.xy$y[1]
   PXY <- FXY/(sum(FXY)*dx*dy)
   PX <- rowSums(PXY)*dy
   PY <- colSums(PXY)*dx
   HXY <- -sum(PXY * log(PXY))*dx*dy
   HX <- -sum(PX * log(PX))*dx
   HY <- -sum(PY * log(PY))*dy

   PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
   PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)

   MI <- sum(PXY * log(PXY/(PX*PY)))*dx*dy
   rho <- cor(x, y)
   SMI <- sign(rho) * MI
   
   IC <- sign(rho) * sqrt(1 - exp(- 2 * MI)) 
   
   NMI <- sign(cor(x, y)) * ((HX + HY)/HXY - 1)  # use peason correlation the get the sign (directionality)

   return(list(MI=MI, SMI=SMI, HXY=HXY, HX=HX, HY=HY, NMI=NMI, IC=IC))
}

MSIG.Gct2Frame <- function(filename = "NULL") { 
#
# Read a gene expression dataset in GCT format and converts it into an R data frame
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

   ds <- read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T, na.strings = "")
   descs <- ds[,1]
   ds <- ds[-1]
   row.names <- row.names(ds)
   names <- names(ds)
   return(list(ds = ds, row.names = row.names, descs = descs, names = names))
}

write.gct.2 <- function(gct.data.frame, descs = "", filename) 
{
    f <- file(filename, "w")
    cat("#1.2", "\n", file = f, append = TRUE, sep = "")
    cat(dim(gct.data.frame)[1], "\t", dim(gct.data.frame)[2], "\n", file = f, append = TRUE, sep = "")
    cat("Name", "\t", file = f, append = TRUE, sep = "")
    cat("Description", file = f, append = TRUE, sep = "")

    colnames <- colnames(gct.data.frame)
    cat("\t", colnames[1], file = f, append = TRUE, sep = "")

    if (length(colnames) > 1) {
       for (j in 2:length(colnames)) {
           cat("\t", colnames[j], file = f, append = TRUE, sep = "")
       }
     }
    cat("\n", file = f, append = TRUE, sep = "\t")

    oldWarn <- options(warn = -1)
    m <- matrix(nrow = dim(gct.data.frame)[1], ncol = dim(gct.data.frame)[2] +  2)
    m[, 1] <- row.names(gct.data.frame)
    if (length(descs) > 1) {
        m[, 2] <- descs
    } else {
        m[, 2] <- row.names(gct.data.frame)
    }
    index <- 3
    for (i in 1:dim(gct.data.frame)[2]) {
        m[, index] <- gct.data.frame[, i]
        index <- index + 1
    }
    write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
    close(f)
    options(warn = 0)

}

OPAM.project.dataset.5 <- function( 
		input.ds,
		output.ds,
		gene.set.databases,
		gene.set.selection  = "ALL",  # "ALL" or list with names of gene sets
		sample.norm.type    = "rank",  # "rank", "log" or "log.rank"
		weight              = 0.25,
		statistic           = "area.under.RES",
		output.score.type   = "ES",  # "ES" or "NES"
		nperm               = 200,  # number of random permutations for NES case
		combine.mode        = "combine.off",  # "combine.off" do not combine *_UP and *_DN versions in 
		# a single score. "combine.replace" combine *_UP and 
		# *_DN versions in a single score that replaces the individual
		# *_UP and *_DN versions. "combine.add" combine *_UP and 
		# *_DN versions in a single score and add it but keeping 
		# the individual *_UP and *_DN versions.
                min.overlap         = 1,
		correl.type  = "rank")             # "rank", "z.score", "symm.rank"
{ #----------------------------------------------------------------------------------------
	
	# Load libraries
        suppressPackageStartupMessages(library(gtools))
        suppressPackageStartupMessages(library(verification))
        suppressPackageStartupMessages(library(RColorBrewer))
	
	# Read input dataset
	
	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
	m <- data.matrix(dataset$ds)
	gene.names <- dataset$row.names
	gene.descs <- dataset$descs
	sample.names <- dataset$names
	Ns <- length(m[1,])
	Ng <- length(m[,1])
	temp <- strsplit(input.ds, split="/") # Extract input file name
	s <- length(temp[[1]])
	input.file.name <- temp[[1]][s]
	temp <- strsplit(input.file.name, split=".gct")
	input.file.prefix <-  temp[[1]][1]
	
	# Sample normalization
	
	if (sample.norm.type == "rank") {
		for (j in 1:Ns) {  # column rank normalization 
			m[,j] <- rank(m[,j], ties.method = "average")
		}
		m <- 10000*m/Ng
	} else if (sample.norm.type == "log.rank") {
		for (j in 1:Ns) {  # column rank normalization 
			m[,j] <- rank(m[,j], ties.method = "average")
		}
		m <- log(10000*m/Ng + exp(1))
	} else if (sample.norm.type == "log") {
		m[m < 1] <- 1
		m <- log(m + exp(1))
	}
	
	# Read gene set databases
	
	max.G <- 0
	max.N <- 0
	for (gsdb in gene.set.databases) {
		GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
		max.G <- max(max.G, max(GSDB$size.G))
		max.N <- max.N +  GSDB$N.gs
	}
	N.gs <- 0
	gs <- matrix("null", nrow=max.N, ncol=max.G)
	gs.names <- vector(length=max.N, mode="character")
	gs.descs <- vector(length=max.N, mode="character")
	size.G <- vector(length=max.N, mode="numeric")
	start <- 1
	for (gsdb in gene.set.databases) {
		GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
		N.gs <- GSDB$N.gs 
		gs.names[start:(start + N.gs - 1)] <- GSDB$gs.names
		gs.descs[start:(start + N.gs - 1)] <- GSDB$gs.desc
		size.G[start:(start + N.gs - 1)] <- GSDB$size.G
		gs[start:(start + N.gs - 1), 1:max(GSDB$size.G)] <- GSDB$gs[1:N.gs, 1:max(GSDB$size.G)]
		start <- start + N.gs
	}
	N.gs <- max.N
	
	# Select desired gene sets
	
	if (gene.set.selection[1] != "ALL") {
		locs <- match(gene.set.selection, gs.names)
#                print(rbind(gene.set.selection, locs))
		N.gs <- sum(!is.na(locs))
		if(N.gs > 1) { 
                  gs <- gs[locs,]
		} else { 
                   gs <- t(as.matrix(gs[locs,]))   # Force vector to matrix if only one gene set specified
                }
		gs.names <- gs.names[locs]
 		gs.descs <- gs.descs[locs]
		size.G <- size.G[locs]
	}

        # Check for redundant gene sets

        tab <- as.data.frame(table(gs.names))
        ind <- order(tab[, "Freq"], decreasing=T)
        tab <- tab[ind,]
#        print(tab[1:10,])
#        print(paste("Total gene sets:", length(gs.names)))
#        print(paste("Unique gene sets:", length(unique(gs.names))))

        # Loop over gene sets
	
	score.matrix <- matrix(0, nrow=N.gs, ncol=Ns)
	for (gs.i in 1:N.gs) {
		#browser()
		gene.set <- gs[gs.i, 1:size.G[gs.i]]
		gene.overlap <- intersect(gene.set, gene.names)
#		print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", length(gene.overlap)))
                if (length(gene.overlap) < min.overlap) { 
			score.matrix[gs.i, ] <- rep(NA, Ns)
			next
		} else {
			gene.set.locs <- match(gene.overlap, gene.set)
			gene.names.locs <- match(gene.overlap, gene.names)
			msig <- m[gene.names.locs,]
			msig.names <- gene.names[gene.names.locs]
			if (output.score.type == "ES") {
				OPAM <- OPAM.Projection.3(data.array = m, gene.names = gene.names, n.cols = Ns, 
						n.rows = Ng, weight = weight, statistic = statistic,
						gene.set = gene.overlap, nperm = 1, correl.type = correl.type)
				score.matrix[gs.i,] <- OPAM$ES.vector
			} else if (output.score.type == "NES") {
				OPAM <- OPAM.Projection.3(data.array = m, gene.names = gene.names, n.cols = Ns, 
						n.rows = Ng, weight = weight, statistic = statistic,
						gene.set = gene.overlap, nperm = nperm, correl.type = correl.type)
  				score.matrix[gs.i,] <- OPAM$NES.vector
			}
		}
	}

        
        locs <- !is.na(score.matrix[,1])
#        print(paste("N.gs before overlap prunning:", N.gs))
        N.gs <- sum(locs)
#        print(paste("N.gs after overlap prunning:", N.gs))
        score.matrix <- score.matrix[locs,]
        gs.names <- gs.names[locs]
        gs.descs <- gs.descs[locs]

	initial.up.entries <- 0
	final.up.entries <- 0
	initial.dn.entries <- 0
	final.dn.entries <- 0
	combined.entries <- 0
	other.entries <- 0
	
	if (combine.mode == "combine.off") {
		score.matrix.2 <- score.matrix
		gs.names.2 <- gs.names
		gs.descs.2 <- gs.descs
	} else if ((combine.mode == "combine.replace") || (combine.mode == "combine.add")) {
		score.matrix.2 <- NULL
		gs.names.2 <- NULL
		gs.descs.2 <- NULL
		k <- 1
		for (i in 1:N.gs) {
			temp <- strsplit(gs.names[i], split="_") 
			body <- paste(temp[[1]][seq(1, length(temp[[1]]) -1)], collapse="_")
			suffix <- tail(temp[[1]], 1)
#			print(paste("i:", i, "gene set:", gs.names[i], "body:", body, "suffix:", suffix))
			if (suffix == "UP") {  # This is an "UP" gene set
				initial.up.entries <- initial.up.entries + 1
				target <- paste(body, "DN", sep="_")
				loc <- match(target, gs.names)            
				if (!is.na(loc)) { # found corresponding "DN" gene set: create combined entry
					score <- score.matrix[i,] - score.matrix[loc,]
					score.matrix.2 <- rbind(score.matrix.2, score)
					gs.names.2 <- c(gs.names.2, body)
					gs.descs.2 <- c(gs.descs.2, paste(gs.descs[i], "combined UP & DN"))
					combined.entries <- combined.entries + 1
					if (combine.mode == "combine.add") {  # also add the "UP entry
						score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
						gs.names.2 <- c(gs.names.2, gs.names[i])
						gs.descs.2 <- c(gs.descs.2, gs.descs[i])
						final.up.entries <- final.up.entries + 1
					}
				} else { # did not find corresponding "DN" gene set: create "UP" entry
					score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
					gs.names.2 <- c(gs.names.2, gs.names[i])
					gs.descs.2 <- c(gs.descs.2, gs.descs[i])
					final.up.entries <- final.up.entries + 1
				}
			} else if (suffix == "DN") { # This is a "DN" gene set
				initial.dn.entries <- initial.dn.entries + 1
				target <- paste(body, "UP", sep="_")
				loc <- match(target, gs.names)            
				if (is.na(loc)) { # did not find corresponding "UP" gene set: create "DN" entry
					score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
					gs.names.2 <- c(gs.names.2, gs.names[i])
					gs.descs.2 <- c(gs.descs.2, gs.descs[i])
					final.dn.entries <- final.dn.entries + 1
				} else { # it found corresponding "UP" gene set
					if (combine.mode == "combine.add") { # create "DN" entry
						score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
						gs.names.2 <- c(gs.names.2, gs.names[i])
						gs.descs.2 <- c(gs.descs.2, gs.descs[i])
						final.dn.entries <- final.dn.entries + 1
					}
				}
			} else { # This is neither "UP nor "DN" gene set: create individual entry
				score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
				gs.names.2 <- c(gs.names.2, gs.names[i])
				gs.descs.2 <- c(gs.descs.2, gs.descs[i])
				other.entries <- other.entries + 1
			}
		} # end for loop over gene sets
#		print(paste("initial.up.entries:", initial.up.entries))
#		print(paste("final.up.entries:", final.up.entries))
#		print(paste("initial.dn.entries:", initial.dn.entries))
#		print(paste("final.dn.entries:", final.dn.entries))
#		print(paste("other.entries:", other.entries))
#		print(paste("combined.entries:", combined.entries))
#		print(paste("total entries:", length(score.matrix.2[,1])))
	}            

        # Check for redundant gene sets

        tab <- as.data.frame(table(gs.names.2))
        ind <- order(tab[, "Freq"], decreasing=T)
        tab <- tab[ind,]
#        print(tab[1:20,])
#        print(paste("Total gene sets:", length(gs.names.2)))
#        print(paste("Unique gene sets:", length(unique(gs.names.2))))
        
	V.GCT <- data.frame(score.matrix.2)
	names(V.GCT) <- sample.names
	row.names(V.GCT) <- gs.names.2
	write.gct.2(gct.data.frame = V.GCT, descs = gs.descs.2, filename = output.ds)  
	
} # end of OPAM.project.dataset.5


OPAM.Projection.3 <- function(
		data.array,
		gene.names,
		n.cols,
		n.rows,
		weight = 0,
		statistic = "Kolmogorov-Smirnov",  # "Kolmogorov-Smirnov", # "Kolmogorov-Smirnov", "Cramer-von-Mises",
		# "Anderson-Darling", "Zhang_A", "Zhang_C", "Zhang_K",
		# "area.under.RES", or "Wilcoxon"
		gene.set,
		nperm = 200,
		correl.type  = "rank")             # "rank", "z.score", "symm.rank"
# Runs a 2-3x faster (2-2.5x for ES statistic and 2.5-3x faster for area.under.ES statsitic)
# version of GSEA.EnrichmentScore.5 internally that avoids overhead from the function call.
{
	
	ES.vector <- vector(length=n.cols)
	NES.vector <- vector(length=n.cols)
	p.val.vector <- vector(length=n.cols)
	correl.vector <- vector(length=n.rows, mode="numeric")
	
# Compute ES score for signatures in each sample
	
#   print("Computing GSEA.....")
	phi <- array(0, c(n.cols, nperm))
	for (sample.index in 1:n.cols) {
		gene.list <- order(data.array[, sample.index], decreasing=T)            
		
		#      print(paste("Computing observed enrichment for UP signature in sample:", sample.index, sep=" ")) 

		gene.set2 <- match(gene.set, gene.names)
		
		if (weight == 0) {
			correl.vector <- rep(1, n.rows)
		} else if (weight > 0) {
			if (correl.type == "rank") {
				correl.vector <- data.array[gene.list, sample.index]
			} else if (correl.type == "symm.rank") {
				correl.vector <- data.array[gene.list, sample.index]
				correl.vector <- ifelse(correl.vector > correl.vector[ceiling(n.rows/2)], 
						correl.vector,
						correl.vector + correl.vector - correl.vector[ceiling(n.rows/2)]) 
			} else if (correl.type == "z.score") {
				x <- data.array[gene.list, sample.index]
				correl.vector <- (x - mean(x))/sd(x)
			}
		}
		### Olga's Additions ###
#		ptm.new = proc.time()
		tag.indicator <- sign(match(gene.list, gene.set2, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
		no.tag.indicator <- 1 - tag.indicator 
		N <- length(gene.list) 
		Nh <- length(gene.set2) 
		Nm <-  N - Nh 
		orig.correl.vector <- correl.vector
		if (weight == 0) correl.vector <- rep(1, N)   # unweighted case
		ind = which(tag.indicator==1)
		correl.vector <- abs(correl.vector[ind])^weight
		
		
		sum.correl = sum(correl.vector)
		up = correl.vector/sum.correl     # "up" represents the peaks in the mountain plot
		gaps = (c(ind-1, N) - c(0, ind))  # gaps between ranked pathway genes
		down = gaps/Nm
		
		RES = cumsum(c(up,up[Nh])-down)
		valleys = RES[1:Nh]-up
		
		max.ES = max(RES)
		min.ES = min(valleys)
		
		if( statistic == "Kolmogorov-Smirnov" ){
			if( max.ES > -min.ES ){
				ES <- signif(max.ES, digits=5)
				arg.ES <- which.max(RES)
			} else{
				ES <- signif(min.ES, digits=5)
				arg.ES <- which.min(RES)
			}
		}
		
		if( statistic == "area.under.RES"){
			if( max.ES > -min.ES ){
				arg.ES <- which.max(RES)
			} else{
				arg.ES <- which.min(RES)
			}
			gaps = gaps+1
			RES = c(valleys,0) * (gaps) + 0.5*( c(0,RES[1:Nh]) - c(valleys,0) ) * (gaps)
			ES = sum(RES)
		}
		GSEA.results = list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator)
#		new.time <<- new.time + (proc.time() - ptm.new)
		### End Olga's Additions ###
		#GSEA.results <- GSEA.EnrichmentScore5(gene.list=gene.list, gene.set=gene.set2,
		#		statistic = statistic, alpha = weight, correl.vector = correl.vector)
		ES.vector[sample.index] <- GSEA.results$ES
		
		if (nperm == 0) {
			NES.vector[sample.index] <- ES.vector[sample.index]
			p.val.vector[sample.index] <- 1
		} else {
			for (r in 1:nperm) {
				reshuffled.gene.labels <- sample(1:n.rows)
				if (weight == 0) {
					correl.vector <- rep(1, n.rows)
				} else if (weight > 0) {
					correl.vector <- data.array[reshuffled.gene.labels, sample.index]
				} 
#				GSEA.results <- GSEA.EnrichmentScore5(gene.list=reshuffled.gene.labels, gene.set=gene.set2,
#						statistic = statistic, alpha = weight, correl.vector = correl.vector)
				### Olga's Additions ###
				tag.indicator <- sign(match(reshuffled.gene.labels, gene.set2, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
				no.tag.indicator <- 1 - tag.indicator 
				N <- length(reshuffled.gene.labels) 
				Nh <- length(gene.set2) 
				Nm <-  N - Nh 
#   orig.correl.vector <- correl.vector
				if (weight == 0) correl.vector <- rep(1, N)   # unweighted case
				ind <- which(tag.indicator==1)
				correl.vector <- abs(correl.vector[ind])^weight   
				
				sum.correl <- sum(correl.vector)
				up = correl.vector/sum.correl
				gaps = (c(ind-1, N) - c(0, ind))
				down = gaps/Nm
				
				RES = cumsum(c(up,up[Nh])-down)
				valleys = RES[1:Nh]-up
				
				max.ES = max(RES)
				min.ES = min(valleys)
				
				if( statistic == "Kolmogorov-Smirnov" ){
					if( max.ES > -min.ES ){
						ES <- signif(max.ES, digits=5)
						arg.ES <- which.max(RES)
					} else{
						ES <- signif(min.ES, digits=5)
						arg.ES <- which.min(RES)
					}
				}
				
				if( statistic == "area.under.RES"){
					if( max.ES > -min.ES ){
						arg.ES <- which.max(RES)
					} else{
						arg.ES <- which.min(RES)
					}
					gaps = gaps+1
					RES = c(valleys,0) * (gaps) + 0.5*( c(0,RES[1:Nh]) - c(valleys,0) ) * (gaps)
					ES = sum(RES)
				}
				
				GSEA.results = list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator)
				### End Olga's Additions ###
				phi[sample.index, r] <- GSEA.results$ES
			}
			if (ES.vector[sample.index] >= 0) {
				pos.phi <- phi[sample.index, phi[sample.index, ] >= 0]
				if (length(pos.phi) == 0) pos.phi <- 0.5
				pos.m <- mean(pos.phi)
				NES.vector[sample.index] <- ES.vector[sample.index]/pos.m
				s <- sum(pos.phi >= ES.vector[sample.index])/length(pos.phi)
				p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
			} else {
				neg.phi <-  phi[sample.index, phi[sample.index, ] < 0]
				if (length(neg.phi) == 0) neg.phi <- 0.5 
				neg.m <- mean(neg.phi)
				NES.vector[sample.index] <- ES.vector[sample.index]/abs(neg.m)
				s <- sum(neg.phi <= ES.vector[sample.index])/length(neg.phi)
				p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
			}
		}
	}
	return(list(ES.vector = ES.vector, NES.vector =  NES.vector, p.val.vector = p.val.vector))
	
} # end of OPAM.Projection.3

Read.GeneSets.db <- function(
		gs.db,
		thres.min = 2,
		thres.max = 2000,
		gene.names = NULL)
{
	
	temp <- readLines(gs.db)
	max.Ng <- length(temp)
	temp.size.G <- vector(length = max.Ng, mode = "numeric") 
	for (i in 1:max.Ng) {
		temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
	}
	max.size.G <- max(temp.size.G)      
	gs <- matrix(rep("null", max.Ng*max.size.G), nrow=max.Ng, ncol= max.size.G)
	temp.names <- vector(length = max.Ng, mode = "character")
	temp.desc <- vector(length = max.Ng, mode = "character")
	gs.count <- 1
	for (i in 1:max.Ng) {
		gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
		gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
		gene.set.name <- gs.line[1] 
		gene.set.desc <- gs.line[2] 
		gene.set.tags <- vector(length = gene.set.size, mode = "character")
		for (j in 1:gene.set.size) {
			gene.set.tags[j] <- gs.line[j + 2]
		}
		if (is.null(gene.names)) {
			existing.set <- rep(TRUE, length(gene.set.tags))
		} else {
			existing.set <- is.element(gene.set.tags, gene.names)
		}
		set.size <- length(existing.set[existing.set == T])
		if ((set.size < thres.min) || (set.size > thres.max)) next
		temp.size.G[gs.count] <- set.size
		gs[gs.count,] <- c(gene.set.tags[existing.set], rep("null", max.size.G - temp.size.G[gs.count]))
		temp.names[gs.count] <- gene.set.name
		temp.desc[gs.count] <- gene.set.desc
		gs.count <- gs.count + 1
	}
	Ng <- gs.count - 1
	gs.names <- vector(length = Ng, mode = "character")
	gs.desc <- vector(length = Ng, mode = "character")
	size.G <- vector(length = Ng, mode = "numeric") 
	
	gs.names <- temp.names[1:Ng]
	gs.desc <- temp.desc[1:Ng]
	size.G <- temp.size.G[1:Ng]
	
	return(list(N.gs = Ng, gs = gs, gs.names = gs.names, gs.desc = gs.desc, size.G = size.G, max.N.gs = max.Ng))
}


SCREENER.v1 <- function(
       ds1,                                # Input dataset with target (CLS or GCT) file
       target.name,                        # Name of target in ds1
       target.combination.op = "max",      # operation to reduce multiple targets to one vector of values
       cond.feature.name = NULL,           # Feature in ds1 to be used as conditional variable for CMI. "TISSUE" uses the tissue type
       ds2,                                # Input feature dataset 
       n.markers = 20,                     # Number of markers for bootstrap, heatmap and mds plot
       n.perm = 3,                         # Number of random permutations
       permutation.test.type = "standard", # balanced  subclass.stratified
       n.boot = 5,                         # Number of bootstrap samples for confidence interval of the n.markers
       seed = 86876,                       # Random number generator seed
       assoc.metric.type = "RNMI",         # Association metric: RNMI, NMI, SMI, AUC.ROC, AUC.REV, DIFF.MEDIANS, DIFF.MEANS, S2N, T.TEST, CORR
       direction = "positive",             # Direction of feature matching
       sort.target = TRUE,                 # Sort columns according to phenotype (heatmap)
       results.file.pdf,                   # PDF output file
       results.file.txt,                   # TXT output file
       results.file.gct = NULL,            # GCT file with top results                 
       sort.columns.inside.classes = T,    # Sort columns in heatmap of top matches inside each target class               
       cluster.top.markers = F,            # Sort rows in hetamap of top matches (F, "each.class" or "both.classes")                   
       consolidate.identical.features = F, # Consolidate identical features: F or "identical" or "similar" 
       cons.features.hamming.thres = 3,    # If consolidate.identical.features = "similar" then consolidate features within this Hamming dist. thres.
       locs.table.file = NULL,             # Table with chromosonal locations per feature (gene)
       save.matched.dataset = F,           # Save target-fetaures matched dataset                                          
       produce.aux.histograms = F,         # Produce histograms of assoc. metric distribution etc.                         
       produce.heat.map = T,               # Produce heatmap                                                               
       produce.mds.plots = T,              # Produce multi-dimensional scaling (mds) plot (Landscape plot)                 
       character.scaling = 1,              # character scaling for heatmap
       mds.plot.type = "smacof",           # mds algorithm                                                                         
       knn = 3,                            # k for knn-based assoc. metrics
       display.cor.coeff.heatmap = F,      # show correlation coeff. in heatmap/text file
       n.grid=25,                          # grid size for kernel-based assoc. metrics                                     
       phen.table = NULL,                  # Table with phenotypes for each sample (optional)
       phen.column = NULL,                 # Column in phen.table containing the relevant phenotype info
       phen.selected = NULL,               # Use only samples of these phenotypes in analysis
       debug.mode = F)                     # Print additional info to diagnose problems

  {

   suppressPackageStartupMessages(library(MASS))
   suppressPackageStartupMessages(library(PerformanceAnalytics))
   suppressPackageStartupMessages(library(parmigene))
   suppressPackageStartupMessages(library(maptools))
   suppressPackageStartupMessages(library(smacof))

   set.seed(seed)

   time1 <- proc.time()

   pdf(file=results.file.pdf, height=11, width=8.5)
   if (is.null(results.file.gct)) {
      results.file.gct = paste(results.file.pdf, ".gct", sep="")
   }
     
   # Read table with HUGO gene symbol vs. chr location
   
   if (!is.null(locs.table.file)) {
      locs.table <- read.table(locs.table.file, header=T, sep="\t", skip=0, colClasses = "character")
    }

   # Read input files 

#   print("Reading features file...")
   
   if (regexpr(pattern=".gct", ds2) != -1) {
      dataset2 <- MSIG.Gct2Frame(filename = ds2)
      m2 <- data.matrix(dataset2$ds)
      feature.names2 <- dataset2$row.names
      sample.names2 <- colnames(m2)
   } else if (regexpr(pattern=".txt", ds2) != -1) {
      df1 <- read.table(ds2, header=T, row.names=1, sep="\t", skip=0)
      df1 <- data.matrix(df1)
      m2 <- t(df1)
      feature.names2 <- row.names(m2)
      sample.names2 <- colnames(m2)
    }
   m2[is.na(m2)] <- 0

   # Filter samples with only the selected phenotypes 

   if (!is.null(phen.selected)) {
#      print(paste("Subselecting samples with phenotype: ", phen.selected))
#      samples.table <- read.table(phen.table, header=T, row.names=1, sep="\t", skip=0)
      samples.table <- read.delim(phen.table, header=T, row.names=1, sep="\t", skip=0)
      table.sample.names <- row.names(samples.table)
      locs1 <- match(colnames(m2), table.sample.names)
      phenotype <- as.character(samples.table[locs1, phen.column])
     # print(paste("phenotype:", phenotype))
      table(phenotype)
      locs2 <- NULL
      for (k in 1:ncol(m2)) {   
         if (!is.na(match(phenotype[k], phen.selected))) {
            locs2 <- c(locs2, k)
         }
      }
#      print(paste("Matching phenotype total number of samples:", length(locs2)))
      m2 <- m2[, locs2]
      sample.names2 <- colnames(m2)
      phenotype <- phenotype[locs2]
#      print(table(phenotype))
#      print(dim(m2))
    }
   
#   print("Reading target file...")
   
   if (regexpr(pattern=".gct", ds1) != -1) {
      dataset1 <- MSIG.Gct2Frame(filename = ds1)
      m1 <- data.matrix(dataset1$ds)
#      m1[is.na(m1)] <- 0          
      feature.names1 <- dataset1$row.names
      sample.names1 <- colnames(m1)

      if (length(target.name) > 1) {  # multiple targets => combine them using target.combination.op 
#         print("multi target")
         target <- apply(m1[target.name,], MARGIN=2, FUN=target.combination.op)
         target.name <- paste(target.name, collapse="__")
#         print(target.name)
      } else { # single target
#         print("single target")        
         target <- m1[target.name,]
#         print(target.name)
      }

     # exclude samples with target == NA

#      print(paste("initial target length:", length(target)))      
      locs <- seq(1, length(target))[!is.na(target)]
      m1 <- m1[,locs]
      target <- target[locs]
      sample.names1 <- sample.names1[locs]

#      print(paste("target length after excluding NAs:", length(target)))      

      overlap <- intersect(sample.names1, sample.names2)
#      print(paste("Size of overlap:", length(overlap)))
      locs1 <- match(overlap, sample.names1)
      locs2 <- match(overlap, sample.names2)
      m1 <- m1[, locs1]
      m2 <- m2[, locs2]
      target <- target[locs1]
      
#      print(paste("final target length:", length(target)))
      
#      print(dim(m1))
#      print(dim(m2))

      if (!is.null(cond.feature.name)) {
         if (cond.feature.name == "TISSUE") {
            tissue.type <- vector(length=ncol(m1), mode="character")
            for (k in 1:ncol(m1)) {
               temp <- strsplit(colnames(m3)[k], split="_") 
               tissue.type[k] <- paste(temp[[1]][2:length(temp[[1]])], collapse="_")
             }
            cond.feature <- match(tissue.type, unique(tissue.type))
         } else {
            cond.feature <- m1[cond.feature.name,]
         }
       }
      
      if (save.matched.dataset == T) {
         write.gct.2(gct.data.frame = m1, descs = row.names(m1), filename = paste(ds1, ".MATCHED.SET.gct", sep=""))
         write.gct.2(gct.data.frame = m2, descs = row.names(m2), filename = paste(ds2, ".MATCHED.SET.gct", sep=""))
#         print(paste(ds1, ".MATCHED.SET.gct", sep=""))
#         print(paste(ds2, ".MATCHED.SET.gct", sep=""))
      }

   } else if (regexpr(pattern=".txt", ds1) != -1) {
      df1 <- read.table(ds1, header=T, row.names=1, sep="\t", skip=0)
#      df1 <- data.matrix(df1)
      m1 <- t(df1)
      m1[is.na(m1)] <- 0
      sample.names1 <- colnames(m1)
      overlap <- intersect(sample.names1, sample.names2)
      locs1 <- match(overlap, sample.names1)
      locs2 <- match(overlap, sample.names2)
      m1 <- m1[, locs1]
      m2 <- m2[, locs2]
      if (length(target.name) > 1) {  # multiple targets => combine them using target.combination.op 
         target <- apply(m1[target.name,], MARGIN=2, FUN=target.combination.op)
         target.name <- paste(target.name, collapse="__")
      } else { # single target
         target <- m1[target.name,]
      }
      classes <- unique(target)
      target <- match(target, classes)
#      print(dim(m1))
#      print(dim(m2))
                        
   } else if (regexpr(pattern=".cls", ds1) != -1) {
      CLS <- MSIG.ReadClsFile(ds1)
      if (is.numeric(CLS$class.list)) {
         target <- as.numeric(CLS$class.list)
       } else {
         target <- match(CLS$class.list, rev(unique(CLS$class.list)))
       }
   }
  target <- as.numeric(target)
   
  if (debug.mode == T) {
#     print("target:")
#     print(target)
#     if (length(table(target)) > 10) print(table(target))
#     print(table(target))
   }
   
  # Filter out features with only 0's

   row.sums <- rowSums(m2)
   locs <- seq(1, nrow(m2))[row.sums == 0]
   if (length(locs) > 1)  m2 <- m2[-locs,]
#   print(dim(m2))

   #  Consolidate identical features

   if (consolidate.identical.features == "identical") {  # This is a very fast way to eliminate perfectly identical features compared
                                                         # with what we do below in "similar"
#      print(paste("Consolidating features..."))
      summary.vectors <- apply(m2, MARGIN=1, FUN=paste, collapse="")
      ind <- order(summary.vectors)
      summary.vectors <- summary.vectors[ind]
      m2 <- m2[ind,]
      taken <- i.count <- rep(0, length(summary.vectors))
      i <- 1
      while (i <= length(summary.vectors)) {
         j <- i + 1
         while ((summary.vectors[i] == summary.vectors[j]) & (j <= length(summary.vectors))) {
           j <- j + 1
         }
        i.count[i] <- j - i
        if (i.count[i] > 1) taken[seq(i + 1, j - 1)] <- 1
        i <- j
      }
      if (sum(i.count) != length(summary.vectors)) stop("ERROR")    # add counts in parenthesis
      i.count <- ifelse(i.count > 1, paste("(", i.count, ")", sep=""), rep(" ", length(i.count)))
      row.names(m2) <- paste(row.names(m2), i.count)
      m2 <- m2[taken == 0,]

   } else if (consolidate.identical.features == "similar") { # this uses the hamming distance to consolidate similar features up to the Hamming dist. threshold 
#   print(paste("Consolidating features..."))
   hamming.matrix <- hamming.distance(m2)
   taken <- rep(0, nrow(m2))
   for (i in 1:nrow(m2)) {
    if (taken[i] == 0) { 
       similar.features <- row.names(m2)[hamming.matrix[i,] <= cons.features.hamming.thres]
       if (length(similar.features) > 1) {
           row.names(m2)[i]  <- paste(row.names(m2)[i], " (", length(similar.features), ")", sep="")  # add counts in brackets
           locs <- match(similar.features, row.names(m2))
           taken[locs] <- 1
           taken[i] <- 0
        }
      }
   }
  m2 <- m2[taken == 0,]
 }
#   print(dim(m2))
   
   # Add location info

   if (!is.null(locs.table.file)) {
#      print(paste("Adding location info..."))
      gene.symbol <- row.names(m2)
      chr <- rep(" ", length(gene.symbol))
      for (i in 1:length(gene.symbol)) {
        temp1 <- strsplit(gene.symbol[i], split="_")
        temp2 <- strsplit(temp1[[1]][1], split="\\.")
        gene.symbol[i] <- ifelse(temp2[[1]][1] == "", temp1[[1]][1], temp2[[1]][1])
        loc <- match(gene.symbol[i], locs.table[,"Approved.Symbol"])
        chr[i] <- ifelse(!is.na(loc), paste("(", locs.table[loc, "Chromosome"], ") ", sep=""), " ")
       }
      row.names(m2)  <- paste(row.names(m2), chr)
#      print(paste("Total unmatched to chromosomal locations:", sum(chr == " "), "out of ", nrow(m2), "features"))
    }
   
   # Add small amount of noise to remaining features

   noise.m <- matrix(rnorm(nrow(m2)*ncol(m2)), nrow=nrow(m2), ncol(m2))
   m2 <- m2 +  10 * .Machine$double.eps * noise.m 

#   for (i in 1:nrow(m2)) {
#      if (i < 10) print(paste(" before:", m2[i,]))
#      m2[i,] <- m2[i,] + 10 * .Machine$double.eps * rnorm(ncol(m2))
#      if (i < 10) print(paste(" after:", m2[i,]))
#    }
   
   if (sort.target == TRUE) {
      ind <- order(target, decreasing=T)
      target <-  target[ind]
      if (!is.null(cond.feature.name)) {
         cond.feature <-  cond.feature[ind]
      }
      m3 <- m2[, ind]
   } else {
      m3 <- m2
    }

   N <- ncol(m3)
   p <- nrow(m3)
   if (direction == "negative") {
      if (length(table(target)) > N*0.5) { # continuous target
         target2 <- -target
      } else {
         target2 <-  1 - target
      }
   } else if (direction == "positive") {
      target2 <- target
   } else {
      stop(paste("Unknown direction:", direction))
   }

   target2 <- target2 + 10 * .Machine$double.eps * rnorm(length(target2))
#   print(paste("Length of target:", length(target2)))
      
   # Find association of target vs. features using selected metric

#   print("Finding association of target vs. features using selected metric...")

   metric <- cor.val <- vector(mode="numeric", length=p)
   metric.rand <- matrix(0, nrow=p, ncol=n.perm)

   target.rand <- matrix(target2, nrow=n.perm, ncol=N, byrow=TRUE)
   if (permutation.test.type == "standard") {    # balanced  subclass.stratified
      for (i in 1:n.perm) target.rand[i,] <- sample(target.rand[i,])
   } else if (permutation.test.type == "balanced") {    # balanced  subclass.stratified

   # we removed this option because of too much noise
     
   } else if (permutation.test.type == "subclass.stratified") {
      subclass.type <- vector(length=ncol(m3), mode="character")
      for (k in 1:ncol(m3)) {
         temp <- strsplit(colnames(m3)[k], split="_") 
         subclass.type[k] <- paste(temp[[1]][2:length(temp[[1]])], collapse="_")
      }
      all.types <- unique(subclass.type)
      for (k in 1: length(all.types)) {
         V2 <- as.matrix(target.rand[, all.types[k] == subclass.type])
#         print(dim(V2))
         if (ncol(V2) > 1) for (i in 1:n.perm) V2[i,] <- sample(V2[i,])
         target.rand[, all.types[k] == subclass.type] <- V2
       }
    } else {
      stop(paste("Unknown permutation test type:", permutation.test.type))
    }

   for (i in 1:p) {
#      if (i %% 100 == 0) print(paste(" feature #", i, " out of ", p))
      feature <- as.numeric(m3[i,])
      metric[i] <- assoc.metric(target = target2, feature = feature, cond.feature = cond.feature, type=assoc.metric.type, knn=knn, n.grid=n.grid)
      cor.val[i] <- cor(target2, feature)
      for (k in 1:n.perm) {
         metric.rand[i, k] <- assoc.metric(target = target.rand[k,], feature = feature, cond.feature = cond.feature, type=assoc.metric.type,
                                           knn=knn, n.grid=n.grid)
       }
    }

   # Sort features

   ind <- order(metric, decreasing=T)
   metric <- metric[ind]
   cor.val <- cor.val[ind]
   metric.rand <- metric.rand[ind,]
   metric.rand.max <- apply(metric.rand, MARGIN=1, FUN=max)
   m3 <- m3[ind,]
   feature.names2 <- row.names(m3) 
 
   n.markers <- ifelse(nrow(m3) < 2*n.markers, floor(nrow(m3)/2), n.markers)
   
   # Compute Confidence Intervals for top features using 0.632 bootstrap

   if (!is.null(n.boot)) {
#      print("computing bootstrap confidence intervals...")
      b.size <- ceiling(0.632*N)
      metric.CI <- matrix(0, nrow=p, ncol=3)
      boots.null <- matrix(0, nrow=2*n.markers, ncol=n.boot)
      for (k in 1:n.boot) {
         locs <- sample(seq(1, N), b.size, replace=T)
         m3.sample <- m3[c(seq(1, n.markers, 1), seq(p, p - n.markers + 1, -1)), locs]
         target2.sample <- target2[locs]
         for (i in 1:(2*n.markers)) {
            feature.sample <- as.numeric(m3.sample[i,])
            boots.null[i, k] <- assoc.metric(target = target2.sample, feature = feature.sample, cond.feature = cond.feature,
                                          type=assoc.metric.type, knn=knn, n.grid=n.grid)
          }
       }
       
       for (i in 1:n.markers) metric.CI[i, ] <- quantile(boots.null[i,], probs = c(0.05, 0.5, 0.95), na.rm = T)
       j <- n.markers + 1
       for (i in seq(p, p - n.markers + 1, -1)) {
          metric.CI[i, ] <- quantile(boots.null[j,], probs = c(0.05, 0.5, 0.95))
          j <- j+1
       }
    } else {
      metric.CI <- matrix(0, nrow=p, ncol=3)
    }
  
   # Make histogram and QQplot

   if (produce.aux.histograms == T) {
      nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), c(1, 1), 1, FALSE)
      h <- hist(metric.rand, breaks=40, col="steelblue", main="Global Null Dist", xlab=paste("Metric:", assoc.metric.type),
                ylab = "P(M)", xlim=range(metric))
      for (i in 1:p) lines(c(metric[i], metric[i]), c(0, -0.025*max(h$counts)), lwd=1, col="black")
      chart.QQPlot(metric, main = paste("QQ Plot for ", assoc.metric.type, " Association Metric"), distribution = 'norm',
                   envelope=0.95, pch=20, cex=0.6, lwd=3)
    }

   # Compute p-values and FDR

#   print("computing p-vals and FDRs...")
   
   p.val1 <- p.val2 <- FDR1 <- FDR1.lower <- FDR2 <- FDR2.lower <- FDR3 <- FDR3.lower <- FWER.p.val <- Bonfe.p.val <- rep(0, p)
   p.val1.tag <- p.val2.tag <- p.val3.tag <- FDR1.tag <- FDR2.tag <- FDR3.tag <- FWER.p.val.tag <- Bonfe.p.val.tag <- rep(0, p)

   for (i in 1:p) {
      p.val1[i] <- sum(metric.rand[i,] > metric[i])/n.perm
      p.val2[i] <- sum(metric.rand > metric[i])/(p*n.perm)
      FWER.p.val[i] <- sum(metric.rand.max > metric[i])/p
   }
   if (assoc.metric.type == "DIFF.MEANS") { # use local p-vals: p.val1
     FDR1 <- p.adjust(p.val1, method = "fdr", n = length(p.val1))
     FDR1.lower <- p.adjust(1 - p.val1, method = "fdr", n = length(p.val1))
     FDR2 <- qvalue(p.val1)$qvalues
     FDR2.lower <- qvalue(1 - p.val2)$qvalues
   } else { # Use p.val2 for the other metrics (have global null distributions)
#     print("computing FDRs using p.adjust...")
     FDR1 <- p.adjust(p.val2, method = "fdr", n = length(p.val2))
     FDR1.lower <- p.adjust(1 - p.val2, method = "fdr", n = length(p.val2))
#     print("computing FDRs using q.values...")     
#     FDR2 <- qvalue(p.val2)$qvalues
#     FDR2.lower <- qvalue(1 - p.val2)$qvalues
     FDR2 <- FDR1    # qvalue crashes sometimes: use FDR1 
     FDR2.lower <- FDR1.lower
   }

   for (i in 1:p) {
      FDR3[i] <- (sum(metric.rand >= metric[i])/(p*n.perm))/(sum(metric >= metric[i])/p)
      FDR3.lower[i] <- (sum(metric.rand <= metric[i])/(p*n.perm))/(sum(metric <= metric[i])/p)
    }
   for (i in 1:p) {
     FDR3[i] <- min(FDR3[i:p]) 
     FDR3.lower[i] <- min(FDR3[1:i])
   }

   lower.side <- rep(0, p)
   for (i in 1:p) {
#      if ((assoc.metric.type == "AUC.ROC") | (assoc.metric.type == "CORR")) {
   if (assoc.metric.type == "AUC.ROC") {
        if (metric[i] < 0.5) lower.side[i] <- 1
      } else { # RNMI, NMI, DIFF.MEDIANS, DIFF.MEANS, S2N, T.TEST
        if (metric[i] < 0) lower.side[i] <- 1
      }
    }
   
   for (i in 1:p) {
      if (lower.side[i] == 1) {
         p.val1[i] <- 1 - p.val1[i]
         p.val2[i] <- 1 - p.val2[i]
         FWER.p.val[i] <- 1 - FWER.p.val[i]
         Bonfe.p.val[i] <- 1 - Bonfe.p.val[i]
         FDR1[i] <- FDR1.lower[i]
         FDR2[i] <- FDR2.lower[i]
         FDR3[i] <- FDR3.lower[i]
       }
      if (p.val1[i] == 0) {
          p.val1[i] <- 1/n.perm
          p.val1.tag[i] <- "<"
      }
      if (p.val2[i] == 0) {
          p.val2[i] <- 1/(p*n.perm)
          p.val2.tag[i] <- "<"
      }
      if (FDR3[i] == 0) {
          FDR3[i] <- (1/(p*n.perm))/(i/p)
          FDR3.tag[i] <- "<"
      }
      if (FWER.p.val[i] == 0) {
          FWER.p.val[i] <- 1/p
          FWER.p.val.tag[i] <- "<"
      }
    }
   for (i in 1:p) Bonfe.p.val[i] <- ifelse(p.val2[i] * p > 1, 1, p.val2[i] * p)

   # Make histograms of p-values and FDRs

   if (produce.aux.histograms == T) {
     
      nf <- layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), 3, 3, byrow=T), c(1, 1, 1), c(1, 1, 1), FALSE)
      h <- hist(p.val1, breaks=40, col="darkblue", main="p-val1", xlab="p-value", ylab = "P(p-val1)")
      h <- hist(p.val2, breaks=40, col="darkblue", main="p-val2", xlab="p-value", ylab = "P(p-va2l)")
      h <- hist(FDR1, breaks=40, col="darkblue", main="FDR1", xlab="FDR (HB)", ylab = "P(FDR1)")
      h <- hist(FDR2, breaks=40, col="darkblue", main="FDR2", xlab="FDR (Storey)", ylab = "P(FDR2)")
      h <- hist(FDR3, breaks=40, col="darkblue", main="FDR3", xlab="FDR (Theor.)", ylab = "P(FDR3)")
      h <- hist(FWER.p.val, breaks=40, col="darkblue", main="FWER p-values", xlab="p-value", ylab = "P(FWER p-val)")
      h <- hist(Bonfe.p.val, breaks=40, col="darkblue", main="Bonferroni p-values", xlab="p-value", ylab = "P(Bonfe p-val)")
    }

   metric.CI <- signif(metric.CI, 4)
   metric.CI[metric.CI == 0] <- "-"                      
   if (display.cor.coeff.heatmap == T) {
      report1 <- cbind(seq(1, p), feature.names2, signif(metric, 3), metric.CI,
                    paste(p.val1.tag, signif(p.val1, 3), sep=""),
                    paste(p.val2.tag, signif(p.val2, 3), sep=""),
                    paste(FDR1.tag, signif(FDR1, 3), sep=""),
                    paste(FDR2.tag, signif(FDR2, 3), sep=""),
                    paste(FDR3.tag, signif(FDR3, 3), sep=""),
                    paste(FWER.p.val.tag, signif(FWER.p.val, 3), sep=""),
                    paste(Bonfe.p.val.tag, signif(Bonfe.p.val, 3), sep=""),
                    signif(cor.val), 3)
       colnames(report1) <- c("Rank", "Feature", assoc.metric.type, "5% CI", "50% CI", "95% CI", "p-val1 (local)", "p-val2 (global)", "FDR1",
                          "FDR2", "FDR3", "FWER p-val", "Bonfe p-val", "Corr. Coef.")
     } else {
      report1 <- cbind(seq(1, p), feature.names2, signif(metric, 3), metric.CI,
                    paste(p.val1.tag, signif(p.val1, 3), sep=""),
                    paste(p.val2.tag, signif(p.val2, 3), sep=""),
                    paste(FDR1.tag, signif(FDR1, 3), sep=""),
                    paste(FDR2.tag, signif(FDR2, 3), sep=""),
                    paste(FDR3.tag, signif(FDR3, 3), sep=""),
                    paste(FWER.p.val.tag, signif(FWER.p.val, 3), sep=""),
                    paste(Bonfe.p.val.tag, signif(Bonfe.p.val, 3), sep=""))
       colnames(report1) <- c("Rank", "Feature", assoc.metric.type, "5% CI", "50% CI", "95% CI", "p-val1 (local)", "p-val2 (global)", "FDR1",
                          "FDR2", "FDR3", "FWER p-val", "Bonfe p-val")
     }
     
#   print(noquote(report1[1:n.markers,]))
#   print(noquote(report1[seq(p, p - n.markers + 1, -1),]))

   write.table(report1, file=results.file.txt, quote=F, col.names = T, row.names = F, append = F, sep="\t")

   metric.sorted <- metric[c(1:n.markers, seq(p - n.markers + 1, p))]
   cor.val.sorted <- cor.val[c(1:n.markers, seq(p - n.markers + 1, p))]
   p.val1.sorted <- p.val1[c(1:n.markers, seq(p - n.markers + 1, p))]
   p.val2.sorted <- p.val2[c(1:n.markers, seq(p - n.markers + 1, p))]
   FDR1.sorted <- FDR1[c(1:n.markers, seq(p - n.markers + 1, p))]

   # Make heatmap of top and bottom features

   if (produce.heat.map == T) {
#      print("making heatmap...")
     
      mycol <- vector(length=512, mode = "numeric")
      for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
      for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
      mycol <- rev(mycol)

      ncolors <- length(mycol)

      suppressPackageStartupMessages(library(RColorBrewer))
      nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(1, 11), FALSE)

      m.V <- m3[c(1:n.markers, seq(p - n.markers + 1, p)),]
      m.V[is.na(m.V)] <- 0
      target.V <- target
      
      if (sort.columns.inside.classes == T & length(table(target.V)) < ncol(m.V)*0.5)  {   # sort columns inside classes if target is not continuous
         num.phen <- length(unique(target.V))
         for (k in unique(target.V)) {
            V3 <- m.V[, target.V == k]
            cl <- target.V[target.V == k]
            c.names <- colnames(m.V)[target.V == k]

#            dist.matrix <- dist(t(V3))
#            HC <- hclust(dist.matrix, method="complete")
#            ind <- HC$order

            dist.matrix <- dist(t(V3))  
            s <- smacofSym(dist.matrix, ndim=1)
            ind <- order(s$conf, decreasing=T)
            
            V3 <- V3[, ind]
            cl <- cl[ind]
            c.names <- c.names[ind]
            target.V[target.V == k] <- cl
            m.V[, target.V == k] <- V3
            colnames(m.V)[target.V == k] <- c.names
         }
      }

      if (cluster.top.markers == "each.class") {
         marker.ind <- c(rep(1, n.markers), rep(2, n.markers))
         num.phen <- length(unique(marker.ind))
         for (k in unique(marker.ind)) {
            V3 <- m.V[marker.ind == k,]
            row.names.V3 <- row.names(m.V)[marker.ind == k]
            metric.sorted.V <- metric.sorted[marker.ind == k]
            cor.val.sorted.V <-cor.val.sorted[marker.ind == k]            
            p.val1.sorted.V <- p.val1.sorted[marker.ind == k]
            p.val2.sorted.V <- p.val2.sorted[marker.ind == k]
            FDR1.sorted.V <- FDR1.sorted[marker.ind == k]

#            dist.matrix <- dist(V3)
#            HC <- hclust(dist.matrix, method="complete")
#            ind <- HC$order
            
            dist.matrix <- dist(V3)  
            s <- smacofSym(dist.matrix, ndim=1)
            ind <- order(s$conf, decreasing=T)

            V3 <- V3[ind,]
            row.names.V3 <- row.names.V3[ind]
            metric.sorted.V <- metric.sorted.V[ind]
            cor.val.sorted.V <- cor.val.sorted.V[ind]            
            p.val1.sorted.V <- p.val1.sorted.V[ind]
            p.val2.sorted.V <- p.val2.sorted.V[ind]
            FDR1.sorted.V <-   FDR1.sorted.V[ind]

            m.V[marker.ind == k,] <- V3
            row.names(m.V)[marker.ind == k] <- row.names.V3            
            metric.sorted[marker.ind == k] <- metric.sorted.V
            cor.val.sorted[marker.ind == k] <- cor.val.sorted.V            
            p.val1.sorted[marker.ind == k] <- p.val1.sorted.V
            p.val2.sorted[marker.ind == k] <- p.val2.sorted.V
            FDR1.sorted[marker.ind == k] <-   FDR1.sorted.V
         }
         if (metric.sorted[1] < 0 ) {
            m.V <- apply(m.V, MARGIN=2, FUN=rev)
            metric.sorted <- rev(metric.sorted )
            cor.val.sorted <- rev(cor.val.sorted )            
            p.val1.sorted <- rev(p.val1.sorted)
            p.val2.sorted <- rev(p.val2.sorted)
            FDR1.sorted <- rev(FDR1.sorted)
         }            
   } else if (cluster.top.markers == "both.classes") {

      dist.matrix <- dist(m.V) 
      s <- smacofSym(dist.matrix, ndim=1)
      ind <- order(s$conf, decreasing=T)

#      dist.matrix <- dist(m.V)
#      HC <- hclust(dist.matrix, method="complete")
#      ind <- HC$order

      m.V <- m.V[ind,]
      metric.sorted <- metric.sorted[ind]
      cor.val.sorted <- cor.val.sorted[ind]      
      p.val1.sorted <- p.val1.sorted[ind]
      p.val2.sorted <- p.val2.sorted[ind]
      FDR1.sorted <- FDR1.sorted[ind]
      if (metric.sorted[1] < 0 ) {
         m.V <- apply(m.V, MARGIN=2, FUN=rev)
         metric.sorted <- rev(metric.sorted )
         cor.val.sorted <- rev(cor.val.sorted )         
         p.val1.sorted <- rev(p.val1.sorted)
         p.val2.sorted <- rev(p.val2.sorted)
         FDR1.sorted <- rev(FDR1.sorted)
      }            
   }

#      V <- target.V
#      max.v <- max(max(V), -min(V))
#      V1 <- ceiling(ncolors * (V - (- max.v))/(1.001*(max.v - (- max.v))))

      cutoff <- 2.5
      x <- as.numeric(target.V)         
      x <- (x - mean(x))/sd(x)         
      ind1 <- which(x > cutoff)
      ind2 <- which(x < -cutoff)
      x[ind1] <- cutoff
      x[ind2] <- -cutoff
      V1 <- ceiling(ncolors * (x + cutoff)/(cutoff*2))
     
      par(mar = c(1, 14, 2, 10))
      image(1:N, 1:1, as.matrix(V1), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
      axis(2, at=1:1, labels=target.name, adj= 0.5, tick=FALSE,
           las = 1, cex=1, cex.axis=0.65*character.scaling, font.axis=1, line=-1)
      if (assoc.metric.type == "DIFF.MEANS") { # use local p-vals: p.val1
          if (display.cor.coeff.heatmap == T) {
            axis(4, at=1:1, labels=paste(assoc.metric.type, "p-val (loc)", "FDR", "cor"), adj= 0.5, tick=FALSE,
                 las = 1, cex=1, cex.axis=0.65*character.scaling, font.axis=1, line=-1)
          } else {
            axis(4, at=1:1, labels=paste(assoc.metric.type, "p-val (loc)", "FDR"), adj= 0.5, tick=FALSE,
                 las = 1, cex=1, cex.axis=0.65*character.scaling, font.axis=1, line=-1)
          }
      } else {
          if (display.cor.coeff.heatmap == T) {
             axis(4, at=1:1, labels=paste(assoc.metric.type, "p-val (glob)", "FDR", "cor"), adj= 0.5, tick=FALSE,
                 las = 1, cex=1, cex.axis=0.65*character.scaling, font.axis=1, line=-1)
          } else {
             axis(4, at=1:1, labels=paste(assoc.metric.type, "p-val (glob)", "FDR"), adj= 0.5, tick=FALSE,
                 las = 1, cex=1, cex.axis=0.65*character.scaling, font.axis=1, line=-1)
          }
      }
      V <- m.V
      for (i in 1:length(V[,1])) {
#           x <- as.numeric(V[i,])
#           V[i,] <- (x - mean(x))/sd(x)
#           max.v <- max(max(V[i,]), -min(V[i,]))
#           V[i,] <- ceiling(ncolors * (V[i,] - (- max.v))/(1.001*(max.v - (- max.v))))

           cutoff <- 2.5
            x <- as.numeric(V[i,])         
            x <- (x - mean(x))/sd(x)         
            ind1 <- which(x > cutoff)
            ind2 <- which(x < -cutoff)
            x[ind1] <- cutoff
            x[ind2] <- -cutoff
            V[i,] <- ceiling(ncolors * (x + cutoff)/(cutoff*2))
      }
      V <- apply(V, MARGIN=2, FUN=rev)
      par(mar = c(8, 14, 1, 10))
      image(1:dim(V)[2], 1:dim(V)[1], t(V), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
      axis(2, at=1:dim(V)[1], labels=row.names(V), adj= 0.5, tick=FALSE,
           las = 1, cex=1, cex.axis=0.65*character.scaling, font.axis=1, line=-1)
      if (assoc.metric.type == "DIFF.MEANS") { # use local p-vals: p.val1
           if (display.cor.coeff.heatmap == T) {
              mi <- paste(signif(metric.sorted, 3), signif(p.val1.sorted, 3), signif(FDR1.sorted, 3), signif(cor.val.sorted, 3), sep="   ")
            } else {
              mi <- paste(signif(metric.sorted, 3), signif(p.val1.sorted, 3), signif(FDR1.sorted, 3), sep="   ")
            }
      } else {
           if (display.cor.coeff.heatmap == T) {
              mi <- paste(signif(metric.sorted, 3), signif(p.val2.sorted, 3), signif(FDR1.sorted, 3), signif(cor.val.sorted, 3), sep="   ")
           } else {
              mi <- paste(signif(metric.sorted, 3), signif(p.val2.sorted, 3), signif(FDR1.sorted, 3), sep="   ")
           }
      }
      axis(4, at=1:dim(V)[1], labels=rev(mi), adj= 0.5, tick=FALSE,
           las = 1, cex=1, cex.axis=0.65*character.scaling, font.axis=1, line=-1)
      axis(1, at=1:dim(V)[2], labels=colnames(V), adj= 0.5, tick=FALSE,
           las = 3, cex=1, cex.axis=0.3*character.scaling, font.axis=1, line=-1)

      # save top marker dataset

       ds <- rbind(target.V, m.V)
       row.names(ds) <- c(target.name, row.names(m.V))
       write.gct.2(gct.data.frame = ds, descs = row.names(ds), filename = results.file.gct)
      
    }

    # Make MDS top features projection

     if (produce.mds.plots == T) {

       #### New circular version

   nf <- layout(matrix(1, 1, byrow=T), 1, 1, FALSE)

   total.points <- n.markers
   V2 <- m3[1:total.points,]
   row.names(V2) <- row.names(m3)[1:total.points]
   metric.matrix <- matrix(0, nrow=nrow(V2), ncol=nrow(V2))
   row.names(metric.matrix)  <- row.names(V2)
   colnames(metric.matrix) <- row.names(V2)
   MI.ref <- metric.sorted[1:total.points]
   for (i in 1:nrow(V2)) {
      for (j in 1:i) {
          metric.matrix[i, j] <- assoc.metric(target = V2[j,], feature = V2[i,], type="RNMI", knn=knn, n.grid=n.grid)
      }
   }

   metric.matrix <- metric.matrix + t(metric.matrix)

   alpha <- 8
   metric.matrix2 <- 1 - ((1/(1+exp(-alpha*metric.matrix))))
   for (i in 1:nrow(metric.matrix2)) metric.matrix2[i, i] <- 0

   smacof.map <- smacofSphere(metric.matrix2, ndim = 2, weightmat = NULL, init = NULL,
                          ties = "primary", verbose = FALSE, modulus = 1, itmax = 1000, eps = 1e-6)
   x0 <- smacof.map$conf[,1]
   y0 <- smacof.map$conf[,2]
   r <- sqrt(x0*x0 + y0*y0)
   radius <-  1 - ((1/(1+exp(-alpha*MI.ref))))
   x <- x0*radius/r
   y <- y0*radius/r
   angles <- atan2(y0, x0)
   par(mar = c(10, 2, 10, 2))
   plot(x, y, pch=20, bty="n", xaxt='n', axes = FALSE, type="n", xlab="", ylab="",
        xlim=1.2*c(-max(radius), max(radius)), ylim=1.2*c(-max(radius), max(radius)))
   line.angle <- seq(0, 2*pi-0.001, 0.001)
   for (i in 1:length(x)) {
      line.max.x <- radius[i] * cos(line.angle)
      line.max.y <- radius[i] * sin(line.angle)
      points(line.max.x, line.max.y, type="l", col="gray80", lwd=1)
      points(c(0, x[i]), c(0, y[i]), type="l", col="gray80", lwd=1)
   }
   line.max.x <- 1.2*max(radius) * cos(line.angle)
   line.max.y <- 1.2*max(radius) * sin(line.angle)
   points(line.max.x, line.max.y, type="l", col="purple", lwd=2)
   points(x, y, pch=21, bg="steelblue", col="darkblue", cex=1.1)
   points(0, 0, pch=20, col="red", cex=2.5)
   text(0, 0, labels=target.name, cex=1, col="red", pos=1)


   d <- density(angles, adjust=0.10, n=4000, from= -2*pi, to=2*pi)
   dx <- d$x[1001:3000]
   dd <- d$y[1001:3000]
   dd[1:1000] <- dd[1:1000] + d$y[3001:4000]
   dd[1001:2000] <- dd[1001:2000] + d$y[1:1000]

   xd <- 1.2*max(radius)*cos(dx)
   yd <- 1.2*max(radius)*sin(dx)
#   ring.colors <- mycol[1 + 511*(dd - min(dd))/(max(dd) - min(dd))]

   ddd <- (dd - min(dd))/(max(dd) - min(dd))
   ring.colors <- rgb(1, 1 - ddd, 1, maxColorValue = 1)

#   ring.colors <- rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
   for (i in 1:length(xd)) points(xd[i], yd[i], type="p", pch=21, cex=2, bg=ring.colors[i], col=ring.colors[i])
   points(x0*1.2*max(radius)/r, y0*1.2*max(radius)/r, pch=21, bg="purple", col="purple", cex=1.1)

   pointLabel(x, y, paste(seq(1, length(x)), ":", labels=colnames(metric.matrix2)), cex=0.70, col="darkblue")


     }
      
   dev.off()

   time2 <- proc.time()
#   print(paste("Total time:", signif(sum((time2 - time1)[1:2]), digits=3), " secs"))

 }

   assoc.metric <- function(target, feature, cond.feature = NULL, type="RNMI", knn=3, n.grid=25) {
     locs <- seq(1, length(feature))[!is.na(feature)]
     if (length(locs) <= 1) return(0)
     feature <- feature[locs]
     target <- target[locs]
     if (type=="TBBW") { # Rescaled normalized mutual information * Diff. means
         target <- target
         feature <- feature 
         x <- split(feature, target)
         m1 <- mean(x[[order(names(x), decreasing=T)[1]]])
         m2 <- mean(x[[order(names(x), decreasing=T)[2]]])
         TBBW <- abs(m1 - m2)*mutual.inf(x = target, y = feature, n.grid=n.grid)$NMI/mutual.inf(x = target, y = target, n.grid=n.grid)$NMI
         return(TBBW)
     } else if (type=="NMI") { # Normalized mutual information (Kernel method)
         return(mutual.inf(x = target, y = feature, n.grid=n.grid)$NMI)
      } else if (type=="RNMI") { # Rescaled normalized mutual information (Kernel method)
         return(mutual.inf(x = target, y = feature, n.grid=n.grid)$NMI/mutual.inf(x = target, y = target, n.grid=n.grid)$NMI)
      } else if (type=="RNMI.A") { # Rescaled normalized mutual information (Adaptive Kernel method)
         adj <- log(1/(abs(cor(target, feature)) + 0.25)) + 0.7768
         delta <- c(adj * bcv(target), adj * bcv(feature))
         return(mutual.inf(x = target, y = feature, n.grid=n.grid, delta = delta)$NMI/
                mutual.inf(x = target, y = target, n.grid=n.grid, delta= c(bcv(target), bcv(target)))$NMI)
      } else if (type=="NMI.knn") { # Rescaled normalized mutual information (knn method)
         return(sign(cor(target, feature))*(knnmi(x = target, y = feature, k=knn, i1_flag = 1)/knnjh(target, feature, k=3, flag = 1)))
      } else if (type=="NMI.LN.kde") {  # a grid-based kernel estimator of normalized Mutual information
         return(NMI.LN.kde(target, feature, n.grid=n.grid, h = c(bcv(target), bcv(feature)))$NMI)
      } else if (type=="NMI.LN.plugin") {  # a plug in estimator of normalized Mutual information
         return(NMI.LN.plugin(target, feature, h = c(bcv(target), bcv(feature))))
      } else if (type=="NMI.LN.hybrid") {  # hybrid NMI estimator
         return(NMI.LN.hybrid(target, feature, k=25, h = c(bcv(target), bcv(feature))))

      } else if (type=="IC") { # Information Coefficient (Kernel method) 
#         return(mutual.inf.v2(x = target, y = feature, n.grid=n.grid, delta= c(bcv(target), bcv(feature)))$IC)
#       return(mutual.inf.v3(x = target, y = feature, n.grid=n.grid, delta= c(bcv(target), bcv(feature)))$IC)
        return(IC.v1(target, feature, n.grid=n.grid))
      } else if (type=="ICR") { # Information Coefficient (Kernel method) of "Ranked" profiles
        target2 <- rank(target)
        feature2 <- rank(feature)           
        return(IC.v1(target2, feature2, n.grid=n.grid))
      } else if (type=="IC2") { # Information Coefficient (Kernel method) 
        return(IC.v2(target, feature, alpha=3, n.grid=n.grid))
      } else if (type=="IC.knn") { # Information Coefficient (knn method, original parmigene, no boundary correction) 
         return(IC.knn(target, feature, k =knn, flag=0))
      } else if (type=="IC.knn.1") { # Information Coefficient (knn method, with boundary correction) 
         return(IC.knn(target, feature, k =knn, flag=1))
      } else if (type=="SMI") { # Signed mutual information (Kernel method) 
         return(mutual.inf.v2(x = target, y = feature, n.grid=n.grid, delta= c(bcv(target), bcv(feature)))$SMI)
      } else if (type=="SMI.LN") { # Signed mutual information (Kernel method) 
         return(mutual.inf.v2(x = target, y = feature, n.grid=n.grid, delta= c(bcv(target), bcv(feature)))$SMI)
      } else if (type=="SMI.LN.kde") {  # a grid-based kernel estimator of signed Mutual information
         return(sign(cor(target, feature))*NMI.LN.kde(target, feature, n.grid=n.grid, h = c(bcv(target), bcv(feature)))$MI)
      } else if (type=="SMI.LN.plugin") {  # a plug in estimator of signed Mutual information
         return(sign(cor(target, feature))*MI.LN.plugin(target, feature, h = c(bcv(target), bcv(feature))))
      } else if (type=="SMI.LN.hybrid") {  # hybrid SMI estimator
         return(sign(cor(target, feature))*MI.LN.hybrid(target, feature, k=25, h = c(bcv(target), bcv(feature))))
      } else if (type=="SMI.knn") { # Signed mutual information (knn method)
         return(sign(cor(target, feature))*(knnmi(x = target, y = feature,  k =knn, i1_flag = 1)))
      } else if (type=="SMI.skew.t") { # Signed mutual information (skew.t method)
         return(sign(cor(target, feature))*(MI.skew.t(x = target, y = feature)))
      } else if (type=="SCMI") { # Signed Conditional mutual information (Kernel method) 
         return(sign(cor(target, feature))*cond.mutual.inf(x = target, y = feature, z=cond.feature, n.grid=n.grid)$CMI)
      } else if (type=="CIC") { # Conditional Information Coefficient (Kernel method) 
         return(cond.mutual.inf(x = target, y = feature, z=cond.feature, n.grid=n.grid)$CIC)
      } else if (type=="DIFF.MEANS") {
         x <- split(feature, target)
         m1 <- mean(x[[order(names(x), decreasing=T)[1]]])
         m2 <- mean(x[[order(names(x), decreasing=T)[2]]])
         return(m1 - m2)
      } else if (type=="AUC.ROC") {
         target <- match(target, sort(unique(target))) - 1
         perf.auc <- roc.area(obs=target, pred= (feature - min(feature))/(max(feature) - min(feature)))
         return(perf.auc$A)
      } else if (type=="AUC.REC") {
         perf.auc <- rec.area(obs=target, pred= (feature - min(feature))/(max(feature) - min(feature)), metric = "squared.error")
         return(perf.auc$A)
      } else if (type=="DIFF.MEDIANS") {
         x <- split(feature, target)
         m1 <- median(x[[order(names(x), decreasing=T)[1]]])
         m2 <- median(x[[order(names(x), decreasing=T)[2]]])
         return(m1 - m2)
      } else if (type=="T.TEST") {
         x <- split(feature, target)
         return(t.test(x=x[[order(names(x), decreasing=T)[1]]], y=x[[order(names(x), decreasing=T)[2]]])$statistic)
      } else if (type=="S2N") {
         x <- split(feature, target)
         m1 <- mean(x[[order(names(x), decreasing=T)[1]]])
         m2 <- mean(x[[order(names(x), decreasing=T)[2]]])
         s1 <- ifelse(length(x[[order(names(x), decreasing=T)[1]]]) > 1, sd(x[[order(names(x), decreasing=T)[1]]]), 0)
         s2 <- ifelse(length(x[[order(names(x), decreasing=T)[2]]]) > 1, sd(x[[order(names(x), decreasing=T)[2]]]), 0)
         s1 <- ifelse(s1 < 0.1*abs(m1), 0.1*abs(m1), s1)
         s2 <- ifelse(s2 < 0.1*abs(m2), 0.1*abs(m2), s2)
         return((m1 - m2)/(s1 + s2 + 0.1))
      } else if (type=="CORR") {           
         return(cor(target, feature))
      } else if (type=="SPEAR") {           
         return(cor(target, feature, method = "spearman"))
      } else {
         stop(paste("Unknow assoc. metric name:", type))
      }
   }

mutual.inf <- function(x, y, n.grid=25, delta = c(bcv(x), bcv(y))) {
    # for definitions of mutual information and the universal metric (NMI) see the 
    # definition of "Mutual Information" in wikipedia and Thomas and Cover's book

   kde2d.xy <- kde2d(x, y, n = n.grid, h = delta)
   X <- kde2d.xy$x
   Y <- kde2d.xy$y
   PXY <- kde2d.xy$z + .Machine$double.eps
   PXY <- PXY/sum(PXY)
   PX <- apply(PXY, MARGIN=1, sum)
   PX <- PX/sum(PX)
   HX <- -sum(PX * log2(PX))
   PX <- matrix(PX, nrow=n.grid, ncol=n.grid)

   PY <- apply(PXY, MARGIN=2, sum)
   PY <- PY/sum(PY)
   HY <- -sum(PY * log2(PY))
   PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)

   MI <- sum(PXY * log2(PXY/(PX*PY)))
   SMI <- sign(cor(x, y)) * MI
   HXY <- - sum(PXY * log2(PXY))

   NMI <- sign(cor(x, y)) * ((HX + HY)/HXY - 1)  # use peason correlation the get the sign (directionality)

   return(list(MI=MI, SMI=SMI, HXY=HXY, HX=HX, HY=HY, NMI=NMI))
}

mutual.inf.v2 <- function(x, y, n.grid=25, delta = c(bcv(x), bcv(y))) {

    # For definitions of mutual information and the universal metric (NMI) see the 
    # definition of "Mutual Information" in wikipedia and Thomas and Cover's book

   rho <- cor(x, y)
   rho2 <- abs(rho)
   delta <- delta*(1 + (-0.75)*rho2)

   kde2d.xy <- kde2d(x, y, n = n.grid, h = delta)
   FXY <- kde2d.xy$z + .Machine$double.eps
   dx <- kde2d.xy$x[2] - kde2d.xy$x[1]
   dy <- kde2d.xy$y[2] - kde2d.xy$y[1]
   PXY <- FXY/(sum(FXY)*dx*dy)
   PX <- rowSums(PXY)*dy
   PY <- colSums(PXY)*dx
   HXY <- -sum(PXY * log(PXY))*dx*dy
   HX <- -sum(PX * log(PX))*dx
   HY <- -sum(PY * log(PY))*dy

   PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
   PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)

   MI <- sum(PXY * log(PXY/(PX*PY)))*dx*dy
   rho <- cor(x, y)
   SMI <- sign(rho) * MI
   
   IC <- sign(rho) * sqrt(1 - exp(- 2 * MI)) 
   
   NMI <- sign(cor(x, y)) * ((HX + HY)/HXY - 1)  # use peason correlation the get the sign (directionality)

   return(list(MI=MI, SMI=SMI, HXY=HXY, HX=HX, HY=HY, NMI=NMI, IC=IC))
}

IC.v1 <-  function(x, y, n.grid=25) {

    # For definitions of mutual information and the universal metric (NMI) see the 
    # definition of "Mutual Information" in wikipedia and Thomas and Cover's book

   x.set <- !is.na(x)
   y.set <- !is.na(y)
   overlap <- x.set & y.set

  # x <- x[overlap] +  0.000000001*runif(length(overlap))
  # y <- y[overlap] +  0.000000001*runif(length(overlap))

   x <- x[overlap] 
   y <- y[overlap] 

   if (length(x) > 2) {
   
      delta = c(bcv(x), bcv(y))
   
      rho <- cor(x, y)
      rho2 <- abs(rho)
      delta <- delta*(1 + (-0.75)*rho2)
      kde2d.xy <- kde2d(x, y, n = n.grid, h = delta)
      FXY <- kde2d.xy$z + .Machine$double.eps
      dx <- kde2d.xy$x[2] - kde2d.xy$x[1]
      dy <- kde2d.xy$y[2] - kde2d.xy$y[1]
      PXY <- FXY/(sum(FXY)*dx*dy)
      PX <- rowSums(PXY)*dy
      PY <- colSums(PXY)*dx
      HXY <- -sum(PXY * log(PXY))*dx*dy
      HX <- -sum(PX * log(PX))*dx
      HY <- -sum(PY * log(PY))*dy
      PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
      PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)
      MI <- sum(PXY * log(PXY/(PX*PY)))*dx*dy
      IC <- sign(rho) * sqrt(1 - exp(- 2 * MI))
      if (is.na(IC)) IC <- 0
   } else {
      IC <- 0
   }

   return(IC)
}

IC.knn <-  function(x, y, k, flag = 0) {

    # Information Coefficient using knn/parmigene
    # k = number of neighbors
    # flag = 0 (original parmigene, no boundary correction), 1 (with boundary correction)

   x.set <- !is.na(x)
   y.set <- !is.na(y)
   overlap <- x.set & y.set

   x <- x[overlap] +  0.000000001*runif(length(overlap))
   y <- y[overlap] +  0.000000001*runif(length(overlap))

   if (length(x) > 2) {
     MI <- knnmi(x = x, y = y, k = k, i1_flag = flag)
     IC <- sign(cor(x, y)) * sqrt(1 - exp(- 2 * MI))
   } else {
      IC <- 0
   }
   return(IC)
 }

IC.v2 <-  function(x, y, n.grid=25, alpha = 1, delta = c(bcv(x), bcv(y))) {

    # For definitions of mutual information and the universal metric (NMI) see the 
    # definition of "Mutual Information" in wikipedia and Thomas and Cover's book

    # This version supports "boosted" MI (alpha parameter)

   x.set <- !is.na(x)
   y.set <- !is.na(y)
   overlap <- x.set & y.set
   x <- x[overlap] +  0.000000001*runif(length(overlap))
   y <- y[overlap] +  0.000000001*runif(length(overlap))

   rho <- cor(x, y)
   rho2 <- abs(rho)
   delta <- delta*(1 + (-0.75)*rho2)
   kde2d.xy <- kde2d(x, y, n = n.grid, h = delta)
   FXY <- kde2d.xy$z + .Machine$double.eps
   dx <- kde2d.xy$x[2] - kde2d.xy$x[1]
   dy <- kde2d.xy$y[2] - kde2d.xy$y[1]
   PXY <- FXY/(sum(FXY)*dx*dy)
   PX <- rowSums(PXY)*dy
   PY <- colSums(PXY)*dx
   HXY <- -sum(PXY * log(PXY))*dx*dy
   HX <- -sum(PX * log(PX))*dx
   HY <- -sum(PY * log(PY))*dy
   PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
   PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)
   MI.alpha <- sum(PXY * (log(PXY/(PX*PY))) ^ alpha) *dx*dy
   IC.alpha <- sign(rho) * sqrt(1 - exp(- 2 * MI.alpha)) 

   return(list(IC.alpha=IC.alpha, HX=HX, HY=HY))
}

 
mutual.inf.v3 <- function(x, y, n.grid=25, delta = c(bcv(x), bcv(y))) {

    # For definitions of mutual information and the universal metric (NMI) see the 
    # definition of "Mutual Information" in wikipedia and Thomas and Cover's book

   rho <- cor(x, y)
   rho2 <- abs(rho)
   delta <- delta*(1 + (-0.75)*rho2)

   kde2d.xy <- kde2d(x, y, n = n.grid, h = delta)
   FXY <- kde2d.xy$z + .Machine$double.eps
   dx <- kde2d.xy$x[2] - kde2d.xy$x[1]
   dy <- kde2d.xy$y[2] - kde2d.xy$y[1]
   PXY <- FXY/(sum(FXY)*dx*dy)
   PX <- rowSums(PXY)*dy
   PY <- colSums(PXY)*dx
   HXY <- -sum(PXY * log(PXY))*dx*dy
   HX <- -sum(PX * log(PX))*dx
   HY <- -sum(PY * log(PY))*dy

   PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
   PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)

   MI <- sum(PXY * log(PXY/(PX*PY)))*dx*dy

   SMI <- sign(rho) * MI
   
   IC <- sign(rho) * sqrt(1 - exp(- 2 * MI)) 
   
   NMI <- sign(cor(x, y)) * ((HX + HY)/HXY - 1)  # use peason correlation the get the sign (directionality)

   return(list(MI=MI, SMI=SMI, HXY=HXY, HX=HX, HY=HY, NMI=NMI, IC=IC))
}


mutual.inf.v4 <- function(x, y, n.grid=25, delta = c(bcv(x), bcv(y))) {

    # For definitions of mutual information and the universal metric (NMI) see the 
    # definition of "Mutual Information" in wikipedia and Thomas and Cover's book

   rho <- cor(x, y)
   rho2 <- ifelse(rho < 0, 0, rho)
   delta <- delta*(1 + (-0.75)*rho2)

   kde2d.xy <- kde2d(x, y, n = n.grid, h = delta)
   FXY <- kde2d.xy$z + .Machine$double.eps
   dx <- kde2d.xy$x[2] - kde2d.xy$x[1]
   dy <- kde2d.xy$y[2] - kde2d.xy$y[1]
   PXY <- FXY/(sum(FXY)*dx*dy)
   PX <- rowSums(PXY)*dy
   PY <- colSums(PXY)*dx
   HXY <- -sum(PXY * log(PXY))*dx*dy
   HX <- -sum(PX * log(PX))*dx
   HY <- -sum(PY * log(PY))*dy

   PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
   PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)

   MI <- sum(PXY * log(PXY/(PX*PY)))*dx*dy

   MI2 <- sum(PXY * (log(PXY/(PX*PY))) ^ 2) *dx*dy

   SMI <- sign(rho) * MI
   
   IC <- sign(rho) * sqrt(1 - exp(- 2 * MI))

   IC2 <- sign(rho) * sqrt(1 - exp(- 2 * MI2)) 
   
   NMI <- sign(cor(x, y)) * ((HX + HY)/HXY - 1)  # use peason correlation the get the sign (directionality)

   return(list(MI=MI, SMI=SMI, HXY=HXY, HX=HX, HY=HY, NMI=NMI, IC=IC, IC2=IC2))
}



qvalue <- function(p=NULL, lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=NULL, robust=FALSE, 
  gui=FALSE, smooth.df = 3, smooth.log.pi0 = FALSE) {
#Input
#=============================================================================
#p: a vector of p-values (only necessary input)
#fdr.level: a level at which to control the FDR (optional)
#lambda: the value of the tuning parameter to estimate pi0 (optional)
#pi0.method: either "smoother" or "bootstrap"; the method for automatically
#           choosing tuning parameter in the estimation of pi0, the proportion
#           of true null hypotheses
#robust: an indicator of whether it is desired to make the estimate more robust
#        for small p-values and a direct finite sample estimate of pFDR (optional)
#gui: A flag to indicate to 'qvalue' that it should communicate with the gui.  ## change by Alan
#     Should not be specified on command line.
#smooth.df: degrees of freedom to use in smoother (optional)
#smooth.log.pi0: should smoothing be done on log scale? (optional)
#
#Output
#=============================================================================
#call: gives the function call
#pi0: an estimate of the proportion of null p-values
#qvalues: a vector of the estimated q-values (the main quantity of interest)
#pvalues: a vector of the original p-values
#significant: if fdr.level is specified, an indicator of whether the q-value
#    fell below fdr.level (taking all such q-values to be significant controls
#    FDR at level fdr.level)

#Set up communication with GUI, if appropriate
#    print(sys.calls())
#    print(sys.frames())

#    if(gui) {
#        idx <- (1:sys.nframe())[as.character(sys.calls()) == "qvalue.gui()"]
#        gui.env <- sys.frames()[[idx]]
#    }

#This is just some pre-processing
    if(is.null(p))  ## change by Alan
      {qvalue.gui(); return("Launching point-and-click...")}
    if(gui & !interactive())  ## change by Alan
      gui = FALSE

    if(min(p)<0 || max(p)>1) {
      if(gui) ## change by Alan:  check for GUI
        eval(expression(postMsg(paste("ERROR: p-values not in valid range.", "\n"))), parent.frame())
      else
#        print("ERROR: p-values not in valid range.")
      return(0)
    }
    if(length(lambda)>1 && length(lambda)<4) {
      if(gui)
        eval(expression(postMsg(paste("ERROR: If length of lambda greater than 1, you need at least 4 values.",
            "\n"))), parent.frame())
      else
#        print("ERROR: If length of lambda greater than 1, you need at least 4 values.")
      return(0)
    }
    if(length(lambda)>1 && (min(lambda) < 0 || max(lambda) >= 1)) { ## change by Alan:  check for valid range for lambda
      if(gui)
        eval(expression(postMsg(paste("ERROR: Lambda must be within [0, 1).", "\n"))), parent.frame())
      else
#        print("ERROR: Lambda must be within [0, 1).")
      return(0)
    }
    m <- length(p)
#These next few functions are the various ways to estimate pi0
    if(length(lambda)==1) {
        if(lambda<0 || lambda>=1) { ## change by Alan:  check for valid range for lambda
          if(gui)
            eval(expression(postMsg(paste("ERROR: Lambda must be within [0, 1).", "\n"))), parent.frame())
          else
#            print("ERROR: Lambda must be within [0, 1).")
          return(0)
        }

        pi0 <- mean(p >= lambda)/(1-lambda)
        pi0 <- min(pi0,1)
    }
    else {
        pi0 <- rep(0,length(lambda))
        for(i in 1:length(lambda)) {
            pi0[i] <- mean(p >= lambda[i])/(1-lambda[i])
        }

        if(pi0.method=="smoother") {
            if(smooth.log.pi0)
              pi0 <- log(pi0)

            spi0 <- smooth.spline(lambda,pi0,df=smooth.df)
            pi0 <- predict(spi0,x=max(lambda))$y

            if(smooth.log.pi0)
              pi0 <- exp(pi0)
            pi0 <- min(pi0,1)
        }
        else if(pi0.method=="bootstrap") {
            minpi0 <- min(pi0)
            mse <- rep(0,length(lambda))
            pi0.boot <- rep(0,length(lambda))
            for(i in 1:100) {
                p.boot <- sample(p,size=m,replace=TRUE)
                for(i in 1:length(lambda)) {
                    pi0.boot[i] <- mean(p.boot>lambda[i])/(1-lambda[i])
                }
                mse <- mse + (pi0.boot-minpi0)^2
            }
            pi0 <- min(pi0[mse==min(mse)])
            pi0 <- min(pi0,1)
        }
        else {  ## change by Alan: check for valid choice of 'pi0.method' (only necessary on command line)
#            print("ERROR: 'pi0.method' must be one of 'smoother' or 'bootstrap'.")
            return(0)
        }
    }
    if(pi0 <= 0) {
      if(gui)
        eval(expression(postMsg(
            paste("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.",
                "\n"))), parent.frame())
      else
#        print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.")
      return(0)
    }
    if(!is.null(fdr.level) && (fdr.level<=0 || fdr.level>1)) {  ## change by Alan:  check for valid fdr.level
      if(gui)
        eval(expression(postMsg(paste("ERROR: 'fdr.level' must be within (0, 1].", "\n"))), parent.frame())
      else
#        print("ERROR: 'fdr.level' must be within (0, 1].")
      return(0)
    }
#The estimated q-values calculated here
    u <- order(p)

    # change by Alan
    # ranking function which returns number of observations less than or equal
    qvalue.rank <- function(x) {
      idx <- sort.list(x)

      fc <- factor(x)
      nl <- length(levels(fc))
      bin <- as.integer(fc)
      tbl <- tabulate(bin)
      cs <- cumsum(tbl)
 
      tbl <- rep(cs, tbl)
      tbl[idx] <- tbl

      return(tbl)
    }

    v <- qvalue.rank(p)
    
    qvalue <- pi0*m*p/v
    if(robust) {
        qvalue <- pi0*m*p/(v*(1-(1-p)^m))
    }
    qvalue[u[m]] <- min(qvalue[u[m]],1)
    for(i in (m-1):1) {
    qvalue[u[i]] <- min(qvalue[u[i]],qvalue[u[i+1]],1)
    }
#The results are returned
    if(!is.null(fdr.level)) {
        retval <- list(call=match.call(), pi0=pi0, qvalues=qvalue, pvalues=p, fdr.level=fdr.level, ## change by Alan
          significant=(qvalue <= fdr.level), lambda=lambda)
    }
    else {
        retval <- list(call=match.call(), pi0=pi0, qvalues=qvalue, pvalues=p, lambda=lambda)
    }
    class(retval) <- "qvalue"
    return(retval)
}

  SE_assoc <- function(x, y) {
           return(2 - sqrt(mean((x - y)^2)))
       }

REVEALER_assess_features.v1 <- function(
      input_dataset,
      target,
      direction              = "positive",
      feature.files,                           
      features,                                
      output.file,
      description            = "",
      sort.by.target         = T,
      character.scaling      = 1.5,
      n.perm                 = 10000,
      create.feature.summary = F,
      feature.combination.op = "max")
 {
   suppressPackageStartupMessages(library(maptools))
   suppressPackageStartupMessages(library(RColorBrewer))

   set.seed(5209761)

   missing.value.color <- "khaki1"
   mycol <- vector(length=512, mode = "numeric")
   for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
   for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
   mycol <- rev(mycol)
   max.cont.color <- 512
   mycol <- c(mycol,
              missing.value.color)                  # Missing feature color

   categ.col <- c("#9DDDD6", # dusty green
                     "#F0A5AB", # dusty red
                     "#9AC7EF", # sky blue
                     "#F970F9", # violet
                     "#FFE1DC", # clay
                     "#FAF2BE", # dusty yellow
                     "#AED4ED", # steel blue
                     "#C6FA60", # green
                     "#D6A3FC", # purple
                     "#FC8962", # red
                     "#F6E370", # orange
                     "#F0F442", # yellow
                     "#F3C7F2", # pink
                     "#D9D9D9", # grey
                     "#FD9B85", # coral
                     "#7FFF00", # chartreuse
                     "#FFB90F", # goldenrod1
                     "#6E8B3D", # darkolivegreen4
                     "#8B8878", # cornsilk4
                     "#7FFFD4") # aquamarine

   cex.size.table <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0.9,   # 1-10 characters
                       0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, # 11-20 characters
                       0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)

   pdf(file=output.file, height=14, width=11)


   n.panels <- length(feature.files) 
   l.panels <- NULL
   for (l in 1:n.panels) l.panels <- c(l.panels, 1.5, length(features[[l]]))
   l.panels[l.panels < 2] <- 1.5
   empty.panel <- 30 - sum(unlist(l.panels))
   l.panels <- c(l.panels, empty.panel)
   n.panels <- length(l.panels)

   nf <- layout(matrix(c(seq(1, n.panels - 1), 0), n.panels, 1, byrow=T), 1, l.panels,  FALSE)
                      
   for (f in 1:length(feature.files)) {   # loop over feature types
#      print(paste("Processing feature file:", feature.files[[f]]))
      dataset <- MSIG.Gct2Frame(filename = input_dataset)
      m.1 <- data.matrix(dataset$ds)
      sample.names.1 <- colnames(m.1)
      Ns.1 <- ncol(m.1)
#      print(paste("Total samples in input file:", Ns.1))

      target.vec <- m.1[target,]

      non.nas <- !is.na(target.vec)

      target.vec <- target.vec[non.nas]
      sample.names.1 <- sample.names.1[non.nas]
      Ns.1 <- length(target.vec)
      m.1 <- m.1[, non.nas]

#      print(dim(m.1))
      
      dataset.2 <- MSIG.Gct2Frame(filename = feature.files[[f]])
      m.2 <- data.matrix(dataset.2$ds)
      dim(m.2)
      row.names(m.2) <- dataset.2$row.names
      Ns.2 <- ncol(m.2)  
      sample.names.2 <- colnames(m.2) <- dataset.2$names
      
      overlap <- intersect(sample.names.1, sample.names.2)
      locs1 <- match(overlap, sample.names.1)
      locs2 <- match(overlap, sample.names.2)
      m.1 <- m.1[, locs1]
      target.vec <- target.vec[locs1]
      m.2 <- m.2[, locs2]
      Ns.1 <- ncol(m.1)
      Ns.2 <- ncol(m.2)
      sample.names.1 <- colnames(m.1)
      sample.names.2 <- colnames(m.2)
      
#      print(paste("feature file overlap with target samples:", ncol(m.2)))
      
      if (sort.by.target == T) {
         if (direction == "positive") {
            ind <- order(target.vec, decreasing=T)
         } else {
            ind <- order(target.vec, decreasing=F)
         }
         target.vec <- target.vec[ind]
         sample.names.1 <- sample.names.1[ind]
         m.1 <- m.1[, ind]
         m.2 <- m.2[, ind]         
         sample.names.2 <- sample.names.2[ind]
      }

    # normalize target
      target.vec.orig <- target.vec
      unique.target.vals <- unique(target.vec)
      n.vals <- length(unique.target.vals)
      if (n.vals >= length(target.vec)*0.5) {    # Continuous value color map        
         cutoff <- 2.5
         x <- target.vec
         x <- (x - mean(x))/sd(x)         
         x[x > cutoff] <- cutoff
         x[x < - cutoff] <- - cutoff      
         x <- ceiling((max.cont.color - 1) * (x + cutoff)/(cutoff*2)) + 1
         target.vec <- x
      }
      if (f == 1) {
          main <- description
      } else {
          main <- ""
      }
      par(mar = c(0, 22, 2, 9))

      target.nchar <- ifelse(nchar(target) > 30, 30, nchar(target))
      cex.axis <- cex.size.table[target.nchar]
#      print(paste("cex.axis:", cex.axis))      
      
      if (n.vals >= length(target.vec)*0.5) {    # Continuous value color map        
          image(1:Ns.1, 1:1, as.matrix(target.vec), zlim = c(0, max.cont.color), col=mycol[1: max.cont.color],
                axes=FALSE, main=main, sub = "", xlab= "", ylab="")
      } else if (n.vals == 2) {  # binary
         image(1:Ns.1, 1:1, as.matrix(target.vec), zlim = range(target.vec), col=brewer.pal(9, "Blues")[3], brewer.pal(9, "Blues")[9],
               axes=FALSE, main=main, sub = "", xlab= "", ylab="")
      } else {  # categorical
         image(1:Ns.1, 1:1, as.matrix(target.vec), zlim = range(target.vec), col=categ.col[1:n.vals],
               axes=FALSE, main=main, sub = "", xlab= "", ylab="")
      }
      axis(2, at=1:1, labels=target, adj= 0.5, tick=FALSE, las = 1, cex=1, cex.axis=cex.axis*character.scaling,
           font.axis=1, line=0, font=2, family="")
      axis(4, at=1:1, labels=paste("   IC       p-val"), adj= 0.5, tick=FALSE, las = 1, cex=1, cex.axis=0.7*character.scaling,
           font.axis=1, line=0, font=2, family="")

      feature.mat <- feature.names <- NULL
      
     for (feat.n in 1:length(features[[f]])) { 
        len <- length(unlist(features[[f]][feat.n]))         
        feature.name <- unlist(features[[f]][feat.n])
#        print(paste("      Feature:", feature.name))
        if (is.na(match(feature.name, row.names(m.2)))) next
        feature.mat <- rbind(feature.mat,  m.2[feature.name,])
        feature.names <- c(feature.names, feature.name)
      }
      feature.mat <- as.matrix(feature.mat)
      row.names(feature.mat) <- feature.names

      if (create.feature.summary == T) {
         summary.feature <- apply(feature.mat, MARGIN=2, FUN=feature.combination.op) + 2
         feature.mat <- rbind(feature.mat, summary.feature)
         row.names(feature.mat) <- c(feature.names, "SUMMARY FEATURE")
     }

      for (i in 1:nrow(feature.mat)) {

           feature.vec <- feature.mat[i,]
           unique.feature.vals <- unique(sort(feature.vec))
           non.NA.vals <- sum(!is.na(feature.vec))
           n.vals <- length(unique.feature.vals)
           if (n.vals > 2) {    # Continuous value color map        
              feature.vals.type <- "continuous"
              cutoff <- 2.5
              x <- feature.vec
              locs.non.na <- !is.na(x)
              x.nonzero <- x[locs.non.na]
              x.nonzero <- (x.nonzero - mean(x.nonzero))/sd(x.nonzero)         
              x.nonzero[x.nonzero > cutoff] <- cutoff
              x.nonzero[x.nonzero < - cutoff] <- - cutoff      
              feature.vec[locs.non.na] <- x.nonzero
              feature.vec2 <- feature.vec
              feature.vec[locs.non.na] <- ceiling((max.cont.color - 2) * (feature.vec[locs.non.na] + cutoff)/(cutoff*2)) + 1
              feature.vec[is.na(x)] <- max.cont.color + 1
              feature.mat[i,] <- feature.vec              
           }
      }
      feature.mat <- as.matrix(feature.mat)

      # compute IC association with target

      IC.vec <- p.val.vec <- stats.vec <- NULL
      sqr_error.vec <- roc.vec <- NULL
      
      for (i in 1:nrow(feature.mat)) {
           feature.vec <- feature.mat[i,]
#           IC <- IC.v1(target.vec, feature.vec)
           IC <- mutual.inf.v2(x = target.vec.orig, y = feature.vec, n.grid=25)$IC

           feature.vec.0.1 <- (feature.vec - min(feature.vec))/(max(feature.vec) - min(feature.vec))           
           target.vec.0.1 <- (target.vec - min(target.vec))/(max(target.vec) - min(target.vec))
           if (direction == "negative") target.vec.0.1 <- 1 - target.vec.0.1
           sqr_error <- SE_assoc(target.vec.0.1, feature.vec.0.1)           
           roc <- roc.area(feature.vec.0.1, target.vec.0.1)$A
           IC <- signif(IC, 3)
           null.IC <- vector(length=n.perm, mode="numeric")
#           for (h in 1:n.perm) null.IC[h] <- IC.v1(feature.vec, sample(target.vec))
       for (h in 1:n.perm) null.IC[h] <- mutual.inf.v2(x = sample(target.vec.orig), y = feature.vec, n.grid=25)$IC
           if (IC >= 0) {
             p.val <- sum(null.IC >= IC)/n.perm
           } else {
             p.val <- sum(null.IC <= IC)/n.perm
           }
           p.val <- signif(p.val, 3)
           if (p.val == 0) {
             p.val <- paste("<", signif(1/n.perm, 3), sep="")
           }

           IC.vec <- c(IC.vec, IC)
           sqr_error.vec <- c(sqr_error.vec, sqr_error)
           roc.vec <- c(roc.vec, roc)
           p.val.vec <- c(p.val.vec, p.val)
           space.chars <- "           "
           IC.char <- nchar(IC)
           pad.char <- substr(space.chars, 1, 10 - IC.char)
           stats.vec <- c(stats.vec, paste(IC, pad.char, p.val, sep=""))
       }

      if (nrow(feature.mat) > 1) {
          V <- apply(feature.mat, MARGIN=2, FUN=rev)
      } else {
          V <- as.matrix(feature.mat)
      }
      
      features.max.nchar <- max(nchar(row.names(V)))
      features.nchar <- ifelse(features.max.nchar > 30, 30, features.max.nchar)
      cex.axis <- cex.size.table[features.nchar]
#      print(paste("cex.axis:", cex.axis))

      par(mar = c(1, 22, 1, 9))
      
     if (n.vals > 2) {    # Continuous value color map        
         image(1:dim(V)[2], 1:dim(V)[1], t(V), zlim = c(0, max.cont.color + 3), col=mycol, axes=FALSE, main=main,
               cex.main=0.8, sub = "", xlab= "", ylab="")
     } else {  # binary
         image(1:dim(V)[2], 1:dim(V)[1], t(V), zlim = c(0, 3), col=c(brewer.pal(9, "Blues")[3], brewer.pal(9, "Blues")[9],
                                                                    brewer.pal(9, "Greys")[2], brewer.pal(9, "Greys")[5]),
               axes=FALSE, main="", cex.main=0.8,  sub = "", xlab= "", ylab="")
     }
     axis(2, at=1:dim(V)[1], labels=row.names(V), adj= 0.5, tick=FALSE, las = 1, cex=1, cex.axis=cex.axis*character.scaling,
           font.axis=1, line=0, font=2, family="")
     axis(4, at=1:dim(V)[1], labels=rev(stats.vec), adj= 0.5, tick=FALSE, las = 1, cex=1, cex.axis=0.7*character.scaling,
           font.axis=1, line=0, font=2, family="")
     }
   dev.off()

#   print(cbind(rev(row.names(V)), IC.vec, sqr_error.vec, roc.vec), )
   
}

args <- commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--ds1", dest="ds1"),
  make_option("--target.name", dest="target.name"),
  make_option("--target.match", dest="target.match"),
  make_option("--ds2", dest="ds2"),
  make_option("--seed.names", dest="seed.names"),            
  make_option("--max.n.iter", dest="max.n.iter"),
#  make_option("--identifier", dest="identifier"),
#  make_option("--n.perm", dest="n.perm"),
  make_option("--n.markers", dest="n.markers"),
  make_option("--locs.table.file", dest="locs.table.file"),
#  make_option("--exclude.features", dest="exclude.features"),
  make_option("--count.thres.low", dest="count.thres.low"),                                
  make_option("--count.thres.high", dest="count.thres.high"),                                
  make_option("--pdf.output.file", dest="pdf.output.file")
)

opt <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE, args=args)
opts <- opt$options

# ex.fea.I <- as.character(opts$exclude.features)
# print(paste("ex.fea I:", ex.fea.I))
# ex.fea.II <- strsplit(ex.fea.I, split=",")
# print(paste("ex.fea II:", ex.fea.II))
# ex.fea.III <- ex.fea.II[[1]]
# print(paste("ex.fea III:", ex.fea.III))


 seed.names.I <- as.character(opts$seed.names)
# print(paste("seed.names I:", seed.names.I))
 seed.names.II <- strsplit(seed.names.I, split=",")
# print(paste("seed.names II:", seed.names.II))
 seed.names.III <- seed.names.II[[1]]
# print(seed.names.III[1])
# print(seed.names.III[2])
# print(seed.names.III[3])

if (seed.names.III[1] == "NULL") {
    seed_names <- NULL
} else {
    seed_names <- seed.names.III
}


REVEALER.v1(ds1 = opts$ds1,
            target.name = opts$target.name,
            target.match = opts$target.match,
            ds2 = opts$ds2,
            seed.names = seed_names,
            max.n.iter = as.numeric(opts$max.n.iter),
#            identifier = opts$identifier,
#            n.perm = as.numeric(opts$n.perm),
            n.markers = as.numeric(opts$n.markers),
            locs.table.file = opts$locs.table.file,
#            exclude.features = opts$exclude.features,
#            exclude.features = ex.fea.III,            
            count.thres.low = as.numeric(opts$count.thres.low),
            count.thres.high = as.numeric(opts$count.thres.high),
            pdf.output.file = opts$pdf.output.file)

#  REVEALER.v1(ds1, target.name, target.match, ds2, seed.names, max.n.iter, identifier, n.perm, n.markers, 
#              locs.table.file, exclude.features, count.thres.low, count.thres.high, pdf.output.file)
# command line:
# <R3.0_script> <libdir>/REVEALER_library.v5C.R --ds1=<ds1> --target.name = <target.name> --target.match = <target.match> --ds2 = <ds2> --seed.names = <seed.names> --max.n.iter = <max.n.iter> --identifier = <identifier> --n.perm = <n.perm> --n.markers = <n.markers> --locs.table.file = <locs.table.file> --excludes.features = <exclude.features> --count.thres.low=<count.thres.low> --count.thres.high=<count.thres.high> --pdf.output.file = <pdf.output.file>
