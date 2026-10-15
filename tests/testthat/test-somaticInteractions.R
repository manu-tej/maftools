# somaticInteractions() builds every 2x2 table from one co-mutation count
# matrix and runs fisher.test() on the 2x2 counts. These tests compare its
# output against a frozen copy of the previous implementation, which ran
# fisher.test(x, y) on the two 0/1 vectors and rebuilt each table with
# apply/paste/table. The results must be identical.

# previous implementation of the statistics, from R/somaticInteractions.R at
# commit 015a4cf (lines 43-67, 71-140 and 152), without the plot
referenceInteractions = function(maf, top = 25, genes = NULL, returnAll = TRUE, pvalue = c(0.05, 0.01),
                                 showCounts = FALSE, countStats = "all", ...){
  # ... absorbs the other plot-only arguments of somaticInteractions()
  if(is.null(genes)){
    genes = getGeneSummary(x = maf)[1:top, Hugo_Symbol]
  }
  om = maftools:::createOncoMatrix(m = maf, g = genes)
  all.tsbs = levels(getSampleSummary(x = maf)[,Tumor_Sample_Barcode])
  mutMat = t(om$numericMatrix)
  missing.tsbs = all.tsbs[!all.tsbs %in% rownames(mutMat)]
  mutMat[mutMat > 0 ] = 1
  if(length(missing.tsbs) > 0){
    missing.tsbs = as.data.frame(matrix(data = 0, nrow = length(missing.tsbs), ncol = ncol(mutMat)),
                                 row.names = missing.tsbs)
    colnames(missing.tsbs) = colnames(mutMat)
    mutMat = rbind(mutMat, missing.tsbs)
  }
  interactions = sapply(1:ncol(mutMat), function(i)
    sapply(1:ncol(mutMat), function(j) {
      f = try(fisher.test(mutMat[, i], mutMat[, j]), silent = TRUE)
      if (class(f) == "try-error"){
        if(all(mutMat[,i] == mutMat[,j])){
          if(colnames(mutMat)[i] != colnames(mutMat)[j]){
            warning("All the samples are in the same direction for the genes ", colnames(mutMat)[i], " and ",  colnames(mutMat)[j], "! Could not perform Fisher test.")
          }
          NA
        }else{
          if(colnames(mutMat)[i] != colnames(mutMat)[j]){
            warning("Contigency table could not created for the genes ", colnames(mutMat)[i], " and ",  colnames(mutMat)[j], "! Could not perform Fisher test.")
          }
          NA
        }
      }else{
        ifelse(f$estimate > 1,-log10(f$p.val), log10(f$p.val))
      }
    }))
  oddsRatio <-
    sapply(1:ncol(mutMat), function(i)
      sapply(1:ncol(mutMat), function(j) {
        f = try(fisher.test(mutMat[, i], mutMat[, j]), silent = TRUE)
        if (class(f) == "try-error")
          NA
        else
          f$estimate
      }))
  rownames(interactions) = colnames(interactions) = rownames(oddsRatio) = colnames(oddsRatio) = colnames(mutMat)
  sigPairs = which(x = 10^-abs(interactions) < 1, arr.ind = TRUE)
  sigPairs2 = which(x = 10^-abs(interactions) >= 1, arr.ind = TRUE)
  sigPairs = rbind(sigPairs, sigPairs2)
  sigPairsTbl = data.table::rbindlist(
    lapply(X = seq_along(1:nrow(sigPairs)), function(i) {
      x = sigPairs[i,]
      g1 = rownames(interactions[x[1], x[2], drop = FALSE])
      g2 = colnames(interactions[x[1], x[2], drop = FALSE])
      tbl = as.data.frame(table(factor(apply(X = mutMat[,c(g1, g2), drop = FALSE], 1, paste, collapse = ""), levels = c("00", "01","11", "10"))))
      combn = data.frame(t(tbl$Freq))
      colnames(combn) = tbl$Var1
      pval = 10^-abs(interactions[x[1], x[2]])
      fest = oddsRatio[x[1], x[2]]
      d = data.table::data.table(gene1 = g1, gene2 = g2, pValue = pval, oddsRatio = fest)
      d = cbind(d, combn)
      d
    }), fill = TRUE)
  sigPairsTbl[, pAdj := p.adjust(pValue, method = 'fdr')]
  sigPairsTbl[is.na(sigPairsTbl)] = 0
  sigPairsTbl$Event = ifelse(test = sigPairsTbl$oddsRatio > 1, yes = "Co_Occurrence", no = "Mutually_Exclusive")
  sigPairsTbl$pair = apply(X = sigPairsTbl[,.(gene1, gene2)], MARGIN = 1, FUN = function(x) paste(sort(unique(x)), collapse = ", "))
  sigPairsTbl[,event_ratio := `01`+`10`]
  sigPairsTbl[,event_ratio := paste0(`11`, '/', event_ratio)]
  sigPairsTblSig = sigPairsTbl[order(as.numeric(pValue))][!duplicated(pair)]
  sigPairsTblSig = sigPairsTblSig[!gene1 == gene2]
  if(nrow(interactions) >= 5 && showCounts && match.arg(countStats, c("all", "sig")) == "all"){
    # the previous plot loop looked every cell up with sigPairsTblSig[pair %in% g12, ...]
    # (lines 243-249), which leaves a data.table auto-index on the returned table
    invisible(sigPairsTblSig[pair %in% rownames(interactions)[1], event_ratio])
  }
  if(!returnAll){
    sigPairsTblSig = sigPairsTblSig[pValue < min(pvalue)]
  }
  sigPairsTblSig
}

# run f, returning its value and the warnings it raised, in order
withWarnings = function(expr){
  w = character()
  v = withCallingHandlers(expr, warning = function(cnd){
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w)
}

expectSameAsReference = function(maf, ...){
  pdf(NULL)
  on.exit(dev.off())
  new = withWarnings(somaticInteractions(maf = maf, ...))
  ref = withWarnings(referenceInteractions(maf = maf, ...))
  expect_identical(new$value, ref$value)
  expect_identical(new$warnings, ref$warnings)
  invisible(new$value)
}

laml.maf = system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml = read.maf(maf = laml.maf, verbose = FALSE)

# a MAF with one mutation per (gene, sample) pair of a 0/1 matrix
mafFromMatrix = function(m){
  idx = which(m == 1, arr.ind = TRUE)
  d = data.table::data.table(Hugo_Symbol = colnames(m)[idx[, 2]], Chromosome = "1",
                             Start_Position = 1000L * idx[, 2], End_Position = 1000L * idx[, 2],
                             Reference_Allele = "A", Tumor_Seq_Allele2 = "T",
                             Variant_Classification = "Missense_Mutation", Variant_Type = "SNP",
                             Tumor_Sample_Barcode = rownames(m)[idx[, 1]])
  read.maf(maf = d, verbose = FALSE)
}

test_that("somaticInteractions matches the previous implementation on TCGA LAML", {
  res = expectSameAsReference(laml, top = 25)
  expect_equal(nrow(res), choose(25, 2))
  expect_true(all(res$`00` + res$`01` + res$`10` + res$`11` == 193L))
  expectSameAsReference(laml, top = 10, returnAll = FALSE)
  expectSameAsReference(laml, genes = c("DNMT3A", "NPM1", "FLT3", "TP53", "IDH1", "IDH2"))
  expectSameAsReference(laml, genes = c("DNMT3A", "NPM1"))
  # plotPadj and showCounts only change the plot
  expectSameAsReference(laml, top = 20, plotPadj = TRUE)
  expectSameAsReference(laml, top = 20, showCounts = TRUE, countStats = "sig", countType = "cooccur")
  expectSameAsReference(laml, top = 20, showCounts = TRUE, countStats = "all", countType = "mutexcl")
  res = expectSameAsReference(laml, top = 25, showCounts = TRUE)
  expect_identical(data.table::indices(res), "pair")
  res = expectSameAsReference(laml, top = 25, showCounts = TRUE, countStats = "sig")
  expect_null(data.table::indices(res))
})

test_that("somaticInteractions matches the previous implementation on simulated data", {
  set.seed(1)
  nS = 400; nG = 40
  freq = 0.35 * seq_len(nG)^-0.45
  burden = rlnorm(nS, 0, 1); burden = burden / mean(burden)
  m = (matrix(runif(nS * nG), nS, nG) < pmin(1, outer(burden, freq))) * 1L
  m[m[, 3] == 1, 4] = 0L  # a mutually exclusive pair
  dimnames(m) = list(sprintf("S%04d", seq_len(nS)), sprintf("GENE%03d", seq_len(nG)))
  maf = mafFromMatrix(m)
  expectSameAsReference(maf, top = 30)
  expectSameAsReference(maf, top = 15, returnAll = FALSE)
})

test_that("somaticInteractions matches the previous implementation for degenerate genes", {
  set.seed(2)
  nS = 150; nG = 10
  m = (matrix(runif(nS * nG), nS, nG) < 0.15) * 1L
  dimnames(m) = list(sprintf("S%04d", seq_len(nS)), sprintf("GENE%03d", seq_len(nG)))
  m = cbind(m, ALLMUT = 1L, ALLMUT2 = 1L)   # mutated in every sample: no Fisher test possible
  m = cbind(m, TWIN001 = m[, "GENE001"])    # identical to another gene
  m = cbind(m, SOFTHYPHEN = m[, "GENE005"])
  colnames(m)[colnames(m) == "SOFTHYPHEN"] = "GENE\u00ad002" # distinct from GENE002 but may collate equal to it
  maf = mafFromMatrix(m)
  res = expectSameAsReference(maf, top = ncol(m))
  expect_false(any(c("ALLMUT", "ALLMUT2") %in% c(res$gene1, res$gene2)))
  expect_true(any(res$gene1 == "TWIN001" | res$gene2 == "TWIN001"))
})
