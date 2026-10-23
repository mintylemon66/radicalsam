# EvoWeaver on the genome-wide set: how DUF4229's pairing with MqnE compares with two null sets.
#   Rscript _scripts/evoweaver/genomewide/run_genomewide.R      # after build_labels.py
#
# Step 8, the last one. Each family is scored against MqnE twice: on presence and absence alone, and
# on presence and absence plus gene organization, using the built-in models the EvoWeaver authors
# trained on KEGG modules. Both rank DUF4229 against the random Pfam families and against the
# families that recur near MqnE, and report an empirical p from its rank in each.
#
# The Ensemble probability saturates at 1, so two genuinely different pairs can both read 1.000. The
# model's score before the logistic step is kept for that reason, along with an unweighted sum of the
# component scores, and a sum without MoransI, which returns exactly 1 in a degenerate case.
#
# Writes ew_scores_profiles.tsv and ew_scores_full.tsv.

suppressPackageStartupMessages({ library(SynExtend); library(DECIPHER) })

dir <- "_data/evoweaver/genomewide/"
tree <- ReadDendrogram(paste0(dir, "taxonomy_tree.nwk"))
# the EvoWeaver authors rescaled their reference trees to a root height of 1
root <- attr(tree, "height")
tree <- dendrapply(tree, function(n) { attr(n, "height") <- attr(n, "height") / root; n })

lab <- read.delim(paste0(dir, "labels.tsv"), stringsAsFactors = FALSE)
fams <- split(lab$label, lab$family)
others <- setdiff(names(fams), "MqnE")
# a two-column table of pairs; in SynExtend 1.24.0 a plain vector of names scores the wrong pairs
pairs <- data.frame(Gene1 = "MqnE", Gene2 = others)
models <- get(data("BuiltInEnsembles", package = "SynExtend"))$KEGG

score <- function(input, model) {
  ew <- EvoWeaver(input, MySpeciesTree = tree)
  set.seed(635L)   # the seed the authors used; the profiling p-values come from permutations
  p <- predict(ew, Method = "Ensemble", Subset = pairs, MySpeciesTree = tree,
               PretrainedModel = "KEGG", Verbose = FALSE)
  p$partner <- ifelse(p$Gene1 == "MqnE", p$Gene2, p$Gene1)
  stopifnot(all(p$Gene1 == "MqnE" | p$Gene2 == "MqnE"), setequal(p$partner, others))
  # the Ensemble probability saturates at 1, so also keep the model's score before the logistic step
  p$EnsembleLogit <- predict(model, p, type = "link")
  stopifnot(isTRUE(all.equal(plogis(p$EnsembleLogit), p$Ensemble, check.attributes = FALSE)))   # same model predict() used
  comp <- setdiff(names(p), c("Gene1", "Gene2", "partner", "Ensemble", "EnsembleLogit"))
  p$ComponentSum <- rowSums(p[comp])                            # unweighted sum, our choice
  p$SumWithoutMoransI <- rowSums(p[setdiff(comp, "MoransI")])   # MoransI returns exactly 1 in a degenerate case
  p
}

# 1. presence/absence only: genome IDs without coordinates (built-in model 1: profiling only)
profiles <- score(lapply(fams, function(x) unique(sub("_.*", "", x))), models[[1]])
# 2. presence/absence plus gene organization: genome_contig_strand_index labels (model 3: profiling + organization)
full <- score(fams, models[[3]])

report <- function(p, title, cols) {
  cat("\n", title, "\n", sep = "")
  duf <- p[p$partner == "DUF4229", ]
  for (g in c("MqnC", "MqnAD", "Fer4")) {
    cat(sprintf("  control MqnE-%-6s Ensemble %.3f, component sum %.3f\n", g,
                p$Ensemble[p$partner == g], p$ComponentSum[p$partner == g]))
  }
  for (set in c("null_", "nbr_")) {
    grp <- p[startsWith(p$partner, set), ]
    cat(sprintf("  DUF4229 against %d %s families\n", nrow(grp),
                if (set == "null_") "random Pfam" else "MqnE-neighborhood"))
    for (col in cols) {
      # rounded so that two scores at the model's ceiling count as a tie
      v <- round(grp[[col]], 6)
      d <- round(duf[[col]], 6)
      k <- sum(v >= d)
      cat(sprintf("    %-18s DUF4229 %8.3f | median %8.3f, max %8.3f | %3d at or above, p = %.4f\n",
                  col, d, median(v), max(v), k, (1 + k) / (1 + length(v))))
    }
  }
}
report(profiles, "EvoWeaver, phylogenetic profiles only", c("Ensemble", "EnsembleLogit", "ComponentSum"))
report(full, "EvoWeaver, phylogenetic profiles + gene organization",
       c("Ensemble", "EnsembleLogit", "ComponentSum", "SumWithoutMoransI", "GeneDistance", "OrientationMI"))

write.table(profiles, paste0(dir, "ew_scores_profiles.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(full, paste0(dir, "ew_scores_full.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
