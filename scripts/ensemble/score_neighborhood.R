suppressPackageStartupMessages({ library(SynExtend); library(DECIPHER); library(glmnet) })

# The fair test. The combiner is fit on the eight benchmark systems exactly as in train_ensemble.R,
# then used to score DUF4229 against the 145 families that recur near MqnE. All of those share the
# closeness to MqnE that found DUF4229, so its rank here is not guaranteed by how it was picked.
#   Rscript _scripts/ensemble/score_neighborhood.R      # after build_neighborhood_labels.py

dir <- "_data/ensemble/"
FEAT <- c("GLMI", "GLDistance", "PAJaccard", "PAOverlap", "GeneDistance", "MoransI", "OrientationMI")

# EvoWeaver component scores for MqnE against every candidate, cached like component_scores.tsv
score_path <- paste0(dir, "neighborhood_component_scores.tsv")
cand <- read.delim(paste0(dir, "neighborhood_candidates.tsv"), stringsAsFactors = FALSE)
if (!file.exists(score_path)) {
  tree <- ReadDendrogram("_data/evoweaver/genomewide/taxonomy_tree.nwk")
  tree <- dendrapply(tree, function(n) { attr(n, "height") <- attr(n, "height") / attr(tree, "height"); n })
  lab <- read.delim(paste0(dir, "labels_futalosine_neighborhood.tsv"), stringsAsFactors = FALSE)
  fams <- split(lab$label, lab$family)
  ew <- EvoWeaver(fams, MySpeciesTree = tree)
  set.seed(635L)
  p <- predict(ew, Method = "Ensemble", Subset = data.frame(Gene1 = "ANCHOR", Gene2 = intersect(cand$candidate, names(fams))),
               MySpeciesTree = tree, PretrainedModel = "KEGG", NoPrediction = TRUE, Verbose = FALSE)
  p$candidate <- ifelse(p$Gene1 == "ANCHOR", p$Gene2, p$Gene1)
  write.table(p, score_path, sep = "\t", quote = FALSE, row.names = FALSE)
}
s <- merge(read.delim(score_path, stringsAsFactors = FALSE), cand, by = "candidate")
s[FEAT] <- lapply(s[FEAT], function(x) { x[is.na(x)] <- 0; x })

# the combiner, fit on the benchmark the same way as the final model in train_ensemble.R
d <- read.delim(paste0(dir, "component_scores.tsv"), stringsAsFactors = FALSE)
d <- merge(d[setdiff(names(d), "is_partner")], read.delim(paste0(dir, "pair_labels.tsv"), stringsAsFactors = FALSE),
           by = c("system", "candidate"))
d[FEAT] <- lapply(d[FEAT], function(x) { x[is.na(x)] <- 0; x })
train <- d[d$system != "futalosine", ]
set.seed(1L)
cv <- cv.glmnet(as.matrix(train[FEAT]), train$is_partner, family = "binomial", alpha = 0, nfolds = 5)
final <- glmnet(as.matrix(train[FEAT]), train$is_partner, family = "binomial", alpha = 0, lambda = cv$lambda.min)

s$model <- as.numeric(predict(final, as.matrix(s[FEAT]), type = "response"))
s <- s[order(-s$model), ]
s$rank <- seq_len(nrow(s))
write.table(s, paste0(dir, "neighborhood_scored.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

nb <- s[s$kind != "known", ]           # DUF4229 and the neighborhood families, without the two known genes
duf <- which(nb$kind == "DUF4229")
cat(sprintf("%d candidates scored: %d neighborhood families, DUF4229, and the 2 known futalosine genes\n",
            nrow(s), sum(s$kind == "neighborhood")))
cat(sprintf("trained model: DUF4229 ranks %d of %d among DUF4229 plus the neighborhood families (p = %.4f)\n",
            duf, nrow(nb), duf / nrow(nb)))
cat("\ntop of the full list:\n")
for (i in seq_len(min(10, nrow(s)))) cat(sprintf("  %2d. %-14s %-28s model %.3f\n", i, s$candidate[i], s$name[i], s$model[i]))
cat("\nDUF4229 rank among DUF4229 plus the neighborhood families, score by score:\n")
for (f in FEAT) cat(sprintf("  %-14s %d\n", f, rank(-nb[[f]], ties.method = "min")[duf]))
