suppressPackageStartupMessages({ library(SynExtend); library(DECIPHER); library(glmnet) })

# The retrained combiner. Same seven EvoWeaver scores and the same ridge fit as train_ensemble.R, but
# the negatives are each system's neighborhood families, the comparison DUF4229 faces near MqnE.
# Judged leave-one-system-out on the eight training systems only. Futalosine is scored once, at the
# end, against its 145 neighborhood families, and the rank is reported whatever it is.
#   Rscript _scripts/ensemble/train_neighbor_model.R      # after build_neighbor_pairs.py

dir <- "_data/ensemble/"
FEAT <- c("GLMI", "GLDistance", "PAJaccard", "PAOverlap", "GeneDistance", "MoransI", "OrientationMI")

# EvoWeaver component scores for every anchor and candidate pair, cached like component_scores.tsv
pairs <- read.delim(paste0(dir, "pair_labels_nb.tsv"), stringsAsFactors = FALSE)
score_path <- paste0(dir, "component_scores_nb.tsv")
if (!file.exists(score_path)) {
  tree <- ReadDendrogram("_data/evoweaver/genomewide/taxonomy_tree.nwk")
  tree <- dendrapply(tree, function(n) { attr(n, "height") <- attr(n, "height") / attr(tree, "height"); n })
  scores <- list()
  for (s in unique(pairs$system)) {
    lab <- read.delim(sprintf("%slabels_nb_%s.tsv", dir, s), stringsAsFactors = FALSE)
    fams <- split(lab$label, lab$family)
    cand <- intersect(pairs$candidate[pairs$system == s], names(fams))
    ew <- EvoWeaver(fams, MySpeciesTree = tree)
    set.seed(635L)
    p <- predict(ew, Method = "Ensemble", Subset = data.frame(Gene1 = "ANCHOR", Gene2 = cand),
                 MySpeciesTree = tree, PretrainedModel = "KEGG", NoPrediction = TRUE, Verbose = FALSE)
    p$system <- s
    p$candidate <- ifelse(p$Gene1 == "ANCHOR", p$Gene2, p$Gene1)
    scores[[s]] <- p
    cat(sprintf("%-14s scored %d candidates\n", s, nrow(p)))
  }
  write.table(do.call(rbind, scores), score_path, sep = "\t", quote = FALSE, row.names = FALSE)
}
d <- merge(read.delim(score_path, stringsAsFactors = FALSE)[c("system", "candidate", FEAT)], pairs, by = c("system", "candidate"))
d[FEAT] <- lapply(d[FEAT], function(x) { x[is.na(x)] <- 0; x })
cat(sprintf("\n%d pairs, %d of them partners, across %d systems\n", nrow(d), sum(d$is_partner), length(unique(d$system))))

auroc <- function(s, y) {
  r <- rank(s); np <- sum(y == 1); nn <- sum(y == 0)
  (sum(r[y == 1]) - np * (np + 1) / 2) / (np * nn)
}
within <- function(df, col) mean(sapply(split(df, df$system), function(x) auroc(x[[col]], x$is_partner)))
fit_ridge <- function(df) {
  x <- as.matrix(df[FEAT]); y <- df$is_partner
  set.seed(1L)
  cv <- cv.glmnet(x, y, family = "binomial", alpha = 0, nfolds = 5)
  glmnet(x, y, family = "binomial", alpha = 0, lambda = cv$lambda.min)
}

# leave one system out
d$loso <- NA_real_
for (s in unique(d$system)) {
  m <- fit_ridge(d[d$system != s, ])
  d$loso[d$system == s] <- as.numeric(predict(m, as.matrix(d[d$system == s, FEAT]), type = "response"))
}

# the first combiner, trained against random families, on the same neighborhood pairs
old <- read.delim(paste0(dir, "component_scores.tsv"), stringsAsFactors = FALSE)
old <- merge(old[setdiff(names(old), "is_partner")], read.delim(paste0(dir, "pair_labels.tsv"), stringsAsFactors = FALSE),
             by = c("system", "candidate"))
old[FEAT] <- lapply(old[FEAT], function(x) { x[is.na(x)] <- 0; x })
d$old <- as.numeric(predict(fit_ridge(old[old$system != "futalosine", ]), as.matrix(d[FEAT]), type = "response"))

cat(sprintf("\nneighborhood benchmark, AUROC pooled / mean within system\n"))
cat(sprintf("  retrained, leave one system out: %.3f / %.3f\n", auroc(d$loso, d$is_partner), within(d, "loso")))
cat(sprintf("  first combiner (random negatives): %.3f / %.3f\n", auroc(d$old, d$is_partner), within(d, "old")))
cat("\nper system, where the known partners rank among that system's candidates (retrained | first):\n")
for (s in unique(d$system)) {
  i <- d$system == s
  r_new <- sort(round(rank(-d$loso[i])[d$is_partner[i] == 1]))
  r_old <- sort(round(rank(-d$old[i])[d$is_partner[i] == 1]))
  cat(sprintf("  %-14s %3d candidates | %-18s | %s\n", s, sum(i), paste(r_new, collapse = ","), paste(r_old, collapse = ",")))
}

final <- fit_ridge(d)
co <- as.matrix(coef(final))
cat("\ncoefficients of the retrained model:\n")
for (n in rownames(co)) cat(sprintf("  %-16s %8.3f\n", n, co[n, 1]))

# applied once to futalosine: DUF4229 against MqnE's 145 neighborhood families
s <- merge(read.delim(paste0(dir, "neighborhood_component_scores.tsv"), stringsAsFactors = FALSE)[c("candidate", FEAT)],
           read.delim(paste0(dir, "neighborhood_candidates.tsv"), stringsAsFactors = FALSE), by = "candidate")
s[FEAT] <- lapply(s[FEAT], function(x) { x[is.na(x)] <- 0; x })
s$model <- as.numeric(predict(final, as.matrix(s[FEAT]), type = "response"))
s <- s[order(-s$model), ]
write.table(s, paste0(dir, "neighborhood_scored_nb.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
nb <- s[s$kind != "known", ]
duf <- which(nb$kind == "DUF4229")
cat(sprintf("\nfutalosine, held out: DUF4229 ranks %d of %d among DUF4229 plus the neighborhood families (p = %.4f)\n",
            duf, nrow(nb), duf / nrow(nb)))
for (k in which(s$kind == "known")) cat(sprintf("  known gene %-10s ranks %d of %d overall\n", s$name[k], k, nrow(s)))
cat("top of the list:\n")
for (i in seq_len(min(10, nrow(s)))) cat(sprintf("  %2d. %-14s %-28s %.3f\n", i, s$candidate[i], s$name[i], s$model[i]))
