suppressPackageStartupMessages({ library(SynExtend); library(DECIPHER); library(glmnet) })

# EvoWeaver supplies the component scores. The model that combines them is fitted here on radical SAM
# accessory proteins, in place of the built-in one trained on KEGG metabolic modules. Each system is held
# out in turn so every prediction comes from a model that never saw it, and futalosine is held out of all
# of it, so scoring MqnE against DUF4229 at the end is a transfer result.
#
# The fit is ridge-penalized. With 28 positives and seven features an unpenalized logistic regression
# separates the classes completely, which sends the coefficients to infinity and makes the fit unusable.

dir <- "_data/ensemble/"
FEAT <- c("GLMI", "GLDistance", "PAJaccard", "PAOverlap", "GeneDistance", "MoransI", "OrientationMI")

score_path <- paste0(dir, "component_scores.tsv")
if (!file.exists(score_path)) {
  tree <- ReadDendrogram("_data/evoweaver/genomewide/taxonomy_tree.nwk")
  tree <- dendrapply(tree, function(n) { attr(n, "height") <- attr(n, "height") / attr(tree, "height"); n })
  pairs <- read.delim(paste0(dir, "pair_labels.tsv"), stringsAsFactors = FALSE)
  scores <- list()
  for (s in unique(pairs$system)) {
    lab <- read.delim(sprintf("%slabels_%s.tsv", dir, s), stringsAsFactors = FALSE)
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
  d <- merge(do.call(rbind, scores), pairs, by = c("system", "candidate"))
  d[FEAT] <- lapply(d[FEAT], function(x) { x[is.na(x)] <- 0; x })
  write.table(d, score_path, sep = "\t", quote = FALSE, row.names = FALSE)
}
d <- read.delim(score_path, stringsAsFactors = FALSE)
# labels are read fresh from pair_labels.tsv on every run, so fixing a label never needs a rescore;
# pairs no longer in the table drop out here
d <- merge(d[setdiff(names(d), "is_partner")], read.delim(paste0(dir, "pair_labels.tsv"), stringsAsFactors = FALSE),
           by = c("system", "candidate"))
d[FEAT] <- lapply(d[FEAT], function(x) { x[is.na(x)] <- 0; x })

train <- d[d$system != "futalosine", ]
held <- d[d$system == "futalosine", ]
cat(sprintf("\n%d pairs scored. Training on %d pairs, %d of them partners, across %d systems.\n",
            nrow(d), nrow(train), sum(train$is_partner), length(unique(train$system))))

auroc <- function(s, y) {
  r <- rank(s); np <- sum(y == 1); nn <- sum(y == 0)
  (sum(r[y == 1]) - np * (np + 1) / 2) / (np * nn)
}
fit_ridge <- function(df) {
  x <- as.matrix(df[FEAT]); y <- df$is_partner
  set.seed(1L)
  cv <- cv.glmnet(x, y, family = "binomial", alpha = 0, nfolds = 5)
  glmnet(x, y, family = "binomial", alpha = 0, lambda = cv$lambda.min)
}

# leave one system out
train$loso <- NA_real_
for (s in unique(train$system)) {
  m <- fit_ridge(train[train$system != s, ])
  train$loso[train$system == s] <- predict(m, as.matrix(train[train$system == s, FEAT]), type = "response")
}
cat(sprintf("\nleave-one-system-out AUROC: %.3f\n", auroc(train$loso, train$is_partner)))

# the built-in KEGG ensemble on the same pairs, as the baseline this has to beat
kegg <- get(data("BuiltInEnsembles", package = "SynExtend"))$KEGG
train$kegg <- predict(kegg[[3]], train, type = "response")
cat(sprintf("built-in KEGG ensemble, same pairs:  %.3f\n", auroc(train$kegg, train$is_partner)))

cat("\nper system, where the known partners rank among that system's candidates:\n")
for (s in unique(train$system)) {
  i <- train$system == s
  ranks <- rank(-train$loso[i])[train$is_partner[i] == 1]
  kranks <- rank(-train$kegg[i])[train$is_partner[i] == 1]
  cat(sprintf("  %-14s %d of %-3d candidates | ours %-14s | KEGG %s\n", s, sum(i), sum(i),
              paste(sort(round(ranks)), collapse = ","), paste(sort(round(kranks)), collapse = ",")))
}

final <- fit_ridge(train)
co <- as.matrix(coef(final))
cat("\ncoefficients of the final model:\n")
for (n in rownames(co)) cat(sprintf("  %-16s %8.3f\n", n, co[n, 1]))

held$ours <- predict(final, as.matrix(held[FEAT]), type = "response")
held$kegg <- predict(kegg[[3]], held, type = "response")
held <- held[order(-held$ours), ]
write.table(held, paste0(dir, "futalosine_scored.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nheld-out futalosine system, scored by a model that never saw it:\n")
for (i in seq_len(nrow(held))) {
  tag <- if (held$candidate[i] == "fam_PF14012") "  <- DUF4229, the candidate" else
         if (held$is_partner[i] == 1) "  (known pathway gene)" else ""
  cat(sprintf("  %2d. %-14s ours %.3f   KEGG %.3f%s\n", i, held$candidate[i], held$ours[i], held$kegg[i], tag))
}
duf <- which(held$candidate == "fam_PF14012")
if (length(duf)) cat(sprintf("\nDUF4229 ranks %d of %d in the held-out system.\n", duf, nrow(held)))
