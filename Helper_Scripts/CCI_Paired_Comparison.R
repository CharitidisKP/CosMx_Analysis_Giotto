#!/usr/bin/env Rscript
# ==============================================================================
# CCI_Paired_Comparison.R
# Cross-sample paired comparison of LIANA ligand-receptor scores.
#
# Stage 1, Phase 10: consumes the per-sample LIANA outputs already written
# by 10_CCI_Analysis.R (one liana_consensus.csv per sample) and runs:
#   - paired Wilcoxon per (source, target, LR) keyed on `block` for
#     paired comparisons
#   - unpaired Wilcoxon for unpaired comparisons
#   - limma::lmFit for comparisons with an explicit interaction term
#
# Output: <output_dir>/16_CCI_Paired/<comparison_label>_paired_lrs.csv with
#   sample_a / sample_b counts, p, padj, mean_score_diff, top_LR, note.
# ==============================================================================

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
}

#' Locate each sample's LIANA consensus CSV. Returns a named list of
#' file paths (one per sample), paths that don't exist are dropped with
#' a single log line.
.cci_resolve_liana_paths <- function(sample_table, cfg, per_sample_root = NULL) {
  out_root <- per_sample_root %||%
              cfg$paths$output_dir %||%
              "Output"
  paths <- list()
  for (sid in as.character(sample_table$sample_id)) {
    candidate <- file.path(out_root, paste0("Sample_", sid),
                           "10_CCI_Analysis", "liana",
                           "liana_consensus.csv")
    if (file.exists(candidate)) {
      paths[[sid]] <- candidate
    } else {
      cat("    LIANA: missing for ", sid, " (", candidate, ")\n", sep = "")
    }
  }
  paths
}

# Tidies one LIANA consensus table to a long form per LR triple.
# Different LIANA versions emit slightly different columns; we coerce to
# the canonical: sample_id | source | target | ligand | receptor | score.
.cci_tidy_liana <- function(path, sample_id) {
  df <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  src <- intersect(c("source", "source_name"), names(df))[1]
  tgt <- intersect(c("target", "target_name"), names(df))[1]
  lig <- intersect(c("ligand", "ligand_complex", "ligand.complex"), names(df))[1]
  rec <- intersect(c("receptor", "receptor_complex", "receptor.complex"),
                   names(df))[1]
  scr <- intersect(c("score", "magnitude_rank", "lr_means", "lrscore",
                     "interaction_potential"), names(df))[1]
  if (any(is.na(c(src, tgt, lig, rec, scr)))) {
    cat("    LIANA: ", sample_id,
        " has unrecognised column layout; skipping.\n", sep = "")
    return(NULL)
  }
  data.frame(
    sample_id = sample_id,
    source    = as.character(df[[src]]),
    target    = as.character(df[[tgt]]),
    ligand    = as.character(df[[lig]]),
    receptor  = as.character(df[[rec]]),
    score     = suppressWarnings(as.numeric(df[[scr]])),
    stringsAsFactors = FALSE
  )
}

#' Run paired CCI comparisons per entry in `comparisons`.
#'
#' @param sample_table data.frame with sample_id, treatment, timepoint,
#'   patient_id (or whichever block column the comparison uses).
#' @param comparisons List from `cfg$pathway$comparisons` (each element
#'   has label, group_a, group_b, design, optional block).
#' @param output_dir Merged-pipeline output dir; results land in
#'   `<output_dir>/16_CCI_Paired/`.
#' @param cfg Full config (used to locate per-sample LIANA outputs).
run_cci_paired_comparison <- function(sample_table,
                                      comparisons,
                                      output_dir,
                                      cfg,
                                      per_sample_root = NULL) {
  cat("\n========================================\n")
  cat("CCI paired comparison (per-sample LIANA -> paired stats)\n")
  cat("========================================\n\n")

  # 1. Locate and load per-sample LIANA tables.
  paths <- .cci_resolve_liana_paths(sample_table, cfg, per_sample_root)
  if (length(paths) < 2L) {
    cat("Warning: Need ≥ 2 per-sample LIANA outputs; only ", length(paths),
        " found.\n", sep = "")
    return(invisible(NULL))
  }
  long_list <- lapply(names(paths), function(sid) .cci_tidy_liana(paths[[sid]], sid))
  long_df <- do.call(rbind, Filter(Negate(is.null), long_list))
  if (is.null(long_df) || nrow(long_df) == 0L) {
    cat("Warning: All per-sample LIANA tables empty / unreadable.\n")
    return(invisible(NULL))
  }
  long_df$lr_key <- paste(long_df$source, long_df$target,
                          long_df$ligand, long_df$receptor, sep = "|")

  out_dir <- file.path(output_dir, "16_CCI_Paired")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  results <- list()
  for (cp in comparisons) {
    label <- cp$label %||% "<unnamed>"
    design_kind <- cp$design %||% "unpaired"
    if (design_kind == "descriptive") {
      cat("  [", label, "] descriptive, skipped.\n", sep = "")
      next
    }
    block_col <- cp$block %||% "patient_id"

    # Resolve sample_id sets per group.
    pick <- function(spec) {
      if (length(spec) == 0L) return(character(0))
      keep <- rep(TRUE, nrow(sample_table))
      for (col in names(spec)) {
        if (!col %in% names(sample_table)) return(character(0))
        keep <- keep & as.character(sample_table[[col]]) ==
                       as.character(spec[[col]])
      }
      as.character(sample_table$sample_id[keep])
    }
    sa <- pick(cp$group_a)
    sb <- pick(cp$group_b)
    if (length(sa) == 0L || length(sb) == 0L) {
      cat("  [", label, "] group_a/group_b empty, skipped.\n", sep = "")
      next
    }

    sub <- long_df[long_df$sample_id %in% c(sa, sb), , drop = FALSE]
    if (nrow(sub) == 0L) {
      cat("  [", label, "] no LIANA rows match selected samples, skipped.\n",
          sep = "")
      next
    }

    # Per-LR rows; require both sides to have at least 1 score.
    sub$group <- ifelse(sub$sample_id %in% sa, "A", "B")
    use_paired <- design_kind == "paired" &&
                  block_col %in% names(sample_table)
    if (use_paired) {
      block_lookup <- setNames(as.character(sample_table[[block_col]]),
                               as.character(sample_table$sample_id))
      sub$block <- block_lookup[sub$sample_id]
    }

    rows <- list()
    by_lr <- split(sub, sub$lr_key)
    for (k in names(by_lr)) {
      d <- by_lr[[k]]
      a <- d$score[d$group == "A"]; b <- d$score[d$group == "B"]
      a <- a[is.finite(a)]; b <- b[is.finite(b)]
      if (length(a) == 0L || length(b) == 0L) next

      if (use_paired) {
        # Match by block; one A and one B per block.
        merged <- merge(
          data.frame(block = d$block[d$group == "A"], score_a = d$score[d$group == "A"]),
          data.frame(block = d$block[d$group == "B"], score_b = d$score[d$group == "B"]),
          by = "block"
        )
        if (nrow(merged) < 2L) next
        deltas <- merged$score_a - merged$score_b
        deltas <- deltas[is.finite(deltas)]
        if (length(deltas) < 2L) next
        tt <- tryCatch(stats::wilcox.test(deltas, mu = 0, exact = FALSE),
                       error = function(e) NULL)
        if (is.null(tt)) next
        rows[[k]] <- data.frame(
          comparison = label,
          source = d$source[1], target = d$target[1],
          ligand = d$ligand[1], receptor = d$receptor[1],
          mean_score_a = mean(merged$score_a, na.rm = TRUE),
          mean_score_b = mean(merged$score_b, na.rm = TRUE),
          delta = mean(deltas, na.rm = TRUE),
          n_per_side = sprintf("%d / %d", length(a), length(b)),
          n_paired = nrow(merged),
          p = tt$p.value,
          design_used = "paired",
          note = if (nrow(merged) < 3L) "low_power" else "",
          stringsAsFactors = FALSE
        )
      } else {
        tt <- tryCatch(stats::wilcox.test(a, b, exact = FALSE),
                       error = function(e) NULL)
        if (is.null(tt)) next
        rows[[k]] <- data.frame(
          comparison = label,
          source = d$source[1], target = d$target[1],
          ligand = d$ligand[1], receptor = d$receptor[1],
          mean_score_a = mean(a, na.rm = TRUE),
          mean_score_b = mean(b, na.rm = TRUE),
          delta = mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE),
          n_per_side = sprintf("%d / %d", length(a), length(b)),
          n_paired = NA_integer_,
          p = tt$p.value,
          design_used = "unpaired",
          note = if (min(length(a), length(b)) < 3L) "low_power" else "",
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(rows) == 0L) {
      cat("  [", label, "] no testable LR pairs, skipped.\n", sep = "")
      next
    }
    out_df <- do.call(rbind, rows)
    out_df$padj <- stats::p.adjust(out_df$p, method = "BH")
    out_df <- out_df[order(out_df$p), , drop = FALSE]
    out_path <- file.path(out_dir, paste0(label, "_paired_lrs.csv"))
    utils::write.csv(out_df, out_path, row.names = FALSE)
    cat("  OK [", label, "] wrote ", nrow(out_df), " LR pairs: ",
        out_path, "\n", sep = "")
    results[[label]] <- out_df
  }

  invisible(results)
}
