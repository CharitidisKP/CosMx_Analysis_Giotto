#!/usr/bin/env Rscript
# ==============================================================================
# 13_Pathway_Analysis.R
# Pathway enrichment (ORA via clusterProfiler::enricher) + GSEA (fgsea) across
# four biologically meaningful cross-sample comparisons, stratified per cell
# type (celltype_*_supervised) and per Leiden cluster. Adds PROGENy pathway
# activity overlays focused on B cells and a multi-page PDF summary.
#
# Optional dependencies are installed via Parameters/Install_CCI_dependencies.R
# (see the `pathway` section comment in Parameters/config.yaml for the list).
#
# INPUT:  merged Giotto object produced by merge_giotto_samples / merge_batch
#         (must contain `treatment`, `timepoint`, a celltype_*_supervised
#         column, and `leiden_clust` in pDataDT).
# OUTPUT: Output/merged/13_Pathway_Analysis/
#           <comparison>/<stratum>/ ... per-stratum CSVs + P1-P4
#           <comparison>/_cross_*.png                ... P5/P6
#           _cross_comparison_*.png                   ... P7/P8
#           <comparison>/leading_edge/                ... P9
#           bcell_focus/                              ... P10-P12
#           _all_results_merged.csv                   ... tidied merged table
# ==============================================================================

current_script_dir_13 <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]),
                                 winslash = "/", mustWork = FALSE)))
  }
  ofiles <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else frame$ofile
  }, character(1))
  ofiles <- ofiles[nzchar(ofiles)]
  if (length(ofiles) > 0) {
    return(dirname(normalizePath(tail(ofiles, 1),
                                 winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

# Source Pipeline_Utils.R (presentation_theme, save_presentation_plot,
# sample_plot_title) if not already loaded.
.pathway_pipeline_utils <- file.path(current_script_dir_13(), "Helper_Scripts",
                                     "Pipeline_Utils.R")
if ((!exists("presentation_theme") ||
     !exists("save_presentation_plot") ||
     !exists("sample_plot_title")) &&
    file.exists(.pathway_pipeline_utils)) {
  source(.pathway_pipeline_utils)
}

# ==============================================================================
# Giotto accessors — mirrored from step 10 to avoid cross-script sourcing and
# to sidestep Giotto's match.call()[[1]] dispatch pitfall (see cerebrum
# 2026-04-19): always call Giotto's namespaced function directly, never via
# a variable-held reference.
# ==============================================================================

.pathway_pdata_dt <- function(gobj) {
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("pDataDT", envir = asNamespace("Giotto"), inherits = FALSE)) {
    return(Giotto::pDataDT(gobj))
  }
  if (requireNamespace("GiottoClass", quietly = TRUE) &&
      exists("pDataDT", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    return(GiottoClass::pDataDT(gobj))
  }
  get("pDataDT", mode = "function")(gobj)
}

.pathway_get_expression <- function(gobj, values = "normalized",
                                    output = "matrix") {
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("getExpression", envir = asNamespace("Giotto"), inherits = FALSE)) {
    return(Giotto::getExpression(gobj, values = values, output = output))
  }
  if (requireNamespace("GiottoClass", quietly = TRUE) &&
      exists("getExpression", envir = asNamespace("GiottoClass"),
             inherits = FALSE)) {
    return(GiottoClass::getExpression(gobj, values = values, output = output))
  }
  get("getExpression", mode = "function")(gobj, values = values,
                                          output = output)
}

# ==============================================================================
# Polygon overlay helpers — adapted from 10_CCI_Analysis.R so step 13 is
# self-contained. Reuses the safe grouping key (cell_ID + geom + part) per
# cerebrum bug-065.
# ==============================================================================

.pathway_extract_polygon_df <- function(gobj) {
  poly_sv <- tryCatch({
    p <- GiottoClass::getPolygonInfo(gobject = gobj, polygon_name = "cell",
                                     return_giottoPolygon = FALSE)
    if (inherits(p, "giottoPolygon")) p@spatVector else p
  }, error = function(e) {
    tryCatch({
      gp <- GiottoClass::getPolygonInfo(gobject = gobj, polygon_name = "cell")
      if (inherits(gp, "giottoPolygon")) gp@spatVector else NULL
    }, error = function(e2) NULL)
  })
  if (is.null(poly_sv)) {
    poly_sv <- tryCatch(gobj@polygon$cell@spatVector, error = function(e) NULL)
  }
  if (is.null(poly_sv)) return(NULL)

  poly_attr <- tryCatch(as.data.frame(poly_sv), error = function(e) NULL)
  if (is.null(poly_attr) || nrow(poly_attr) == 0) return(NULL)
  poly_coords <- tryCatch(terra::geom(poly_sv, df = TRUE), error = function(e) NULL)
  if (is.null(poly_coords) || nrow(poly_coords) == 0) return(NULL)

  id_col <- intersect(c("poly_ID", "cell_ID", "id"), names(poly_attr))
  if (length(id_col) == 0) return(NULL)
  poly_coords$cell_ID <- poly_attr[[id_col[1]]][poly_coords$geom]
  poly_coords$poly_group <- paste(poly_coords$cell_ID, poly_coords$geom,
                                  poly_coords$part, sep = "_")
  poly_coords
}

.pathway_resolve_bcell_ids <- function(gobj, celltype_col,
                                       focus_celltype = "^B cell$") {
  if (is.null(celltype_col) || !nzchar(celltype_col)) return(character(0))
  meta <- as.data.frame(.pathway_pdata_dt(gobj))
  if (!celltype_col %in% names(meta)) return(character(0))
  hit <- grepl(focus_celltype, as.character(meta[[celltype_col]]),
               ignore.case = TRUE)
  unique(as.character(meta$cell_ID[hit]))
}

# ==============================================================================
# Pathway-name cleaning — produces "<collection>[:subcollection] - <Description>"
# with underscores/slashes/periods stripped, title-cased on the description
# half only (collection + sub-collection codes stay uppercase).
# ==============================================================================

.pathway_title_case <- function(x) {
  x <- tolower(as.character(x))
  parts <- strsplit(x, " ", fixed = TRUE)
  out <- vapply(parts, function(words) {
    words <- words[nzchar(words)]
    if (length(words) == 0) return("")
    keep_upper <- c("dna", "rna", "mrna", "ncrna", "nk", "tca", "tgf", "tnf",
                    "jak", "stat", "nfkb", "nf", "mhc", "gpcr", "cdk", "er",
                    "atp", "gtp", "mtor", "pi3k", "map", "mapk")
    cap <- ifelse(words %in% keep_upper,
                  toupper(words),
                  paste0(toupper(substr(words, 1, 1)),
                         substr(words, 2, nchar(words))))
    paste(cap, collapse = " ")
  }, character(1))
  out
}

.pathway_clean_name <- function(id, collection = NA, subcollection = NA,
                                description = NA) {
  id <- as.character(id)
  description <- ifelse(is.na(description) | !nzchar(description), id,
                        as.character(description))
  collection <- ifelse(is.na(collection) | !nzchar(collection), "",
                       as.character(collection))
  subcollection <- ifelse(is.na(subcollection) | !nzchar(subcollection), "",
                          as.character(subcollection))
  # Strip MSigDB collection prefix from the description if still present.
  # Common prefixes: HALLMARK_, KEGG_, REACTOME_, WP_, BIOCARTA_, PID_,
  # GOBP_, GOMF_, GOCC_, GO_, C7_, IMMUNE_, etc. Only the first one matches.
  desc_clean <- gsub("[_/\\.]+", " ", description)
  desc_clean <- gsub("^(HALLMARK|KEGG|REACTOME|WP|WIKIPATHWAYS|BIOCARTA|PID|GOBP|GOMF|GOCC|GO)\\s+", "",
                     desc_clean, ignore.case = TRUE)
  desc_clean <- gsub("\\s+", " ", trimws(desc_clean))
  desc_clean <- .pathway_title_case(desc_clean)

  prefix <- collection
  if (nzchar(subcollection)) prefix <- paste0(prefix, ":", subcollection)
  if (!nzchar(prefix)) return(desc_clean)
  paste0(prefix, " - ", desc_clean)
}

# ==============================================================================
# MSigDB gene-set loader — returns per-database "by_id" list (pathway_id →
# gene SYMBOL vector) plus a $meta data frame with collection / subcollection
# / description for each pathway ID. Uses the msigdbr CRAN package.
# ==============================================================================

.pathway_load_genesets <- function(cfg_pathway) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("msigdbr is required. Install via Rscript Parameters/Install_CCI_dependencies.R msigdbr")
  }
  species <- cfg_pathway$species %||% "Homo sapiens"
  if (is.null(`%||%`)) {
    `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
  }

  dbs <- cfg_pathway$databases
  if (is.null(dbs) || length(dbs) == 0) {
    stop("pathway.databases must list at least one MSigDB collection")
  }

  out <- list()
  for (spec in dbs) {
    coll <- spec$collection
    subs <- spec$subcollections
    if (is.null(subs) || length(subs) == 0) {
      subs <- list(NULL)   # pull the whole collection
    }
    for (sub in subs) {
      df <- tryCatch(
        if (is.null(sub)) {
          msigdbr::msigdbr(species = species, category = coll)
        } else {
          msigdbr::msigdbr(species = species, category = coll,
                           subcategory = sub)
        },
        error = function(e) {
          cat("⚠ msigdbr(", coll,
              if (!is.null(sub)) paste0(",", sub) else "",
              ") failed: ", conditionMessage(e), "\n", sep = "")
          NULL
        }
      )
      if (is.null(df) || nrow(df) == 0) next

      id_col <- if ("gs_name" %in% names(df)) "gs_name" else "standard_name"
      sym_col <- if ("gene_symbol" %in% names(df)) "gene_symbol" else "human_gene_symbol"
      if (!id_col %in% names(df) || !sym_col %in% names(df)) {
        cat("⚠ Unexpected msigdbr output for ", coll, ":",
            sub %||% "", " — skipping\n", sep = "")
        next
      }

      db_key <- if (is.null(sub)) coll else paste0(coll, ":", sub)
      by_id <- split(as.character(df[[sym_col]]), as.character(df[[id_col]]))
      meta <- unique(data.frame(
        pathway       = as.character(df[[id_col]]),
        collection    = coll,
        subcollection = sub %||% "",
        description   = as.character(df[[id_col]]),
        stringsAsFactors = FALSE
      ))
      meta$clean_name <- vapply(seq_len(nrow(meta)), function(i) {
        .pathway_clean_name(meta$pathway[i], meta$collection[i],
                            meta$subcollection[i], meta$description[i])
      }, character(1))
      out[[db_key]] <- list(by_id = by_id, meta = meta)
      cat("  ✓ Loaded ", length(by_id), " gene sets for ", db_key, "\n",
          sep = "")
    }
  }
  if (length(out) == 0) {
    stop("No MSigDB gene sets could be loaded. Check species + connectivity.")
  }
  out
}

# ==============================================================================
# Resolve comparison groups — returns list(label=..., cell_ids_a=..., cell_ids_b=...)
# ==============================================================================

.pathway_match_meta_filter <- function(meta, filter_spec) {
  # Return logical vector of length nrow(meta) where ALL keys in filter_spec match.
  if (is.null(filter_spec) || length(filter_spec) == 0) {
    return(rep(TRUE, nrow(meta)))
  }
  hits <- rep(TRUE, nrow(meta))
  for (key in names(filter_spec)) {
    want <- as.character(filter_spec[[key]])
    if (!key %in% names(meta)) {
      cat("  ⚠ Comparison filter key '", key,
          "' not found in cell metadata; treating as mismatch\n", sep = "")
      return(rep(FALSE, nrow(meta)))
    }
    have <- as.character(meta[[key]])
    hits <- hits & (have %in% want)
  }
  hits
}

.pathway_resolve_groups <- function(cfg_pathway, cell_meta) {
  comps <- cfg_pathway$comparisons
  if (is.null(comps) || length(comps) == 0) {
    stop("pathway.comparisons is empty — nothing to analyse.")
  }
  out <- vector("list", length(comps))
  for (i in seq_along(comps)) {
    comp <- comps[[i]]
    label <- comp$label %||% sprintf("comparison_%02d", i)
    mask_a <- .pathway_match_meta_filter(cell_meta, comp$group_a)
    mask_b <- .pathway_match_meta_filter(cell_meta, comp$group_b)
    ids_a <- as.character(cell_meta$cell_ID[mask_a])
    ids_b <- as.character(cell_meta$cell_ID[mask_b])
    out[[i]] <- list(
      label    = label,
      group_a  = comp$group_a,
      group_b  = comp$group_b,
      ids_a    = ids_a,
      ids_b    = ids_b,
      n_a      = length(ids_a),
      n_b      = length(ids_b)
    )
  }
  out
}

# ==============================================================================
# DE — presto::wilcoxauc (very fast Wilcoxon on sparse matrix). Returns a
# tidy data frame suitable for both GSEA ranking and ORA filtering.
# ==============================================================================

.pathway_run_de <- function(expr_mat, ids_a, ids_b) {
  if (!requireNamespace("presto", quietly = TRUE)) {
    stop("presto is required. Install with remotes::install_github('immunogenomics/presto').")
  }
  both_ids <- c(ids_a, ids_b)
  both_ids <- intersect(both_ids, colnames(expr_mat))
  ids_a <- intersect(ids_a, both_ids)
  ids_b <- intersect(ids_b, both_ids)
  if (length(ids_a) < 3L || length(ids_b) < 3L) {
    return(NULL)
  }
  y <- ifelse(both_ids %in% ids_a, "A",
              ifelse(both_ids %in% ids_b, "B", NA))
  keep <- !is.na(y)
  y <- y[keep]; both_ids <- both_ids[keep]
  sub_expr <- expr_mat[, both_ids, drop = FALSE]

  wx <- tryCatch(
    presto::wilcoxauc(sub_expr, y = y),
    error = function(e) {
      cat("  ⚠ presto::wilcoxauc error: ", conditionMessage(e), "\n", sep = "")
      NULL
    }
  )
  if (is.null(wx) || nrow(wx) == 0) return(NULL)
  group_a_rows <- wx[wx$group == "A", , drop = FALSE]
  if (nrow(group_a_rows) == 0) return(NULL)
  de <- data.frame(
    gene    = as.character(group_a_rows$feature),
    logFC   = suppressWarnings(as.numeric(group_a_rows$logFC)),
    auc     = suppressWarnings(as.numeric(group_a_rows$auc)),
    pval    = suppressWarnings(as.numeric(group_a_rows$pval)),
    padj    = suppressWarnings(as.numeric(group_a_rows$padj)),
    pct_a   = suppressWarnings(as.numeric(group_a_rows$pct_in)),
    pct_b   = suppressWarnings(as.numeric(group_a_rows$pct_out)),
    stringsAsFactors = FALSE
  )
  de$pval[!is.finite(de$pval)] <- 1
  de$padj[!is.finite(de$padj)] <- 1
  de$logFC[!is.finite(de$logFC)] <- 0
  de$stat <- sign(de$logFC) * -log10(pmax(de$pval, 1e-300))
  de <- de[order(-abs(de$stat)), , drop = FALSE]
  de
}

# ==============================================================================
# GSEA — run fgsea on each database's gene sets. Returns a long-format tibble
# with collection/subcollection/clean_name columns joined in.
# ==============================================================================

.pathway_build_rank_vec <- function(de) {
  ranks <- de$stat
  names(ranks) <- de$gene
  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- ranks[is.finite(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)
  ranks
}

.pathway_run_gsea <- function(de, genesets_by_db, cfg_pathway) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("fgsea is required. Install via BiocManager::install('fgsea').")
  }
  min_size <- cfg_pathway$gsea$min_size %||% 10
  max_size <- cfg_pathway$gsea$max_size %||% 500
  nperm_seed <- cfg_pathway$gsea$nperm_seed %||% 42

  ranks <- .pathway_build_rank_vec(de)
  if (length(ranks) < min_size) return(NULL)

  set.seed(nperm_seed)
  out_list <- list()
  for (db_key in names(genesets_by_db)) {
    entry <- genesets_by_db[[db_key]]
    res <- tryCatch(
      suppressWarnings(fgsea::fgsea(
        pathways  = entry$by_id,
        stats     = ranks,
        minSize   = min_size,
        maxSize   = max_size,
        eps       = 0
      )),
      error = function(e) {
        cat("  ⚠ fgsea error for ", db_key, ": ", conditionMessage(e),
            "\n", sep = "")
        NULL
      }
    )
    if (is.null(res) || nrow(res) == 0) next
    res_df <- as.data.frame(res)
    # fgsea returns `leadingEdge` as a list-column — collapse to comma-string.
    if ("leadingEdge" %in% names(res_df)) {
      res_df$leadingEdge <- vapply(
        res_df$leadingEdge,
        function(x) if (length(x) == 0) "" else paste(x, collapse = ","),
        character(1)
      )
    }
    meta <- entry$meta
    joined <- merge(
      res_df, meta,
      by.x = "pathway", by.y = "pathway",
      all.x = TRUE, sort = FALSE
    )
    out_list[[db_key]] <- joined
  }
  if (length(out_list) == 0) return(NULL)
  do.call(rbind, lapply(out_list, function(d) {
    d[, c("pathway", "pval", "padj", "ES", "NES", "size",
          "leadingEdge", "collection", "subcollection", "clean_name")]
  }))
}

# ==============================================================================
# ORA — clusterProfiler::enricher with TERM2GENE per database. Runs once on
# up-regulated genes (logFC > threshold & padj < threshold) and once on
# down-regulated (same thresholds, opposite sign), then merges with a
# "direction" column.
# ==============================================================================

.pathway_run_ora <- function(de, genesets_by_db, cfg_pathway) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler is required. Install via BiocManager::install('clusterProfiler').")
  }
  fc_thr   <- cfg_pathway$ora$logfc_threshold %||% 0.5
  padj_thr <- cfg_pathway$ora$padj_threshold %||% 0.05
  universe <- unique(as.character(de$gene))
  up_genes <- de$gene[de$logFC >= fc_thr & de$padj <= padj_thr]
  dn_genes <- de$gene[de$logFC <= -fc_thr & de$padj <= padj_thr]

  run_one <- function(genes, direction) {
    if (length(genes) < 5) return(NULL)
    out <- list()
    for (db_key in names(genesets_by_db)) {
      entry <- genesets_by_db[[db_key]]
      t2g <- do.call(rbind, lapply(names(entry$by_id), function(p) {
        data.frame(term = p, gene = entry$by_id[[p]], stringsAsFactors = FALSE)
      }))
      res <- tryCatch(
        clusterProfiler::enricher(
          gene          = unique(genes),
          TERM2GENE     = t2g,
          universe      = universe,
          minGSSize     = 5,
          maxGSSize     = 1000,
          pvalueCutoff  = 1,
          qvalueCutoff  = 1
        ),
        error = function(e) {
          cat("  ⚠ enricher error for ", db_key, " (", direction, "): ",
              conditionMessage(e), "\n", sep = "")
          NULL
        }
      )
      if (is.null(res)) next
      df <- as.data.frame(res)
      if (nrow(df) == 0) next
      df$direction <- direction
      meta <- entry$meta
      df <- merge(df, meta[, c("pathway", "collection", "subcollection",
                               "clean_name")],
                  by.x = "ID", by.y = "pathway", all.x = TRUE, sort = FALSE)
      out[[db_key]] <- df
    }
    if (length(out) == 0) return(NULL)
    do.call(rbind, out)
  }

  up_df <- run_one(up_genes, "up")
  dn_df <- run_one(dn_genes, "down")
  combined <- rbind(up_df, dn_df)
  if (is.null(combined) || nrow(combined) == 0) return(NULL)
  keep_cols <- intersect(
    c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust",
      "qvalue", "Count", "geneID", "direction", "collection", "subcollection",
      "clean_name"),
    names(combined)
  )
  combined[, keep_cols, drop = FALSE]
}

# ==============================================================================
# PROGENy activity via decoupleR (preferred) or direct progeny fallback.
# Returns a cell x 14 matrix of activity scores.
# ==============================================================================

.pathway_run_progeny <- function(expr_mat, top_n_genes = 500) {
  if (requireNamespace("decoupleR", quietly = TRUE) &&
      requireNamespace("progeny", quietly = TRUE)) {
    model <- tryCatch(
      decoupleR::get_progeny(organism = "human", top = top_n_genes),
      error = function(e) NULL
    )
    if (!is.null(model) && nrow(model) > 0) {
      res <- tryCatch(
        decoupleR::run_mlm(
          mat = as.matrix(expr_mat),
          net = model,
          .source = "source",
          .target = "target",
          .mor    = "weight",
          minsize = 5
        ),
        error = function(e) {
          cat("  ⚠ decoupleR::run_mlm failed: ", conditionMessage(e),
              "\n", sep = "")
          NULL
        }
      )
      if (!is.null(res) && nrow(res) > 0) {
        wide <- tryCatch({
          # res is long: condition, source, score, p_value, statistic
          # Pivot score by source x condition (condition = cell_id)
          d <- as.data.frame(res)
          d <- d[d$statistic == "mlm", , drop = FALSE]
          mat <- matrix(
            NA_real_,
            nrow = length(unique(d$condition)),
            ncol = length(unique(d$source)),
            dimnames = list(sort(unique(d$condition)),
                            sort(unique(d$source)))
          )
          for (i in seq_len(nrow(d))) {
            mat[d$condition[i], d$source[i]] <- d$score[i]
          }
          mat
        }, error = function(e) NULL)
        if (!is.null(wide)) return(wide)
      }
    }
  }
  if (requireNamespace("progeny", quietly = TRUE)) {
    scores <- tryCatch(
      progeny::progeny(as.matrix(expr_mat), scale = TRUE,
                       organism = "Human", top = top_n_genes, perm = 1),
      error = function(e) {
        cat("  ⚠ progeny::progeny failed: ", conditionMessage(e),
            "\n", sep = "")
        NULL
      }
    )
    if (!is.null(scores)) return(scores)
  }
  NULL
}

# ==============================================================================
# Plot helpers — each returns the saved file path (or NULL). All PNG writes
# go through save_presentation_plot() (falls back to direct grDevices::png /
# ggsave as needed, see cerebrum bug-060).
# ==============================================================================

.pathway_safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) "unnamed" else x
}

.pathway_plot_gsea_top_bar <- function(gsea_df, outfile, title,
                                       top_n = 20, width = 9, height = 7) {
  df <- gsea_df[!is.na(gsea_df$NES), , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  df <- df[order(-abs(df$NES)), , drop = FALSE]
  df <- head(df, top_n)
  df$label <- ifelse(!is.na(df$clean_name) & nzchar(df$clean_name),
                     df$clean_name, df$pathway)
  df$label <- factor(df$label, levels = df$label[order(df$NES)])
  p <- ggplot2::ggplot(df, ggplot2::aes(x = label, y = NES, fill = NES)) +
    ggplot2::geom_col(colour = "grey30", linewidth = 0.15) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_gradient2(low = "steelblue", mid = "grey85",
                                  high = "firebrick", midpoint = 0) +
    ggplot2::labs(title = title, x = NULL, y = "NES",
                  fill = "NES") +
    ggplot2::theme_minimal(base_size = 10)
  if (exists("save_presentation_plot")) {
    save_presentation_plot(p, outfile, width = width, height = height,
                           dpi = 150)
  } else {
    ggplot2::ggsave(outfile, p, width = width, height = height, dpi = 150)
  }
  outfile
}

.pathway_plot_ora_top_bar <- function(ora_df, outfile, title,
                                      top_n = 20, width = 9, height = 7) {
  if (is.null(ora_df) || nrow(ora_df) == 0) return(NULL)
  df <- ora_df
  df$padj <- suppressWarnings(as.numeric(df$p.adjust))
  df <- df[!is.na(df$padj), , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  df <- df[order(df$padj), , drop = FALSE]
  df <- head(df, top_n)
  df$nlog10 <- -log10(pmax(df$padj, 1e-300))
  df$count <- suppressWarnings(as.numeric(df$Count))
  df$label <- ifelse(!is.na(df$clean_name) & nzchar(df$clean_name),
                     df$clean_name, df$ID)
  df$label <- factor(df$label, levels = df$label[order(df$nlog10)])
  p <- ggplot2::ggplot(df, ggplot2::aes(x = label, y = nlog10, fill = count)) +
    ggplot2::geom_col(colour = "grey30", linewidth = 0.15) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_viridis_c(option = "D") +
    ggplot2::labs(title = title, x = NULL,
                  y = "-log10(p.adjust)", fill = "gene count") +
    ggplot2::facet_wrap(~ direction, ncol = 1, scales = "free_y") +
    ggplot2::theme_minimal(base_size = 10)
  if (exists("save_presentation_plot")) {
    save_presentation_plot(p, outfile, width = width, height = height,
                           dpi = 150)
  } else {
    ggplot2::ggsave(outfile, p, width = width, height = height, dpi = 150)
  }
  outfile
}

.pathway_plot_gsea_volcano <- function(gsea_df, outfile, title,
                                       top_n = 15, width = 8, height = 7) {
  df <- gsea_df[!is.na(gsea_df$NES) & !is.na(gsea_df$padj), , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  df$nlog10 <- -log10(pmax(df$padj, 1e-300))
  df$label <- ifelse(!is.na(df$clean_name) & nzchar(df$clean_name),
                     df$clean_name, df$pathway)
  df$sig <- df$padj <= 0.05
  top <- df[order(-abs(df$NES) * df$nlog10), , drop = FALSE]
  top <- head(top, top_n)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = NES, y = nlog10, colour = sig)) +
    ggplot2::geom_point(size = 1.2, alpha = 0.7) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "firebrick",
                                            `FALSE` = "grey60"),
                                 name = "padj<=0.05") +
    ggplot2::labs(title = title, x = "NES", y = "-log10(padj)")
  if (requireNamespace("ggrepel", quietly = TRUE) && nrow(top) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = top, ggplot2::aes(label = label), size = 2.5,
      max.overlaps = 30, colour = "grey20"
    )
  }
  p <- p + ggplot2::theme_minimal(base_size = 10)
  if (exists("save_presentation_plot")) {
    save_presentation_plot(p, outfile, width = width, height = height,
                           dpi = 150)
  } else {
    ggplot2::ggsave(outfile, p, width = width, height = height, dpi = 150)
  }
  outfile
}

.pathway_plot_dotplot_by_db <- function(gsea_df, outfile, title,
                                        top_per_db = 5, width = 10,
                                        height = 8) {
  df <- gsea_df[!is.na(gsea_df$NES) & !is.na(gsea_df$padj), , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  df$nlog10 <- -log10(pmax(df$padj, 1e-300))
  df$db <- ifelse(nzchar(df$subcollection),
                  paste0(df$collection, ":", df$subcollection),
                  df$collection)
  top <- do.call(rbind, lapply(split(df, df$db), function(g) {
    head(g[order(-abs(g$NES)), , drop = FALSE], top_per_db)
  }))
  if (is.null(top) || nrow(top) == 0) return(NULL)
  top$label <- ifelse(!is.na(top$clean_name) & nzchar(top$clean_name),
                      top$clean_name, top$pathway)
  top$label <- factor(top$label,
                      levels = unique(top$label[order(top$db, top$NES)]))
  p <- ggplot2::ggplot(top, ggplot2::aes(x = NES, y = label,
                                         size = nlog10, colour = NES)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient2(low = "steelblue", mid = "grey85",
                                    high = "firebrick", midpoint = 0) +
    ggplot2::scale_size_continuous(name = "-log10(padj)") +
    ggplot2::facet_wrap(~ db, scales = "free_y", ncol = 2) +
    ggplot2::labs(title = title, x = "NES", y = NULL) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
  if (exists("save_presentation_plot")) {
    save_presentation_plot(p, outfile, width = width, height = height,
                           dpi = 150)
  } else {
    ggplot2::ggsave(outfile, p, width = width, height = height, dpi = 150)
  }
  outfile
}

.pathway_plot_cross_stratum_heatmap <- function(gsea_long, outfile, title,
                                                top_n = 30, width = 12,
                                                height = 9) {
  # gsea_long: rows across multiple strata for a single comparison, with
  # stratum column labelling rows.
  df <- gsea_long
  df <- df[!is.na(df$NES) & !is.na(df$padj), , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  df$signed_nlog10 <- sign(df$NES) * -log10(pmax(df$padj, 1e-300))
  df$label <- ifelse(!is.na(df$clean_name) & nzchar(df$clean_name),
                     df$clean_name, df$pathway)
  # Pick top pathways by max |signed_nlog10| across strata.
  agg <- aggregate(abs(df$signed_nlog10),
                   by = list(pathway = df$label), FUN = max)
  agg <- agg[order(-agg$x), , drop = FALSE]
  top_paths <- head(agg$pathway, top_n)
  df2 <- df[df$label %in% top_paths, , drop = FALSE]
  if (nrow(df2) == 0) return(NULL)
  df2$label <- factor(df2$label, levels = rev(top_paths))
  p <- ggplot2::ggplot(df2, ggplot2::aes(x = stratum, y = label,
                                         fill = signed_nlog10)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.15) +
    ggplot2::scale_fill_gradient2(
      low = "steelblue", mid = "white", high = "firebrick",
      midpoint = 0, name = "sign(NES) * -log10(padj)"
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1))
  if (exists("save_presentation_plot")) {
    save_presentation_plot(p, outfile, width = width, height = height,
                           dpi = 150)
  } else {
    ggplot2::ggsave(outfile, p, width = width, height = height, dpi = 150)
  }
  outfile
}

.pathway_plot_cross_comparison_heatmap <- function(gsea_by_comp, stratum_label,
                                                   outfile, width = 9,
                                                   height = 10) {
  # gsea_by_comp: named list (comparison -> gsea_df for the given stratum)
  all_rows <- do.call(rbind, lapply(names(gsea_by_comp), function(nm) {
    d <- gsea_by_comp[[nm]]
    if (is.null(d) || nrow(d) == 0) return(NULL)
    d$comparison <- nm
    d
  }))
  if (is.null(all_rows) || nrow(all_rows) == 0) return(NULL)
  sig_paths <- unique(all_rows$pathway[!is.na(all_rows$padj) &
                                       all_rows$padj <= 0.05])
  if (length(sig_paths) == 0) return(NULL)
  all_rows <- all_rows[all_rows$pathway %in% sig_paths, , drop = FALSE]
  all_rows$signed_nlog10 <- sign(all_rows$NES) *
    -log10(pmax(all_rows$padj, 1e-300))
  all_rows$label <- ifelse(!is.na(all_rows$clean_name) &
                             nzchar(all_rows$clean_name),
                           all_rows$clean_name, all_rows$pathway)
  all_rows$label <- factor(all_rows$label,
                           levels = unique(all_rows$label[order(
                             -abs(all_rows$signed_nlog10))]))
  p <- ggplot2::ggplot(all_rows, ggplot2::aes(x = comparison, y = label,
                                              fill = signed_nlog10)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.15) +
    ggplot2::scale_fill_gradient2(low = "steelblue", mid = "white",
                                  high = "firebrick", midpoint = 0) +
    ggplot2::labs(title = paste0("Cross-comparison pathway heatmap — ",
                                 stratum_label),
                  x = NULL, y = NULL,
                  fill = "sign(NES) * -log10(padj)") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  if (exists("save_presentation_plot")) {
    save_presentation_plot(p, outfile, width = width, height = height,
                           dpi = 150)
  } else {
    ggplot2::ggsave(outfile, p, width = width, height = height, dpi = 150)
  }
  outfile
}

.pathway_plot_upset <- function(gsea_by_comp, outfile, width = 9, height = 6) {
  if (!requireNamespace("UpSetR", quietly = TRUE)) return(NULL)
  sets <- lapply(gsea_by_comp, function(d) {
    if (is.null(d) || nrow(d) == 0) return(character(0))
    sig <- d[!is.na(d$padj) & d$padj <= 0.05, "pathway"]
    as.character(sig)
  })
  sets <- sets[lengths(sets) > 0]
  if (length(sets) < 2) return(NULL)
  ok <- tryCatch({
    grDevices::png(outfile, width = width, height = height, units = "in",
                   res = 150, bg = "white", type = "cairo")
    on.exit(grDevices::dev.off(), add = TRUE)
    print(UpSetR::upset(UpSetR::fromList(sets), order.by = "freq",
                        nsets = length(sets)))
    TRUE
  }, error = function(e) {
    try(grDevices::dev.off(), silent = TRUE)
    cat("  ⚠ upset plot error: ", conditionMessage(e), "\n", sep = "")
    FALSE
  })
  if (isTRUE(ok)) outfile else NULL
}

.pathway_plot_leading_edge <- function(de, pathway_id, leading_edge,
                                       expr_mat, cell_meta,
                                       sample_col, group_col,
                                       outfile, title, width = 8, height = 6) {
  genes <- strsplit(leading_edge, ",", fixed = TRUE)[[1]]
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) < 3) return(NULL)
  ids <- intersect(colnames(expr_mat), as.character(cell_meta$cell_ID))
  if (length(ids) == 0) return(NULL)
  sub_expr <- as.matrix(expr_mat[genes, ids, drop = FALSE])
  meta_sub <- cell_meta[match(ids, cell_meta$cell_ID), , drop = FALSE]
  sample_ids <- as.character(meta_sub[[sample_col]])
  agg <- matrix(NA_real_, nrow = length(genes),
                ncol = length(unique(sample_ids)),
                dimnames = list(genes, sort(unique(sample_ids))))
  for (s in colnames(agg)) {
    cells_s <- which(sample_ids == s)
    if (length(cells_s) == 0) next
    agg[, s] <- rowMeans(sub_expr[, cells_s, drop = FALSE], na.rm = TRUE)
  }
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    cat("  ⚠ pheatmap missing; skipping leading-edge plot for ",
        pathway_id, "\n", sep = "")
    return(NULL)
  }
  # Optional column annotation by group_col
  anno <- NULL
  if (!is.null(group_col) && group_col %in% names(meta_sub)) {
    g_map <- tapply(as.character(meta_sub[[group_col]]), sample_ids,
                    function(v) names(sort(table(v), decreasing = TRUE))[1])
    anno <- data.frame(group = unname(g_map[colnames(agg)]),
                       row.names = colnames(agg),
                       stringsAsFactors = FALSE)
  }
  ok <- tryCatch({
    grDevices::png(outfile, width = width, height = height, units = "in",
                   res = 150, bg = "white", type = "cairo")
    on.exit(grDevices::dev.off(), add = TRUE)
    pheatmap::pheatmap(agg, main = title, scale = "row",
                       annotation_col = anno, silent = FALSE,
                       color = grDevices::colorRampPalette(
                         c("steelblue", "white", "firebrick"))(100))
    TRUE
  }, error = function(e) {
    try(grDevices::dev.off(), silent = TRUE)
    cat("  ⚠ leading-edge heatmap error: ", conditionMessage(e), "\n",
        sep = "")
    FALSE
  })
  if (isTRUE(ok)) outfile else NULL
}

.pathway_plot_progeny_polygons <- function(poly_df, progeny_mat,
                                           bcell_ids, outfile,
                                           ncol_grid = 4, width = 14,
                                           height = 10) {
  if (is.null(poly_df) || is.null(progeny_mat)) return(NULL)
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    cat("  ⚠ patchwork missing; skipping PROGENy polygon grid\n")
    return(NULL)
  }
  pathways <- colnames(progeny_mat)
  poly_df$.is_bcell <- poly_df$cell_ID %in% bcell_ids
  panels <- lapply(pathways, function(p) {
    vals <- stats::setNames(progeny_mat[, p], rownames(progeny_mat))
    poly_df2 <- poly_df
    poly_df2$.val <- unname(vals[poly_df2$cell_ID])
    bcell_df <- poly_df2[poly_df2$.is_bcell, , drop = FALSE]
    ggplot2::ggplot() +
      ggplot2::geom_polygon(
        data = poly_df2,
        mapping = ggplot2::aes(x = x, y = y, group = poly_group, fill = .val),
        colour = NA, linewidth = 0
      ) +
      ggplot2::geom_polygon(
        data = bcell_df,
        mapping = ggplot2::aes(x = x, y = y, group = poly_group),
        fill = NA, colour = "black", linewidth = 0.15
      ) +
      ggplot2::scale_fill_gradient2(low = "steelblue", mid = "white",
                                    high = "firebrick", midpoint = 0,
                                    na.value = "grey90") +
      ggplot2::coord_equal() +
      ggplot2::labs(title = p, x = NULL, y = NULL) +
      ggplot2::theme_void(base_size = 8) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 9, face = "bold"),
        legend.position = "none"
      )
  })
  combined <- patchwork::wrap_plots(panels, ncol = ncol_grid) +
    patchwork::plot_annotation(
      title = "PROGENy pathway activity — B cells outlined"
    )
  ok <- tryCatch({
    ggplot2::ggsave(outfile, combined, width = width, height = height,
                    dpi = 150)
    TRUE
  }, error = function(e) {
    cat("  ⚠ PROGENy polygon save failed: ", conditionMessage(e),
        "\n", sep = "")
    FALSE
  })
  if (isTRUE(ok)) outfile else NULL
}

.pathway_plot_progeny_bcell_by_comparison <- function(progeny_mat, cell_meta,
                                                      comparisons,
                                                      bcell_ids,
                                                      outfile,
                                                      width = 12, height = 8) {
  bcell_ids <- intersect(bcell_ids, rownames(progeny_mat))
  if (length(bcell_ids) == 0) return(NULL)
  meta <- cell_meta[match(bcell_ids, cell_meta$cell_ID), , drop = FALSE]
  long <- list()
  for (comp in comparisons) {
    mask_a <- .pathway_match_meta_filter(meta, comp$group_a)
    mask_b <- .pathway_match_meta_filter(meta, comp$group_b)
    ids_a <- bcell_ids[mask_a]
    ids_b <- bcell_ids[mask_b]
    for (pw in colnames(progeny_mat)) {
      if (length(ids_a) >= 3) long[[length(long) + 1]] <- data.frame(
        comparison = comp$label, pathway = pw, arm = "A",
        score = progeny_mat[ids_a, pw], stringsAsFactors = FALSE
      )
      if (length(ids_b) >= 3) long[[length(long) + 1]] <- data.frame(
        comparison = comp$label, pathway = pw, arm = "B",
        score = progeny_mat[ids_b, pw], stringsAsFactors = FALSE
      )
    }
  }
  if (length(long) == 0) return(NULL)
  df <- do.call(rbind, long)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = arm, y = score, fill = arm)) +
    ggplot2::geom_boxplot(outlier.size = 0.3, alpha = 0.75) +
    ggplot2::facet_grid(pathway ~ comparison, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(A = "firebrick", B = "steelblue")) +
    ggplot2::labs(title = "PROGENy activity — B cells, by comparison",
                  x = NULL, y = "activity z-score") +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(legend.position = "none",
                   strip.text.y = ggplot2::element_text(angle = 0, size = 7))
  if (exists("save_presentation_plot")) {
    save_presentation_plot(p, outfile, width = width, height = height,
                           dpi = 150)
  } else {
    ggplot2::ggsave(outfile, p, width = width, height = height, dpi = 150)
  }
  outfile
}

.pathway_build_pdf <- function(image_paths, out_pdf, title) {
  image_paths <- image_paths[file.exists(image_paths)]
  if (length(image_paths) == 0) return(NULL)
  # Simplest portable approach: use `magick` if available.
  if (requireNamespace("magick", quietly = TRUE)) {
    ok <- tryCatch({
      imgs <- lapply(image_paths, magick::image_read)
      combined <- Reduce(c, imgs)
      magick::image_write(combined, out_pdf, format = "pdf")
      TRUE
    }, error = function(e) {
      cat("  ⚠ magick PDF build failed: ", conditionMessage(e), "\n",
          sep = "")
      FALSE
    })
    if (isTRUE(ok)) return(out_pdf)
  }
  # Fallback: one-page-per-image PDF via R's grDevices::pdf().
  ok <- tryCatch({
    grDevices::pdf(out_pdf, width = 11, height = 8.5, onefile = TRUE,
                   title = title)
    on.exit(grDevices::dev.off(), add = TRUE)
    for (pth in image_paths) {
      img <- tryCatch(grDevices::as.raster(png::readPNG(pth)),
                      error = function(e) NULL)
      if (is.null(img)) next
      graphics::plot.new()
      graphics::grid.raster(img)
    }
    TRUE
  }, error = function(e) {
    try(grDevices::dev.off(), silent = TRUE)
    cat("  ⚠ PDF fallback failed: ", conditionMessage(e), "\n",
        sep = "")
    FALSE
  })
  if (isTRUE(ok)) out_pdf else NULL
}

# ==============================================================================
# Metadata sanity guard + column resolution
# ==============================================================================

.pathway_resolve_celltype_column <- function(meta, explicit = NULL) {
  if (!is.null(explicit) && nzchar(explicit) && explicit %in% names(meta)) {
    return(explicit)
  }
  candidates <- grep("^celltype_.*_supervised$", names(meta), value = TRUE)
  if (length(candidates) > 0) return(candidates[1])
  if ("celltype" %in% names(meta)) return("celltype")
  NA_character_
}

.pathway_ensure_metadata <- function(gobj, sample_table = NULL,
                                     required = c("treatment", "timepoint")) {
  meta <- as.data.frame(.pathway_pdata_dt(gobj))
  missing <- setdiff(required, names(meta))
  if (length(missing) == 0) return(meta)
  if (is.null(sample_table) || !"sample_id" %in% names(meta)) {
    stop("Required metadata columns missing and no sample_table provided: ",
         paste(missing, collapse = ", "))
  }
  st <- as.data.frame(sample_table)
  need <- intersect(missing, names(st))
  if (length(need) == 0) {
    stop("sample_table has no matching columns for: ",
         paste(missing, collapse = ", "))
  }
  add <- st[match(as.character(meta$sample_id), as.character(st$sample_id)),
            need, drop = FALSE]
  meta <- cbind(meta, add)
  cat("  ℹ Injected ", length(need),
      " sample-level column(s) into cell metadata: ",
      paste(need, collapse = ", "), "\n", sep = "")
  meta
}

# ==============================================================================
# Orchestrator
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

run_pathway_analysis <- function(gobj,
                                 sample_id = "merged",
                                 cfg,
                                 output_dir,
                                 sample_table = NULL,
                                 celltype_column = NULL,
                                 leiden_column   = "leiden_clust",
                                 focus_celltype  = "^B cell$") {
  cat("\n========================================\n")
  cat("STEP 13: Pathway Enrichment + GSEA\n")
  cat("Run label:", sample_id, "\n")
  cat("========================================\n\n")

  cfg_p <- cfg$pathway
  if (is.null(cfg_p) || !isTRUE(cfg_p$enabled)) {
    cat("⚠ pathway.enabled is FALSE or the pathway: block is missing. ",
        "Nothing to run.\n", sep = "")
    return(invisible(list(status = "disabled")))
  }

  out_root <- file.path(output_dir, "13_Pathway_Analysis")
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

  if (is.character(gobj)) {
    cat("Loading Giotto object from:", gobj, "\n")
    gobj <- loadGiotto(gobj)
    cat("✓ Loaded\n")
  }

  meta <- .pathway_ensure_metadata(
    gobj, sample_table = sample_table,
    required = c("treatment", "timepoint", "sample_id")
  )
  if (is.null(celltype_column) || !nzchar(celltype_column)) {
    celltype_column <- cfg_p$celltype_column %||%
      .pathway_resolve_celltype_column(meta)
  }
  if (is.na(celltype_column) || !celltype_column %in% names(meta)) {
    stop("Could not resolve a supervised cell-type column (tried ",
         "celltype_*_supervised). Set cfg$pathway$celltype_column ",
         "or annotate first.")
  }
  leiden_column <- leiden_column %||% cfg_p$leiden_column %||% "leiden_clust"
  if (!leiden_column %in% names(meta)) {
    cat("⚠ Leiden column '", leiden_column,
        "' not in metadata — will skip cluster stratification\n", sep = "")
    leiden_column <- NA_character_
  }
  cat("Cell-type column :", celltype_column, "\n")
  cat("Leiden column    :", leiden_column, "\n")
  cat("Focus cell type  :", focus_celltype %||% "<none>", "\n\n")

  # Load MSigDB gene sets once.
  cat("Loading MSigDB gene sets...\n")
  genesets_by_db <- .pathway_load_genesets(cfg_p)

  # Resolve all 4 comparison ID vectors.
  comparisons <- .pathway_resolve_groups(cfg_p, meta)
  valid_comparisons <- Filter(function(cp) {
    min_cells <- cfg_p$min_cells_per_stratum %||% 30
    ok <- cp$n_a >= min_cells && cp$n_b >= min_cells
    if (!ok) {
      cat("⚠ Comparison '", cp$label,
          "' skipped: groups too small (A=", cp$n_a,
          ", B=", cp$n_b, ")\n", sep = "")
    }
    ok
  }, comparisons)
  if (length(valid_comparisons) == 0) {
    stop("No comparisons passed the min_cells_per_stratum filter.")
  }

  # Pre-extract normalized expression for the merged object — reused for
  # all comparisons × strata.
  cat("Extracting normalized expression matrix...\n")
  expr_mat <- .pathway_get_expression(gobj, values = "normalized",
                                      output = "matrix")
  expr_mat <- as.matrix(expr_mat)
  cat("✓ Expression matrix: ", paste(dim(expr_mat), collapse = " x "),
      "\n", sep = "")

  # Run per comparison / per stratification / per stratum.
  all_results <- list()
  bcell_plot_paths <- character(0)
  # For cross-comparison summaries: collect per-(comparison, stratum_type,
  # stratum) GSEA data frames keyed on comparison.
  cross_comp_store <- list()   # list per (stratum_type, stratum) of list(comparison -> df)

  focus_snake <- "B_cell"
  bcell_ids_all <- .pathway_resolve_bcell_ids(
    gobj, celltype_column,
    focus_celltype = focus_celltype %||% "^B cell$"
  )
  cat("B cells in merged object :", length(bcell_ids_all), "\n\n")

  for (cp in valid_comparisons) {
    cat("\n--- Comparison: ", cp$label, " (A=", cp$n_a, ", B=", cp$n_b,
        ") ---\n", sep = "")
    comp_dir <- file.path(out_root, cp$label)
    dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

    # For each stratification, iterate over strata.
    strata_types <- cfg_p$stratifications %||% c("celltype", "leiden")
    cross_stratum_gsea <- list()     # list(stratum_type -> list(stratum -> gsea_df))
    for (st_type in strata_types) {
      st_col <- switch(st_type,
                       celltype = celltype_column,
                       leiden   = leiden_column,
                       NA_character_)
      if (is.na(st_col) || !st_col %in% names(meta)) {
        cat("  - ", st_type, ": column not available; skipping\n", sep = "")
        next
      }
      st_vec <- as.character(meta[[st_col]])
      st_levels <- sort(unique(st_vec[!is.na(st_vec) & nzchar(st_vec)]))
      cross_stratum_gsea[[st_type]] <- list()
      for (lvl in st_levels) {
        lvl_cells <- as.character(meta$cell_ID[st_vec == lvl])
        ids_a_s <- intersect(cp$ids_a, lvl_cells)
        ids_b_s <- intersect(cp$ids_b, lvl_cells)
        min_cells <- cfg_p$min_cells_per_stratum %||% 30
        if (length(ids_a_s) < min_cells || length(ids_b_s) < min_cells) next

        stratum_tag <- paste0(st_type, "_", .pathway_safe_name(lvl))
        stratum_dir <- file.path(comp_dir, stratum_tag)
        dir.create(stratum_dir, recursive = TRUE, showWarnings = FALSE)

        cat("    [", stratum_tag, "] A=", length(ids_a_s),
            " B=", length(ids_b_s), "\n", sep = "")

        de <- .pathway_run_de(expr_mat, ids_a_s, ids_b_s)
        if (is.null(de) || nrow(de) < (cfg_p$min_genes_tested %||% 50)) {
          cat("      skipped (too few genes)\n")
          next
        }
        utils::write.csv(
          de,
          file.path(stratum_dir, paste0(stratum_tag, "_de_ranks.csv")),
          row.names = FALSE
        )

        gsea_df <- .pathway_run_gsea(de, genesets_by_db, cfg_p)
        if (!is.null(gsea_df) && nrow(gsea_df) > 0) {
          utils::write.csv(
            gsea_df,
            file.path(stratum_dir, paste0(stratum_tag, "_gsea_all.csv")),
            row.names = FALSE
          )
          gsea_df$comparison   <- cp$label
          gsea_df$stratum_type <- st_type
          gsea_df$stratum      <- lvl
          all_results[[length(all_results) + 1]] <- transform_gsea_to_tidy(gsea_df)
          cross_stratum_gsea[[st_type]][[lvl]] <- gsea_df

          title_main <- sprintf("%s | %s", cp$label, stratum_tag)
          p1 <- .pathway_plot_gsea_top_bar(
            gsea_df,
            file.path(stratum_dir,
                      paste0(stratum_tag, "_gsea_top20_bar.png")),
            title_main)
          p3 <- .pathway_plot_gsea_volcano(
            gsea_df,
            file.path(stratum_dir,
                      paste0(stratum_tag, "_gsea_volcano.png")),
            title_main)
          p4 <- .pathway_plot_dotplot_by_db(
            gsea_df,
            file.path(stratum_dir,
                      paste0(stratum_tag, "_dotplot_by_db.png")),
            title_main)

          # Leading-edge heatmaps (top N GSEA hits, guarded by cap)
          cap <- cfg_p$leading_edge_top_n %||% 5
          if (cap > 0 && !is.null(de)) {
            le_dir <- file.path(comp_dir, "leading_edge")
            dir.create(le_dir, recursive = TRUE, showWarnings = FALSE)
            df_top <- gsea_df[!is.na(gsea_df$padj), , drop = FALSE]
            df_top <- df_top[order(df_top$padj), , drop = FALSE]
            df_top <- head(df_top, cap)
            for (r in seq_len(nrow(df_top))) {
              pid <- df_top$pathway[r]
              .pathway_plot_leading_edge(
                de, pid, df_top$leadingEdge[r], expr_mat, meta,
                sample_col = "sample_id",
                group_col = celltype_column,
                outfile = file.path(le_dir,
                                    paste0(stratum_tag, "_",
                                           .pathway_safe_name(pid),
                                           ".png")),
                title = df_top$clean_name[r]
              )
            }
          }

          # Collect for B-cell-focused PDF
          if (st_type == "celltype" &&
              grepl(focus_celltype, lvl, ignore.case = TRUE)) {
            if (!is.null(p1)) bcell_plot_paths <- c(bcell_plot_paths, p1)
            if (!is.null(p3)) bcell_plot_paths <- c(bcell_plot_paths, p3)
          }
        }

        ora_df <- .pathway_run_ora(de, genesets_by_db, cfg_p)
        if (!is.null(ora_df) && nrow(ora_df) > 0) {
          utils::write.csv(
            ora_df,
            file.path(stratum_dir, paste0(stratum_tag, "_ora_all.csv")),
            row.names = FALSE
          )
          ora_tidy <- ora_df
          ora_tidy$comparison   <- cp$label
          ora_tidy$stratum_type <- st_type
          ora_tidy$stratum      <- lvl
          all_results[[length(all_results) + 1]] <- transform_ora_to_tidy(ora_tidy)

          title_main <- sprintf("%s | %s", cp$label, stratum_tag)
          p2 <- .pathway_plot_ora_top_bar(
            ora_df,
            file.path(stratum_dir,
                      paste0(stratum_tag, "_ora_top20_bar.png")),
            title_main)
          if (st_type == "celltype" &&
              grepl(focus_celltype, lvl, ignore.case = TRUE) &&
              !is.null(p2)) {
            bcell_plot_paths <- c(bcell_plot_paths, p2)
          }
        }
      }
    }

    # P5 / P6 — per-comparison cross-stratum heatmaps
    for (st_type in names(cross_stratum_gsea)) {
      store <- cross_stratum_gsea[[st_type]]
      if (length(store) < 2) next
      long <- do.call(rbind, lapply(names(store), function(nm) {
        d <- store[[nm]]; d$stratum <- nm; d
      }))
      outfile <- file.path(comp_dir, paste0("_cross_", st_type, "_heatmap.png"))
      .pathway_plot_cross_stratum_heatmap(
        long, outfile,
        title = sprintf("%s — top pathways across %s strata",
                        cp$label, st_type)
      )
    }

    # Stash per-comparison GSEA for P7 later
    for (st_type in names(cross_stratum_gsea)) {
      for (lvl in names(cross_stratum_gsea[[st_type]])) {
        key <- paste0(st_type, "::", lvl)
        cross_comp_store[[key]] <- cross_comp_store[[key]] %||% list()
        cross_comp_store[[key]][[cp$label]] <-
          cross_stratum_gsea[[st_type]][[lvl]]
      }
    }

    cleanup_memory(verbose = FALSE)
  }

  # ---------- P7 + P8 : cross-comparison summaries ----------
  cat("\n--- Cross-comparison summaries ---\n")
  for (key in names(cross_comp_store)) {
    parts <- strsplit(key, "::", fixed = TRUE)[[1]]
    st_type <- parts[1]; lvl <- parts[2]
    store <- cross_comp_store[[key]]
    if (length(store) < 2) next
    outfile <- file.path(
      out_root,
      paste0("_cross_comparison_pathway_heatmap_",
             .pathway_safe_name(st_type), "_",
             .pathway_safe_name(lvl), ".png")
    )
    .pathway_plot_cross_comparison_heatmap(store, lvl, outfile)
  }
  # P8 — union upset across the 4 comparisons, aggregated across ALL strata
  # (pool any significant pathway from any stratum into the comparison's set).
  per_comp_paths <- list()
  for (key in names(cross_comp_store)) {
    store <- cross_comp_store[[key]]
    for (cp_name in names(store)) {
      d <- store[[cp_name]]
      if (is.null(d) || nrow(d) == 0) next
      sig <- d$pathway[!is.na(d$padj) & d$padj <= 0.05]
      if (length(sig) == 0) next
      per_comp_paths[[cp_name]] <- unique(c(per_comp_paths[[cp_name]],
                                            sig))
    }
  }
  if (length(per_comp_paths) >= 2) {
    .pathway_plot_upset(
      lapply(per_comp_paths, function(x)
        data.frame(pathway = x, padj = 0.01)),
      file.path(out_root, "_cross_comparison_upset.png")
    )
  }

  # ---------- PROGENy (P10 / P11) ----------
  progeny_cfg <- cfg_p$progeny %||% list(enabled = FALSE)
  if (isTRUE(progeny_cfg$enabled)) {
    bcell_dir <- file.path(out_root, "bcell_focus")
    dir.create(bcell_dir, recursive = TRUE, showWarnings = FALSE)
    cat("\n--- PROGENy pathway activity ---\n")
    progeny_mat <- .pathway_run_progeny(
      expr_mat,
      top_n_genes = progeny_cfg$top_n_genes %||% 500
    )
    if (!is.null(progeny_mat) && nrow(progeny_mat) > 0) {
      bcell_rows <- intersect(bcell_ids_all, rownames(progeny_mat))
      if (length(bcell_rows) > 0) {
        df <- as.data.frame(progeny_mat[bcell_rows, , drop = FALSE])
        df$cell_ID <- rownames(df)
        m_sub <- meta[match(bcell_rows, meta$cell_ID), , drop = FALSE]
        extra <- intersect(c("sample_id", "treatment", "timepoint"),
                           names(m_sub))
        df <- cbind(df, m_sub[, extra, drop = FALSE])
        utils::write.csv(
          df, file.path(bcell_dir, "progeny_bcell_activity.csv"),
          row.names = FALSE
        )

        poly_df <- .pathway_extract_polygon_df(gobj)
        if (!is.null(poly_df)) {
          p10 <- .pathway_plot_progeny_polygons(
            poly_df, progeny_mat, bcell_rows,
            outfile = file.path(bcell_dir,
                                "progeny_bcell_activity_polygons.png")
          )
          if (!is.null(p10)) bcell_plot_paths <- c(bcell_plot_paths, p10)
        }
        p11 <- .pathway_plot_progeny_bcell_by_comparison(
          progeny_mat, meta, valid_comparisons, bcell_rows,
          outfile = file.path(bcell_dir,
                              "progeny_bcell_activity_by_comparison.png")
        )
        if (!is.null(p11)) bcell_plot_paths <- c(bcell_plot_paths, p11)
      }
    } else {
      cat("  ⚠ PROGENy scores unavailable (no decoupleR/progeny?); ",
          "skipping P10/P11\n", sep = "")
    }
  }

  # ---------- Merged results CSV ----------
  if (length(all_results) > 0) {
    merged_df <- do.call(rbind, all_results)
    utils::write.csv(merged_df, file.path(out_root, "_all_results_merged.csv"),
                     row.names = FALSE)
    cat("✓ Wrote merged results: ",
        nrow(merged_df), " rows → _all_results_merged.csv\n", sep = "")
  } else {
    cat("⚠ No enrichment results produced — merged CSV skipped\n")
  }

  # ---------- P12 : B-cell PDF ----------
  if (isTRUE(cfg_p$pdf_summary) && length(bcell_plot_paths) > 0) {
    pdf_out <- file.path(out_root, "bcell_focus", "bcell_pathway_summary.pdf")
    dir.create(dirname(pdf_out), recursive = TRUE, showWarnings = FALSE)
    .pathway_build_pdf(bcell_plot_paths, pdf_out,
                       title = "B-cell pathway summary")
    cat("✓ Wrote B-cell PDF: ", pdf_out, "\n", sep = "")
  }

  cat("\n=== Step 13 summary ===\n")
  for (cp in valid_comparisons) {
    cat(sprintf("  ✓ %-28s A=%d B=%d\n", cp$label, cp$n_a, cp$n_b))
  }
  cat("\n✓ STEP 13 complete for ", sample_id, "\n\n", sep = "")
  invisible(list(status = "ok", results = all_results,
                 bcell_plots = bcell_plot_paths))
}

# ==============================================================================
# Tidy helpers — normalise GSEA/ORA outputs into a common "long" schema for
# the merged results CSV.
# ==============================================================================

transform_gsea_to_tidy <- function(gsea_df) {
  data.frame(
    comparison     = gsea_df$comparison,
    stratum_type   = gsea_df$stratum_type,
    stratum        = gsea_df$stratum,
    test           = "GSEA",
    collection     = gsea_df$collection,
    subcollection  = gsea_df$subcollection,
    pathway        = gsea_df$pathway,
    clean_name     = gsea_df$clean_name,
    stat           = gsea_df$NES,
    pval           = gsea_df$pval,
    padj           = gsea_df$padj,
    n_genes        = gsea_df$size,
    direction      = ifelse(is.na(gsea_df$NES), NA_character_,
                            ifelse(gsea_df$NES >= 0, "up", "down")),
    leading_edge   = gsea_df$leadingEdge,
    stringsAsFactors = FALSE
  )
}

transform_ora_to_tidy <- function(ora_df) {
  enrich_stat <- suppressWarnings(-log10(pmax(
    as.numeric(ora_df$p.adjust), 1e-300)))
  data.frame(
    comparison     = ora_df$comparison,
    stratum_type   = ora_df$stratum_type,
    stratum        = ora_df$stratum,
    test           = "ORA",
    collection     = ora_df$collection,
    subcollection  = ora_df$subcollection,
    pathway        = ora_df$ID,
    clean_name     = ora_df$clean_name,
    stat           = enrich_stat,
    pval           = suppressWarnings(as.numeric(ora_df$pvalue)),
    padj           = suppressWarnings(as.numeric(ora_df$p.adjust)),
    n_genes        = suppressWarnings(as.numeric(ora_df$Count)),
    direction      = ora_df$direction,
    leading_edge   = ora_df$geneID,
    stringsAsFactors = FALSE
  )
}
