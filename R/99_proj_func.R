
collapse_to_genes <- function(tt, feature_tbl) {
  tt <- tt |>
    tibble::as_tibble()
  
  feature_tbl <- feature_tbl |>
    tibble::as_tibble()
  
  tt |>
    dplyr::select(
      ProbeID,
      logFC,
      t,
      P.Value,
      adj.P.Val
    ) |>
    dplyr::left_join(
      feature_tbl |>
        dplyr::select(
          ProbeID,
          Gene.Symbol = `Gene.Symbol`,
          ENTREZ_GENE_ID
        ),
      by = "ProbeID"
    ) |>
    dplyr::filter(
      !is.na(Gene.Symbol),
      Gene.Symbol != ""
    ) |>
    dplyr::group_by(Gene.Symbol) |>
    dplyr::slice_min(
      order_by = adj.P.Val,
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::ungroup() |>
    dplyr::rename(
      symbol = Gene.Symbol,
      entrez = ENTREZ_GENE_ID
    )
}


pick_seeds <- function(
    tt_gene,
    fdr = 0.01,
    lfc = 1.5,
    top_n = NULL
) {
  seeds <- tt_gene |>
    dplyr::filter(
      adj.P.Val < fdr,
      abs(logFC) >= lfc
    ) |>
    dplyr::arrange(adj.P.Val) |>
    dplyr::pull(symbol) |>
    unique()
  
  if (!is.null(top_n)) {
    seeds <- head(seeds, top_n)
  }
  
  seeds
}

get_hub_logical <- function(g, q = 0.90) {
  deg <- igraph::degree(g)
  thr <- stats::quantile(
    deg,
    probs = q,
    na.rm = TRUE
  )
  deg >= thr
}

plot_net_modules_hubs <- function(
    g,
    tag,
    output_dir,
    lay = NULL
) {
  if (is.null(lay)) {
    set.seed(1)
    lay <- igraph::layout_with_fr(g)
  }
  
  deg <- igraph::degree(g)
  deg_min <- min(deg)
  deg_max <- max(deg)
  
  if (deg_max > deg_min) {
    igraph::V(g)$size <- 4 +
      14 * (deg - deg_min) / (deg_max - deg_min)
  } else {
    igraph::V(g)$size <- 8
  }
  
  modules <- sort(unique(igraph::V(g)$module))
  pal <- grDevices::rainbow(length(modules))
  
  igraph::V(g)$color <- pal[
    match(igraph::V(g)$module, modules)
  ]
  
  hubs <- get_hub_logical(g, q = 0.90)
  igraph::V(g)$frame.color <- ifelse(
    hubs,
    "black",
    NA
  )
  
  ew <- igraph::E(g)$weight
  ew_min <- min(ew)
  ew_max <- max(ew)
  
  if (ew_max > ew_min) {
    igraph::E(g)$width <- 0.5 +
      2.5 * (ew - ew_min) / (ew_max - ew_min)
  } else {
    igraph::E(g)$width <- 1
  }
  
  igraph::E(g)$color <- "grey80"
  
  out_file <- file.path(
    output_dir,
    paste0(
      tag,
      "_network_modules.png"
    )
  )
  
  grDevices::png(
    out_file,
    width = 1400,
    height = 1000,
    res = 130,
    type = "cairo"
  )
  
  plot(
    g,
    layout = lay,
    vertex.label = NA,
    main = paste0(
      tag,
      " Reactome co-membership network (modules & hubs)"
    )
  )
  
  graphics::legend(
    "topleft",
    legend = paste("Module", modules),
    col = pal,
    pch = 16,
    bty = "n",
    cex = 0.8
  )
  
  graphics::legend(
    "bottomleft",
    legend = c("Hubs (top 10% degree)"),
    pch = 21,
    pt.bg = "white",
    col = "black",
    bty = "n",
    cex = 0.8
  )
  
  grDevices::dev.off()
}

enrich_module <- function(
    genes_symbol,
    universe_symbol,
    tag,
    mod_id,
    output_dir
) {
  if (length(genes_symbol) < 5) {
    return(invisible(NULL))
  }
  
  sig_m <- tt_gene |>
    filter(
      symbol %in% genes_symbol,
      !is.na(entrez)
    )
  
  uni_m <- tt_gene |>
    filter(
      symbol %in% universe_symbol,
      !is.na(entrez)
    )
  
  gene <- sig_m$entrez
  universe <- uni_m$entrez
  
  enricher(
    gene = gene,
    universe = universe,
    TERM2GENE = reactome_t2g,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
}

