update_rctd <- function(res_RCDT_df) {

  # NOTE!! for singlet results from RCTD, we only used "first_type"!
  res_RCDT_df[res_RCDT_df$spot_class == "singlet", "second_type"] <- res_RCDT_df[res_RCDT_df$spot_class == "singlet", "first_type"]

  # NOTE!! RCTD results provide first_type and second_type information, But we need just cell type information in the doublet. The cell type ordering is not important.
  update_res <- pbapply::pblapply(seq_len(nrow(res_RCDT_df)), function(i) {
    sel <- res_RCDT_df[i, ]

    if (sel$spot_class == "singlet") {
      sel
    } else {
      celltype1 <- sel$first_type %>% as.character()
      celltype2 <- sel$second_type %>% as.character()
      celltype <- c(celltype1, celltype2) %>% sort()

      sel$first_type <- celltype[1]
      sel$second_type <- celltype[2]
      sel
    }
  })
  update_res <- do.call(rbind.data.frame, update_res)
  return(update_res)
}

summary_rctd <- function(myRCTD) {

  # counts and coord of ST
  coords <- myRCTD@spatialRNA@coords

  # results of RCTD
  rawres_RCDT_df <- myRCTD@results$results_df[, 1:3]
  all.equal(rownames(coords), rownames(rawres_RCDT_df))
  # rawres_RCDT_df %>% head

  # update RCTD resutls
  res_RCDT_df <- update_rctd(rawres_RCDT_df)
  res_RCDT_df %>% head
  res_RCDT_df %>% tail

  # geterate data.frame for CellNeighborEX analysis
  res_RCDT_df <- data.frame(
    spot_class = res_RCDT_df$spot_class,
    barcode = rownames(res_RCDT_df),
    first_type = res_RCDT_df$first_type,
    second_type = res_RCDT_df$second_type,
    celltype1 = res_RCDT_df$first_type,
    celltype2 = res_RCDT_df$second_type
  )
  # res_RCDT_df %>% head

  all.equal(levels(res_RCDT_df$first_type), levels(res_RCDT_df$second_type))
  celltypes <- as.character(sort(unique(c(res_RCDT_df$first_type, res_RCDT_df$second_type))))

  # make factor for celltypes
  levels(factor(celltypes))
  levels(factor(celltypes, labels = seq_len(length(celltypes))))

  res_RCDT_df$first_type <- factor(as.character(res_RCDT_df$first_type), levels = levels(factor(celltypes)), labels = seq_len(length(celltypes)))
  res_RCDT_df$second_type <- factor(as.character(res_RCDT_df$second_type), levels = levels(factor(celltypes)), labels = seq_len(length(celltypes)))
  # res_RCDT_df$first_type %>% levels
  # res_RCDT_df$second_type %>% levels

  res_RCDT_df$celltype1 <- factor(as.character(res_RCDT_df$celltype1), levels = levels(factor(celltypes)))
  res_RCDT_df$celltype2 <- factor(as.character(res_RCDT_df$celltype2), levels = levels(factor(celltypes)))
  # res_RCDT_df$celltype1 %>% levels
  # res_RCDT_df$celltype2 %>% levels
  # res_RCDT_df %>% View

  # cell type proportion
  res_RCDT_prop <- as.data.frame(as.matrix(myRCTD@results$weights_doublet))
  res_RCDT_prop %>% head
  colnames(res_RCDT_prop) <- c("prop1", "prop2")
  res_RCDT_prop %>% head()
  all.equal(res_RCDT_df$barcode, rownames(res_RCDT_prop))
  all.equal(rownames(coords), rownames(res_RCDT_prop))
  summ_rctd <- cbind(res_RCDT_df, coords, res_RCDT_prop)

  print(head(summ_rctd))

  return(summ_rctd)
}

summary_st <- function(myRCTD, barcode) {
  counts <- myRCTD@spatialRNA@counts[, barcode]
  dat_seurat <- Seurat::CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
  dat_seurat <- Seurat::NormalizeData(dat_seurat, verbose = FALSE)

  df_log_data <- dat_seurat[["RNA"]]$data %>%
    as.matrix() %>%
    as.data.frame()

  df_gene_name <- rownames(df_log_data)
  df_cell_id <- colnames(df_log_data)

  return(list(df_log_data = df_log_data, df_cell_id = df_cell_id, df_gene_name = df_gene_name))
}
