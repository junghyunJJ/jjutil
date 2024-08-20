ll <- function(x){
  return(length(x))
}

summry_gprofiler2 <- function(enrich, plot = FALSE) {

  raw_res_enrich <- enrich$result %>%
    mutate(GeneRatio = paste0(intersection_size, "/", query_size)) %>%
    mutate(BgRatio = paste0(term_size, "/", effective_domain_size)) %>%
    mutate(n_GeneRatio = intersection_size / query_size)

  res_enrich <- raw_res_enrich %>%
    select(query, p_value, source, term_name, term_id, GeneRatio, BgRatio, n_GeneRatio, intersection_size, intersection) %>%
    mutate(term = paste0(term_name, "\n(", term_id, ")")) %>%
    select(term, query, p_value, source, everything())

  if (plot) {
    gprofiler2::gostplot(enrich, capped = TRUE, interactive = TRUE)
  }
  return(res_enrich)
}

updatejjutil <- function(x){
  devtools::install_github("junghyunJJ/jjutil")
}

mouse_to_human <- function(mouse_ids, input = "ENSEMBL", target = "SYMBOL") {
  mouse_geneid <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, mouse_ids, "ENTREZID", input)
  mapped <- AnnotationDbi::select(Orthology.eg.db::Orthology.eg.db, mouse_geneid, "Homo_sapiens", "Mus_musculus")
  human_target <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, as.character(mapped[, 2]), target, "ENTREZID")
  return(data.frame(mouse_input = mouse_ids, mapped, mouse_target = human_target[, 2]))
}

human_to_mouse <- function(human_ids, input = "ENSEMBL", target = "SYMBOL") {
  human_geneid <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, human_ids, "ENTREZID", input)
  mapped <- AnnotationDbi::select(Orthology.eg.db::Orthology.eg.db, human_geneid, "Mus_musculus", "Homo_sapiens")
  mouse_target <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, as.character(mapped[, 2]), target, "ENTREZID")
  return(data.frame(human_input = human_ids, mapped, mouse_target = mouse_target[, 2]))
}

# jjutil::mouse_to_human(mouse_ids = c("Mapk10", "Syt5", "Kank4", "Crb2"), input = "SYMBOL", target = "ENSEMBL")
# jjutil::human_to_mouse(c("ENSG00000109339", "ENSG00000129990", "ENSG00000132854", "ENSG00000148204"))


# hENSEMBL_to_mSYMBOL <- function(human_ids) {
#
#   # human_geneid <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, human_ids, "ENTREZID", input)
#   anno_human <- inner_join(
#     AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL),
#     AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egENSEMBL),
#     by = "gene_id"
#   )
#
#   human_geneid <- anno_human %>%
#     filter(ensembl_id %in% human_ids) %>%
#     pull(gene_id) %>%
#     unique
#
#   mapped <- AnnotationDbi::select(Orthology.eg.db::Orthology.eg.db, human_geneid, "Mus_musculus", "Homo_sapiens")
#   mapped <- mapped[!is.na(mapped$Mus_musculus), ]
#
#   # mouse_target <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, as.character(mapped[, 2]), target, "ENTREZID")
#   anno_mouse <- inner_join(
#     AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egSYMBOL),
#     AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egENSEMBL),
#     by = "gene_id"
#   )
#
#   mouse_target <- anno_mouse %>%
#     filter(gene_id %in% mapped$Mus_musculus) %>%
#     pull(symbol) %>%
#     unique
#
#   return(list(human_input = human_ids, mapped, mouse_target = mouse_target))
# }



hENSEMBL_to_mSYMBOL <- function(human_ids, removena = FALSE) {

  # human_geneid <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, human_ids, "ENTREZID", input)
  anno_human <- inner_join(
    AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL),
    AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egENSEMBL),
    by = "gene_id"
  )

  f_anno_human <- anno_human[match(human_ids, anno_human$ensembl_id), ]
  f_anno_human$ensembl_id <- human_ids
  colnames(f_anno_human) <- paste0("human_", colnames(f_anno_human))

  mapped <- AnnotationDbi::select(Orthology.eg.db::Orthology.eg.db, f_anno_human$human_gene_id, "Mus_musculus", "Homo_sapiens")

  # anno_inputid <- cbind(f_anno_human, mapped)
  anno_inputid <- cbind(f_anno_human, mouse_gene_id = mapped$Mus_musculus)

  # mouse_target <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, as.character(mapped[, 2]), target, "ENTREZID")
  anno_mouse <- AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egSYMBOL)
  f_anno_mouse <- anno_mouse[match(anno_inputid$mouse_gene_id, anno_mouse$gene_id), ]
  colnames(f_anno_mouse) <- paste0("mouse_", colnames(f_anno_mouse))
  # final_mouse <- cbind(anno_inputid, f_anno_mouse)
  final_mouse <- cbind(anno_inputid, mouse_symbol = f_anno_mouse$mouse_symbol)

  if (removena) {
    final_mouse <- final_mouse[apply(is.na(final_mouse), 1, sum) == 0, ]
  }


  return(final_mouse)
}

hSYMBOL_to_mSYMBOL <- function(human_ids, removena = FALSE) {

  # human_geneid <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, human_ids, "ENTREZID", input)
  anno_human <- inner_join(
    AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL),
    AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egENSEMBL),
    by = "gene_id"
  )

  f_anno_human <- anno_human[match(human_ids, anno_human$symbol), ]
  f_anno_human$ensembl_id <- human_ids
  colnames(f_anno_human) <- paste0("human_", colnames(f_anno_human))

  mapped <- AnnotationDbi::select(Orthology.eg.db::Orthology.eg.db, f_anno_human$human_gene_id, "Mus_musculus", "Homo_sapiens")

  # anno_inputid <- cbind(f_anno_human, mapped)
  anno_inputid <- cbind(f_anno_human, mouse_gene_id = mapped$Mus_musculus)

  # mouse_target <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, as.character(mapped[, 2]), target, "ENTREZID")
  anno_mouse <- AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egSYMBOL)
  f_anno_mouse <- anno_mouse[match(anno_inputid$mouse_gene_id, anno_mouse$gene_id), ]
  colnames(f_anno_mouse) <- paste0("mouse_", colnames(f_anno_mouse))
  # final_mouse <- cbind(anno_inputid, f_anno_mouse)
  final_mouse <- cbind(anno_inputid, mouse_symbol = f_anno_mouse$mouse_symbol)

  if (removena) {
    final_mouse <- final_mouse[apply(is.na(final_mouse), 1, sum) == 0, ]
  }

  return(final_mouse)
}
