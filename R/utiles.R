ll <- function(x){
  return(length(x))
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


hENSEMBL_to_mSYMBOL <- function(human_ids) {

  # human_geneid <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, human_ids, "ENTREZID", input)
  anno_human <- inner_join(
    AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL),
    AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egENSEMBL),
    by = "gene_id"
  )

  human_geneid <- anno_human %>%
    filter(ensembl_id %in% human_ids) %>%
    pull(gene_id) %>%
    unique

  mapped <- AnnotationDbi::select(Orthology.eg.db::Orthology.eg.db, human_geneid, "Mus_musculus", "Homo_sapiens")
  mapped <- mapped[!is.na(mapped$Mus_musculus), ]

  # mouse_target <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, as.character(mapped[, 2]), target, "ENTREZID")
  anno_mouse <- inner_join(
    AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egSYMBOL),
    AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egENSEMBL),
    by = "gene_id"
  )

  mouse_target <- anno_mouse %>%
    filter(gene_id %in% mapped$Mus_musculus) %>%
    pull(symbol) %>%
    unique

  return(list(human_input = human_ids, mapped, mouse_target = mouse_target))
}
