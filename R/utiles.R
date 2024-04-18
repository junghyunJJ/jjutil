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



