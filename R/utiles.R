ll <- function(x){
  return(length(x))
}

updatejjutil <- function(x){
  devtools::install_github("junghyunJJ/jjutil")
}

mouse_to_human <- function(mouse_ids, input = "ENSEMBL", target = "SYMBOL", horg = org.Hs.eg.db, morg = org.Mm.eg.db, orth = Orthology.eg.db) {
  mouse_geneid <- AnnotationDbi::mapIds(morg, mouse_ids, "ENTREZID", input)
  mapped <- AnnotationDbi::select(orth, mouse_geneid, "Homo_sapiens", "Mus_musculus")
  human_target <- AnnotationDbi::select(horg, as.character(mapped[, 2]), target, "ENTREZID")
  return(data.frame(mouse_input = mouse_ids, mapped, mouse_target = human_target[, 2]))
}

human_to_mouse <- function(human_ids, input = "ENSEMBL", target = "SYMBOL", horg = org.Hs.eg.db, morg = org.Mm.eg.db, orth = Orthology.eg.db) {
  human_geneid <- AnnotationDbi::mapIds(horg, human_ids, "ENTREZID", input)
  mapped <- AnnotationDbi::select(orth, human_geneid, "Mus_musculus", "Homo_sapiens")
  mouse_target <- AnnotationDbi::select(morg, as.character(mapped[, 2]), target, "ENTREZID")
  return(data.frame(human_input = human_ids, mapped, mouse_target = mouse_target[, 2]))
}
