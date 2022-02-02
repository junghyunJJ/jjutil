jj <- function() { message("Hello, my name JJ!") }


# tmp_anno1 <- inner_join(AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egCHR),
#                         AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL), by = "gene_id")
# tmp_anno2 <- inner_join(tmp_anno1,
#                         AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egENSEMBL), by = "gene_id")
# anno_human <- inner_join(tmp_anno2,
#                          AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egGENENAME), by = "gene_id") %>%
#   arrange(chromosome, symbol) %>% as_tibble()
#
# save(anno_human, file = "data/anno_human.RData")




