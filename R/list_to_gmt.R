# read.gmtfunctions in qusage r package
read.gmt <- function(file){
  if (!grepl("\\.gmt$", file)[1]) {
    stop("Pathway information must be a .gmt file")
  }
  geneSetDB = readLines(file)
  geneSetDB = strsplit(geneSetDB, "\t")
  names(geneSetDB) = sapply(geneSetDB, "[", 1)
  geneSetDB = lapply(geneSetDB, "[", -1:-2)
  geneSetDB = lapply(geneSetDB, function(x) {
    x[which(x != "")]
  })
  return(geneSetDB)
}

# write.gmt <- function(object, fname){
#   if (class(object) != "list") stop("object should be of class 'list'")
#   if(file.exists(fname)) unlink(fname)
#   for (iElement in 1:length(object)){
#     write.table(t(c(make.names(rep(names(object)[iElement],2)),object[[iElement]])),
#                 sep="\t",quote=FALSE,
#                 file=fname,append=TRUE,col.names=FALSE,row.names=FALSE)
#   }
# }

list_to_gmt <- function(listgeneset = listgeneset, file){
  if(is.null(file)) stop("'quote' must be 'TRUE', 'FALSE' or numeric")
  gs_len <- length(listgeneset)
  gs_maxlen <- max(unlist(lapply(listgeneset,length)))
  raw_gmt <- data.frame(matrix(nrow = gs_maxlen,ncol = gs_len))

  for(i in 1:gs_len) raw_gmt[,i] <- c(listgeneset[[i]],rep(NA,(gs_maxlen - length(listgeneset[[i]]))))

  colnames(raw_gmt) <- names(listgeneset)
  save_gmt <- data.frame("NA",t(raw_gmt))
  write.table(save_gmt,file = paste0(file,".gmt"),row.names = T,col.names = F,quote = F,sep = "\t",na = "")
  cat(paste0(gs_len," gene set were saved (",paste0(file,".gmt"),")\n"))
}

list_to_starfysh <- function(listgeneset = listgeneset, file){
  if(is.null(file)) stop("'quote' must be 'TRUE', 'FALSE' or numeric")
  gs_len <- length(listgeneset)
  gs_maxlen <- max(unlist(lapply(listgeneset,length)))
  raw_gmt <- data.frame(matrix(nrow = gs_maxlen,ncol = gs_len))

  for(i in 1:gs_len) raw_gmt[,i] <- c(listgeneset[[i]],rep(NA,(gs_maxlen - length(listgeneset[[i]]))))

  colnames(raw_gmt) <- names(listgeneset)
  # save_gmt <- data.frame("NA",t(raw_gmt))
  write.csv(raw_gmt,file = paste0(file,".csv"),na = "", row.names = FALSE, quote = FALSE)
  cat(paste0(gs_len," gene set were saved (",paste0(file,".csv"),")\n"))
}

list_to_DAVID <- function(listgeneset = listgeneset, file){
  if(is.null(file)) stop("'quote' must be 'TRUE', 'FALSE' or numeric")
  gs_len <- length(listgeneset)
  gs_maxlen <- max(unlist(lapply(listgeneset,length)))
  raw_gmt <- data.frame(matrix(nrow = gs_maxlen,ncol = gs_len))

  for(i in 1:gs_len) raw_gmt[,i] <- c(listgeneset[[i]],rep(NA,(gs_maxlen - length(listgeneset[[i]]))))

  colnames(raw_gmt) <- names(listgeneset)
  # save_gmt <- data.frame("NA",t(raw_gmt))
  write.table(raw_gmt, file = paste0(file,".txt"), na = "", sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste0(gs_len," gene set were saved (",paste0(file,".txt"),")\n"))
}


h <- function(x) return(x[1:5,1:5])
hh <- function(x) return(x[1:10,1:10])

e <- function(x) return(x[(nrow(x)-5):nrow(x),(ncol(x)-5):ncol(x)])

dd <- function(x) return(dim(x))
ll <- function(x) return(length(x))

# GCSC
# https://github.com/ksiewert/GCSC

# gene_universe <- readLines("~/tools/GCSC/gene_universe.txt")
# save(gene_universe, file = "data/gene_universe.RData")
list_to_GCSC <- function(glist, inputtype=c("gene_id","symbol", "ensembl_id"), gene_universe=gene_universe, th=10, file){

  cat(length(glist),"gene set\n")
  if(is.na(inputtype) | sum(inputtype %in%c("gene_id","symbol", "ensembl_id")) !=1 ){
    message("Please set input type among gene_id, symbol, and ensembl_id\n")
  }
  # load annotation
  anno <- inner_join(AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL),
                     AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egENSEMBL), by = "gene_id")

  if(inputtype == "symbol"){
    cat("sym -> ensembl_id\n")
    ensembl_list <- lapply(glist, function(x){
      anno %>% filter(anno$symbol %in% x) %>% dplyr::select(ensembl_id) %>% pull
    })
  }else if(inputtype == "gene_id"){
    cat("gene_id -> ensembl_id\n")
    ensembl_list <- lapply(glist, function(x){
      anno %>% filter(anno$gene_id %in% x) %>% dplyr::select(ensembl_id) %>% pull
    })
  }else if(inputtype == "ensembl_id"){
    # just reshape
    ensembl_list <- glist
  }

  # Transformation -  GCSC format
  save_dat <- data.frame(gene_universe)
  raw_res <- lapply(ensembl_list, function(xx){
    data.frame(t(as.numeric(!is.na(match(gene_universe,xx)))))
  }) %>% rbindlist

  colnames(raw_res) <- gene_universe
  res <- cbind(names(glist),data.frame(raw_res))
  colnames(res)[1] <- ""

  # cat(apply(res[,-1], 1, sum),"\n")

  idx <- which(apply(res[,-1], 1, sum) < th)
  if(length(idx) > 0){
    cat(length(idx),"genesets filtered (gene set threshold = 10)\n")
    res <- res[-idx,]
  }
  # cat(apply(res[,-1], 1, sum),"\n")
  write.csv(res,paste0(file,".cvs"), quote = F, row.names = F)
  cat(paste0((length(glist) - length(idx))," gene set were saved (",paste0(file,".cvs"),")\n"))
}


list_to_GCSC_GTEXv8 <- function(glist, gene_universe=gene_universe, th=10, out){

  cat(length(glist),"gene set\n")
  anno <- jjutil::anno.GTExv8.gencode.v26.annotation

  cat("sym -> ensembl_id\n")
  ensembl_list <- lapply(glist, function(x){
    anno %>% filter(anno$symbol %in% x) %>% dplyr::select(ID) %>% pull
  })

  # Transformation -  GCSC format
  save_dat <- data.frame(gene_universe)
  raw_res <- lapply(ensembl_list, function(xx){
    data.frame(t(as.numeric(!is.na(match(gene_universe,xx)))))
  }) %>% rbindlist

  colnames(raw_res) <- gene_universe
  res <- cbind(names(glist),data.frame(raw_res))
  colnames(res)[1] <- ""

  # cat(apply(res[,-1], 1, sum),"\n")

  idx <- which(apply(res[,-1], 1, sum) < th)
  if(length(idx) > 0){
    cat(length(idx),"genesets filtered (gene set threshold = 10)\n")
    res <- res[-idx,]
  }
  # cat(apply(res[,-1], 1, sum),"\n")
  write.csv(res,paste0(out,".cvs"), quote = F, row.names = F)
  cat(paste0(length(glist)," gene set were saved (",paste0(out,".cvs"),")\n"))
}



list_to_mesc <- function(listgeneset, file) {

  list.ensemble <- lapply(listgeneset, function(sel){
    anno.GTExv8.gencode.v26.annotation %>%
      filter(symbol %in% sel) %>%
      dplyr::select(ID) %>%
      unique %>%
      pull
  })

  if (is.null(file))
    stop("'quote' must be 'TRUE', 'FALSE' or numeric")

  gs_len <- length(list.ensemble)
  gs_maxlen <- max(unlist(lapply(list.ensemble, length)))
  raw_gmt <- data.frame(matrix(nrow = gs_maxlen, ncol = gs_len))
  for (i in 1:gs_len) raw_gmt[, i] <- c(list.ensemble[[i]], rep(NA, (gs_maxlen - length(list.ensemble[[i]]))))
  colnames(raw_gmt) <- names(list.ensemble)
  save_gmt <- data.frame(t(raw_gmt))
  write.table(save_gmt, file = paste0(file, ".mesc"), row.names = T,
              col.names = F, quote = F, sep = "\t", na = "")
  cat(paste0(gs_len, " gene set were saved (", paste0(file,".mesc"), ")\n"))
}

