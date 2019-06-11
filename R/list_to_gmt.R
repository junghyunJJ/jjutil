
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

h <- function(x) return(x[1:5,1:5])

e <- function(x) return(x[(nrow(x)-5):nrow(x),(ncol(x)-5):ncol(x)])

d <- function(x) return(dim(x))
