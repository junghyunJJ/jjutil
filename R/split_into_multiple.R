
split_into_multiple <- function(data,pattern){
  for(i in 1:nrow(data)){
    if( i == 1){
      if (!require(plyr)) install.packages("plyr") 
      if (!require(stringr)) install.packages("stringr") 
      if (!require(tidyverse)) install.packages("tidyverse") 
      
      
      progress.bar <- plyr::create_progress_bar("text")
      progress.bar$init(nrow(data))
      t_data <- c()
    }
    imsi <- stringr::str_split_fixed(as.character(data[i,]), pattern, n = Inf) %>% 
      t 
    
    for(j in 1:ncol(imsi)){
      if(sum(imsi[,j] == "") > 0){
        imsi[,j] <- imsi[1,j]
      }  
    }
    
    t_data <- rbind(t_data,imsi)
    progress.bar$step()
  }
  colnames(t_data) <- colnames(data)
  t_data <- as.tibble(t_data)
  return(t_data)
}







