library(tidyverse)
# qplot(1,1) + 
#   annotation_compass('testN') + 
#   annotation_compass('testE','E') + 
#   annotation_compass('testSW','SW') + 
#   annotation_compass('testW','W')

annotation_compass <- function(label,
                               position = c('N','NE','E','SE','S','SW','W','NW'),
                               padding = grid::unit(c(0.5,0.5),"line"), ...){
  position <- match.arg(position)
  x <- switch (position,
               N = 0.5,
               NE = 1,
               E = 1,
               SE = 1,
               S = 0.5, 
               SW = 0,
               W = 0, 
               NW = 0
  )
  y <- switch (position,
               N = 1,
               NE = 1,
               E = 0.5,
               SE = 0,
               S = 0, 
               SW = 0,
               W = 0.5, 
               NW = 1
  )
  hjust <- switch (position,
                   N = 0.5,
                   NE = 1,
                   E = 1,
                   SE = 1,
                   S = 0.5, 
                   SW = 0,
                   W = 0, 
                   NW = 0
  )
  vjust <- switch (position,
                   N = 1,
                   NE = 1,
                   E = 0.5,
                   SE = 0,
                   S = 0, 
                   SW = 0,
                   W = 0.5, 
                   NW = 1
  )
  f1 <- switch (position,
                N = 0,
                NE = -1,
                E = -1,
                SE = -1,
                S = 0, 
                SW = 1,
                W = 1, 
                NW = 1
  )
  f2 <- switch (position,
                N = -1,
                NE = -1,
                E = 0,
                SE = 1,
                S = 1, 
                SW = 1,
                W = 0, 
                NW = -1
  )
  annotation_custom(grid::textGrob(label, 
                                   x=grid::unit(x,"npc") + f1*padding[1] , 
                                   y=grid::unit(y,"npc") + f2*padding[2],
                                   hjust=hjust,vjust=vjust, ...))
}



inflation <- function(p, method = c("regression","median")){
  out <- list()
  
  if(method == "median"){
    chisq <- qchisq(1-p,1)
    lambda<-median(chisq)/qchisq(0.5,1)
    out$estimate <- lambda # lambda 
    out$se <- NA
    out
  }else{
    data <- qchisq(p, 1, lower.tail = FALSE)
    data <- sort(data)
    ppoi <- ppoints(data) #Generates the sequence of probability points
    ppoi <- sort(qchisq(ppoi, df = 1, lower.tail = FALSE))
    s <- summary(lm(data ~ 0 + ppoi))$coeff
    out$estimate <- s[1, 1] # lambda 
    out$se <- s[1, 2]
    out
    # median method
    # out$estimate <- median(data, na.rm = TRUE)/qchisq(0.5, 1)
  }
}

gg_qqplot <- function(p, lamda, ci = 0.95, title=NULL){
  n  <- length(p)
  df <- data.frame(
    observed = -log10(sort(p)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  gg_df <- ggplot(df) +
    theme_bw() +
    geom_point(aes(expected, observed), shape = 1, size = 2) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2) +
    geom_line(aes(expected, clower), linetype = 2) +
    xlab(log10Pe) +
    ylab(log10Po)+
    theme(axis.text = element_text(color = "black")) +
    annotation_compass(paste0("\u03BB"," = ",round(lamda,2)))
  
  if(is.null(title)){
    gg_df
  }else{
    gg_df + ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5))
  }
}



gg_manh<-function(dat, sigth = 5e-08, X_chrm = FALSE){
  sig <- dat %>% filter(P < sigth) %>% dplyr::select(SNP) %>% pull
  
  raw_ggdat <- dat %>% group_by(CHR) %>% 
    count %>% 
    ungroup %>% 
    dplyr::rename(chr_len=n) %>% 
    mutate(tot=cumsum(chr_len)-chr_len) 
  
  if(X_chrm == TRUE){
    n_CHR <- 23
  }else{
    n_CHR <- 22
  }
  
  ggdat <- dat %>% 
    mutate(bp = unlist(lapply(1:n_CHR, function(i){seq(1,raw_ggdat$chr_len[i])}))) %>% 
    left_join(raw_ggdat, ., by=c("CHR"="CHR")) %>%
    dplyr::select(-chr_len) %>% 
    arrange(CHR, bp) %>%
    mutate( bpcum=bp+tot) %>% 
    mutate( is_highlight=ifelse(SNP %in% sig, "yes", "no"))
  
  axisdf <- ggdat %>% dplyr::select(CHR, bpcum) %>% group_by(CHR) %>% summarize(center=( max(bpcum) + min(bpcum) ) / 2 ,.groups = 'drop')
  maxp <- max(-log10(dat$P))
  nbreaks <- ifelse(maxp < 10, 1,2)
  
  gg <- ggplot(ggdat, aes(x=bpcum, y=-log10(P))) +
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=0.5) +
    geom_hline(yintercept = -log10(1e-05), color = "orange",size=0.5, linetype = "dashed") +
    geom_hline(yintercept = -log10(sigth), color = "red",size=0.5, linetype = "dashed") +
    theme_bw() +
    geom_point(data=subset(ggdat, is_highlight=="yes"), color="red", size=2) +
    scale_color_manual(values = rep(c("black", "grey60"), n_CHR )) +
    scale_x_continuous(expand = c(0.01,0.01),label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0,0),limits=c(0,maxp + (maxp * 0.1)), breaks = seq(0,maxp+(maxp * 0.1), nbreaks)) +
    xlab("") +
    theme(
      axis.text = element_text(colour = "black"),
      legend.position="none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  return(gg)
}

# ggrepel::geom_text_repel(data=anno_sig_dat,aes(label = gene), size=3,
#                          box.padding = unit(0.35, "lines"),
#                          point.padding = unit(0.3, "lines")) 