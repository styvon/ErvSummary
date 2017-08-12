#!/usr/bin/Rscript

cat("Loading packages...\n")
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(tidyr))
suppressPackageStartupMessages(require(reshape2))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(RColorBrewer))

weighted.sd <- function(x, wt){
  wm <- weighted.mean(x,wt,na.rm=T)
  out <- sum(wt * (x - wm)* (x - wm), na.rm=T)
}

plotGen <- function(filedir, nmer, fign=1, choice=1){
  df <- read.table(filedir, header=T, fill=T)
  pat1 <- c(0:(nmer-1))
  pat2 <- c((nmer-1):0)
  patterns <- paste0(pat1, pat2)
  
  sum_motif <- df %>%
    select(pattern, nMotifs) %>%
    group_by(pattern) %>%
    summarise(sum=sum(as.numeric(nMotifs)))
  step <- (4^(nrow(sum_motif) -1)) # number of subtypes per mutation type
  for(i in seq_along(sum_motif$sum)) {
    df$wt[(step * (i-1)+1) : (step * i) ] = df$nMotifs[(step * (i-1)+1) : (step * i) ]/ sum_motif$sum[i]
  }
  
  summary_all <-  df %>%
    group_by(pattern) %>%
    summarise(ll = sum(nERVs * log(ERV_rel_rate), na.rm=T) ) %>%
    mutate(aic = 2*step*6 - 2*ll)
  # sd = weighted.sd(ERV_rel_rate, wt) too small returned NA_real
  
  summary_grouped <-  df %>%
    group_by(pattern, mtype) %>%
    summarise(ll = sum(nERVs * log(ERV_rel_rate), na.rm=T) ) %>%
    mutate(aic = 2*step - 2*ll)
  # sd = weighted.sd(ERV_rel_rate, wt) too small for some subtypes
  
  aic_grouped <- summary_grouped %>%
    group_by(mtype) %>%
    mutate(.min = min(aic), .rank = rank(aic),pattern_num = factor(pattern,labels= patterns ) ) %>%
    ungroup() %>%
    arrange(mtype, aic) %>%
    mutate(.r = row_number(), .diff = aic - .min)
  
  if(choice==1){
    plot <- ggplot(aic_grouped,aes(x = .r, y = aic, color = aic)) + 
      scale_fill_gradient2() +
      geom_point(stat="identity", size=5) +
      geom_text(aes(label = formatC(aic, format="e", digits=3)), vjust = 2) +
      ggtitle(paste0("Fig.", fign, " Comparison of AIC by subtype (", nmer, "-mer)")) +
      facet_wrap(~mtype,shrink=FALSE, scales = "free", ncol=3) +
      scale_x_continuous(  
        breaks = aic_grouped$.r,     
        labels = aic_grouped$pattern_num
      )
  }

  if(choice ==2){
    plot <- ggplot(aic_grouped,aes(x = .r, y = .diff, fill = aic)) + 
      scale_fill_gradient2() +
      geom_bar(stat="identity") + 
      ggtitle(paste0("Fig.", fign, " Comparison of AIC difference by subtype (", nmer, "-mer)")) +
      facet_wrap(~mtype,shrink=FALSE, scales = "free", ncol=6) +
      scale_x_continuous(  
        breaks = aic_grouped$.r,     
        labels = aic_grouped$pattern_num
      )
  }
  
  if(choice ==3){
    suppressMessages( 
      aic_hist <- aic_grouped %>%
      dcast(pattern_num+mtype+aic+.diff ~ .rank, length) %>%
      select(-mtype, -aic, -.diff) %>%
      group_by(pattern_num) %>%
      summarise_all(sum)
    )
    
    plot <- suppressWarnings(
      ggplot(aic_grouped, aes(x = as.factor(.rank))) +
      geom_histogram(stat="count") +
      facet_wrap(~pattern_num,shrink=FALSE, ncol=3) +
      ggtitle(paste0("Fig.", fign, " Distributions of AIC ranks (", nmer, "-mer)"))
    )
  }
  
  plot
  
  
}

