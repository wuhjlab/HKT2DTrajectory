
# This code is for directivity test for significant pairs.
# df_stat is the output of calculating the statistics of disease pairs.

df_statistics<-df_stat%>%filter(p_bonferroni<0.05)
dat_dia<- dat %>%
  group_by(id) %>%
  filter(any(icd_main %in% df_statistics$exposure)) %>%
  ungroup()

unique_pairs <- df_statistics %>%
  select(exposure, outcome) %>%
  distinct()%>%split(.,seq(nrow(.)))

dat_date<-dat%>%select(icd_main,date,id)

res_directivity<-mclapply(unique_pairs,FUN=function(x){
  dat_dire <- dat_date %>%
    group_by(id) %>%
    filter(any(icd_main == x$exposure) & any(icd_main == x$outcome)) %>%
    ungroup() %>%
    filter(icd_main %in% c(x$exposure, x$outcome)) %>%
    distinct() %>%
    group_by(id, icd_main) %>%
    filter(date == min(date)) %>%
    ungroup() %>%
    arrange(id, date) %>%
    group_by(id) %>%
    summarise(order = ifelse(date[1]==date[2], "same", paste0(first(icd_main), "_", last(icd_main)))) %>%
    ungroup() %>%
    count(order)%>%
    mutate(#directivity=ifelse(order==paste0(x$exposure,"_",x$outcome),T,F),
      directivity_pct=n/sum(n))
  
  dat_dire$p_directivity<-NA
  for(i in 1:nrow(dat_dire)){
    dat_dire$p_directivity[i]<-binom.test(dat_dire$n[i], sum(dat_dire$n),0.5)$p.value
    
  }
  
  dat_dire$directivity<-ifelse(dat_dire$directivity_pct>0.5&dat_dire$p_directivity<0.001,T,F)
  return(dat_dire)
},mc.cores = thread)

res_directivity<-bind_rows(res_directivity)%>%
  separate(order,c("exposure",'outcome'),"\\_")
res_directivity<-left_join(df_statistics,res_directivity%>%
                             select(exposure,outcome,p_directivity,directivity_pct,directivity),by=c('exposure','outcome'))

