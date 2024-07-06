## This code is for calculating the statistics of disease pairs.

## input
# df_patients: A dataframe containing patient information, including patient IDs, sex, age, diagnosis dates, and the primary diagnosis ICD codes.
# dat: A dataframe containing the complete dataset of patient records, including IDs, diagnosis dates, and ICD codes.
# N: The number of control patients to match for each exposed patient.

## output
# outcome: The ICD code of the disease.
# expo_freq: The frequency of the disease in the exposed group.
# con_freq_avg: The average frequency of the disease in the control group.
# p: The proportion of control groups where the disease frequency in the exposed group is higher.
# p_BD: The p-value from the binomial test.
# p_BD_bonferroni: The Bonferroni-adjusted p-value for multiple comparisons.
# RR: The relative risk, calculated as the frequency of the disease in the exposed group divided by the average frequency in the control group.
# exposure: The primary diagnosis ICD code for the exposed group being analyzed.



diseases_summary_count <- function(df_patients, dat, N, seed = NULL, match_paralle = F) {

  year_range <- function(date) {
    year_range <- c(as.numeric(str_sub(date, 1, 4)) - 2, as.numeric(str_sub(date, 1, 4)) + 2)
    year_range <- c(paste0(year_range[1], str_sub(date, 5, 6)), paste0(year_range[2], str_sub(date, 5, 6)))
    year_range <- sapply(year_range, as.numeric)
    year_range <- data.frame(from = year_range[1], to = year_range[2])
    return(year_range)
  }
  expo_id <- unique(df_patients$id)
  yr <- rbindlist(pblapply(df_patients$date, year_range))
  df_patients$from <- yr$from
  df_patients$to <- yr$to
  expo <- dat %>% filter(id %in% unique(df_patients$id)) %>% left_join(df_patients %>% select(id, exposed_date = date), by = 'id') %>% filter(date > exposed_date) %>% left_join(df_patients %>% select(id, from, to), by = 'id')
  if (nrow(expo) == 0) return()
  count_expo <- expo %>% group_by(icd_main) %>% summarise(freq = n()) %>% mutate(group = "exposed")
  match_dat <- df_patients %>% select(id, sex, age, from, to, dm_years) %>% distinct()

  match_controls <- function(y, dat, expo_id, N, seed) {
    dat_match <- dat[dat$sex == y$sex & (dat$age - 2 <= y$age & dat$age + 2 >= y$age) & (dat$dm_years - 2 <= y$dm_years & dat$dm_years + 2 >= y$dm_years), ] %>% filter(between(date, y$from, y$to)) %>% select(id, date, icd_main) %>% distinct() %>% filter(!id %in% expo_id) %>% arrange(id, date) %>% group_by(id) %>% filter(!duplicated(id))
    id <- data.frame(id = dat_match$id)
    match_con <- id %>% left_join(dat %>% select(id, date, icd_main), by = 'id') %>% left_join(dat_match %>% select(id, exposed_date = date), by = 'id') %>% filter(date > exposed_date) %>% distinct()
    if (nrow(match_con) == 0) return()
    id_con <- data.frame(id = sample(unique(match_con$id), N, replace = T), group = 1:N)
    count_con <- id_con %>% left_join(match_con %>% select(id, icd_main), by = 'id', relationship = "many-to-many") %>% select(group, icd_main)
    return(count_con)
  }
  match_dat <- split(match_dat, seq(nrow(match_dat)))
  dat_con <- pblapply(match_dat, FUN = function(x) match_controls(x, dat, expo_id, N, seed))
  if (length(dat_con) == 0) return() else dat_con <- compact(dat_con)
  if (sum(unlist(lapply(dat_con, nrow))) > 1e+9) {
    part_size <- ceiling(length(dat_con) / 2)
    l <- 1:2
    names(l) <- 1:2
    count_con1 <- mclapply(1:2, function(i) {
      if (i != 384) {
        count_con <- rbindlist(dat_con[part_size * (i - 1) + 1:part_size * 1]) %>% group_by(group, icd_main) %>% summarise(freq = n())
      } else {
        count_con <- rbindlist(dat_con[part_size * (i - 1) + 1:sn]) %>% group_by(group, icd_main) %>% summarise(freq = n())
      }
      return(count_con)
    }, mc.cores = 384, mc.preschedule = T)
    count_con <- merge(as.data.table(count_con1[[1]]), as.data.table(count_con1[[2]]), by = c("group", "icd_main"), suffixes = c(".dt1", ".dt2"))
    count_con[, freq := freq.dt1 + freq.dt2]
    count_con <- count_con[, .(group, icd_main, freq)]
  } else {
    count_con <- rbindlist(compact(dat_con)) %>% group_by(group, icd_main) %>% summarise(freq = n())
  }
  icd_e <- df_patients$icd_main1[1]
  all_icds <- unique(count_expo$icd_main)
  count_con <- count_con %>% ungroup() %>% group_by(group) %>% tidyr::complete(icd_main = all_icds, fill = list(freq = 0)) %>% ungroup()
  count_con_avg <- count_con %>% mutate(pr = freq / N) %>% group_by(icd_main) %>% summarise(con_avg = mean(freq), pr = sum(pr) / length(unique(df_patients$id)), sum_count = sum(freq))
  count_expo <- count_expo %>% left_join(count_con_avg, by = 'icd_main')
  num_test <- length(unique(df_patients$id))
  p_BD <- unlist(pblapply(split(count_expo, count_expo$icd_main), FUN = function(x) {
    print(x)
    binom.test(x$freq[1], num_test, x$pr[1], alternative = "greater")$p.value
  }))
  count_con_bd <- data.frame(p_BD = p_BD, icd_main = names(p_BD))
  count_con_bd$p_BD_bonferroni <- p.adjust(count_con_bd$p_BD)
  count_con_p <- count_con %>% left_join(count_expo %>% select(icd_main, freq_expo = freq), by = 'icd_main') %>% mutate(P_logit = ifelse(freq_expo < freq, T, F)) %>% group_by(icd_main) %>% summarise(p = sum(P_logit) / N)
  RR <- data.frame(icd_main = all_icds) %>% left_join(count_expo %>% select(icd_main, expo_freq = freq), by = 'icd_main') %>% left_join(count_con_avg %>% select(icd_main, con_freq_avg = con_avg), by = 'icd_main') %>% left_join(count_con_p %>% select(icd_main, p), by = 'icd_main') %>% left_join(count_con_bd, by = 'icd_main') %>% mutate(RR = expo_freq / con_freq_avg) %>% rename('outcome' = icd_main) %>% mutate(exposure = df_patients$icd_main[1])
  return(RR)
}

dat <- as.data.frame(dat)
if (!is.null(disease_icd)) {
  dat_cases <- dat %>% mutate(date = as.numeric(date)) %>% filter(icd_main %in% disease_icd) %>% mutate(id1 = id, icd_main1 = icd_main) %>% group_by(id, icd_main) %>% filter(date == min(date)) %>% distinct() %>% ungroup() %>% group_by(icd_main) %>% nest()
} else {
  dat_cases <- dat %>% mutate(date = as.numeric(date)) %>% mutate(id1 = id, icd_main1 = icd_main) %>% group_by(id, icd_main) %>% filter(date == min(date)) %>% distinct() %>% ungroup() %>% group_by(icd_main) %>% nest()
}

ll <- dat_cases$icd_main
dat_cases <- dat_cases$data
names(dat_cases) <- ll
dat_cases <- mclapply(dat_cases, function(x) as.data.frame(x), mc.cores = 384)
df_stat <- mclapply(dat_cases, function(x) {
  diseases_summary_count(x, dat = dat, N = N)
}, mc.cores = 384)

df_stat <- as.data.frame(bind_rows(df_stat))
df_stat$p <- ifelse(df_stat$p < 1e-300, 1e-300, df_stat$p)
df_stat$p_BD <- ifelse(df_stat$p_BD < 1e-300, 1e-300, df_stat$p_BD)
if (is.null(disease_icd)) {
  df_stat$p_bonferroni <- p.adjust(df_stat$p, 'bonferroni')
  df_stat$p_BD_bonferroni <- p.adjust(df_stat$p_BD, 'bonferroni')
} else {
  dat_stat <- df_stat %>% filter(outcome %in% c('D', disease_icd))
  dat_stat$p_BD_bonferroni <- p.adjust(dat_stat$p_BD, 'bonferroni')
  dat_stat$p_bonferroni <- p.adjust(dat_stat$p, 'bonferroni')
}





