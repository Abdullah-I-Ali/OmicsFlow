df <- data.frame(
  hgnc_symbol = paste0('GENE', 1:5000),
  chromosome_name = sample(1:22, 5000, replace=TRUE),
  start_position = sample(1:10000000, 5000, replace=TRUE)
)
df$end_position <- df$start_position + 1000
saveRDS(df, 'fake_cache.rds')
