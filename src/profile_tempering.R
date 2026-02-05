
shrink_random_profiles_with_seed <- function(data, test_size = .1){
  data_degraded <- data
  data_indexes <- unique(data$INDEX)
  n_random <- round(length(data_indexes)*test_size)
  set.seed(123)
  random_profiles <- sample(data_indexes, n_random)
  for(index in random_profiles){
    set.seed(222 + index)
    upper_cut <- sample(3:10, 1)
    set.seed(111 + index)
    lower_cut <- sample(40:48, 1)
    data_degraded <- data_degraded %>% filter(!(INDEX == index & (DEPTH < upper_cut | DEPTH > lower_cut)))
  }
  return(data_degraded)
}
