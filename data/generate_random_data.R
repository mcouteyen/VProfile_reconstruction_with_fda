

my_data <- data.frame()

args <- seq(1, 50, .5)
n_args <- length(args)

for(i in 1:30){
  min <- sample(c(13:16, 20:22), 1)
  val <- runif(n_args, min, min + 0.01)
  my_data <- rbind(my_data, data.frame(INDEX = i,
                                       DEPTH = args,
                                       TEMPERATURE = val))
}

for(i in 31:65){
  val <- sort(runif(n_args, 16, 19), decreasing = TRUE)
  my_data <- rbind(my_data, data.frame(INDEX = i,
                                       DEPTH = args,
                                       TEMPERATURE = val))
}

for(i in 66:100){
  val <- sort(runif(n_args, 16, 25), decreasing = TRUE)
  my_data <- rbind(my_data, data.frame(INDEX = i,
                                       DEPTH = args,
                                       TEMPERATURE = val))
}

my_data$TEMPERATURE <- round(my_data$TEMPERATURE, 5)

head(my_data)

source("src/profile_tempering.R")
sparse_data <- shrink_random_profiles_with_seed(my_data)

write.csv(sparse_data, "data/example_data.csv", row.names = FALSE)


