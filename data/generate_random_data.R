

somlit_data <- read.csv("C:/Users/mcouteyen/Dropbox/Thèse/publicgit_Couteyen_Carpaye_thesis/data/fpca/somlit_ctd_marseille.csv")
fdobj <- create_fdobj(somlit_data, "PROFONDEUR", "TEMPERATURE", lambda = bspline_smooth_lambda, K = 50)
new_data <- eval.fd(1:50, fdobj) %>% as.data.frame()
new_data <- cbind(data.frame(DEPTH = 1:50), new_data)
new_data <- pivot_longer(new_data, 2:ncol(new_data), values_to = "TEMPERATURE", names_to = "INDEX") %>% as.data.frame()
indexes <- unlist(strsplit(new_data$INDEX, " "))
new_data$INDEX <- parse_integer(indexes [seq(2, length(indexes), 2)])
new_data <- new_data %>% arrange(INDEX, DEPTH)
new_data <- new_data %>% filter(INDEX < 251)

source("src/profile_tempering.R")
sparse_data <- shrink_random_profiles_with_seed(new_data)

write.csv(sparse_data, "data/example_data.csv", row.names = FALSE)


