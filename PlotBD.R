fixed_birth <- 4
fixed_death <- 1

birth_rate <- seq(1.1, 10, 0.1)
b_likelyhood <- lapply(birth_rate, bd.likelyhood.function, fixed_death, finches)
b_log_likelyhood <- lapply(birth_rate, bd.log.likelyhood.function, fixed_death, finches)

plot(birth_rate, b_likelyhood, xlab = "Birth Rate", ylab = "Likelihood")
title(main ="Likelihood of Birth-death model vs Birth Rate", sub = paste("Death rate of", fixed_death))

plot(birth_rate, b_log_likelyhood, xlab = "Birth Rate", ylab = "Log Likelihood")
title("Log Likelihood of Birth-death model vs Birth Rate", sub = paste("Death rate of", fixed_death))

death_rate <- seq(0.05, 3.95, 0.05)
d_log_likelyhood <- sapply(death_rate, function(x) bd.log.likelyhood.function(fixed_birth, x, finches))
d_likelyhood <- sapply(death_rate, function(x) bd.likelyhood.function(fixed_birth, x, finches))

plot(death_rate, d_likelyhood, xlab = "Death Rate", ylab = "Likelihood")
title(main ="Likelihood of Birth-death model vs Death Rate", sub = paste("Birth rate of", fixed_birth))

plot(death_rate, d_log_likelyhood, xlab = "Death Rate", ylab = "Log Likelihood")
title(main ="Log Likelihood of Birth-death model vs Death Rate", sub = paste("Birth rate of", fixed_birth))