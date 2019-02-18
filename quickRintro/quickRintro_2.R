# QuickRIntro - objects in R

a <- sqrt (25)
b <- round (x = 23.4567, digits = 1)

d <- "apple"  # ' or " both ok
e <- FALSE  # TRUE vs FALSE, T vs F, not: true vs false, True vs False

A <- "banana"

var1 <- 1
var.1 <- 2
var_1 <- 3
Var_1 <- 4

# not allowed: 1_var - cannot start with number; var-1

vec_num <- c(2.3, 3.4, 4.5, 5.6)
vec_char <- c('apple', 'carrot', 'lemon')
vec_log <- c(TRUE, FALSE, TRUE)

length (vec_num)
length (vec_char)

vec_num [1:2]
vec_num [c(TRUE, TRUE, FALSE, FALSE)]
vec_char[vec_log]

mat_num <- matrix (1:80, nrow = 10, ncol = 8, byrow = TRUE)
dat_num <- as.data.frame (mat_num)

dim (cars)
names (cars)  # also colnames
rownames (cars)

cars [, "speed"]
cars$speed

list_1 <- list (com_1 = vec_num, com_2 = mat_num, com_3 = list (com_3.1 = vec_char, com_3.2 = cars))

list_1$com_1[2:4]
list_1[[1]][2:4]

list_1$com_3$com_3.2$dist


LETTERS
list_1[[3]][[2]][,"dist"]
