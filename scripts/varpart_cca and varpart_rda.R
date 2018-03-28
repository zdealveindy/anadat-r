# does dbRDA with Chi-square distance give the same results as CCA?

vltava.spe <- read.delim ('http://www.davidzeleny.net/anadat-r/data-download/vltava-spe.txt', row.names = 1)
vltava.env <- read.delim ('http://www.davidzeleny.net/anadat-r/data-download/vltava-env.txt')

# varpart na tb-RDA
library (vegan)
vltava.spe.hell <- decostand (sqrt (vltava.spe), method = 'hell')
env1 <- vltava.env[, c('pH', 'SOILDPT')]
env2 <- vltava.env[, c('LITHIC', 'SKELETIC', 'CAMBISOL', 'FLUVISOL')]

vp.rda.1 <- varpart (vltava.spe.hell, ~ pH + SOILDPT, ~ LITHIC + SKELETIC + CAMBISOL + FLUVISOL, data = vltava.env)
plot (vp.rda.1)


rda.env1.env2 <- rda (vltava.spe.hell ~ pH + SOILDPT + LITHIC + SKELETIC + CAMBISOL + FLUVISOL, data = vltava.env)
rda.env1_env2 <- rda (vltava.spe.hell ~ pH + SOILDPT + Condition (LITHIC + SKELETIC + CAMBISOL + FLUVISOL), data = vltava.env)
rda.env2_env1 <- rda (vltava.spe.hell ~ LITHIC + SKELETIC + CAMBISOL + FLUVISOL + Condition (pH + SOILDPT), data = vltava.env)

ABC <- RsquareAdj (rda.env1.env2)$adj.r.squared
A <- RsquareAdj (rda.env1_env2)$adj.r.squared
C <- RsquareAdj (rda.env2_env1)$adj.r.squared
B <- ABC - A - C

A
B
C


# CCA (without Hellinger transformation)

cca.env1.env2 <- cca (vltava.spe ~ pH + SOILDPT + LITHIC + SKELETIC + CAMBISOL + FLUVISOL, data = vltava.env)
cca.env1_env2 <- cca (vltava.spe ~ pH + SOILDPT + Condition (LITHIC + SKELETIC + CAMBISOL + FLUVISOL), data = vltava.env)
cca.env2_env1 <- cca (vltava.spe ~ LITHIC + SKELETIC + CAMBISOL + FLUVISOL + Condition (pH + SOILDPT), data = vltava.env)

anova (cca.env1.env2)
anova (cca.env1_env2)
anova (cca.env2_env1)

ABC <- RsquareAdj (cca.env1.env2)$adj.r.squared
A <- RsquareAdj (cca.env1_env2)$adj.r.squared
C <- RsquareAdj (cca.env2_env1)$adj.r.squared
B <- ABC - A - C

A
B
C


# general function - two data matrices
varpart_rda <- function (Y, X1, X2)
{
  require (vegan)
  rda.env1.env2 <- rda (Y, Y = cbind (X1, X2))
  rda.env1_env2 <- rda (Y, Y = X1, Z = X2)
  rda.env2_env1 <- rda (Y, Y = X2, Z = X1)
  
  ABC <- RsquareAdj (rda.env1.env2)$adj.r.squared
  A <- RsquareAdj (rda.env1_env2)$adj.r.squared
  C <- RsquareAdj (rda.env2_env1)$adj.r.squared
  B <- ABC - A - C
  
  P_ABC <- anova (rda.env1.env2)$`Pr(>F)`[1]
  P_A <- anova (rda.env1_env2)$`Pr(>F)`[1]
  P_C <- anova (rda.env2_env1)$`Pr(>F)`[1]
  
 list (ABC = ABC, A = A, B = B, C = C, P_ABC = P_ABC, P_A = P_A, P_C = P_C)
}

vp.rda <- varpart_rda (vltava.spe.hell, X1 = vltava.env[,c('pH', 'SOILDPT')], X2 = vltava.env[,c('LITHIC', 'SKELETIC', 'CAMBISOL', 'FLUVISOL')])

varpart_cca <- function (Y, X1, X2, ...)
{
  require (vegan)
  cca.env1.env2 <- cca (Y, Y = cbind (X1, X2))
  cca.env1_env2 <- cca (Y, Y = X1, Z = X2)
  cca.env2_env1 <- cca (Y, Y = X2, Z = X1)
  
  ABC <- RsquareAdj (cca.env1.env2)$adj.r.squared
  A <- RsquareAdj (cca.env1_env2)$adj.r.squared
  C <- RsquareAdj (cca.env2_env1)$adj.r.squared
  B <- ABC - A - C
  
  P_ABC <- anova (cca.env1.env2, ...)$`Pr(>F)`[1]
  P_A <- anova (cca.env1_env2, ...)$`Pr(>F)`[1]
  P_C <- anova (cca.env2_env1, ...)$`Pr(>F)`[1]
  
  out <- list (ABC = ABC, A = A, B = B, C = C, P_ABC = P_ABC, P_A = P_A, P_C = P_C)
  out$indfract <- matrix (ncol = 3, nrow = 4)
  out$indfract[,3] <- c(A, B, C, 1-(A+B+C))
  out$nsets <- 2
  class (out) <- "varpart_rda_cca"
  return (out)
}

print.varpart_rda_cca <- function (x)
{
  output <- data.frame (Adj.R.squared = c(x$A, x$B, x$C, 1 - (x$A+x$B+x$C)), P.value = c(x$P_A, '', x$P_C, ''))
  rownames (output) <- c("[a] = X1|X2", "[b]", "[c] = X2|X1", "[d] = Residuals")
  print (output)
}

plot.varpart_rda_cca <- function (x, ...)
{
  vegan:::plot.varpart234 (x, ...)
}

vp.cca <- varpart_cca (vltava.spe.hell, X1 = vltava.env[,c('pH', 'SOILDPT')], X2 = vltava.env[,c('LITHIC', 'SKELETIC', 'CAMBISOL', 'FLUVISOL')])

plot (vp.cca)
varpart_cca (vltava.spe.hell, X1 = vltava.env[,c('pH', 'SOILDPT')], X2 = vltava.env[,c('LITHIC', 'SKELETIC', 'CAMBISOL', 'FLUVISOL')])
