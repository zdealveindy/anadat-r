# Importing the data via clipboard
#vltava.spe <- read.delim (file = "clipboard", row.names = 1)
#vltava.spe.t <- t (vltava.spe)

# Importing the data from the file
setwd ('C:/Users/Zeleny/Dropbox/anadatr/data')
dir ()
vltava.spe <- t(read.delim (file = 'vltava.spe.txt', row.names = 1))

# DCA on vltava.spe data
install.packages ("vegan")
library (vegan)
decorana (veg =  log1p(vltava.spe))

# now calculating CA on vltava.spe
CA <- cca (X = log1p (vltava.spe))
ordiplot (CA, display = 'sites', type = 't')


# After this comes Exercise 1 and 3 from Unconstrained ordinations > CA+DCA > Exercises
# (Exercise 2 is homework)