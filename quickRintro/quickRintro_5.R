# Installing packages in R

# install.packages("vegan")

library (vegan)

data (dune)

specnumber (dune)

detach("package:vegan", unload=TRUE)  # detached the vegan package from Global Environment

vegan::specnumber (dune)
