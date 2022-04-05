
#install rrtools
remotes::install_github("benmarwick/rrtools")

#install development package version for ENMeval
devtools::install_github("jamiemkass/ENMeval")

#install development package version for kuenm
devtools::install_github("marlonecobos/kuenm", force = TRUE)

#make compendium
rrtools::use_compendium("../spermwhale", open = FALSE)

#load all functions in R
devtools::load_all()

#document functions
devtools::document()

#utiliser les pipes dans les fonctions
usethis::use_pipe()

#to automatically add/remove dependencies
rcompendium::add_dependencies(here::here())


