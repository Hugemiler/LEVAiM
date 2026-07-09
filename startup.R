######
# Pre-
######

install.packages(c(
  "usethis",
  "devtools",
  "roxygen2",
  "testthat",
  "cli"
))

usethis::use_mit_license("Guilherme Fahur Bottino")
usethis::use_readme_rmd()
usethis::use_testthat()
usethis::use_news_md()

usethis::use_package("cli")

######
# Post-
######

devtools::document()
devtools::load_all()
devtools::test()
