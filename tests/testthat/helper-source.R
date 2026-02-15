# Source all R files for testing without package installation
r_files <- list.files(
  file.path(dirname(dirname(getwd())), "R"),
  pattern = "\\.R$", full.names = TRUE
)
# If running from testthat dir, try project root
if (length(r_files) == 0) {
  r_files <- list.files("../../R", pattern = "\\.R$", full.names = TRUE)
}
# Also try from test_dir context
if (length(r_files) == 0) {
  r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
}
for (f in r_files) {
  source(f)
}
