#shamelessly stollen from: https://raw.githubusercontent.com/r-lib/usethis/82f3e2ec17113bd0126bb9c47dd322c143377411/R/style.R 

bullet <- function(lines, bullet) {
  lines <- paste0(bullet, " ", lines)
  cat_line(lines)
}

failed <- function(...) {
  bullet(paste0(...), bullet = crayon::red(clisymbols::symbol$bullet))
}

passed <- function(...) {
  bullet(paste0(...), bullet = crayon::green(clisymbols::symbol$tick))
}

cat_line <- function(...) {
  cat(..., "\n", sep = "")
}