#shamelessly stollen from: https://raw.githubusercontent.com/r-lib/usethis/82f3e2ec17113bd0126bb9c47dd322c143377411/R/style.R 

# bullet <- function(lines, bullet) {
#   lines <- paste0(bullet, " ", lines)
#   cat_line(lines)
# }
# 
# failed <- function() {
#   crayon::red(clisymbols::symbol$cross)
# }
# 
# warned <- function() {
#   crayon::yellow(clisymbols::symbol$square_small_filled)
# }
# 
# passed <- function() {
#   crayon::green(clisymbols::symbol$tick)
# }
# 
# blue_arrow <- function() {
#   crayon::cyan(clisymbols::symbol$arrow_right)
# }

cat_line <- function(...) {
  cat(..., "\n", sep = "")
}
