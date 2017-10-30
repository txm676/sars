
#' @export

plot.coleman <- function(x, ps = 3, sh = 21, ls1 = 1.2, ls2 = 1.2, as1 = 14, as2 = 13,
                         c1 = "darkblue", c2 = "black", c3 = "#5e5e5e",
                         xl = "Relative area (log transformed)", yl = "Species richness", th = NULL){

  df <- data.frame(x$Predicted_values, x$Standard_deviation,
                   x$Relative_areas, x$Species_richness)
  colnames(df) <- c("pv", "sd", "ra", "os")
  ggplot2::ggplot(data = df) + ggplot2::geom_point(ggplot2::aes(x = log(ra),
                                                                y = os), size = ps, col = c1, shape = sh) +
    th + ggplot2::geom_line(ggplot2::aes(x = log(ra), y = pv), col = c2, size = ls1) +
    ggplot2::geom_line(ggplot2::aes(x = log(ra), y = (pv + sd)), col = c3, size = ls2) +
    ggplot2::geom_line(ggplot2::aes(x = log(ra), y = (pv - sd)), col = c3, size = ls2) +
    ggplot2::xlab(xl) + ggplot2::ylab(yl) +
    ggplot2::theme(axis.title = ggplot2::element_text(size = as1), axis.text =
                     ggplot2::element_text(size = as2)) 
}

