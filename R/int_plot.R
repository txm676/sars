
##internal plot function (Not to be exported)

int_plot <- function(x, title, sh1 = 16, s1 = 2, s2 = 1, s3 = 14, s4 = 13, s5 = 15,
                     c1 = "darkred", c2 = "black", xl = "Area", yl = "Species richness",
                     p1 = 0)
{
  df <- x$data
  xx <- df$A
  yy <- df$S
  ff <- x$calculated

  g <- ggplot2::ggplot(data = df) +
    ggplot2::geom_point(ggplot2::aes(x = xx, y = yy), size = s1, col = c1, shape = sh1) +
    ggplot2::geom_line(ggplot2::aes(x = xx, y = ff), size = s2, col = c2) +
    ggplot2::xlab(xl) + ggplot2::ylab(yl) +
    ggplot2::ggtitle(title) + ggplot2::theme_bw() +
    ggplot2::theme(axis.title = ggplot2::element_text(size = s3), axis.text =
                     ggplot2::element_text(size = s4)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = p1, size = s5))

  return(g)
}
