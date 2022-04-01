#!/usr/bin/env Rscript

library(ggplot2)
library(getopt)

command <- matrix(c(
  'input', 'i', 1, 'character', 'input file containing corresponding relationships bewteen postion and phastCons score',
  'title', 't', 1, 'character', 'title of output figure',
  'xlab', 'x', 1, 'character', 'label of x-axis',
  'ylab', 'y', 1, 'character', 'label of y-axis',
  'color', 'c', 1, 'character', 'color of line',
  'low', 'a', 1, 'character', 'color of low score, suppressed when specify color',
  'high', 'z', 1, 'character', 'color of high score, suppressed when specify color',
  'output', 'o', 1, 'character', 'output PDF file',
  'help', 'h', 0, 'logical', 'show this help message and exit'
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

if (!is.null(opts$help)) {
  cat(getopt(command, usage = TRUE))
  q(status=2)
}
if (is.null(opts$input) || is.null(opts$output)) {
  cat(getopt(command, usage = TRUE))
  q(status=3)
}
if (is.null(opts$title) || is.null(opts$xlab) || is.null(opts$ylab)) {
  cat(getopt(command, usage = TRUE))
  q(status=4)
}
if (is.null(opts$low)) {opts$low = "cyan"}
if (is.null(opts$high)) {opts$high = "magenta"}

df <- data.frame(read.delim(opts$input, header = FALSE, stringsAsFactors=FALSE))
colnames(df) <- c("position", "score")

if (!is.null(opts$color)) {
  p <- ggplot(df, aes(x = position, y = score)) +
    geom_line(size = 0.75, colour = paste("#", opts$color, sep = ""))
} else {
  p <- ggplot(df, aes(x = position, y = score, colour = score)) +
    scale_colour_gradient(low = opts$low, high = opts$high) +
    geom_line(size = 0.75)
}

# p <- p + scale_x_continuous(breaks = seq(0, nrow(df), 10 ** (nchar(nrow(df)) - 1))) +
#   theme_bw()
p <- p + theme_bw()
p <- p + ggtitle(opts$title) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
p <- p + xlab(opts$xlab) +
  theme(axis.title.x = element_text(size = 16, face = "bold"))
p <- p + ylab(opts$ylab) +
  theme(axis.title.y = element_text(size = 16, face = "bold"))
p <- p + theme(panel.grid.minor = element_line(linetype = "dashed"))

# ggsave(opts$output, width = 18, height = 10, units = "in", dpi = 600)
ggsave(opts$output, device = "pdf", width = 18, height = 10, units = "in")
