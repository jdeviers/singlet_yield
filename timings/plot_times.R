#!/usr/bin/Rscript

library(ggplot2)

df <- read.table("timings_ISCA_32threads_zoom.dat", sep = "", header = FALSE)

names(df) <- c('N','m1','m2s','m2p')

p <- ggplot(df, aes(N)) +
        geom_line(aes(y = m1,  colour = "m1"))  + 
        geom_line(aes(y = m2s, colour = "m2s")) + 
        geom_line(aes(y = m2p, colour = "m2p")) +
        ggtitle("Algorithms runtimes with 16cores/32 threads, on ISCA node") +
        theme(plot.title = element_text(hjust = 0.5)) +
        xlab("Matrix size (NxN)") + ylab("Runtime (s)")

ggsave("timings.pdf", plot = p)
