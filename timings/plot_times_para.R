#!/usr/bin/Rscript

library(ggplot2)

df <- read.table("timings_GPU.dat", sep = "", header = FALSE)

names(df) <- c('N','m2p')

p <- ggplot(df, aes(N)) +
        geom_line(aes(y = m2p, colour = "m2p")) +
        ggtitle("Algorithms runtimes with 32cores/64 threads, on PC.") +
        theme(plot.title = element_text(hjust = 0.5)) +
        xlab("Matrix size (NxN)") + ylab("Runtime (s)")

ggsave("timings_GPU.pdf", plot = p)
