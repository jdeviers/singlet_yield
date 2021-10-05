#!/usr/bin/Rscript

library(ggplot2)

df <- read.table("timings_ISCA_16v32.dat", sep = "", header = TRUE)

names(df) <- c('N','cores16','cores32')

p <- ggplot(df, aes(N)) +
        geom_line(aes(y = cores16, colour = "16 threads")) +
        geom_line(aes(y = cores32, colour = "32 threads")) +
        ggtitle("Algorithms runtimes with 16 vs. 32 threads.") +
        theme(plot.title = element_text(hjust = 0.5)) +
        xlab("Matrix size (NxN)") + ylab("Runtime (s)")

ggsave("timings_ISCA_16v32.pdf", plot = p)
