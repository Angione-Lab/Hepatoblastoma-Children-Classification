library(matrixStats)
library(tidyr)
library(ggplot2)
library(dplyr)
library(purrr)
library(ggpubr)
library(reshape2)

input_files <- "Health status\\health_status_no_sig_average_pathway_weights.csv"
out_put_files <- "Health status\\health_status_no_sig_average_pathway_weigthbarplot.pdf"

vals<- c("Control" = "#4A9AAA", "Tumour" = "red")

br = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100)
lim = c(0, 1100)
percent_threshold = 0.80
y_title <- expression(paste("Mean ", italic("flux rate"), ""), font.lab=2)
pathway.cor <- read.csv(file = file_input)
colnames(pathway.cor) <- c("Pathway", "Mean", "Status")
pathway.cor <- pathway.cor[Reduce(`&`, lapply(pathway.cor, function(x) !(is.na(x)|x==""))),]  # delete empty/Nan entries

p <- pathway.cor %>% filter(Mean > 0.0) %>%            # in order to disregard the pathways with 0 weight
  filter(quantile(Mean, percent_threshold)<Mean) %>%
  ggplot(aes(x = Pathway, y = Mean), fill = Status, color = as.factor(Status)) +
  geom_bar(aes(fill = Status), stat = "identity", alpha = 0.6) +
  theme_bw() + 
  scale_y_continuous(breaks = br, limits = lim) +         # only positive value
  xlab("") + ylab(y_title) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        plot.margin=unit(c(0,0,-5,0),"mm"), # to remove white space below
        axis.ticks = element_blank(),
        text = element_text(size=20),   # y-axis
        axis.text.x = element_text(hjust = 1, vjust = 1.0, angle = 90, size = 20),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position="right", legend.text=element_text(size = 20)) +
  scale_fill_manual(values = vals)

ggsave(filename = file_output, width = 9, height = 8, dpi=600, device="pdf")
p
dev.off()
p
