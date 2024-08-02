# mk simple coverage plot

require(tidyverse)
require("BrattleExtras")

args <- commandArgs(trailingOnly=TRUE)
aav_name = args[1]
run = args[2]
cov_file = args[3]
outfile = args[4]

covDf = read_tsv(cov_file)

(covDf %>% 
  ggplot() + geom_line(aes(x=pos, y=depth, colour=sample_name)) +
		   theme_bw() + ggtitle(aav_name) + ggsubtitle(run) + xlim(0, 4000)) %>% ggsave(filename=outfile, width = 30, height = 25, units = "cm", plot=.)
