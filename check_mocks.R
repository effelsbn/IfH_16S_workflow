library(tidyverse)

load("DADA2_pipeline_outputs/Out_B8779_2022-07-26_11-31-50/DADA2_out_seqtab.RData")
load("DADA2_pipeline_outputs/Out_B8779_2022-07-26_11-31-50/taxa.RData")
zymo <- read.delim("../16S_Strohstall_Kiel/input/zymo_mock_genus.tsv")
cols26 <- read_lines("../16S_Strohstall_Kiel/input/26colors.txt")

OTU <- t(st.nochim)

controls2 <- merge(taxa, OTU, by = 0, all = TRUE) %>% 
  select(Genus, WDHzymoII) %>% 
  filter(WDHzymoII > 0) %>% 
  mutate(rel_abund = (WDHzymoII/sum(WDHzymoII)*100)) %>% 
  group_by(Genus) %>% 
  summarize(ctrl = sum(rel_abund))

comptab2 <- merge(zymo, controls2, by = "Genus", all = TRUE)

comptab_stack2 <- arrange(comptab, desc(Zymo)) %>%
  pivot_longer(cols = !Genus, names_to = "sample")

comptab_stack2$Genus <- factor(comptab_stack2$Genus, levels = unique(comptab2$Genus))

ggplot(comptab_stack, aes(fill=Genus, y=value, x=sample)) + 
  geom_bar(position="stack", stat="identity", color = "black") +
  scale_fill_manual(values = cols26) +
  ylab("relative abundance") +
  xlab(NULL) +
  scale_x_discrete(labels = c("theoretical_composition", "B8779")) +
  theme_classic()

ggsave("mockplot2.png")
