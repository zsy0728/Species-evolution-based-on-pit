rm(list = ls())
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(org.Hs.eg.db)
library(tidyverse)

out_dir <- "Y:/sda/results/Evolution_cortical_shape/homo_with_related_species/enrich_results/"

# 定义一个函数，封装富集分析及绘图流程
run_enrichGO_and_plot <- function(gene_file, prefix) {
  gene_df <- read.csv(file.path(out_dir, gene_file))
  gene_symbols <- gene_df$entrez_id
  
  eG <- enrichGO(gene = gene_symbols,
                 OrgDb = org.Hs.eg.db,
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.01,
                 ont = "all",
                 readable = TRUE)
  
  # 保存富集结果
  write.table(eG, file = file.path(out_dir, paste0(prefix, "_eG.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # 读取富集结果，计算Enrichment Factor
  eGo <- read.table(file.path(out_dir, paste0(prefix, "_eG.txt")),
                    header = TRUE, sep = "\t", quote = "")
  eGo <- separate(eGo, GeneRatio, into = c("GR1", "GR2"), sep = "/")
  eGo <- separate(eGo, BgRatio, into = c("BR1", "BR2"), sep = "/")
  eGo <- mutate(eGo, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2)) /
                  (as.numeric(BR1)/as.numeric(BR2)))
  
  # 取BP、MF、CC分别前10条
  eGoBP <- eGo %>% filter(ONTOLOGY == "BP") %>% slice_head(n = 10)
  eGoCC <- eGo %>% filter(ONTOLOGY == "CC") %>% slice_head(n = 10)
  eGoMF <- eGo %>% filter(ONTOLOGY == "MF") %>% slice_head(n = 10)
  eGo10 <- bind_rows(eGoBP, eGoMF, eGoCC)
  
  # 画气泡图
  p <- ggplot(eGo10, aes(enrichment_factor, Description)) + 
    geom_point(aes(size = Count, color = -log10(pvalue), shape = ONTOLOGY)) +
    scale_color_gradient(low = "green", high = "red") + 
    labs(color = expression(-log[10](p_value)), size = "Count", shape = "Ontology",
         x = "Enrichment Factor", y = "GO term", title = paste0(prefix, " GO enrichment")) + 
    theme_bw() +
    facet_wrap(~ ONTOLOGY, ncol = 1, scales = "free")
  
  return(p) 
}
# 保存图片的函数，可以自定义比例
save_enrichment_plot <- function(plot_obj, prefix, width = 10, height = 6, dpi = 300) {
  ggsave(filename = file.path(out_dir, paste0(prefix, "_GO_enrichment.pdf")),
         plot = plot_obj, width = width, height = height, dpi = dpi)
}
# 运行shared区域分析
shared_plot <- run_enrichGO_and_plot("shared_top_2000_geneID.csv", "shared")
# 保存shared图片，可自定义宽高比例
save_enrichment_plot(shared_plot, "shared", width = 8, height = 8)

# 运行unique区域分析
unique_plot <- run_enrichGO_and_plot("unique_top_2000_geneID.csv", "unique")
# 保存unique图片，可自定义宽高比例
save_enrichment_plot(unique_plot, "unique", width = 7.2, height = 8)    