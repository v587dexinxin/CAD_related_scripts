library(tidyverse)

# 1. 读取 FPKM 矩阵
fpkm <- read.table(
  "H:/work/Postdoctoral/GWAS疾病位点检测/results/HiRPC/RNA-seq/reads_count/PCA/union_all_reads_count.csv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

# 2. 设置目标基因
genes_show <- c(
  "KLF2", "KLF4", "NOS3", "THBD",
  "NQO1", "HMOX1",
  "ICAM1", "SELE", "CCL2", "CXCL8", "EDN1"
)

# 3. 提取目标基因
genes_found <- intersect(genes_show, rownames(fpkm))
fpkm_sub <- fpkm[genes_found, ]

# 4. 转成长表
fpkm_long <- fpkm_sub %>%
  rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "FPKM"
  )

# 5. 根据样本名添加分组信息
fpkm_long <- fpkm_long %>%
  mutate(
    Group = case_when(
      grepl("WT", Sample, ignore.case = TRUE) ~ "WT",
      grepl("LSS", Sample, ignore.case = TRUE) ~ "LSS",
      grepl("OSS", Sample, ignore.case = TRUE) ~ "OSS",
      TRUE ~ "Other"
    )
  )

# 6. 计算每组均值和标准差
fpkm_summary <- fpkm_long %>%
  group_by(Gene, Group) %>%
  summarise(
    mean_FPKM = mean(FPKM, na.rm = TRUE),
    sd_FPKM = sd(FPKM, na.rm = TRUE),
    .groups = "drop"
  )

# 7. 固定基因和分组顺序
fpkm_summary$Gene <- factor(fpkm_summary$Gene, levels = genes_show)
fpkm_summary$Group <- factor(fpkm_summary$Group, levels = c("WT", "LSS", "OSS"))

fpkm_long$Gene <- factor(fpkm_long$Gene, levels = genes_show)
fpkm_long$Group <- factor(fpkm_long$Group, levels = c("WT", "LSS", "OSS"))

# 8. 画柱状图 + replicate 点
p <- ggplot(fpkm_summary, aes(x = Group, y = mean_FPKM, fill = Group)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.25) +
  geom_errorbar(
    aes(ymin = mean_FPKM - sd_FPKM, ymax = mean_FPKM + sd_FPKM),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_jitter(
    data = fpkm_long,
    aes(x = Group, y = FPKM),
    width = 0.12,
    size = 1.8,
    color = "black",
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(
    WT  = "#4C78A8",
    LSS = "#59A14F",
    OSS = "#E15759"
  )) +
  labs(
    title = "FPKM expression of selected marker genes",
    x = NULL,
    y = "FPKM"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(p)

ggsave(
  "selected_genes_FPKM_barplot.pdf",
  p,
  width = 10,
  height = 7
)

ggsave(
  "selected_genes_FPKM_barplot.png",
  p,
  width = 10,
  height = 7,
  dpi = 300
)