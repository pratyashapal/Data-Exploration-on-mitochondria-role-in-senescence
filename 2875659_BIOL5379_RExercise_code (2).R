# LOADING NECESSARY LIBRARIES ####
library(ggplot2)
library(reshape2)
library(ggrepel)
library(amap)
library(ggridges)
library(clusterProfiler)
library(org.Hs.eg.db)
library(devEMF)
library(dplyr)

# READING THE TABLES ####
em = read.table("C:/Users/Pratyasha/OneDrive/Desktop/data ex/Dataset and Description-20240116/EM.csv", header = TRUE, sep = '\t', row.names = 1)
annot = read.table("C:/Users/Pratyasha/OneDrive/Desktop/data ex/Dataset and Description-20240116/Human_Background_GRCh38.p13.csv", header = TRUE, sep = '\t', row.names = 1)
ss = read.table("C:/Users/Pratyasha/OneDrive/Desktop/data ex/Dataset and Description-20240116/sample_sheet.csv", header = TRUE, sep = '\t', row.names = 1)
de_s_p = read.table("C:/Users/Pratyasha/OneDrive/Desktop/data ex/Dataset and Description-20240116/DE_Senes_vs_Prolif.csv", header = TRUE, sep = '\t', row.names = 1)
de_mt_p = read.table("C:/Users/Pratyasha/OneDrive/Desktop/data ex/Dataset and Description-20240116/DE_Senes_MtD_vs_Prolif.csv", header = TRUE, sep = '\t', row.names = 1)
de_mt_s = read.table("C:/Users/Pratyasha/OneDrive/Desktop/data ex/Dataset and Description-20240116/DE_Senes_MtD_vs_Senes.csv", header = TRUE, sep = '\t', row.names = 1)

# MERGING EM AND ANNOTATION FILES ####
master_temp = merge(em, annot, by.x = 0, by.y = 0)
rownames(master_temp) = master_temp[,1]

# GETTING EM_SYMBOLS AND EM_SCALED ####
rownames(master_temp) = master_temp[,"SYMBOL"]
em_symbol = master_temp[,as.vector(rownames(ss))]
em_symbol = na.omit(em_symbol)
em_scaled = na.omit(data.frame(t(scale(t(em_symbol)))))


# DENSITY PLOT ####
#melting the em symbol
em_symbol.m = melt(em_symbol)


em_symbol.m$group = ifelse(grepl("Prolif", em_symbol.m$variable), "Prolif",
                           ifelse(grepl("Senes_MtD", em_symbol.m$variable), "Senes_MtD",
                                  ifelse(grepl("Senes", em_symbol.m$variable), "Senes",  NA)))



color_palette = c("olivedrab4", "violetred", "steelblue")  

density_plot = ggplot(em_symbol.m, aes(x = log10(value), fill = group)) +
  geom_density(colour = "black", size = 0.5, alpha = 0.2) +
  facet_wrap(~variable, ncol = 3) +
  scale_fill_manual(values = color_palette) +
  labs(title = "Density Plot" , x = "log10(Expression Value)", y = "Density") +
  theme_grey() +  # Here's the corrected placement
  theme(
    text = element_text(family = "TT Times New Roman"),  # Set font family
    plot.title = element_text(size = 18, face = "bold"),  # Adjust title size and style
    axis.title = element_text(size = 14),  # Adjust axis label size
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 12),  # Adjust legend text size
    legend.position = "bottom",  # Position legend
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Adjust panel border
    panel.spacing = unit(1, "lines")
  )

density = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/density.emf", density_plot,3 ,4)

# PCA PLOT ####

numeric_matrix = as.matrix(sapply(em_scaled, as.numeric))
pca = prcomp(t(numeric_matrix))
pca_coordinates = data.frame(pca$x)


vars = apply(pca$x, 2, var)
prop_x = round(vars["PC1"] / sum(vars),4) * 100
prop_y = round(vars["PC2"] / sum(vars),4) * 100
x_axis_label = paste("PC1 ", " (",prop_x, "%)",sep="")
y_axis_label = paste("PC2 ", " (",prop_y, "%)",sep="")

ss$SAMPLE_GROUP = as.vector(ss$SAMPLE_GROUP)
ss$sample = row.names(ss)

group_centroids = pca_coordinates %>%
  group_by(ss$SAMPLE_GROUP) %>%
  summarize(PC1_centroid = mean(PC1), PC2_centroid = mean(PC2))

# PCA Plot final
pca_plot = ggplot() +
  geom_point(data = pca_coordinates, aes(x = PC1, y = PC2, color = ss$SAMPLE_GROUP), size = 3, alpha = 0.8) +
  geom_point(data = group_centroids, aes(x = PC1_centroid, y = PC2_centroid), color = "black", size = 10, alpha = 0.1) +
  
  scale_color_manual(values = c("olivedrab4", "steelblue", "violetred")) +
  geom_text_repel(data = pca_coordinates, aes(x = PC1, y = PC2, label = ss$sample), size = 4, box.padding = 0.5, segment.color = "Transparent") +
  labs(title = "PCA Plot", x = x_axis_label, y = y_axis_label) +
  theme_grey() +  # Use minimal theme for cleaner look
  theme(
    text = element_text(family = "TT Times New Roman"),  
    plot.title = element_text(size = 18, face = "bold"),  
    axis.title = element_text(size = 14),  
    legend.title = element_blank(),  
    legend.text = element_text(size = 12),  
    legend.position = "bottom",  
    panel.border = element_rect(colour = "black", fill = NA, size = 1), 
    panel.spacing = unit(1, "lines")
  )

pca_plot
pca = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/pca.emf", pca_plot,5 ,6)


# SENESCENT VS PROLIFERATING CELLS ANALYSIS ####
# PARSING THE DE TABLE
master_de_sp = parse_de(de_s_p, master_temp)

# GETTING EXPRESSION TABLE FOR JUST THE SIG GENES
em_scaled_sig_sp = sig_genes_exp(master_de_sp,em_scaled)

# VOLCANO PLOT 
volcano_plot_sp = volcano_plot(master_de_sp, "SP Volcano Plot")
print(volcano_plot_sp)
vol_sp = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/volcano_plot_sp.emf", volcano_plot_sp,6 ,7)

# MULTIGENE BOXPLOT
candidate_genes_sp = head(master_de_sp$SYMBOL,15)
box_sp = faceted_boxplot(candidate_genes_sp,em_scaled,ss, 15, "SP Boxplot for Top 15 Sig Genes")
print(box_sp)
box_plot_sp = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/box_plot_sp.emf", box_sp,4,6)

# MA PLOT
ma_plot_sp = ma_plot(master_de_sp, "MA Plot for SP")
print(ma_plot_sp)
ma_sp = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/ma_plot_sp.emf", ma_plot_sp,3 ,4)

# HEATMAP
em_scaled_sp = em_scaled_sig_sp[,c(1:6)]
heat_sp = heatmap(em_scaled_sp[1:100,])
print(heat_sp)
heat_sp = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/heat_sp.emf", heat_sp,5 ,4)

# gsea
gea_sp = gse_analysis(master_de_sp)
ridgeplot(gea_sp)

# ora
results_sp = ora_up_down(master_de_sp)
up_results_sp = results_sp$ora_result_up
down_results_sp = results_sp$ora_result_down
bar_down_sp = barplot(down_results_sp, showCategory=10)
bar_down_sp
barplot_down_sp = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/bar_down_sp.emf", bar_down_sp,4 ,5)
barplot_up_sp = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/bar_up_sp.emf", bar_up_sp,5 ,7)
bar_up_sp = barplot(up_results_sp, showCategory=10)
cnetplot(up_results_sp, showCategory=2)

# SENESCENT VS MITOCHONDRIA DEPLETED CELLS ANALYSIS ####

# PARSING THE DE TABLE
master_de_mts = parse_de(de_mt_s, master_temp)

# GETTING EXPRESSION TABLE FOR JUST THE SIG GENES
em_scaled_sig_mts = sig_genes_exp(master_de_mts,em_scaled)

# VOLCANO PLOT
volcano_plot_mts = volcano_plot(master_de_mts, " MTS Volcano Plot")
print(volcano_plot_mts)
vol_mts = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/volcano_plot_mts.emf", volcano_plot_mts,6 ,7)


# MULTIGENE BOXPLOT
candidate_genes_mts = head(master_de_mts$SYMBOL,15)
box_mts = faceted_boxplot(candidate_genes_mts,em_scaled,ss, 15,"MTS Boxplot for Top 15 Sig Genes")
print(box_mts)
boxplot_mtp = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/boxplot_mts.emf", box_mts,4 ,6)

# MA PLOT
ma_plot_mts = ma_plot(master_de_mts, "MA Plot for MTS")
print(ma_plot_mts)
ma_mts = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/ma_plot_mts.emf", ma_plot_mts,3 ,4)

# HEATMAP
em_scaled_mts = em_scaled_sig_mts[,c(4:9)]
heat_mts = heatmap(em_scaled_mts[1:100,])
print(heat_mts)
heat_mts = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/heat_mts.emf", heat_mts,5 ,4)

# gsea
gea = gse_analysis(master_de_mts)

ggp = ridgeplot(gea) + theme(axis.text.y = element_text(size = 6))
ggp

#1 gene
rip = single_gene_plot("RIPOR3",em_scaled,ss,geom_boxplot, box)
print(rip)

# ora
results_mts = ora_up_down(master_de_mts)
up_results_mts = results_mts$ora_result_up
down_results_mts = results_mts$ora_result_down
bar_up_mts = barplot(up_results_mts, showCategory=10)
bar_down_mts = barplot(down_results_mts, showCategory=10)
barplot_down_mts = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/bar_down_mts.emf", bar_down_mts,4 ,4)
barplot_up_mts = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/bar_up_mts.emf", bar_up_mts,4 ,4)


# PROLIFERATING VS MITOCHONDRAIL DEPLETED CELLS ANALYSIS ####

# PARSING THE DE TABLE
master_de_mtp = parse_de(de_mt_p, master_temp)

# GETTING EXPRESSION TABLE FOR JUST THE SIG GENES
em_scaled_sig_mtp = sig_genes_exp(master_de_mtp,em_scaled)

# VOLCANO PLOT
volcano_plot_mtp = volcano_plot(master_de_mtp, "MTP Volcano Plot")
print(volcano_plot_mtp)
vol_mtp = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/volcano_plot_mtp.emf", volcano_plot_mtp,6 ,7)

# MULTIGENE BOXPLOT
candidate_genes = head(master_de_mtp$SYMBOL,15)
box_mtp = faceted_boxplot(candidate_genes,em_scaled,ss, 15, " MTP Boxplot")
print(box)
boxplot_mtp = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/boxplot_mtp.emf", box_mtp,4 ,6)

# MA PLOT
ma_plot_mtp = ma_plot(master_de_mtp, "MA Plot for MTP")
print(ma_plot_mtp)
ma_mtp = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/ma_plot_mtp.emf", ma_plot_mtp,3 ,4)

# HEATMAP
em_scaled_mtp = em_scaled_sig_mtp[,-c(4:6)]
heat_mtp = heatmap(em_scaled_mtp)
print(heat_mtp)
heat_mtp = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/heat_mtp.emf", heat_mtp,5 ,4)


# gsea
gea = gse_analysis(master_de_mtp)
gea_mtp = ridgeplot(gea)
g_mtp = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/g1_mtp.emf", g_mtp,7 ,6)

# ora
results_mtp = ora_up_down(master_de_mtp)
up_results_mtp = results_mtp$ora_result_up
down_results_mtp = results_mtp$ora_result_down
bar_down_mtp = barplot(down_results_mtp, showCategory=10)
bar_up_mtp = barplot(up_results_mtp, showCategory=10)
barplot_down_mtp = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/bar_down_mts.emf", bar_down_mtp,4 ,5)
barplot_up_mtp = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/bar_up_mts.emf", bar_up_mtp,4 ,5)


# SIGNATURE ANALYSIS ####
master = merge(de_mt_p, de_mt_s, by.x = 0 , by.y = 0, suffixes = c(".mtp",".mts"))
master = merge(master, de_s_p, by.x = 1, by.y = 0, suffixes = c(".merged", ".sp"))
master = merge(master, master_temp, by.x = 1, by.y = 1)

master$sig.mts = factor(master$p.adj.mts < 0.001 & abs(master$log2fold.mts) > 2)
master$sig.mtp = factor(master$p.adj.mtp < 0.001 & abs(master$log2fold.mtp) > 2)
master$sig = factor(master$p.adj < 0.001 & abs(master$log2fold) > 2)
rownames(master) = master[,"SYMBOL"]

sig_all = subset(master, sig.mts == TRUE | sig.mtp == TRUE | sig == TRUE )
vec = rownames(sig_all)
sig_all_expression = em_scaled[vec,]


sig1 = row.names(subset(master,(sig.mts == TRUE & log2fold.mts > 2) & (sig == TRUE & log2fold > 2)))

sig1_expression = em_scaled[sig1,]

sig2 = row.names(subset(master,(sig.mts == TRUE & log2fold.mts < -2) & (sig == TRUE & log2fold < -2)))

sig2_expression = em_scaled[sig2,]

sig3 = row.names(subset(master,(sig.mts == TRUE & log2fold.mts < -2 ) & (sig == TRUE & log2fold > 2)))

sig3_expression = em_scaled[sig3,]

sig4 = row.names(subset(master,(sig.mts == TRUE & log2fold.mts > 2 ) & (sig == TRUE & log2fold < -2)))

sig4_expression = em_scaled[sig4,]


# CLUSTERED HEATMAP ###
all_heatmap = heatmap(sig_all_expression)
all_heatmap
heat_all = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/heat_all.emf", all_heatmap,5 ,4)


sig1_box = faceted_boxplot(sig1,em_scaled,ss,17)
sig1_box


# getting the cnetplot
cet1 = ora_general(sig1, 1,1)
cet1 = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/cnet1.emf", cet1,5 ,7)

cet2 = ora_general(sig2, 1,1)
cet2 = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/cnet2.emf", cet2,5 ,7)

cet3 = ora_general(sig3, 1,1)
cet3 = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/cnet23.emf", cet3,7 ,7)

cet4 = ora_general(sig4, 1,1)
cet4 = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/cnet4.emf", cet4,5 ,7)


#fold vs fold

#fold_vs = fold_plot(master,log2fold.mts,log2fold.mtp)
ggplot(master, aes(x= log2fold, y = log2fold.mts)) + geom_point()

cor.test(master$log2fold, master$log2fold.mtp)
ggplot(master, aes(x = log2fold.mtp, y = log2fold)) +
  geom_point(aes(colour = 'a'), size = 0.7) +
  geom_point(data = master$sig, aes(colour = 'b'), size = 0.7) + 
  geom_point(data = master$sig.mtp, aes(colour = 'c'), size = 0.7) +
  labs(x = "log2 fold change (Condition 1)", y = "log2 fold change (Condition 2)", title = "Fold Change vs. Fold Change Plot") +
  scale_color_manual(values = c("red", "grey" ,"blue"), guide = FALSE) +  # Specify colors for significant and non-significant points
  theme_minimal()

mtp_sp = ggplot(master, aes(x = log2fold.mtp, y = log2fold)) +
  geom_point(data = master, aes(colour = ifelse(as.logical(sig) & as.logical(sig.mtp), "Both", ifelse(as.logical(sig), "Sig", ifelse(as.logical(sig.mtp), "MTP", "Other")))), size = 0.7) +
  labs(x = "log2 fold change Prolif vs Sen_MTD", y = "log2 fold change Prolif vs Senes", title = "Fold Change vs. Fold Change Plot") +
  scale_color_manual(values = c("Sig" = "violetred", "MTP" = "lightseagreen", "Both" = "orange", "Other" = "black"),guides(color = ""),
                     labels = c("SIG_SP","SIG_MTP","SIG_BOTH","SIG_NONE")) +  # Specify colors for significant and non-significant points
  theme_grey() +
  theme(
    text = element_text(family = "Times New Roman"),  # Set font family
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title size and style
    axis.title = element_text(size = 12),  # Adjust axis label size
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 12),  # Adjust legend text size
    legend.position = "right",  # Position legend
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Adjust panel border
    panel.spacing = unit(1, "lines")
  )


mtp_sp_fold = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/mtp_sp.emf", mtp_sp,4 ,5)

mts_mtp = ggplot(master, aes(x = log2fold.mts, y = log2fold.mts)) +
  geom_point(data = master, aes(colour = ifelse(as.logical(sig) & as.logical(sig.mtp), "Both", ifelse(as.logical(sig), "Sig", ifelse(as.logical(sig.mtp), "MTP", "Other")))), size = 0.7) +
  labs(x = "log2 fold change Senes vs Sen_MTD", y = "log2 fold change Prolif vs Senes_MTD", title = "Fold Change vs. Fold Change Plot") +
  scale_color_manual(values = c("Sig" = "violetred", "MTP" = "lightseagreen", "Both" = "orange", "Other" = "black"),guides(color = ""),
                     labels = c("SIG_SP","SIG_MTP","SIG_BOTH","SIG_NONE")) +  # Specify colors for significant and non-significant points
  theme_grey()+
  theme(
    text = element_text(family = "Times New Roman"),  # Set font family
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title size and style
    axis.title = element_text(size = 12),  # Adjust axis label size
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 12),  # Adjust legend text size
    legend.position = "right",  # Position legend
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Adjust panel border
    panel.spacing = unit(1, "lines")
  )

mtp_mts_fold = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/mts_mtp.emf", mts_mtp,4,5)

#metagene violin plot ####
meta1 = create_metagene_plot(em_scaled, sig1, ss,"Signature 1 Violin Plot")
meta1
meta2 = create_metagene_plot(em_scaled, sig2, ss,"Signature 2 Violin Plot")
meta2
meta3 = create_metagene_plot(em_scaled, sig3, ss,"Signature 3 Violin Plot")
meta3
meta4 = create_metagene_plot(em_scaled, sig4, ss,"Signature 4 Violin Plot")
meta4
metasig1 = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/metasig1.emf", meta1,2 ,4)
metasig2 = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/metasig2.emf", meta2,2 ,4)
metasig3 = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/metasig3.emf", meta3,2 ,4)
metasig4 = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/metasig4.emf", meta4,2 ,4)

#VENN ####
mts = nrow(subset(master,sig.mts == TRUE))
mtp = nrow(subset(master,sig.mtp == TRUE))
sp = nrow(subset(master,sig == TRUE))
overlap_mts_sp = nrow(subset(master, sig.mts == TRUE & sig == TRUE))
overlap_sp_mtp = nrow(subset(master, sig.mtp == TRUE & sig == TRUE))
overlap_mts_mtp = nrow(subset(master, sig.mts == TRUE & sig.mtp == TRUE))
overlap =  nrow(subset(master, sig.mts == TRUE & sig.mtp == TRUE & sig == TRUE))

mts = rownames(master,sig.mts == TRUE)
mtp = rownames(master,sig.mtp == TRUE)
sp = rownames(master,sig == TRUE)
copy(mts)


# RUG ####
groups_data = as.matrix(as.numeric(as.factor(ss$SAMPLE_GROUP)))
groups_data = melt(groups_data)

colours = c("olivedrab4", "steelblue", "violetred")
rug = ggplot(groups_data, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(linetype = "blank") +
  scale_fill_gradientn(colours = colours) +
  labs(x = "",
       y = "") +
  theme_classic()
rug
rugs = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/rug.emf", rug,3,5)

# bar plot for up and down reg genes ####
plot_data = data.frame(
  Group = rep(c("SP", "MTS", "MTP"), each = 2),
  Regulation = rep(c("Up", "Down"), times = 3),
  Count = c(862,429,466,1106,1158,1021)
)

barupdown = ggplot(plot_data, aes(x = Group, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Up and Down-regulated Genes",
       x = "Groups",
       y = "Count") +
  scale_fill_manual(values = c("Up" = "violetred", "Down" = "steelblue")) +
  theme_grey()+
  theme(
    text = element_text(family = "Times New Roman"),  # Set font family
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title size and style
    axis.title = element_text(size = 12),  # Adjust axis label size
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 12),  # Adjust legend text size
    legend.position = "right",  # Position legend
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Adjust panel border
    panel.spacing = unit(1, "lines")
  )
barupdownplot = save_image("C:/Users/Pratyasha/OneDrive/Desktop/data ex/barupdown.emf", barupdown,3,4)
