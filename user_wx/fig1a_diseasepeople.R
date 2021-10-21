fig1a <- readRDS("~/ses-1/user_wx/fig1a_withdiseasedpeople.rds")
res = fig1a
res$pval2 = -log10(res$p)
res$pval2[which(res$pval2>-log10(0.0001))] = -log10(0.0001)
res$treatment = as.character(res$treatment)
res$gene_set_name = as.character(res$gene_set_name)
res$treatment = factor(res$treatment, levels = c("SES Composite", "Education","Income","Occupation","Subjective Social Status"))
res$gene_set_name = factor(res$gene_set_name, levels = c("Rheumatoid Arthritis","Hypertension","Diabetes","Depression","CVD","COPD","CKD","Asthma","Aortic Aneurysm", "Alzheimers","1KI"))
colnames(res) = c(colnames(res)[1:3], "Cluster", colnames(res)[5])
res = res[order(res$Cluster),]
res$Instance = 0
for (i in 1:length(res$treatment)) {
res$Instance[i] = length(which(res$treatment==res$treatment[i] & res$gene_set_name==res$gene_set_name[i]))
}
res$gene_set_name = factor(res$gene_set_name, levels = rev(levels(res$gene_set_name)))
p = ggplot(res[res$Instance==1,], aes(x = treatment,y  = gene_set_name, color = Cluster, fill = Cluster,  size = pval2)) +
geom_point(alpha  = 0.5,shape = 21, stroke = 1) +
geom_point(data = res[res$Instance==2 & !duplicated(res[,c(1,2)]),], inherit.aes = T, position = position_nudge(y = +0.15), alpha = 0.5, stroke = 1,show.legend = F) +
geom_point(data = res[res$Instance==2 & duplicated(res[,c(1,2)]),], inherit.aes = T, position = position_nudge(y = -0.15), alpha = 0.5, stroke = 1,show.legend = F) +
theme_bw(base_size = 12) +
scale_color_manual(values = c("#01665e","#ff7f0e"), name = "1KI Genes") +
scale_fill_manual(values = c("#01665e","#ff7f0e"), name = "1KI Genes") +
scale_size_continuous(range = c(1,12),
name = "Adjusted\n p-value",
limits = c(-log10(0.05), -log10(0.0001)),
breaks = c(-log10(0.05),-log10(0.01),-log10(0.001), -log10(0.0001)),
labels = c("p<0.05", "p<0.01","p<0.001","p<0.0001")) +
xlab("SES Indicators") + ylab("mRNA Signatures") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
theme(axis.title = element_text(size = 11, face = "bold", family = "Calibri")) +
theme(axis.text.y = element_text(size=10, face = "bold", family = "Calibri")) +
theme(axis.text.x = element_text(size=10,face = "bold" ,family = "Calibri", angle = 30,hjust = 1)) +
theme(legend.text=element_text(size=10, family = "Calibri")) +
theme(legend.title = element_text(size = 10, face = "bold", family = "Calibri")) +
theme(legend.position="right") +
scale_y_discrete(limits = rev)  +
guides(color = guide_legend(override.aes = list(size = 7))) +
guides(fill = guide_legend(override.aes = list(size = 7))) +
guides(shape = guide_legend(override.aes = list(size = 7)))



tiff(str_c(here(), "1A_diseasepeople.tiff"), units="px", width=(6*600), height=(6*750), res=600)
p
dev.off()
