# Load necessary libraries
library(EnhancedVolcano)
library(ggplot2)

# Read the data from the file
Diff_v_Undiff <- read.delim("Diff_v_Undiff.txt", header = TRUE, sep = "\t")

# Flip the logFC values
Diff_v_Undiff$logFC <- -Diff_v_Undiff$logFC

# Define specific genes to be labeled
specific_genes <- c("MAGEL2", "NDN", "SNRPN", "SNHG14", "UBE3A", "PWAR1", "IPW", "SNORD116-1", "SNORD116-29", 
                    "KIF18B", "IFITM3", "PBK", "H3C8", "AURKB", "H3C3", "BIRC5", "RRM2", "GJA1", "KIF20A", 
                    "H3C2", "LDLRAP1", "UBE2C", "PIMREG", "CA7", "FOXM1", "CCND1", "MKI67", "SLC47A1", "TRG-AS1", 
                    "STMN2", "PARM1", "P2RX3", "SYT4", "SLC17A6", "CNR1", "ELAVL4", "DCC", "NCAN", "CELF3", 
                    "SEZ6L", "SST", "SHISA6", "SCRT1", "EPHA5", "GRIA2", "SYNPR", "DCX", "KCNH7", "SYT13", 
                    "ALCAM", "MAP2", "RTN1", "NCAM1", "CNTN2", "H1-5", "AKAP6", "KIF5A", "SCD5", "ROBO2", 
                    "NRG1", "SEMA6D", "RMC1", "MDGA1", "H2AC11", "ASS1", "TNC", "H2BC18", "NCAPD2", "CCNB1", 
                    "SMC4", "SUSD2", "HMGA2", "CENPF", "AURKB", "BIRC5", "RRM2", "GJA1", "KIF20A", "LDLRAP1", 
                    "UBE2C", "PIMREG", "CA7", "FOXM1", "CCND1", "SLC47A1", "TRG-AS1", "STMN2", "P2RX3", "SYT4", 
                    "SLC17A6", "CNR1", "ELAVL4", "DCC", "NCAN", "CELF3", "SEZ6L", "SST", "SHISA6", "SCRT1", 
                    "EPHA5", "GRIA2", "SYNPR", "DCX", "KCNH7", "SYT13")

# Create the Volcano Plot
volcano_plot <- EnhancedVolcano(Diff_v_Undiff,
                                lab = Diff_v_Undiff$Gene,
                                x = 'logFC',
                                y = 'adj.P.Val',
                                xlim = c(-17, 17),
                                ylim = c(0, 15),
                                title = '',
                                pCutoff = 0.05,
                                FCcutoff = 1.5,
                                pointSize = 2.0,
                                labSize = 4.04,
                                col = c('forestgreen', 'grey30', 'royalblue', 'red2'),
                                legendPosition = 'top',
                                legendLabSize = 12,
                                legendIconSize = 4.0,
                                selectLab = specific_genes,
                                labFace = 'bold')

# Add custom annotations for upregulated and downregulated labels
volcano_plot <- volcano_plot +
  annotate("text", x = 12, y = 13, label = "Upregulated in Neurons", size = 5, fontface = "bold", color = "red", hjust = 0.8) +
  annotate("text", x = -12, y = 13, label = "Downregulated in Neurons", size = 5, fontface = "bold", color = "darkblue", hjust = 0.1)

# Adjusting the y-axis label
volcano_plot <- volcano_plot + ylab(expression("-Log"[10]*" Adjusted P-value"))

# Function to format axis labels
axis_format_x <- function(x) {
  ifelse(x == 0, "", as.character(x))  # Remove '0' label
}

# Function to format y-axis labels
axis_format_y <- function(x) {
  ifelse(x == 1.3, sprintf("%.1f", x), as.character(as.integer(x)))
}

# Adjusting the number of ticks and labels on both axes
volcano_plot <- volcano_plot +
  scale_x_continuous(breaks = c(-12, -8, -4, -1.5, 0, 1.5, 4, 8, 12, 16), labels = axis_format_x) +
  scale_y_continuous(breaks = c(1.3, 3, 6, 9, 12), labels = axis_format_y, limits = c(0, 13))

# Adjust the x-axis label position (Fold Change)
volcano_plot <- volcano_plot + theme(axis.title.x = element_text(vjust = 1, hjust = 0.4))

# Save the plot
ggsave("volcano_plot.pdf",
       plot = volcano_plot,
       width = 8,
       height = 10)
