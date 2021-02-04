
# Import and tidy data ---------------------------------------------------

#install.packages("tidyverse")
library(tidyverse)


path <- file.path(getwd(), "data", "counts.xander.UI.txt")

counts.x <- read.table(file = path,
                     sep = '\t',
                     header = FALSE)
#sep = 't\' tells R that the file is tab-delimited. 
#The second line looks "back" in another folder (up in the folder tree).




counts.x <- tail(counts.x, -4)
#Returns the first or last parts of a vector, matrix, table, data frame or function. WHY -4???

counts.x <- counts.x[ , c(1, seq(4, 24, 4))]
## Super cool: selects only column 3 alignments, starting with 4 and alternating every 4!

#seq creates a sequence from one number to another, in steps shown by argument 3.

row_names <- counts.x[1]
#only 1 gene list needed
row_names <- str_sub(row_names$V1, 6, -1) #tidyverse

#str_sub: substring, (input, start, finish)

rownames(counts.x) <- row_names

##create columns
counts.x <- counts.x[ , 2:7]
#note on 2:7, 7 columns including gene name.  
##why is the head at 6 instead of 5? Wouldn't that miss the first gene? 
#these steps that the row names are the first line of the counts matrix (line 13)


col_names <- c("SRR3119158", "SRR3119165", "SRR3119172", "SRR3119178",
               "SRR3119184", "SRR3119190")

# Data frame, rows = observation, column = variables
colnames(counts.x) <- col_names

#note change to code!
#path.met <- file.path(getwd(), "data", "xander_metadata.csv")
#metadata <- read.table(file = path.met,
#                       sep = '\t',
   #                    header = TRUE)
#subset_meta <- metadata[metadata$name == "Unilateral_Incompatability" |
      #                    metadata$name == "Interspecific_Compatability", ]

#counts.x <- counts.x %>% select(subset_meta$run) #tidyverse


colnames(counts.x) <- c("Unilateral_Incompatability_rep_1",
                      "Unilateral_Incompatability_rep_2",
                      "Unilateral_Incompatability_rep_3",
                      "Interspecific_Compatability_rep_1",
                      "Interspecific_Compatability_rep_2",
                      "Interspecific_Compatability_rep_3")

coldata <- data.frame(row.names = colnames(counts.x),
                      condition = factor(c(rep("pen_lyc", 3), 
                                           rep("lyc_pen", 3))))
coldata$condition <- factor(coldata$condition, levels = c("pen_lyc", "lyc_pen"))

# DESeq Object  -----------------------------------------------------------

##Need to download DESeq2, DESeq did not work with my version of R
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("DESeq2")

library(DESeq2)

dds.x <- DESeqDataSetFromMatrix(countData = counts.x,
                              colData = coldata,
                              design = ~ condition)

# Calculate differential expression ---------------------------------------
dds.x <- DESeq(dds.x)

res <- results(dds.x)

res
summary(res)

#Shrink log fold change to make more accurate with genes with small counts 
#normalizes data

##need to install apeglm first:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("apeglm")

resultsNames(dds.x)

resLFC <- lfcShrink(dds.x,
                    coef = "condition_lyc_pen_vs_pen_lyc",
                    type = "apeglm")
####NOTE: this switched the conditions!!! Should i be concerned?

#LFC log fold change

# VISUALIZATION -----------------------------------------------------------

##MA Plot, for showing log fold change vs mean of normalized counts

plotMA(res, ylim=c(-2,2))  #original not adjusted lfc
plotMA(resLFC, ylim=c(-2,2)) #with shrunken lfc

##PCA
#Transform data with variance stabilizing transformation to remove dependence on 
#variance on the mean.  Useful for visualization and downsream analysis
#Basic PCA remove returnData = True

vsd <- vst(dds.x, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition, shape = condition)) +
geom_point(aes(color = condition, shape = condition), size = 5) +
labs(title = "S. lyc. pistil with S. pen. pollen (IC)/S. pen. pistil with S. lyc pollen
     PCA", 
     x = paste0("PC1: ",percentVar[1],"% variance"), 
     y = paste0("PC2: ",percentVar[2],"% variance")) +
scale_color_manual(values = c("#0072B2", "#009E73")) +
theme_bw() +
theme(axis.title = element_text(size = 26, face = 'bold'),
      axis.text = element_text(size = 22, face = 'bold', color = 'black'),
      axis.text.x = element_text(size = 20, face = 'bold', color = 'black'),
      plot.title = element_text(size = 28, face = 'bold', margin = margin(0, 0, 10, 0)),
      axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
      panel.border = element_blank(),
      axis.line = element_line(size = 1, color = 'black'),
      axis.ticks = element_line(size = 1, color = 'black'),
      axis.ticks.length = unit(8, 'pt'),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 20, face = 'bold'),
      legend.position = 'right',
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(size = 26, face = 'bold'))
ggsave(filename = './PCA_IC_vs_UI_S_lyc',
       device = 'png',
       width = 9,
       height = 7.5,
       dpi = 400,
       units = 'in')
