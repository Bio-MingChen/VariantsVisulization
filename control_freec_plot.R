##提供如下功能
#根据control-freec绘制各条染色体的倍数图
#提供log2转换
#提供指定区间的CNV图

library(optparse)
option_list <- list(
  make_option(c("-f","--file"),type="character",
              help="control-freec生成的倍数文件"),
  make_option(c('-b','--bed'),type="character",
              help="指定绘制区间的bed文件")
)
opt <- parse_args(OptionParser(option_list=option_list))
if(opt$bed)
ratio_data <- read.table('chr22.bam_ratio.txt',header=TRUE,
                     sep="\t",stringsAsFactors = FALSE)
#物种染色体倍数
ploidy = 2

head(ratio_data)
clean_data <- ratio_data[ratio_data$CopyNumber != -1 & ratio_data$Ratio != -1,]
head(clean_data)
head(clean_data$Ratio[clean_data$Ratio < 0])
max_limit <- max(clean_data$CopyNumber) + 1
library(ggplot2)
library(RColorBrewer)
#设置科学记数法位数
opt <- options()
options(scipen = 100)

ggplot(data=clean_data) + 
  geom_point(aes(x=Start,y=Ratio*ploidy,color=factor(CopyNumber)),stat="identity") +
  scale_y_continuous(breaks =seq(0,6),labels=seq(0,6),limits=c(0,max_limit)) +
  scale_color_discrete('CopyNumber')
  # scale_colour_manual(values=c('black','green','red')) +
  theme_bw() + 
  labs(title = "Copy Number Variants Visulization",
       x=sprintf("Position: Chr%s", 22),
       y="Normalized Copy Number Profile"
       ) +
  theme(
    plot.title = element_text(size=20,hjust=0.5,face='bold'),
    #legend.position = 'none'
    )
?guides()
#还原科学计数法位数
optionstheme_bw()
