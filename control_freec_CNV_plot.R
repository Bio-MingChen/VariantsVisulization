library(optparse)
library(ggplot2)
option_list <- list(
  make_option(c('-f','--file'),type="character",help="输入文件"),
  make_option(c('-b',"--bed"),type="character",help="bed文件"),
  make_option(c('-c',"--chromosome"),
              help="需要绘制的染色体号,逗号分隔的字符串"),
  make_option(c('-l',"--log2"),action='store_true',default=FALSE,
              help="是否对数据进行log2处理"),
  make_option(c('-g',"--basic_graph"),action='store_true',default=FALSE,
              help="使用基础绘图系统绘制，默认使用ggplot2绘制"),
  make_option(c("-p","--ploidy"),type="integer",default=2,
              help="物种倍数"),
  make_option(c("-s","--stretch_window"),type="integer",help=''),
  make_option(c("-t","--stretch_threshold"),type="integer",help=''),
  make_option(c('-x','--prefix'),default='foo',help='输出文件的前缀，默认为空')
                
)
#-s/--stretch_window
#将窗口展示为一条线，而不是默认情况下的点，仅在指定了bed文件时使用,
#在选择区间较短时使用可以增加可读性，但是在窗口较大时会大幅增加运行时间，
#默认不使用，使用示例：-s 2500"

#-t/--stretch_threthold
#当bed文件低于某一个数量的行数时进行延展

opt_parser = OptionParser(usage="%prog [options]",option_list)
opt <- parse_args(opt_parser)
print(opt)
if (is.null(opt$file)) {
  print_help(opt_parser)
  stop("--file argument is required")
}

if (is.null(opt$bed)) {
  bed_file <- NULL
} else {
  bed_file <- read.table(opt$bed,sep="\t",header=FALSE,
                         stringsAsFactors = FALSE)
  bed_file <- bed_file[,1:3]
  colnames(bed_file) <- c('Chr','Start','End')
}

if (is.null(opt$chromosome)) {
  chrs = NULL
} else {
  chrs = unlist(strsplit(opt$chromosome,split=",")) #character
}

infile <- read.table(opt$file,sep="\t",header=TRUE,
                     stringsAsFactors = FALSE)

#Ensuring chromosome is character
infile$Chromosome <- as.character(infile$Chromosome)

#去除缺失数据
omit_condition <- which(infile$CopyNumber > 0 &
                          infile$Ratio > 0)
infile <- infile[omit_condition,]
infile$tag[infile$CopyNumber==opt$ploidy] <- sprintf('%s',opt$ploidy)
infile$color[infile$CopyNumber==opt$ploidy] <- 'green'
infile$tag[infile$CopyNumber>opt$ploidy] <- sprintf('>%s',opt$ploidy)
infile$color[infile$CopyNumber>opt$ploidy] <- 'red'
infile$tag[infile$CopyNumber<opt$ploidy] <- sprintf('<%s',opt$ploidy)
infile$color[infile$CopyNumber<opt$ploidy] <- 'blue'
#log2 handling
if(opt$log2 == TRUE) {
  #data handler
  infile$Ratio <- log2(infile$Ratio)
  infile$CopyNumber <- log2(infile$CopyNumber/opt$ploidy)
  #ylab
  ylab = "Normalized Copy Number Profile(log2)"
} else {
  infile$Ratio <- infile$Ratio * opt$ploidy
  #ylab
  ylab = "Normalized Copy Number Profile"
}
## scipen
options(scipen=100)
##plot function
basic_func <- function(clean_data,ofile,xlab,ylab){
  #basic graph plot
  png(filename = ofile, width = 15, height = 12,
      units = "cm", bg = "white", res = 300,type='cairo')
  plot(clean_data$Start,clean_data$Ratio,xlab = xlab,ylab = ylab,pch = ".",col = colors()[88])
  tt <- which(clean_data$tag==sprintf('>%s',opt$ploidy) )
  points(clean_data$Start[tt],clean_data$Ratio[tt],pch = ".",col = colors()[136])
  tt <- which(clean_data$tag==sprintf('<%s',opt$ploidy))
  points(clean_data$Start[tt],clean_data$Ratio[tt],pch = ".",col = colors()[461])
  points(clean_data$Start,clean_data$CopyNumber, pch = ".", col = colors()[24],cex=4)
  title(main="Copy Number Variants Visulization")
  print(sprintf('Plotting %s completed',ofile))
  dev.off()
}

gg_func <- function(clean_data,ofile,xlab,ylab,mn,mx) {

  levels <- c(sprintf('%s',opt$ploidy),sprintf('<%s',opt$ploidy),sprintf('>%s',opt$ploidy))
  gg_break <- factor(unique(clean_data$tag),levels=levels)
  color_order <- c('green','blue','red')
  gg_color <- intersect(color_order,unique(clean_data$color))
  print(gg_break)
  print(gg_color)
  ggplot(data=clean_data) + 
    geom_point(aes(x=Start,y=Ratio,color=tag),stat="identity",shape=15) +
    geom_point(aes(x=Start,y=CopyNumber),color="black",shape=15) +
    scale_y_continuous(breaks =seq(mn,mx),labels=seq(mn,mx),limits=c(mn,mx)) +
    #scale_color_discrete('CopyNumber') +
    scale_colour_manual('CopyNumber',breaks=gg_break,
                        labels=gg_break,values=gg_color ) +
    theme_bw() + 
    labs(title = "Copy Number Variants Visulization",
         x=xlab,
         y=ylab
    ) +
    theme(
      plot.title = element_text(size=20,hjust=0.5,face='bold'),
      legend.position = c(0.9,0.9)
    )
  ggsave(ofile,width=10,height=8,type='cairo')
  print(sprintf('Plotting %s completed',ofile))
}
# bed stretch function
bed_stretch <- function(clean_data,window) {
  new_df <- NULL
  for(i in 1:nrow(clean_data)){
    new_i_df <- NULL
    for (j in 1:window) {
      if(is.null(new_i_df)){
        new_i_df <- clean_data[i,]
      } else {
        new_i_df[j,] <- clean_data[i,]
        new_i_df[j,]$Start <- clean_data[i,]$Start + j - 1
      }
    }
    new_df <- rbind(new_df,new_i_df)
  }
  return(new_df)
}
#Get bed data 
if (!is.null(bed_file)) {
  for(i in 1:nrow(bed_file)) {
    add_region = (bed_file[i,]$End - bed_file[i,]$Start)*1/6
    new_start = bed_file[i,]$Start - add_region
    new_end = bed_file[i,]$End + add_region
    condition = which(infile$Chromosome == bed_file[i,]$Chr &
            infile$Start >= new_start &
            infile$Start <= new_end)

    bed_data <- infile[condition,]
    #y axis range
    mx <- round(max(bed_data$Ratio) + 1)
    mn <- round(min(bed_data$Ratio) - 1)
    #plot function
    xlab = sprintf("Position: Chr%s:%s-%s", bed_file[i,]$Chr,
                   bed_file[i,]$Start,bed_file[i,]$End)
    ofile = sprintf("%s.chr%s_%s_%s.png",opt$prefix,bed_file[i,]$Chr,
                    bed_file[i,]$Start,bed_file[i,]$End)
    #print(bed_data)
    #print(ofile)
    print(!is.null(opt$stretch_window))
    print(!is.null(opt$stretch_threshold))
    if(!is.null(opt$stretch_window) & !is.null(opt$stretch_threshold)) {
      print(nrow(bed_data))
      if(nrow(bed_data)<=opt$stretch_threshold){
        bed_data <- bed_stretch(bed_data,opt$stretch_window)
      }
    } else if(!is.null(opt$stretch_window) & is.null(opt$stretch_threshold)) {
        bed_data <- bed_stretch(bed_data,opt$stretch_window)
    }
    if(opt$basic_graph==FALSE){
      gg_func(bed_data,ofile,xlab,ylab,mn,mx)
    } else {
      basic_func(bed_data,ofile,xlab,ylab)
    }
  }
}

if (!is.null(chrs)) {
  for(chr in chrs) {
    chr_infile <- infile[infile$Chromosome==chr,]
    #plot function
    xlab = sprintf("Position: Chr%s", chr)
    ofile = sprintf("%s.chr%s.png",opt$prefix,chr)
    #y axis range
    mx <- round(max(chr_infile$Ratio) + 1)
    mn <- round(min(chr_infile$Ratio) - 1)
    if(opt$basic_graph==FALSE){
      gg_func(chr_infile,ofile,xlab,ylab,mn,mx)
    } else {
      basic_func(chr_infile,ofile,xlab,ylab)
    }
  }
}

if (is.null(bed_file)& is.null(chrs)) {
  for(chr in unique(infile$Chromosome)) {
    chr_infile <- infile[infile$Chromosome==chr,]
    #plot function
    xlab = sprintf("Position: Chr%s", chr)
    ofile = sprintf("%s.chr%s.png",opt$prefix,chr)
    #y axis range
    mx <- round(max(chr_infile$Ratio) + 1)
    mn <- round(min(chr_infile$Ratio) - 1)
    if(opt$basic_graph==FALSE){
      gg_func(chr_infile,ofile,xlab,ylab,mn,mx)
    } else {
      basic_func(chr_infile,ofile,xlab,ylab)
    }
  }
}
