library(optparse)

option_list <- list(
  make_option(c('-f','--file'),type="character",help="输入文件"),
  make_option(c('-b',"--bed"),type="character",help="bed文件"),
  make_option(c('-chr',"--chromosome"),
              help="需要绘制的染色体号,逗号分隔的字符串"),
  makr_option(c('-l',"--log2"),action='store_ture',default=FALSE,
              help="是否对数据进行log2处理")
)
opt_parser = OptionParser(usage="%prog [options]",option_list)
opt <- parse_args(opt_parser)

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
omit_condition <- which(infile$CopyNumber != -1 &
                          infile$Ratio != -1)
infile <- infile[omit_condition,]

#Get bed data 
if (!is.null(bed_file)) {
  for(i in 1:nrow(bed_file)) {
    condition = which(infile$Chromosome == bed_file[i,]$Chr &
            infile$Start >= bed_file[i,]$Start &
            infile$Start <= bed_file[i,]$End)
      
    bed_file <- infile[condition,]
    #plot function
  }
}

if (!is.null(chrs)) {
  for(chr in chrs) {
    chr_infile <- infile[infile$Chromosome==chr,]
    #plot function
  }
}

if (is.null(bed_file)& is.null(chrs)) {
  for(chr in unique(infile$Chromosome)) {
    chr_infile <- infile[infile$Chromosome==chr,]
    #plot function
  }
}