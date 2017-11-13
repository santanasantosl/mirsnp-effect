#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-m", "--mirna"), type="character", default=NULL, 
              help="mirna name", metavar="character"),

  make_option(c("-d", "--distribution_dir"), type="character", default=NULL,
              help="folder where distribution files are stored", metavar="character"),

  make_option(c("-r", "--reference"), type="character", default=NULL,
              help="reference allele binding energy", metavar="character"),
  
  make_option(c("-a", "--alternative"), type="character", default=NULL, 
              help="alternative allele binding energy", metavar="character"),
  
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (length(opt) != 6){
  print_help(opt_parser)
  stop("Please provide all necessary (--mirna, --mirna_length, --reference, and --alternative)\n", call.=FALSE)
} # end if

## Get variables
mirna.name=opt$mirna
dist.dir=opt$distribution_dir
reference.energy=opt$reference
alternative.energy=opt$alternative
out.file=opt$out

## Load RDA for background and ratio distributions
#distibution.directory="distribution_rda"
#background.binding.dist = read.table(file.path(distibution.directory, paste("background.dist.", gsub("-", "_", mirna.name), ".txt", sep="")),header=TRUE)
#ratio.binding.dist = read.table(file.path(distibution.directory, paste("ratio.dist.", gsub("-", "_", mirna.name), ".txt", sep="")), header=TRUE)
#difference.binding.dist = read.table(file.path(distibution.directory, paste("difference.dist.", gsub("-", "_", mirna.name), ".txt", sep="")), header=TRUE)

#distibution.directory=paste("background_dist/", mirna.length, "mer", sep="")
background.binding.dist = read.table(file.path(dist.dir, paste("background.dist.", gsub("-", "_", mirna.name), ".txt", sep="")),header=TRUE)
ratio.binding.dist = read.table(file.path(dist.dir, paste("ratio.dist.", gsub("-", "_", mirna.name), ".txt", sep="")), header=TRUE)
#difference.binding.dist = read.table(file.path(distibution.directory, paste("difference.dist.", gsub("-", "_", mirna.name), ".txt", sep="")), header=TRUE)

#distibution.directory="distribution_rda"
#background.binding.dist = read.table(file.path(distibution.directory, paste("expanded.background.dist.", gsub("-", "_", mirna.name), ".txt", sep="")),header=TRUE)
#ratio.binding.dist = read.table(file.path(distibution.directory, paste("expanded.ratio.dist.", gsub("-", "_", mirna.name), ".txt", sep="")), header=TRUE)
#difference.binding.dist = read.table(file.path(distibution.directory, paste("expanded.difference.dist.", gsub("-", "_", mirna.name), ".txt", sep="")), header=TRUE)

ecdf.fn = ecdf(background.binding.dist[,1])
ecdf.ratio.fn = ecdf(ratio.binding.dist[,1])


### Get p-value from background distribution
pvalue.reference = ecdf.fn(reference.energy)
pvalue.alternative = ecdf.fn(alternative.energy)

### Get p-value for the ratio
log.pvalue.ratio = abs(log(pvalue.reference/pvalue.alternative))
#pvalue.log.ratio = 1
pvalue.log.ratio = 1 - ecdf.ratio.fn(log.pvalue.ratio)

### Get p-value for the difference

energy.difference = abs(as.numeric(reference.energy) - as.numeric(alternative.energy))
#pvalue.energy.difference = 1 - ecdf.difference.fn(energy.difference)
pvalue.energy.difference = 1

### Format results
results.df = data.frame(ENERGY.REFERENCE=reference.energy,
                        PVALUE.REFERENCE=pvalue.reference, 
                        ENERGY.ALTERNATICE=alternative.energy,
                        PVALUE.ALTERNATIVE=pvalue.alternative,
                        LOG.PVALUE.RATIO=log.pvalue.ratio,
                        PVALUE.OF.LOG.RATIO=pvalue.log.ratio,
                        ENERGY.DIFFERENCE=energy.difference,
                        PVALUE.OF.ENERGY.DIFFERENCE=pvalue.energy.difference)

## Print results and save it to a file
#print(results.df)
write.table(results.df, out.file, row.names=FALSE, quote=FALSE, sep="\t")