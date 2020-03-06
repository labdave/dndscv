library(optparse)
option_list <- list ( 
  make_option(c("-o", "--outputdir"), dest = "out_dir", type = "character", help = "output dir name"),
  make_option(c("-r", "--refpath"), dest = "ref_path", type = "character", help = "reference file path"),
  make_option (c("-f","--filelist"),default="blah.txt", help="comma separated list of files (default %default)")
  
)

parser <-OptionParser(option_list=option_list)
arguments <- parse_args (parser, positional_arguments=TRUE)
opt <- arguments$options
args <- arguments$args

#myfilelist <- strsplit(opt$filelist, ",")
myfilelist <- unlist(strsplit(opt$filelist, ","))

output.prefix <- opt$out_dir
supp_file = opt$ref_path

all_vcfs = data.frame()
# Loop through recoded-vcfs converting to dndscv format
for (i in myfilelist) {
  print(paste("Reading", i, "..."))
  x = read.table(i, sep = "\t",header=T, stringsAsFactors = F)
  # Filter based on read-count
  rc = names(x)[(ncol(x)-3)] # (assumes read-count column is 4th from last)
  x = x[which(x[,rc] >= 0.3),] # WARNING: hard-coded!!!
  # Filter based on the mutation being exonic
  x = x[x$ExonicFunc.refGene != ".",]
  # Apply QC filer
  quality.variants.moderate <- !((as.numeric(x$QD)<2 & !is.na(x$QD)) | 
                                   (as.numeric(x$FS)>100 & !is.na(x$FS)) | 
                                   (as.numeric(x$SOR)>3 & !is.na(x$SOR)) | 
                                   (as.numeric(x$MQ)<40 & !is.na(x$MQ)) |
                                   (as.numeric(x$MQRankSum)<(-2.5) & !is.na(x$MQRankSum)) | 
                                   (as.numeric(x$ReadPosRankSum)<(-8) & !is.na(x$ReadPosRankSum)))
  x = x[quality.variants.moderate,]
  # Check if any mutations passed filters
  if (nrow(x) > 0) {
    # Add sample-name column
    sn = substr(rc,2,nchar(rc)) # get sample name (assumes R added an "X" to the column)
    x[,"sampleID"] = sn
    # Create output table
    print(paste("Saving", nrow(x), "variants!"))
    outx = x[,c("sampleID","CHROM","POS","REF","ALT")]
    all_vcfs = rbind(all_vcfs,outx)
  } else {
    print("No variants passed filters.")
  }
  
}
# Format and save input table
names(all_vcfs) = c("sampleID","chr","pos","ref","mut")
write.table(all_vcfs,paste(output.prefix, "input.tsv",sep="."),sep="\t",quote = F,row.names = F,col.names = T)

# Load driver-finding library
library("dndscv")

# Run driver detection
#dndsout = dndscv(all_vcfs,refdb = "~/tmp/RefCDS_human_GRCh38.p12.rda", cv=NULL)
dndsout = dndscv(all_vcfs,refdb = supp_file, cv=NULL)

# Save Outputs
# Main table
write.table(dndsout$sel_cv,paste(output.prefix,"tsv",sep="."),sep="\t",quote=F,col.names = T,row.names = F)
print(output.prefix)
# RDS object
saveRDS(dndsout, file = paste(output.prefix, "RDS",sep="."))