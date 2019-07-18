args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2)
{
  stop("Usage: Rscript script <kfile> <outputname>", call.=FALSE)
}
if (!file.exists(args[1]) | !file.exists("rnaseq_config.yaml"))
{
  cat("\n Input file(s) not found!")
  stop("Usage: Rscript script <kfile> <outputname>", call.=FALSE)
}
suppressPackageStartupMessages({
library(tximport)
library(readr)
library(yaml)
})

kfile = args[1]
input = read.table( kfile,
                      sep="\t",header=T,stringsAsFactors=F)
files = input$path
names(files) = input$names

species = gsub("kfile_","",basename(kfile))
species = gsub(".txt","",species)
config = read_yaml("rnaseq_config.yaml")

tfile = config[[paste0(species,"_tx2gene")]]
txdb = read.table(tfile,sep="\t",header=T,
                        stringsAsFactors = FALSE)
cat("\ntxdb:\n")
print(head(txdb))

outputname = args[2]
cat("\noutputname:",outputname)
df = tximport(files,
             type = "kallisto", tx2gene = txdb,
             ignoreTxVersion = TRUE,
             countsFromAbundance = "no")
#colnames(df$counts) = outputname
#colnames(df$abundance) = outputname
mainDir = getwd()
subDir="abundance"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)

write.table(df$counts,file=paste0("abundance/kallisto_",outputname,"_counts.tsv"),
            sep="\t",quote=F,col.names=NA)
write.table(df$abundance,file=paste0("abundance/kallisto_",outputname,"_abundances.tsv"),
            sep="\t",quote=F,col.names=NA)
