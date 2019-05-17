library(dplyr)
source('genpred_library.R')
#get_data()


map_pheno <- function(pheno_raw, meta_data){
  
 
  pheno_raw <- mutate(pheno_raw, NPPC_code = sapply(pheno_raw[['id']], function(x) 
    as.numeric(strsplit(as.character(x), '\\-')[[1]][2])))
  
  pheno_raw <- merge(pheno_raw, meta_data, by.x = 'NPPC_code', by.y = 'NPPC.code',
                       all.y = FALSE)
  pheno_raw <- pheno_raw[seq(1, nrow(pheno_raw), by = 2),]
  pheno_raw$W8 <- gsub(' ', '_', pheno_raw[['W8']])
  i = which(colnames(pheno_raw) == 'W8')
  colnames(pheno_raw)[i] <- 'ind'
  return(pheno_raw)
  

}

meta_data <- read.csv('W8_BARCODE.csv', header = TRUE)
pheno_raw <- read.csv('growth_model_traits.csv', header = TRUE)
pheno_data <- map_pheno(pheno_raw, meta_data)
write.table(pheno_data[, c(22, 3:20)], file = 'MAGIC/magic_pheno_w8_dl.csv',
            row.names = FALSE, sep = ',')



#load('w8data')
#names(d)
# trait_data <- d$traitdata
# file <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex/DOex.zip"
# DOex <- read_cross2(file)

cl <- makeCluster(detectCores() -1)
registerDoParallel(cl)

magicex <- read_cross2('MAGIC/magic_w8_dl.json')

pr <- calc_genoprob(magicex, error_prob=0.002, cores=2)
save(pr, file = 'pr.dat')
apr <- genoprob_to_alleleprob(pr)
k <- calc_kinship(apr, "loco")
break


pdf('plots.pdf')
par(mfrow = c(4,2))
for(i in 1:ncol(magicex$pheno)){
  
  out <- scan1(apr, magicex$pheno[,i])
  plot(out, magicex$gmap, main = colnames(magicex$pheno)[i],
       las = 2, altcol="green4", gap=0)
  abline(h=3, col="red", lwd = 1.5, lty=2)
  
}
dev.off()


