BINARY_DHS_MATRIX <- '/lab/work/porchard/sn-muscle-project/data/meuleman/dat_bin_FDR01_hg38.RData'
METADATA <- '/lab/work/porchard/sn-muscle-project/data/meuleman/DHS_Index_and_Vocabulary_metadata.tsv'
DHS_INDEX <- '/lab/work/porchard/sn-muscle-project/data/meuleman/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz'
dhs_index <- read.table(gzfile(DHS_INDEX), head=T, sep='\t')
dhs_index$peak <- paste(dhs_index$seqname, paste(dhs_index$start, dhs_index$end, sep='-'), sep=':')
indices <- dhs_index$identifier

load(BINARY_DHS_MATRIX) # object: dat_bin

sample_info <- read.table(METADATA, head=T, as.is=T, sep='\t', quote = '', comment.char = '')
sample_info$library_name <- paste(sample_info$Biosample.name, sample_info$Altius.Biosample.ID, sep='.')
stopifnot(all(sample_info$library_name==colnames(dat_bin)))

INDEX_OF_INTEREST <- '17.712748' # NOG SNP peak
INDEX_OF_INTEREST <- '17.712715' # other nog SNP
tmp <- with(dhs_index[dhs_index$identifier==INDEX_OF_INTEREST,], glue('{seqname}:{start}-{end}'))
tmp <- dat_bin[rownames(dat_bin)==tmp,]
number_samples_with_peak <- sum(tmp)
sample_info$has_peak <- sample_info$library_name %in% names(tmp[tmp==T])
sample_info %>%
  dplyr::group_by(Organ, has_peak) %>%
  dplyr::summarize(count=n()) %>%
  as.data.frame()
