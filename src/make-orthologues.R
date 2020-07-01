library(RMySQL)
library(glue)
library(dplyr)
library(tidyr)


get_one2one_orthologues <- function(species_1, species_2) {
  # connect to database
  host <- 'ensembldb.ensembl.org'
  user <- 'anonymous'
  port <- 3306
  database <- 'ensembl_compara_75'
  
  con <- dbConnect(RMySQL::MySQL(), dbname = database, host=host, user=user, port=port)
  
  
  # determine the genome IDs for our species
  species_genome_db_ids <- dbGetQuery(con, glue('SELECT * FROM genome_db WHERE name = "{species_1}" OR name = "{species_2}"'))
  stopifnot(nrow(species_genome_db_ids)==2)
  taxon_id_1 <- species_genome_db_ids$taxon_id[species_genome_db_ids$name==species_1]
  taxon_id_2 <- species_genome_db_ids$taxon_id[species_genome_db_ids$name==species_2]
  
  # determine the species set/species set ID that corresponds to this relationship
  # get all the species sets
  species_sets <- dbGetQuery(con, 'SELECT * FROM species_set')
  
  # now narrow down to the one containing our species of interest
  tmp <- species_sets %>% dplyr::mutate(is_one_of_our_species = genome_db_id %in% species_genome_db_ids$genome_db_id) %>%
    dplyr::group_by(species_set_id) %>% 
    dplyr::summarize(number_total_species=n(), number_our_species=sum(is_one_of_our_species)) %>%
    dplyr::filter(number_total_species == 2 & number_our_species == 2)
  
  # verify that we found one and only one such set
  stopifnot(nrow(tmp)==1)
  
  our_species_set <- tmp$species_set_id
  
  
  # we'll want the Ensembl orthologues
  tmp <- dbGetQuery(con, 'SELECT * FROM method_link WHERE type="ENSEMBL_ORTHOLOGUES"')
  stopifnot(nrow(tmp)==1)
  method_link_id <- tmp$method_link_id
  
  # now get the method link species set id
  cmd <- glue("SELECT * FROM method_link_species_set WHERE method_link_id = {method_link_id} AND species_set_id = {our_species_set}")
  tmp <- dbGetQuery(con, cmd)
  stopifnot(nrow(tmp)==1)
  method_link_species_id <- tmp$method_link_species_set_id
  
  
  # now grab the homologies

human_rat <- get_one2one_orthologues('homo_sapiens', 'rattus_norvegicus')

human_rat_gene_name <- human_rat %>%
  dplyr::select(-ensembl_gene_id) %>%
  tidyr::spread(key = species, value = gene_name) %>%
  dplyr::select(-homology_id) %>%
  unique() %>%
  dplyr::rename(hg19=homo_sapiens, rn6=rattus_norvegicus) %>%
  dplyr::filter(!is.na(hg19)) %>%
  dplyr::filter(!is.na(rn6))
counts_rat <- table(human_rat_gene_name$rn6)
counts_human <- table(human_rat_gene_name$hg19)
human_rat_gene_name <- human_rat_gene_name[human_rat_gene_name$hg19 %in% names(counts_human[counts_human==1]),]
human_rat_gene_name <- human_rat_gene_name[human_rat_gene_name$rn6 %in% names(counts_rat[counts_rat==1]),]


human_rat_ensembl <- human_rat %>%
  dplyr::select(-gene_name) %>%
  tidyr::spread(key = species, value = ensembl_gene_id) %>%
  dplyr::select(-homology_id) %>%
  unique() %>%
  dplyr::rename(hg19=homo_sapiens, rn6=rattus_norvegicus)

write.table(human_rat_gene_name, file = 'human-rat.gene_name.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = T)
write.table(human_rat_ensembl, file = 'human-rat.ensembl.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = T)
