library(tidyverse)
library(jsonlite)
library(GenomicRanges)
options(stringsAsFactors = F)

import_tidal <- function(x) {
  z <- paste0(x,"/*result/*Inserts_Annotated.txt")
  map_df(Sys.glob(z), read_tsv) %>%
  mutate(strand = as.character(sign(TE_coord_end - TE_coord_start))) %>%
  mutate(strand=ifelse(strand > 0 ,'+','-')) %>%
  dplyr::select(Chr,Chr_coord_5p,Chr_coord_3p,strand,TE) %>%
  distinct() %>%
  dplyr::rename(start = Chr_coord_5p, end=Chr_coord_3p)
}

dgrp_ins <- import_tidal(snakemake@input[['dgrp']])

mods <- jsonlite::read_json(snakemake@input[['geps']])

lookup <- read_tsv(snakemake@input[['lookup']]) %>%
  dplyr::select(merged_te,Flybase_name,Blast_candidate) %>%
  mutate(Flybase_name = ifelse(is.na(Flybase_name),Blast_candidate,Flybase_name)) %>%
  dplyr::select(merged_te,Flybase_name) %>%
  distinct()

top_mods <- mods %>%
  enframe() %>% unnest(cols =c(value)) %>%
  mutate_at(vars('value'),unlist) %>%
  filter(!str_detect(value,"FBgn")) %>% group_by(name) %>%
  tally() %>%
  arrange(desc(n))

top_mod <- unlist(mods[[top_mods$name[1]]]) %>%
  .[!str_detect(.,"FBgn")]

top_mod_df <- top_mod %>% tibble(merged_te = .) %>%
  left_join(lookup)

# get sizes of chroms
sizes <- read_tsv(snakemake@input[['sizes']],col_names = c('Chr','Chr_len'))

# sum autosomes, no chrom4 because ancient sex chrom
sizes <- sizes %>%
  filter(Chr!='chr4') %>%
  mutate(Chr_type = ifelse(!Chr %in% c("chrX",'chrY'),'autosome',Chr)) %>%
  group_by(Chr_type) %>%
  summarize(Chr_len = sum(Chr_len))

# TART/A/B/C appear to be inconsistently reconciled by TIDAL annotation.
# I will make sure that my top-te mod remains correct by manually changing
# TART-A's REpbase name to TART-A in the lookup.

lookup <- lookup %>%
  mutate(Flybase_name = ifelse(merged_te == 'TART-A','TART-A',Flybase_name))

lookup <- lookup %>%
  mutate(Flybase_name = ifelse(merged_te == 'TART-B','TART-B',Flybase_name))

lookup <- lookup %>%
  mutate(Flybase_name = ifelse(merged_te == 'TART-C','TART-C',Flybase_name))

#rename of stalker3t to stalker3: https://github.com/bergmanlab/transposons/blob/master/current/transposon_sequence_set.readme.txt

lookup <- lookup %>%
  mutate(Flybase_name = ifelse(merged_te == 'Stalker3','Stalker3T',Flybase_name))

lookup <- lookup %>%
  mutate(Flybase_name = ifelse(merged_te == 'TLD2','TLD2_LTR',Flybase_name))

ins <- ins %>%
  mutate(TE = ifelse(TE=="jockey","Jockey",TE))

# only include discoverable TEs (ie TEs that remain in our scRNA dataset)
ins <- ins %>% filter(TE %in% lookup$Flybase_name)


# only consider DGRP
# do not consider 4, which may have sex chrom origins
sex_auto_ratio_df <- ins %>%
  filter(Chr!='chr4')

# count insertions
sex_auto_ratio_df <- group_by(sex_auto_ratio_df, TE, Chr) %>%
  tally() %>%
  ungroup() %>%
  complete(Chr,TE,fill=list(n=0)) %>%
  mutate_at(vars(Chr),as.character)

# normalize to length
sex_auto_ratio_df <- mutate(sex_auto_ratio_df, Chr_type = ifelse(!Chr %in% c("chrX",'chrY'),'autosome',Chr)) %>%
  group_by(TE,Chr_type) %>%
  summarize(n=sum(n)) %>%
  left_join(sizes) %>%
  mutate(ins_per_mb = n/(Chr_len/1e6)) %>%
  ungroup() %>%
  dplyr::select(-Chr_len,-n) %>%
  spread(Chr_type,ins_per_mb)

sex_auto_ratio_df2 <- sex_auto_ratio_df %>%
  #filter(chrY > 0) %>%
  filter((autosome) > 0 & (chrX) > 0) %>%
  mutate(chrX_Auto_ratio = (chrX)/(autosome)) %>%
  dplyr::select(TE,chrX_Auto_ratio) %>%
  mutate(is_top_mod = TE %in% top_mod_df$Flybase_name) %>%
  gather(comparison,ratio,-TE, -is_top_mod)

sex_auto_ratio_df2 <- sex_auto_ratio_df2 %>%
  dplyr::rename(GEP = 'is_top_mod') %>%
  mutate(GEP = ifelse(GEP,'TEP','other')) %>%
  filter(comparison=='chrX_Auto_ratio') %>%
  mutate(group = snakemake@params[['group']])

write_csv(sex_auto_ratio_df2, snakemake@output[[1]])
