################################################
## function to compile overlap with GWAS catalog
## at different ld thresholds

map.gwas.catalog <- function(gwas.catalogue, ukb.snps, res.var, ld.proxies, ld.thr) {
  
  # Subset for the rigt proxies
  sub.proxies <- ld.proxies %<>% filter(R2 >= ld.thr)
  
  # Loop over dataframe and, 
    # 1. Identify all proxies, 
    # 2. Map to GWAS Catalogue findings
    # 3. Reduce redundancy
  
  res.gwas <- bettermc::mclapply(1:nrow(res.var), function(i) {
    
    # Right variant
    leadvariant <- res.var[[i, 'rsid']]
    
    # LD proxies for the right region and phenotype
    sub.proxies %>% 
      filter(leadvar == leadvariant) %>% 
      pull(proxy) -> proxies
    
    # Add the lead variant too. # Lready in there
    # snps <- c(leadvariant, proxies)
    snps <- proxies # proxies contains the lead variant here
    
    gwas.catalogue %>% filter(SNPS %in% snps)  -> sub.gwas.catalogue
    
    # Also possibly merge on chromosome and position? e.g. markerName
    ukbb.snps %>% filter(rsid %in% snps) %>% mutate(snp.id = paste0(chromosome, ':', position)) %>% pull(snp.id) -> tmp
    gwas.catalogue %>% filter(snp.id %in% tmp) -> sub.gwas.catalogue2
    
    sub.gwas.catalogue <- rbind(sub.gwas.catalogue, sub.gwas.catalogue2)
    
    if (nrow(sub.gwas.catalogue) > 0) {
      info              <- data.table(leadvariant           = leadvariant,
                                      rsid.gwas             = paste(sort(unique(sub.gwas.catalogue$SNPS)), collapse = "||"),
                                      trait_reported        = paste(sort(unique(sub.gwas.catalogue$TRAIT)), collapse = "||"),
                                      mapped_trait          = paste(sort(unique(sub.gwas.catalogue$MAPPED_TRAIT)), collapse = "||"),
                                      source_gwas           = paste(sort(unique(sub.gwas.catalogue$PUBMEDID)), collapse = "||"),
                                      num_reported          = nrow(sub.gwas.catalogue))
      
    } else {
      info              <- data.table(leadvariant           = leadvariant,
                                      rsid.gwas             = "",
                                      trait_reported        = "",
                                      mapped_trait          = "",
                                      source_gwas           = "",
                                      num_reported          = 0)
      
    }
    
  }, mc.cores = 3) %>% do.call(rbind, .) %>% return()
  
}
