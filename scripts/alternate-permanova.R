ait <- ps %>%
  # transform to relative abundance
  tax_transform("compositional", keep_counts = FALSE) %>%
  dist_calc("aitchison")


# test beta dispersion
ait %>% dist_bdisp(variables = "Male.Female") %>% bdisp_get() # p=0.261

# test with PERMANOVA
mod1 <- ait %>%
  dist_permanova(
    seed = 81299,
    variables = "Male.Female + Run + Batch + affiliation + Conventional.Organic",
    n_perms = 9999
  )

mod1

