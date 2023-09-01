ait <- ps2 %>%
  # transform to relative abundance
  tax_transform("compositional", keep_counts = FALSE) %>%
  dist_calc("aitchison")


# test beta dispersion
ait %>% dist_bdisp(variables = "Male.Female") %>% bdisp_get() # p=2.24e-05

# test with PERMANOVA
mod1 <- ait %>%
  dist_permanova(
    seed = 81299,
    variables = "Male.Female + Run + Male.Female*Run",
    n_perms = 9999
  )

mod1

