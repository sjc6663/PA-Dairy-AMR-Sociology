# Descriptive Information from Metadata
# April 11, 2023 - SAB

# load phyloseq object
ps <- readRDS("data/ransom/decontam-ps.RDS")

subset_samples(ps,
              Conventional.Organic == "Conventional")

subset_samples(ps,
               Conventional.Organic == "Organic")

subset_samples(ps,
               Male.Female == "Male")

subset_samples(ps,
               Male.Female == "Female")

subset_samples(ps,
               Herd.Size == "0-50")

subset_samples(ps,
               Herd.Size == "50-100")

subset_samples(ps,
               Herd.Size == "100-150")

subset_samples(ps,
               Herd.Size == "150-200")

subset_samples(ps,
               Herd.Size == "200+")

subset_samples(ps,
               Cultural.Language.Barriers == "Yes")

subset_samples(ps,
               Cultural.Language.Barriers == "No")

subset_samples(ps,
               Formal.Team.Meetings.Frequency == "Never")

subset_samples(ps,
               Formal.Team.Meetings.Frequency == "1 or 2 times/year")

subset_samples(ps,
               Formal.Team.Meetings.Frequency == "Quarterly")

subset_samples(ps,
               Formal.Team.Meetings.Frequency == "Once a month")

subset_samples(ps,
               Formal.Team.Meetings.Frequency == "At least twice a month")
