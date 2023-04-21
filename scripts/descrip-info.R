# Descriptive Information from Metadata
# April 11, 2023 - SAB

# load phyloseq object
ps <- readRDS("data/full-run/decontam-ps.RDS")

cows <- subset_samples(ps,
              Group == "Cows")

calves <- subset_samples(ps,
                         Group == "Calves")

conv_cow <- subset_samples(conv,
               Group == "Cows")

subset_samples(conv_cow,
               Male.Female == "Male")

subset_samples(conv,
               Herd.Size == "0-50")

subset_samples(conv,
               Herd.Size == "50-100")

subset_samples(conv,
               Herd.Size == "100-150")

subset_samples(conv,
               Herd.Size == "150-200")

subset_samples(conv,
               Herd.Size == "200+")

subset_samples(conv,
               Cultural.Language.Barriers == "Yes")

subset_samples(conv,
               Cultural.Language.Barriers == "No")

subset_samples(conv,
               Formal.Team.Meetings.Frequency == "Never")

subset_samples(conv,
               Formal.Team.Meetings.Frequency == "1 or 2 times/year")

subset_samples(conv,
               Formal.Team.Meetings.Frequency == "Quarterly")

subset_samples(conv,
               Formal.Team.Meetings.Frequency == "Once a month")

subset_samples(conv,
               Formal.Team.Meetings.Frequency == "At least twice a month")
