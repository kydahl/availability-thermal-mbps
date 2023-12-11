library(tidyverse)


# Mordecai 2013 -----------------------------------------------------------
data.Mordecai2013 <- read.csv("data/raw/Mordecai_2013/Mordecai_2013_supp_data.csv") %>% 
  select(trait.name, ref, mosquito_species, pathogen) %>% 
  unique() %>% 
  mutate(source = "Mordecai 2013")


# Mordecai 2017 -----------------------------------------------------------
data.Mordecai2017.Aegypti <- read.csv("data/raw/Mordecai_2017/aegyptiDENVmodelTempData_2016-03-30.csv", header = TRUE) %>%
  # Exclude the Focks & Barrera 2006 data because they're from a model
  filter(ref != "Focks_Barrera_2006_Research&TrainingTropicalDis_Geneva_Paper") %>%
  # Exclude the Rohani et al 2009 data because it has "unrealistically long lifespans"
  filter(ref != "Rohani_et_al_2009_SEJTropMedPH") %>%
  # Exclude outliers from Beserra 2009 data
  filter(!(ref == "Beserra_2009" & trait.name == "EFD" & trait > 15)) %>% 
  mutate(mosquito_species = "Aedes aegypti") %>%
  mutate(pathogen = ifelse(trait.name %in% c("b", "c", "PDR", "EIP"), case_when(
    ref == "Alto&Bettinardi_2013_AJTMH" ~ "DENV",
    ref == "Carrington_et_al_2013_PNTD" ~ "DENV",
    ref == "Davis_1932_AmJEpidemiology" ~ "YFV",
    ref == "Focks_et_al_1995_AJTMH" ~ "DENV",
    ref == "Lambrects_et_al_2011_PNAS" ~ "DENV",
    ref == "McLeah_et_al_1975_MosquitoNews" ~ "DENV",
    ref == "McLean_et_al_1974_CanJMicobiol" ~ "DENV",
    ref == "Watts_et_al_1987_AJTMH" ~ "DENV"
  ), NA)) %>% 
  select(trait.name, ref, mosquito_species, pathogen) %>% 
  unique()

data.Mordecai2017.Albopictus <- read.csv("data/raw/Mordecai_2017/albopictusCHIKVmodelTempData_2016-03-26.csv", header = TRUE) %>%
  # The starved mosquitoes had much shorter survival than all other data, so remove them
  filter(trait2 %in% c("sugar-fed", NA)) %>%
  # This study only considered Aedes aegypti
  filter(ref != "Lambrects_et_al_2011_PNAS") %>% 
  mutate(mosquito_species = "Aedes albopictus") %>%
  mutate(pathogen = ifelse(trait.name %in% c("b", "c", "PDR", "EIP"), case_when(
    ref == "Alto&Bettinardi_2013_AJTMH" ~ "DENV",
    ref == "Westbrook_et_al_2010_VecBorn&ZooDis" ~ "CHIKV",
    ref == "Westbrook_Thesis_2010" ~ "CHIKV",
    ref == "Xiao_et_al_2014_Arch Virol" ~ "DENV"
  ), NA)) %>% 
  select(trait.name, ref, mosquito_species, pathogen) %>% 
  unique()

data.Mordecai2017.AedesSpp <- read.csv("data/raw/Mordecai_2017/Aedes_prior_data.csv", header = TRUE) %>%
  mutate(mosquito_species = "Aedes spp.") %>%
  dplyr::select(-mosquito.species) %>% 
  mutate(pathogen = ifelse(trait.name %in% c("b", "c", "PDR", "EIP"), "other", NA))  %>% 
  select(trait.name, ref, mosquito_species, pathogen) %>% 
  unique()


data.Mordecai2017.PDRaddl <- read.csv("data/raw/Mordecai_2017/EIP_priors_2015-12-04.csv", header = TRUE) %>%
  rename(pathogen = virus, mosquito_species = mosquito) %>% 
  select(trait.name, ref, mosquito_species, pathogen) %>% 
  unique()

data.Mordecai2017 <- rbind(data.Mordecai2017.AedesSpp, 
                           data.Mordecai2017.Aegypti, 
                           data.Mordecai2017.Albopictus, 
                           data.Mordecai2017.PDRaddl) %>% 
  unique() %>% 
  mutate(source = "Mordecai 2017")


# Shocket 2018 ------------------------------------------------------------
data.EFD.proc <- read.csv("data/raw/Shocket_2018/RRVTraitData.csv", header = TRUE) %>%
  # remove McDonald_1980_AusJofEco EFD data which will be dealt with separately
  filter(ref == "McDonald_1980_AusJofEco" & trait.name == "EFD") %>% 
  mutate(trait2 = as.double(trait2)) %>% 
  group_by(T) %>%
  mutate(trait = sum(trait)/max(trait2)) 

data.Shocket2018 <- read.csv("data/raw/Shocket_2018/RRVTraitData.csv", header = TRUE) %>%
  # remove EFD data from Williams_2011_AusJEntrom which did not count actual eggs laid 
  filter(!(ref == "Williams_2011_AusJEntrom" & trait.name == "EFD")) %>% 
  # remove McDonald_1980_AusJofEco EFD data which will be dealt with separately
  filter(!(ref == "McDonald_1980_AusJofEco" & trait.name == "EFD")) %>%
  # add in processed McDonald_1980_AusJofEco EFD data
  rbind(data.EFD.proc) %>% 
  mutate(mosquito_species = case_when(
    host.code == "Ovig" ~ "Ochlerotatus vigilax",
    host.code == "Cann" ~ "Culex annulirostris",
    host.code == "Anot" ~ "Aedes notoscriptus",
    host.code == "Ocam" ~ "Ochlerotatus camptorhynchus"
  )) %>%
  # Rename pEA traits with the actual trait measured
  mutate(trait.name = ifelse(trait.name == "pEA", notes, trait.name)) %>% 
  mutate(lead_author = "Shocket") %>%
  mutate(year = "2018") %>% 
  select(trait.name, ref, mosquito_species, pathogen = paras.code) %>% 
  unique() %>% 
  mutate(source = "Shocket 2018")

# Tesla 2018 --------------------------------------------------------------

data.Tesla2018 <- read.csv("data/raw/Tesla_2018/zikv_traits.csv", header = TRUE) %>%
  # NB: MDR, pEA, EFD, and a are same as Mordecai 2017
  filter(!trait.name %in% c("MDR", "pEA", "EFD", "a")) %>% 
  mutate(mosquito_species = "Aedes aegypti") %>%
  mutate(infection_status = stringr::word(rep, 2, 2, sep = "-")) %>% 
  # Remove extremely long lifespans for infected mosquitoes
  filter(trait.name != "lf" | infection_status != "inf") %>%
  # Although Tesla et al kept track of infection status, they combined lifespan
  # data across infected and uninfected mosquitoes
  mutate(pathogen = case_when(
    trait.name %in% c("bc", "EIR") ~ "ZIKV",
    (infection_status == "inf" & trait.name != "lf") ~ "ZIKV",
    TRUE ~ NA
  )) %>%
  mutate(ref = "Tesla 2018") %>%
  select(trait.name, ref, mosquito_species, pathogen) %>% 
  unique()  %>% 
  mutate(source = "Tesla 2018")


# Shocket 2020 ------------------------------------------------------------

data.Shocket2020 <- read.csv("data/raw/Shocket_2020/TraitData_a.csv", header = TRUE) %>%
  rename(trait.name = Trait.Name) %>%
  dplyr::select(ref = Citation, trait.name, host.code, paras.code) %>% 
  rbind(read.csv("data/raw/Shocket_2020/TraitData_bc.csv", header = TRUE) %>%
          dplyr::select(ref = Citation, trait.name, host.code, paras.code)) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_EFD.csv", header = TRUE) %>%
          rename(trait.name = Trait.Name) %>%
          dplyr::select(ref = Citation, trait.name, host.code, paras.code)) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_EV.csv", header = TRUE) %>%
          # Following Shocket et al 2020: Remove Oda data - it shows no temp-dep, and briere fit doesn't work with it (T0 = 0)
          filter(Citation != "Oda_1980_TropMed") %>% 
          dplyr::select(ref = Citation, trait.name, host.code, paras.code)) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_lf.csv", header = TRUE) %>%
          dplyr::select(ref = Citation, trait.name, host.code, paras.code)) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_MDR.csv", header = TRUE) %>%
          dplyr::select(ref = Citation, trait.name, host.code, paras.code)) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_PDR.csv", header = TRUE) %>%
          rename(trait.name = Trait.Name) %>%
          dplyr::select(ref = Citation, trait.name, host.code, paras.code)) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_pLA.csv", header = TRUE) %>%
          dplyr::select(ref = Citation, trait.name, host.code, paras.code)) %>%
  rename(pathogen = paras.code) %>%
  mutate(mosquito_species = case_when(
    host.code == "Cpip" ~ "Culex pipiens",
    host.code == "Cqui" ~ "Culex quinquefasciatus",
    host.code == "Ctar" ~ "Culex tarsalis",
    host.code == "Cuni" ~ "Culex univittatus",
    host.code == "Cthe" ~ "Culex theileri",
    host.code == "Atri" ~ "Aedes triseriatus",
    host.code == "Atae" ~ "Aedes taeniorhynchus",
    host.code == "Avex" ~ "Aedes vexans",
    host.code == "Cmel" ~ "Culiseta melanura",
    host.code == "Cmol" ~ "Culex pipiens molestus",
    host.code == "Cpal" ~ "Culex pipiens pallens",
    host.code == "Cres" ~ "Culex restuans",
    host.code == "Ador" ~ "Aedes dorsalis",
    host.code == "Anig" ~ "Aedes nigromaculis",
    host.code == "Asol" ~ "Aedes sollicitans",
    host.code == "Asal" ~ "Aedes salinarius",
    TRUE ~ "missing code"
  )) %>%
  filter(mosquito_species != "missing code") %>%
  select(trait.name, ref, mosquito_species, pathogen) %>% 
  unique()  %>% 
  mutate(source = "Shocket 2020")


# New fecundity data ------------------------------------------------------

data.Anga_fecundity <- read.csv("data/raw/Anopheles_fecundity.csv", header = TRUE) %>%
  mutate(lead_author = "Dahlin") %>%
  mutate(year = "2023") %>%
  mutate(mosquito_species = case_when(
    host.code == "Anga" ~ "Anopheles gambiae",
    host.code == "Ansu" ~ "Anopheles superpictus"
  )) %>%
  mutate(pathogen = NA,
         source = "Dahlin 2023") %>% 
  dplyr::select(ref, trait.name, mosquito_species, pathogen, source) %>% 
  unique()


# Villena 2022 ------------------------------------------------------------
data.Villena2022 <- read.csv("data/raw/Villena_2022/traits.csv", header = TRUE) %>%
  # make mosquito names consistent with other data
  mutate(mosquito_species = case_when(
    specie %in% c("An. stephensi") ~ "Anopheles stephensi",
    specie %in% c("An. gambiae") ~ "Anopheles gambiae",
    specie %in% c("An. culicifacies") ~ "Anopheles culicifacies",
    specie %in% c("An. maculipennies") ~ "Anopheles maculipennies",
    specie %in% c("An. quadrimaculatus", "An. Quadrimaculatus") ~ "Anopheles quadrimaculatus",
    specie %in% c("An. albimanus") ~ "Anopheles albimanus",
    specie %in% c("An. arabiensis", "An. Arabiensis") ~ "Anopheles arabiensis",
    specie %in% c("An. funestus") ~ "Anopheles funestus",
    specie %in% c("An. Pseudopunctipennis") ~ "Anopheles pseudopunctipennis"
  )) %>% 
  # make parasite names consistent with other data
  mutate(pathogen = case_when(
    parasite %in% c("P. falciparum") ~ "Plasmodium falciparum",
    parasite %in% c("P. malariae") ~ "Plasmodium malariae",
    parasite %in% c("P. vivax") ~ "Plasmodium vivax",
    parasite %in% c("P. yoelii") ~ "Plasmodium yoelii",
    parasite %in% c("P. berghei") ~ "Plasmodium berghei",
    TRUE ~ "none"
  )) %>% 
  # make trait names consistent with other data
  mutate(trait.name = case_when(
    trait.name == "pdr" ~ "PDR",
    trait.name == "efd" ~ "EFD",
    trait.name == "bc.succ" ~ "bc",
    trait.name == "mdr" ~ "MDR",
    TRUE ~ trait.name
  )) %>% 
  # Turn vector competence into a probability
  mutate(trait = ifelse(trait.name == "bc", trait/100, trait)) %>% 
  # Add lead author and year of publication
  mutate(source = "Villena 2022") %>%
  # Restrict to relevant columns
  dplyr::select(trait.name, mosquito_species, pathogen, ref, source) %>% 
  unique()


# Combine all data --------------------------------------------------------

All_data <- rbind(data.Mordecai2013, data.Mordecai2017, data.Shocket2018, 
                  data.Tesla2018, data.Shocket2020, data.Anga_fecundity, 
                  data.Villena2022) %>% 
  unique() %>% 
  filter(!is.na(ref)) %>% 
  group_by(ref, mosquito_species, pathogen, source) %>% 
  summarise(traits = toString(unique(trait.name)))

write_csv(All_data, "data/raw/raw_references.csv")


# Clean up references again -----------------------------------------------


new_refs <- read_csv("data/clean/references.csv", show_col_types = FALSE) %>% 
  filter(!is.na(`Mosquito species`)) %>% 
  filter(!is.na(`Citation key`)) %>% 
  # Separate out trait names again
  tidyr::separate_rows(`Trait names`, sep=", ") %>% 
  arrange(`Mosquito species`, `Parasite`, trait = `Trait names`, `Citation key`) %>% 
  # Combine equivalent traits
  mutate(trait = case_when(
    `Trait names` %in% c("GCR", "GCD", "a", "") ~ "a",
    `Trait names` %in% c("b") ~ "b",
    `Trait names` %in% c("bc") ~ "bc",
    `Trait names` %in% c("c") ~ "c",
    `Trait names` %in% c("e2a", "Egg viability", "pEA", "pLA", "pRH") ~ "e2a",
    `Trait names` %in% c("EFD", "EFGC", "EFOC", "EPR", "TFD") ~ "EFD",
    `Trait names` %in% c("1/mu", "lf", "mu", "p/days", "prop.dead") ~ "lf",
    `Trait names` %in% c("1/MDR", "MDR", "mdr") ~ "MDR",
    `Trait names` %in% c("nLR") ~ "nLR",
    `Trait names` %in% c("EIP", "EIR", "PDR") ~ "PDR",
    `Trait names` %in% c("pO") ~ "pO"
  )) %>% 
  group_by(`Mosquito species`, `Parasite`) %>% 
  mutate(traits = toString(unique(trait))) %>% 
  select(-trait, -`Trait names`) %>% 
  unique() %>% 
  group_by(`Mosquito species`, `Parasite`, traits) %>% 
  # Make citation key "nice" for input into LaTeX
  summarise(refs = toString(unique(`Citation key`))) %>% 
  mutate(LaTeX_cite = paste0("\\cite{", refs, "}")) %>% 
  select(-refs)

write_csv(new_refs, "data/clean/latex_references.csv")

# Just get references for each genera
genera_refs <- read_csv("data/clean/references.csv", show_col_types = FALSE) %>% 
  filter(!is.na(`Mosquito species`)) %>% 
  filter(!is.na(`Citation key`)) %>% 
  select(`Mosquito species`, `Citation key`) %>% 
  separate(`Mosquito species`, into = "Genera", sep = " ") %>% 
  group_by(Genera) %>% 
  # Make citation key "nice" for input into LaTeX
  summarise(refs = toString(unique(`Citation key`))) %>% 
  mutate(LaTeX_cite = paste0("\\citep{", refs, "}")) %>% 
  select(-refs)

write_csv(genera_refs, "data/clean/genera_references.csv")
  
