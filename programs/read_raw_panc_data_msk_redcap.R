# --------------------------------------------------------------------------
# Read in raw data via API:
#
# The purpose of this document is to read in PANC data from the MSK REDCap via API.
# --------------------------------------------------------------------------

library(tidyverse)
library(redcapAPI)

# # Note: need to define REDCap API token for GENIE BPC PANC cohort below in order to run script
# genie_bpc_panc <- "xxx"



# -------------------------------------- Read in data 
read_data_panc_list = REDCapR::redcap_read(
  redcap_uri = "https://redcap.mskcc.org/api/",
  token = genie_bpc_panc,
  raw_or_label = "label",
  export_checkbox_label = TRUE
)

read_data_panc = as.data.frame(read_data_panc_list$data)
# nrow(read_data_panc)
# ncol(read_data_panc)
# table(read_data_panc$redcap_repeat_instrument)


# -------------------------------------- Set up data

raw_data_panc <- read_data_panc %>% 
  mutate(cohort = "PANC",
         institution = "MSK",
         redcap_repeat_instrument = case_when(
           redcap_repeat_instrument == "Ca Directed Drugs" ~ "ca_drugs",
           redcap_repeat_instrument == "Cancer Diagnosis" ~ "ca_dx",
           redcap_repeat_instrument == "Cancer Panel Test" ~ "cpt",
           redcap_repeat_instrument == "Prissmm Imaging" ~ "prissmm_i",
           redcap_repeat_instrument == "Prissmm Med Onc Assessment" ~ "prissmm_md",
           redcap_repeat_instrument == "Prissmm Pathology" ~ "prissmm_path",
           redcap_repeat_instrument == "Ca Directed Radtx" ~ "ca_rt",
           redcap_repeat_instrument == "Prissmm Tumor Marker" ~ "prissmm_tm",
           TRUE ~ "pt_char"
         ),
         dup_rec = case_when(
           word(record_id, 2, sep = "_") == 2 ~ 1,
           TRUE ~ 0
         )
  ) %>% 
  # remove duplicate records and cases not meeting incl criteria
  filter(dup_rec != 1)

# table(raw_data_panc$redcap_repeat_instrument, useNA = "always")



# Split rows into separate datasets by structure
# Dataset of everything that's 1 rec/pt and then separate datasets for other structures with repeating records per patient
raw_data_split_panc <- raw_data_panc %>%
  split(.$redcap_repeat_instrument)


# names of each instrument / names of dfs in the list
names(raw_data_split_panc)



# -------------------------------------- Subset on variables in the data

# need to get list of variables in the data dictionary by form
# and have to remove variables that aren't in the data

# want variables record_id, redcap_repeat_instance, and group (deliverable/non-deliverable) on each dataset
rec_id_panc <- tibble(field_name = c(rep("record_id", length(unique(raw_data_panc$redcap_repeat_instrument))),
                                rep("cohort", length(unique(raw_data_panc$redcap_repeat_instrument))),
                                rep("institution", length(unique(raw_data_panc$redcap_repeat_instrument))),
                                rep("redcap_repeat_instance", length(unique(raw_data_panc$redcap_repeat_instrument))),
                                rep("group",
                                    length(unique(raw_data_panc$redcap_repeat_instrument)))),
                 form_name = rep(unique(raw_data_panc$redcap_repeat_instrument), 5))

complete_vars_panc = raw_data_panc %>% 
  pivot_longer(cols = contains("complete"), names_to = "field_name", values_drop_na = TRUE) %>% 
  select(redcap_repeat_instrument, field_name) %>% 
  group_by(redcap_repeat_instrument, field_name) %>% 
  unique() %>% 
  rename(form_name = redcap_repeat_instrument) %>% 
  tibble()

# get meta-data (data dictionary) from REDCap
dictionary_panc <- bind_rows(
  REDCapR::redcap_metadata_read(redcap_uri = "https://redcap.mskcc.org/api/",
                     token = genie_bpc_panc)$data,
  rec_id_panc, 
  complete_vars_panc
) %>%
  mutate(instrument = case_when(
    form_name == "ca_directed_drugs" ~ "ca_drugs",
    form_name == "cancer_diagnosis" ~ "ca_dx",
    form_name == "cancer_panel_test" ~ "cpt",
    form_name == "prissmm_imaging" ~ "prissmm_i",
    form_name == "prissmm_med_onc_assessment" ~ "prissmm_md",
    form_name == "prissmm_pathology" ~ "prissmm_path",
    form_name %in% c("prissmm_tumor_marker", "prissmm_tm") ~ "prissmm_tm",
    form_name == "ca_directed_radtx" ~ "ca_rt",
    form_name == "patient_characteristics" ~ "pt_char",
    form_name == "curation_completion" ~ "pt_char", 
    form_name == "curation_initiation_eligibility" ~ "pt_char",
    form_name == "quality_assurance" ~ "pt_char",
    TRUE ~ form_name)) %>% 
  # subset on available forms until more data is added
  filter(instrument %in% c("pt_char", "ca_dx", "ca_drugs", "ca_rt",
                           "prissmm_path", "prissmm_i", "prissmm_md", "prissmm_tm", "cpt"))

# get names of variables in the data
# have to account for checkbox variables
# e.g. symptoms in data dictionary vs symptoms___1 - symptoms___16 in redcap data
redcap_var_names_panc <- tibble(
  var_name = names(raw_data_panc),
  field_name = word(names(raw_data_panc), 1, sep = "___")
)

# remove variables that are in data dictionary but are not in data
dictionary_subset_panc <- inner_join(redcap_var_names_panc,
                                dictionary_panc,
                                by = c("field_name")
) 



# -------------------------------------- Split into separate datasets

# keep only variables of interest on the dataset
# get data dictionary (restricted to available variables) as a list for each dataset
dictionary_split_panc <- dictionary_subset_panc %>%
  # filter(instrument != "active_treatment") %>% 
  split(.$instrument) %>%
  map(., select, var_name) %>%
  map(., as_vector) %>%
  # need to unname the vector so that the vector consists only of the variable names
  # otherwise each variable is field_name1, field_name2, etc.
  map(., unname)

# check order before using mapply
names(raw_data_split_panc) == names(dictionary_split_panc)

# map over each dataset to keep only the variables belonging to that dataset
new_split_panc <- mapply(select, raw_data_split_panc, dictionary_split_panc)

for (i in 1:length(names(new_split_panc))) {
  names(new_split_panc)[i] = paste0(names(new_split_panc)[i], "_panc")
}


# -------------------------------------- Save data

list2env(new_split_panc, envir = .GlobalEnv)

save(
  pt_char_panc,
  ca_dx_panc,
  ca_drugs_panc, 
  prissmm_i_panc, 
  prissmm_path_panc, 
  prissmm_md_panc,
  prissmm_tm_panc, 
  ca_rt_panc,
  cpt_panc,
  file = here::here("data", "msk_raw_redcap_data.rdata")
)


