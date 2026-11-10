# reference: https://www.karissawhiting.com/oncokbR/
# Date      Note
# 26JUN2026 Use cbioportal data as default and add parameter synapse vs. cbioportal
##################################################################

# remotes::install_github('karissawhiting/oncokbR')
library(oncokbR)


####################################################### #
# Load mutations, SV, CNA data --------------------------
####################################################### #

if (source_cbioportal) {

  cbioportalR::set_cbioportal_db("https://genie.cbioportal.org/api")

  mutations_extended_panc <- cbioportalR::get_mutations_by_study("panc_genie_bpc")
  sv_panc                 <- cbioportalR::get_structural_variants_by_study("panc_genie_bpc")
  cna_panc                <- cbioportalR::get_cna_by_study("panc_genie_bpc")

}

SAMPLE_ID_COL <- if (use_synapse) "sample_id" else "sampleId"

####################################################### #
# MUTATIONS ---------------------------------------------
####################################################### #

if (!use_synapse) {
  mutations_extended_tumor_type <- mutations_extended_panc %>%
    ## merge on oncotree code from clinical data to get tumor_type variable in example
    left_join(
      .,
      panc_v1_0$PANC_v1.0$cpt %>%
        select(cpt_genie_sample_id, cpt_oncotree_code) %>%
        rename(tumor_type = cpt_oncotree_code),
      by = setNames("cpt_genie_sample_id", SAMPLE_ID_COL)
    ) 
} else {
  ## code below is only needed for data from Synapse (not needed for cbioportal data)
  ## filter out variant classifications that are not included in the consequence_map
  ## -- silent, rna, 5'utr, 3'utr
  mutations_extended_tumor_type <- mutations_extended_tumor_type %>%
    filter(
      str_to_lower(Variant_Classification) %in% c(consequence_map$variant_classification),
      !is.na(HGVSp_Short),
      !is.na(Protein_position)
    )
}

## annotate with oncokbR
## - default is annotate by protein change
oncokb_mutations_extended_panc <- annotate_mutations(
  mutations = mutations_extended_tumor_type
  )


####################################################### #
# FUSIONS/SV --------------------------------------------
####################################################### #

## call variant_class unknown if need to use Synapse genomic data
sv_tumor_type <- sv_panc %>%
  left_join(
    .,
    panc_v1_0$PANC_v1.0$cpt %>%
      select(cpt_genie_sample_id, cpt_oncotree_code) %>%
      rename(tumor_type = cpt_oncotree_code),
    by = setNames("cpt_genie_sample_id", SAMPLE_ID_COL)
  )

## annotate with oncokbR
oncokb_sv_panc <- annotate_sv(sv = sv_tumor_type)


####################################################### #
# CNA ---------------------------------------------------
####################################################### #

## merge on oncotree code
cna_tumor_type <- cna_panc %>%
  left_join(
    .,
    panc_v1_0$PANC_v1.0$cpt %>%
      select(cpt_genie_sample_id, cpt_oncotree_code) %>%
      rename(tumor_type = cpt_oncotree_code),
    by = setNames("cpt_genie_sample_id", SAMPLE_ID_COL)
  )

# annotate with oncokbR
oncokb_cna_panc <- annotate_cna(cna_tumor_type)


# Save results ------------------------------------------

save(
  oncokb_mutations_extended_panc,
  oncokb_sv_panc,
  oncokb_cna_panc,
  file = here::here(out.dir, "panc_v1_0_oncokb_annotation.rdata")
)
