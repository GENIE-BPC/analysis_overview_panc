## adapted from Yufei's code: sample-oncoprint-with-annotations 1.R

library(ComplexHeatmap)
library(tidyverse)

out.dir <- here::here("data", "v1_0_public")
fig.dir <- here::here("documents", "figures")

## load oncokb annotations and panel info
load(here::here(out.dir, "panc_v1_0_oncokb_annotation.rdata"))
load(here::here(out.dir, "panc_sample_panels_and_gen_dat_summary.rdata"))


## 6/2026: column (sample) order for the oncoprint
# NOTE: this need to be manually generated via the GENIE cbioportal:
# 1) visit the GENIE cbio webpage for the BPC cohort (https://genie.cbioportal.org/ or https://genie-private.cbioportal.org/)
# 2) query by samples (GENIE-XX-XXX-XXX) and genes (see below top 25 genes), this directs to the summary page
# 3) on the "OncoPrint" tab, click "View" and choose "Events per sample", and click "Download" -> "Sample order".
sample_order_pre <-
  read.table(here::here(out.dir, "panc_v1_0_sample_order.txt"))[,1]
sample_order <- sample_order_pre[sample_order_pre %in% sample_panels$sample_id]

## top 25 most frequently altered genes
top_25_genes <- c(
  "KRAS",
  "TP53",
  "CDKN2A",
  "SMAD4",
  "CDKN2B",
  "ARID1A",
  "ATM",
  "KDM6A",
  "RNF43",
  "BRCA2",
  "AKT2",
  "MYC",
  "PIK3CA",
  "MTAP",
  "KMT2D",
  "GNAS",
  "RB1",
  "SMARCA4",
  "MAP2K4",
  "PBRM1",
  "ARID2",
  "BRCA1",
  "ERBB2",
  "STK11",
  "CCNE1",
  "KMT2C",
  "U2AF1"
)

# ComplexHeatmap version 2.16.0
# If you use it in published research, please cite either one:
# - Gu, Z. Complex Heatmap Visualization. iMeta 2022.
# - Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional
# genomic data. Bioinformatics 2016.



# 1) Prepare Oncoprint Matrix  ---------------------------------------------

## SNV - Select gene, mutation type --------------------
## NOTE: mut is a MAF formatted data frame

mut_df <- oncokb_mutations_extended_panc %>%
  filter(oncokb_oncogenic %in% c("Oncogenic", "Likely Oncogenic"),
         sample_id %in% c(sample_panels %>% select(sample_id) %>% unlist()),
         hugo_symbol %in% top_25_genes) %>%
  select(ID = sample_id,
         Gene = hugo_symbol,
         SNV_TYPE=variant_classification,
         chr = chromosome,
         pos = start_position,
         ref = reference_allele,
         ## 7/9/2025 (SB): is `allele` the correct column? Yufei's column was called `variantAllele`
         alt = allele
         ) %>%
  filter(chr != 'MT') %>%
  mutate(
    SNV_TYPE = case_when(
      SNV_TYPE == 'Missense_Mutation' ~ 'MISSENSE',
      SNV_TYPE == 'Nonsense_Mutation' ~ 'TRUNC',
      grepl("Frame_Shift", SNV_TYPE) ~ 'TRUNC',
      grepl("In_Frame", SNV_TYPE) ~ 'INFRAME',
      grepl("splice", tolower(SNV_TYPE)) ~ 'SPLICE',
      TRUE ~ 'OTHER'
    )
  ) %>%
  select(ID, Gene, SNV_TYPE) %>%
  distinct()

mut_df_wide <- mut_df %>%
  group_by(ID, Gene) %>%
  mutate(snv_index = paste0("SNV_TYPE", row_number())) %>%
  pivot_wider(names_from = snv_index, values_from = SNV_TYPE)


## CNA - Select genes, alteration types -----------------

cna_df <- oncokb_cna_panc %>%
  filter(oncokb_oncogenic %in% c("Oncogenic", "Likely Oncogenic"),
         sample_id %in% c(sample_panels %>% select(sample_id) %>% unlist()),
         hugo_symbol %in% top_25_genes) %>%
  select(ID = sample_id,
         Gene = hugo_symbol,
         CNA_TYPE = alteration) %>%
  mutate(CNA_TYPE = case_when(
    CNA_TYPE == 'AMPLIFICATION'  ~ 'AMP',
    CNA_TYPE == 'DELETION' ~ 'HOMDEL',
    TRUE ~ CNA_TYPE
  ))


## SV - (code by SB) -----------------------------

sv_df <- oncokb_sv_panc %>%
  filter(oncokb_oncogenic %in% c("Oncogenic", "Likely Oncogenic"),
         sample_id %in% c(sample_panels %>% select(sample_id) %>% unlist()),
          (grepl(paste0(top_25_genes, collapse = "|"), site_1_hugo_symbol)) |
          (grepl(paste0(top_25_genes, collapse = "|"), site_2_hugo_symbol))
         ) %>%
  select(ID = sample_id,
         Gene1 = site_1_hugo_symbol,
         Gene2 = site_2_hugo_symbol) %>%
  pivot_longer(cols = c(starts_with("Gene")),
               names_to = "Gene_number",
               values_to = "Gene") %>%
  select(-Gene_number) %>%
  distinct() %>%
  mutate(SV_TYPE = "SV")


## Fuse all alteration types------------------------------

full_mat <- plyr::join_all(list(mut_df_wide, cna_df, sv_df),
                           type = "full",
                           by = c('ID','Gene'))
full_mat[is.na(full_mat)] <- ''

full_mat$Alteration <- full_mat$SNV_TYPE1
for (i in 1:nrow(full_mat)) {
  # check SNV alteration type
  snv1 <- full_mat[i, 'SNV_TYPE1']
  snv2 <- full_mat[i, 'SNV_TYPE2']
  s <- ifelse(nchar(snv1)==0 || nchar(snv2)==0, '', ';')

  comb <- paste(snv1, snv2, sep = s)

  # check CNV amp / homdel
  amp_del <- full_mat[i, 'CNA_TYPE']
  s <- ifelse(nchar(comb)==0 || nchar(amp_del)==0, '', ';')

  comb <- paste(comb, amp_del, sep = s)

  ## SB: add structural variants
  sv <- full_mat[i, 'SV_TYPE']
  s <- ifelse(nchar(comb)==0 || nchar(sv)==0, '', ';')

  comb <- paste(comb, sv, sep = s)

  full_mat[i, "Alteration"] <- comb
}
### Example "full_mat" rows:
# ID                  Gene    SNV_TYPE1 SNV_TYPE2 CNA_TYPE Alteration
# <chr>               <chr>   <chr>     <chr>     <chr>    <chr>
# 1 P-0033093-T01-IM6 TP53    MISSENSE  ""        ""       MISSENSE
# 2 P-0033093-T01-IM6 EGFR    INFRAME   ""        "AMP"    INFRAME;AMP
# 3 P-0033093-T01-IM6 SMAD4   TRUNC     ""        ""       TRUNC
# 4 P-0033093-T01-IM6 PIK3C3  MISSENSE  ""        ""       MISSENSE



## Create final data frame in oncoprint format -------------

genes <- top_25_genes #unique(full_mat$Gene)
samples <- sample_panels %>%
  select(sample_id) %>%
  unlist() #unique(full_mat$ID)

full_mat_final <- data.frame(matrix('', nrow=length(genes), ncol=length(samples)))
rownames(full_mat_final) <- genes
colnames(full_mat_final) <- samples
# dim(full_mat_final)

## Encode NAs for samples without gene coverage so the % can be calculated with the correct denominator
gen_dat_summary_top_genes <- gen_dat_summary %>%
  rename_with(~ gsub("_oncokb", "", .x)) %>%
  select(sample_id, top_25_genes) %>%
  column_to_rownames('sample_id') %>%
  t() %>% data.frame(check.names = F)


for (g in rownames(full_mat_final)) {
  for (s in colnames(full_mat_final)) {

    if (is.na(gen_dat_summary_top_genes[g, s])) {
      alt <- 'NA'

    } else {
      alt <- full_mat[full_mat$ID == s &
                        full_mat$Gene == g, 'Alteration']
      alt <- ifelse(length(alt)==0, '', alt)
    }
    full_mat_final[g, s] <- alt
  }
}

# 10-27-2026 check percent tumor without any mutation in selected genes
is_all_na <- apply(full_mat_final, MARGIN = 2, function(x) all(x==''|x=='NA'))
round(sum(is_all_na) / length(is_all_na) * 100, 3)
write.csv(full_mat_final[, is_all_na] %>% colnames(), file = here::here("data", "v1_0_public", "ids_no_alteration.txt"),
          quote = F, row.names = F)


# 2) Oncoprint settings -----------------------------------------------------------

## Color mappings -------------

col = c(
  "AMP" = "red2",
  "HOMDEL" = 'blue3',
  # "GAIN" = "salmon",
  # "LOSS" = 'steelblue',
  "SV" = "purple4",
  "MISSENSE" = "#008000",
  "TRUNC" = "black",
  "INFRAME" = "red4",
  "SPLICE" = "orange",
  "OTHER" = "gold",
  "NA" = "white"
)



## Aesthetics ------------------
# NOTE: can use alter_graphic() instead of function()

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w*0.98, h*0.98, gp = gpar(fill = "lightgray", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["AMP"], col = NA))
  },
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["HOMDEL"], col = NA))
  },
  # GAIN = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["GAIN"], col = NA))
  # },
  # LOSS = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["LOSS"], col = NA))
  # },
  SV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, gp = gpar(fill = col["SV"], col = NA))
  },
  MISSENSE = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["MISSENSE"], col = NA))
  },
  TRUNC = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["TRUNC"], col = NA))
  },
  INFRAME = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["INFRAME"], col = NA))
  },
  SPLICE = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["SPLICE"], col = NA))
  },
  OTHER = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["OTHER"], col = NA))
  },
  # `NA` = function(x, y, w, h) {
  #   grid.rect(x, y, w, h, gp = gpar(fill = col["NA"], col = 'white'))
  # }
  `NA` = function(x, y, w, h) {
    grid.rect(x, y, w, h, gp = gpar(fill = col["NA"], col = NA))

    # draw thin horizontal line
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.08, gp = gpar(fill = 'lightgrey', col = NA))
  }
)

## Column annotation ---------------------------

profiled_cna_yes_no <- ifelse(
  grepl('UHN', sample_order), 'No', 'Yes')

profiled_sv_yes_no <- ifelse(
  (!grepl('UHN', sample_order)) |
   grepl('GENIE-UHN-OCT271391-ARC1|GENIE-UHN-OCT296767-ARC1|GENIE-UHN-OCT420107-ARC1', sample_order),
  'Yes', 'No')

names(profiled_cna_yes_no) <- names(profiled_sv_yes_no) <- sample_order


ha <- HeatmapAnnotation(
  `Profiled for CNA`= anno_simple(profiled_cna_yes_no,
                                  col = c('No' = "white", 'Yes' = "black"),
                                  border = T),
  `Profiled for SV` = anno_simple(profiled_sv_yes_no,
                                  col = c('No' = "white", 'Yes' = "black"),
                                  border = T),
  height = unit(.8, 'cm'),
  # annotation_height = unit(c(8,8), "mm"),
  gap = unit(1, 'mm'),
  show_annotation_name = T,
  annotation_name_side = 'left',
  annotation_name_gp   = gpar(fontsize = 8),
  show_legend = FALSE
)


## Row annotations -------------------------------------------

## calc correct percentages excluding samples w/out gene coverage
any_na <- apply(full_mat_final, 1, function(x) any(x=='NA'))
profiled_n <- rowSums(full_mat_final != 'NA')
altered_n  <- rowSums(full_mat_final != "" & full_mat_final != 'NA')
pct <- sort(round(100 * altered_n / profiled_n), decreasing = T)
gene_order <- names(pct)

ra <- rowAnnotation(
  Frequency = anno_text(
    sapply(gene_order, function(g) ifelse(any_na[[g]], sprintf("%.0f%%*", pct[[g]]), sprintf("%.0f%%", pct[[g]]))),
    just = "left",
    gp = gpar(fontsize = 8))
  )



## Legends -----------------

lgd1 <- Legend(title = "Genomic alterations",
               labels = c("Amplification",
                          "Deep Deletion",
                          # "Gain of function",
                          # "Loss of function",
                          'Structural Variant',
                          'Missense',
                          'Truncating',
                          'Inframe InDel',
                          "Splice",
                          'Other Mutation'),
               graphics = alter_fun[c("AMP", "HOMDEL", #"GAIN", "LOSS",
                                      "SV", "MISSENSE",
                                      "TRUNC", "INFRAME","SPLICE", "OTHER")],
               labels_gp = gpar(fontsize=8),
               title_gp = gpar(fontsize=10), #title_position = "topcenter",
               nrow = 3)

# lgd2 <- Legend(labels = "Gene not profiled", type = "points", title = NULL,
#                background = "white",  pch = "*",
#                legend_gp = gpar(col = 'grey30'),
#                labels_gp = gpar(fontsize=8))
lgd2 <- Legend(labels = "Gene not profiled", graphics = alter_fun["NA"],
               title = NULL,
               labels_gp = gpar(fontsize=8))

lgd_cna <- Legend(
  title = "Profiled for Copy-Number Alterations (CNA)",
  at = c("No","Yes"),
  legend_gp = gpar(fill = c("white","black"), col='black', lwd=0.1),
  labels_gp = gpar(fontsize=8),
  title_gp = gpar(fontsize=8), #title_position  = "topcenter",
  ncol = 2, border = 'black'
)

lgd_sv <- Legend(
  title = "Profiled for Structural Variants (SV)",
  at = c("No","Yes"),
  legend_gp = gpar(fill = c("white","black")),
  labels_gp = gpar(fontsize=8),
  title_gp = gpar(fontsize=8),# title_position  = "topcenter",
  ncol = 2, border = 'black'
)

lgd_col <- packLegend(lgd_cna, lgd_sv, direction = "vertical", row_gap = unit(2.5, "mm"))

heatmap_legend <- packLegend(lgd1, lgd_col, lgd2,
                             max_width = unit(18, "cm"),
                             direction = "horizontal",
                             column_gap = unit(5, "mm")
  )



# 3) Draw oncoprint ----------------------------------------------------

# * Main oncoprint -----------------------------------

p <- oncoPrint(
  mat = full_mat_final[gene_order, sample_order],
  alter_fun = alter_fun,
  alter_fun_is_vectorized = T,
  col = col,

  top_annotation = ha,
  right_annotation = ra,
  left_annotation = NULL,

  row_split = factor(gene_order, levels = gene_order, labels = gene_order),
  row_title = NULL,
  row_order = gene_order,
  show_row_names = T,
  row_names_gp = gpar(fontsize = 8),
  row_names_side = "left",

  column_order = sample_order,  # E.g., order by TMB (low to high)
  show_column_names = F, # NOTE: change this to TRUE to show sample IDs at bottom

  show_heatmap_legend = F, # Hide default alteration legend and use customized ones
  show_pct = F,            # Hide default (wrong numerator due to NAs)
  width = unit(18,'cm'),
  height = unit(10, 'cm')
)


## [NEW 6/26/2026] ---------------------

# * Group by stage ---------------------------------------

## create named vector for stage (must align with sample_order)
## note 6/26/2026: stage annotation applied only to samples acquired within
## 60 days time window from dx (derived var `tt_dx_proc_days` <= 60 days)
stage_annotation <- (sample_panels %>%
  mutate(stage_filtered = case_when(
    tt_dx_proc_days <= 60 ~ stage,
    TRUE ~ "Sample acquired > 60 days from dx"
  )) %>% pull(stage_filtered))[match(sample_order, sample_panels$sample_id)]

names(stage_annotation) <- sample_order


## row annotation object
ha_w_stage <- HeatmapAnnotation(
  `Stage at Diagnosis (dx)` = anno_simple(as.character(stage_annotation),
                                     col = c("Stage I-III"="#4DAF4A", "Stage IV"="#E41A1C",
                                             "Sample acquired > 60 days from dx"="grey"),
                                     border = T),
  `Profiled for CNA` = anno_simple(profiled_cna_yes_no,
                                   col = c('No' = "white", 'Yes' = "black"),
                                   border = T),
  `Profiled for SV`  = anno_simple(profiled_sv_yes_no,
                                   col = c('No' = "white", 'Yes' = "black"),
                                   border = T),
  height = unit(1.2, 'cm'),
  gap = unit(1, 'mm'),
  show_annotation_name  = T,
  annotation_name_side  = 'left',
  annotation_name_gp    = gpar(fontsize = 8),
  show_legend = FALSE
)

# ## legend
# lgd_stage <- Legend(
#   title  = "Stage at Diagnosis",
#   at     = c("Stage I-III", "Stage IV"),
#   legend_gp  = gpar(fill = stage_col),
#   labels_gp  = gpar(fontsize = 8),
#   title_gp   = gpar(fontsize = 8),
#   ncol   = 2,
#   border = 'black'
# )
# lgd_col1 <- packLegend(lgd_stage, lgd_cna, lgd_sv, direction = "vertical", row_gap = unit(2, "mm"))
# heatmap_legend1 <- packLegend(lgd1, lgd_col, lgd2,
#                              max_width = unit(18, "cm"),
#                              direction = "horizontal",
#                              column_gap = unit(5, "mm"))



## oncoprint: by stage

# column splitting factor
column_split_stage <- factor(
  stage_annotation[sample_order],
  levels = c("Stage I-III", "Stage IV","Sample acquired > 60 days from dx")
)

p.by.stage <- oncoPrint(
  mat = full_mat_final[gene_order, sample_order],
  alter_fun = alter_fun,
  alter_fun_is_vectorized = T,
  col = col,

  top_annotation = ha_w_stage,
  right_annotation = ra,
  left_annotation = NULL,

  row_split  = factor(gene_order, levels = gene_order, labels = gene_order),
  row_title  = NULL,
  row_order  = gene_order,
  show_row_names  = T,
  row_names_gp    = gpar(fontsize = 8),
  row_names_side  = "left",

  column_order = sample_order,
  column_split = column_split_stage,
  column_gap   = unit(3, "mm"),
  column_title_gp = gpar(fontsize = 9),

  show_column_names    = F,
  show_heatmap_legend  = F,
  show_pct             = F,
  width  = unit(18, 'cm'),
  height = unit(10, 'cm')
)



# * Group by prior chemo ---------------------------------------

## create named vector for group label (must align with sample_order)
chemo_annotation <- (sample_panels %>% mutate(chemo = case_when(
    any_chemo == 1 ~ "Yes",
    any_chemo == 0 ~ "No",
    TRUE ~ "Unknown")) %>% pull(chemo))[match(sample_order, sample_panels$sample_id)]
names(chemo_annotation) <- sample_order


## row annotation object
ha_w_chemo <- HeatmapAnnotation(
  `Prior Chemotherapy` = anno_simple(as.character(chemo_annotation),
                                     col = c("Yes"='#8fc33e', 'No'="#0181E8", 'Unknown'='grey'),
                                     border = T),
  `Profiled for CNA` = anno_simple(profiled_cna_yes_no,
                                   col = c('No' = "white", 'Yes' = "black"),
                                   border = T),
  `Profiled for SV`  = anno_simple(profiled_sv_yes_no,
                                   col = c('No' = "white", 'Yes' = "black"),
                                   border = T),
  height = unit(1.2, 'cm'),
  gap = unit(1, 'mm'),
  show_annotation_name  = T,
  annotation_name_side  = 'left',
  annotation_name_gp    = gpar(fontsize = 8),
  show_legend = FALSE
)


## Oncooprint: by prior chemo status

column_split_chemo <- factor(chemo_annotation[sample_order],
                             levels = c("Yes", "No", "Unknown"))

p.by.chemo <- oncoPrint(
  mat = full_mat_final[gene_order, sample_order],
  alter_fun = alter_fun,
  alter_fun_is_vectorized = T,
  col = col,

  top_annotation = ha_w_chemo,
  right_annotation = ra,
  left_annotation = NULL,

  row_split  = factor(gene_order, levels = gene_order, labels = gene_order),
  row_title  = NULL,
  row_order  = gene_order,
  show_row_names  = T,
  row_names_gp    = gpar(fontsize = 8),
  row_names_side  = "left",

  column_order = sample_order,
  column_split = column_split_chemo,
  column_gap   = unit(3, "mm"),
  column_title_gp = gpar(fontsize = 9),

  show_column_names    = F,
  show_heatmap_legend  = F,
  show_pct             = F,
  width  = unit(18, 'cm'),
  height = unit(10, 'cm')
)



# 4) Save output  -------------------

# * Save to PDF ----

pdf(here::here(fig.dir, paste0("Oncoprint_", Sys.Date(), "-test.pdf")), width = 12, height = 8)
draw(p,
     heatmap_legend_side = "bottom",
     heatmap_legend_list = heatmap_legend,
     merge_legend        = T,
     padding = unit(c(6, 6, 14, 6), "mm")
  )
dev.off()

pdf(here::here(fig.dir, paste0("Oncoprint_by_stage_", Sys.Date(), ".pdf")), width = 12, height = 8)
draw(p.by.stage,
     heatmap_legend_side = "bottom",
     heatmap_legend_list = heatmap_legend,
     merge_legend        = T,
     padding = unit(c(6, 6, 14, 6), "mm")
  )
dev.off()

pdf(here::here(fig.dir, paste0("Oncoprint_by_chemo_", Sys.Date(), ".pdf")), width = 12, height = 8)
draw(p.by.chemo,
     heatmap_legend_side = "bottom",
     heatmap_legend_list = heatmap_legend,
     merge_legend        = T,
     padding = unit(c(6, 6, 14, 6), "mm")
)
dev.off()


# * save to JPG -----

jpeg(file = here::here(fig.dir, paste0("Oncoprint_", Sys.Date(), ".jpg")), width = 4000, height = 3000, res = 450)
draw(p,
     heatmap_legend_side = "bottom",
     heatmap_legend_list = heatmap_legend,
     merge_legend        = T,
     padding = unit(c(6, 6, 14, 6), "mm")
)
dev.off()

jpeg(here::here(fig.dir, paste0("Oncoprint_by_stage_", Sys.Date(), ".jpg")), width = 4000, height = 3000, res = 450)
draw(p.by.stage,
     heatmap_legend_side = "bottom",
     heatmap_legend_list = heatmap_legend,
     merge_legend        = T,
     padding = unit(c(6, 6, 14, 6), "mm")
)
dev.off()

jpeg(here::here(fig.dir, paste0("Oncoprint_by_chemo_", Sys.Date(), ".jpg")), width = 4000, height = 3000, res = 450)
draw(p.by.chemo,
     heatmap_legend_side = "bottom",
     heatmap_legend_list = heatmap_legend,
     merge_legend        = T,
     padding = unit(c(6, 6, 14, 6), "mm")
)
dev.off()


# ## 4) Export to ppt
#
# library(officer)
# library(rvg)
# library(grid)
#
# save_plots <- list(p)
#
# plot_output <- read_pptx()
#
# plot_output <- add_slide(plot_output, layout = "Title and Content", master = "Office Theme")
#
# plot_output <- ph_with(
#   plot_output,
#   value = dml(code = {
#     draw(
#       save_plots[[1]],
#       heatmap_legend_side = "bottom",
#       heatmap_legend_list = list(heatmap_legend),  # IMPORTANT: list()
#       merge_legends = TRUE,                        # correct arg name
#       padding = unit(c(6, 6, 14, 6), "mm")
#     )
#   }),
#   location = ph_location_fullsize()
# )
#
#
# print(
#   plot_output,
#   target = here::here("documents/figures",
#                       paste0("panc_oncoprint_", format(Sys.Date(), "%Y-%m-%d"), ".pptx"))
#   )


# library(officer)
# library(rvg)
#
# oncoprint <- draw(p,
#                   heatmap_legend_side = "bottom",
#                   heatmap_legend_list = heatmap_legend,
#                   merge_legend        = T,
#                   padding = unit(c(6, 6, 14, 6), "mm")
# )
#
# save_plots = c(list(oncoprint))
#
# # Create the ppt
# plot_output = read_pptx()
#
# # Create a loop to run through all plots
# for(i in 1:length(save_plots)){
#   # Create a new page
#   plot_output = add_slide(plot_output, 'Title and Content', 'Office Theme')
#
#   # Plot
#   # save_plots[[i]]
#
#   plot_output = ph_with(plot_output,
#                         # dml(ggobj = last_plot()),
#                         value = dml(code = print(save_plots[[i]], newpage = FALSE)),
#                         location = ph_location_fullsize()
#   )
#
# }
#
# print(plot_output, target = here::here(paste0("panc_oncoprint_", Sys.time() %>% format("%Y-%m-%d"), ".pptx")))


