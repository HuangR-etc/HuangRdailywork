library(dplyr)
library(stringr)
library(ape)

## -----------------------------
## Cricetidae sensitivity subsets:
## species membership table for Supplementary Table S2
## -----------------------------

cricetidae_pool_dir <- "/home/galileo-group/galileouser05/huangr/selected_subsets"
cricetidae_tree_file <- "/home/galileo-group/galileouser05/huangr/selected_subsets/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre"
cricetidae_taxonomy_file <- "/home/galileo-group/galileouser05/huangr/selected_subsets/taxonomy_mamPhy_5911species_toPublish.csv"
cricetidae_output <- "/home/galileo-group/galileouser05/huangr/selected_subsets/Supplementary_Table_S2_Cricetidae_membership.csv"
cricetidae_formatted_output <- "/home/galileo-group/galileouser05/huangr/selected_subsets/Supplementary_Table_S2_Cricetidae_membership_formatted.csv"

cricetidae_subset_order <- c("C512", "C256", "C128", "C64", "C32")
cricetidae_pool_files <- file.path(
  cricetidae_pool_dir,
  paste0("Cricetidae_", cricetidae_subset_order, "_candidate_pool.csv")
)

names(cricetidae_pool_files) <- cricetidae_subset_order

missing_pool_files <- cricetidae_pool_files[!file.exists(cricetidae_pool_files)]
if (length(missing_pool_files) > 0) {
  stop(
    "Missing Cricetidae candidate pool files: ",
    paste(missing_pool_files, collapse = ", ")
  )
}

cricetidae_pools <- lapply(cricetidae_pool_files, function(pool_file) {
  pool_tbl <- read.csv(pool_file, stringsAsFactors = FALSE)

  if (!"Species" %in% names(pool_tbl)) {
    stop("Candidate pool file has no Species column: ", pool_file)
  }

  unique(pool_tbl$Species[!is.na(pool_tbl$Species) & pool_tbl$Species != ""])
})

cricetidae_tree <- read.nexus(cricetidae_tree_file)
cricetidae_tree <- drop.tip(
  cricetidae_tree,
  intersect(
    c(
      "_Anolis_carolinensis",
      "Nyctinomops_kalinowskii_MOLOSSIDAE_CHIROPTERA"
    ),
    cricetidae_tree$tip.label
  )
)

cricetidae_tax <- read.csv(cricetidae_taxonomy_file, stringsAsFactors = FALSE)
cricetidae_tax <- cricetidae_tax[
  toupper(cricetidae_tax[["fam"]]) == "CRICETIDAE" &
    cricetidae_tax[["extinct."]] == 0,
]

if ("samp" %in% names(cricetidae_tax)) {
  cricetidae_tax <- cricetidae_tax[cricetidae_tax$samp == "sampled", ]
}

cricetidae_all_tips <- sort(intersect(
  cricetidae_tax$tiplabel,
  cricetidae_tree$tip.label
))

not_in_full_cricetidae <- setdiff(
  unique(unlist(cricetidae_pools)),
  cricetidae_all_tips
)

if (length(not_in_full_cricetidae) > 0) {
  stop(
    "Candidate pool species not found in the full Cricetidae tree/taxonomy pool. Example species: ",
    paste(head(not_in_full_cricetidae, 10), collapse = ", ")
  )
}

for (i in seq_len(length(cricetidae_subset_order) - 1)) {
  parent_subset <- cricetidae_subset_order[i]
  child_subset <- cricetidae_subset_order[i + 1]
  not_nested <- setdiff(
    cricetidae_pools[[child_subset]],
    cricetidae_pools[[parent_subset]]
  )

  if (length(not_nested) > 0) {
    stop(
      child_subset,
      " is not nested within ",
      parent_subset,
      ". Example species: ",
      paste(head(not_nested, 10), collapse = ", ")
    )
  }
}

cricetidae_members_long <- bind_rows(
  lapply(names(cricetidae_pools), function(subset_name) {
    data.frame(
      Subset = subset_name,
      tip_label = cricetidae_pools[[subset_name]],
      stringsAsFactors = FALSE
    )
  })
) %>%
  mutate(
    Latin.name = tip_label %>%
      str_remove("_[A-Z]+_[A-Z]+$") %>%
      str_replace_all("_", " ")
  ) %>%
  distinct(Latin.name, Subset)

cricetidae_species <- sort(str_replace_all(
  str_remove(cricetidae_all_tips, "_[A-Z]+_[A-Z]+$"),
  "_",
  " "
))

cricetidae_membership <- data.frame(
  Latin.name = cricetidae_species,
  stringsAsFactors = FALSE
)

for (subset_name in cricetidae_subset_order) {
  species_in_subset <- cricetidae_members_long$Latin.name[
    cricetidae_members_long$Subset == subset_name
  ]

  cricetidae_membership[[subset_name]] <- ifelse(
    cricetidae_membership$Latin.name %in% species_in_subset,
    "Yes",
    "No"
  )
}

cricetidae_membership_out <- cricetidae_membership
names(cricetidae_membership_out)[1] <- "Latin name"

write.csv(
  cricetidae_membership_out,
  cricetidae_output,
  row.names = FALSE
)

cricetidae_formatted <- rbind(
  data.frame(
    `Latin name` = "Supplementary Table S2. Species membership in the Cricetidae nested empirical candidate pools.",
    C512 = "",
    C256 = "",
    C128 = "",
    C64 = "",
    C32 = "",
    check.names = FALSE
  ),
  data.frame(
    `Latin name` = "",
    C512 = "",
    C256 = "",
    C128 = "",
    C64 = "",
    C32 = "",
    check.names = FALSE
  ),
  cricetidae_membership_out
)

write.csv(
  cricetidae_formatted,
  cricetidae_formatted_output,
  row.names = FALSE,
  na = ""
)

## -----------------------------
## Carnivora empirical candidate pool:
## species membership table for Supplementary Table S1
## -----------------------------

selected_subsets_dir <- "/home/galileo-group/galileouser05/huangr/selected_subsets"
carnivora_tree_file <- file.path(
  selected_subsets_dir,
  "MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre"
)
carnivora_taxonomy_file <- file.path(
  selected_subsets_dir,
  "taxonomy_mamPhy_5911species_toPublish.csv"
)
carnivora_mdd_file <- file.path(
  selected_subsets_dir,
  "MDD_v2.0_6759species.csv"
)
carnivora_selected_file <- "/home/galileo-group/galileouser05/huangr/2026_04_distance/results/carnivora/carnivora_selected_species_s20.csv"
carnivora_output <- file.path(
  selected_subsets_dir,
  "Carnivora_clustered_dispersed_neither.csv"
)
carnivora_unmatched_na_output <- file.path(
  selected_subsets_dir,
  "Carnivora_clustered_dispersed_neither_unmatched_as_NA.csv"
)
carnivora_formatted_output <- file.path(
  selected_subsets_dir,
  "Supplementary_Table_S1_Carnivora_membership_formatted.csv"
)
carnivora_manual_synonym_output <- file.path(
  selected_subsets_dir,
  "Carnivora_MDD_manual_synonym_crosswalk.csv"
)

tree <- read.nexus(carnivora_tree_file)
tree <- drop.tip(
  tree,
  intersect(
    c(
      "_Anolis_carolinensis",
      "Nyctinomops_kalinowskii_MOLOSSIDAE_CHIROPTERA"
    ),
    tree$tip.label
  )
)

tax <- read.csv(carnivora_taxonomy_file, stringsAsFactors = FALSE)
mdd <- read.csv(carnivora_mdd_file, stringsAsFactors = FALSE)
carnivora_selected <- read.csv(carnivora_selected_file, stringsAsFactors = FALSE)

tax_carn <- tax[
  toupper(tax[["ord"]]) == "CARNIVORA" &
    tax[["extinct."]] == 0,
]

if ("samp" %in% names(tax_carn)) {
  tax_carn <- tax_carn[tax_carn$samp == "sampled", ]
}

candidate_tips <- sort(intersect(tax_carn[["tiplabel"]], tree$tip.label))
carn_tree <- drop.tip(tree, setdiff(tree$tip.label, candidate_tips))

if (Ntip(carn_tree) != 261) {
  warning("Carnivora candidate pool has ", Ntip(carn_tree), " species, not 261.")
}

tree_species <- data.frame(
  tip_label = carn_tree$tip.label,
  stringsAsFactors = FALSE
) %>%
  mutate(
    sciName = str_remove(tip_label, "_[A-Z]+_CARNIVORA$")
  ) %>%
  left_join(
    tax_carn %>%
      select(tiplabel, tax_family = fam),
    by = c("tip_label" = "tiplabel")
  )

selected_clean <- bind_rows(
  carnivora_selected %>%
    transmute(
      `Subset membership` = "Dispersed",
      tip_label = Dispersed
    ),
  carnivora_selected %>%
    transmute(
      `Subset membership` = "Clustered",
      tip_label = Clustered
    )
) %>%
  filter(!is.na(tip_label), tip_label != "") %>%
  mutate(
    sciName = str_remove(tip_label, "_[A-Z]+_CARNIVORA$")
  ) %>%
  distinct(sciName, `Subset membership`)

duplicated_selected <- selected_clean %>%
  count(sciName) %>%
  filter(n > 1)

if (nrow(duplicated_selected) > 0) {
  stop(
    "Some Carnivora species occur in more than one selected subset: ",
    paste(duplicated_selected$sciName, collapse = ", ")
  )
}

selected_not_in_pool <- setdiff(selected_clean$sciName, tree_species$sciName)
if (length(selected_not_in_pool) > 0) {
  stop(
    "Selected Carnivora species not found in the 261-species pool: ",
    paste(selected_not_in_pool, collapse = ", ")
  )
}

mdd_carn <- mdd %>%
  filter(order == "Carnivora") %>%
  select(
    sciName,
    mainCommonName,
    family,
    CMW_sciName,
    MSW3_sciName
  )

mdd_lookup <- bind_rows(
  mdd_carn %>%
    transmute(matchName = sciName, mainCommonName, family),
  mdd_carn %>%
    filter(!is.na(CMW_sciName), CMW_sciName != "NA", CMW_sciName != "") %>%
    transmute(matchName = CMW_sciName, mainCommonName, family),
  mdd_carn %>%
    filter(!is.na(MSW3_sciName), MSW3_sciName != "NA", MSW3_sciName != "") %>%
    transmute(matchName = MSW3_sciName, mainCommonName, family)
) %>%
  distinct(matchName, .keep_all = TRUE)

## Manual MDD synonym crosswalk.
## These names are retained as Latin names in the phylogeny/taxonomy input,
## but MDD v2.0 treats them as synonyms, older combinations, or replaced
## generic combinations of the accepted names listed in mdd_sciName. This
## crosswalk is used only to retrieve common names and family assignments from
## MDD; it does not rename the species in the final supplementary table.
carnivora_manual_synonyms <- data.frame(
  matchName = c(
    "Canis_anthus",
    "Herpestes_edwardsii",
    "Lutra_maculicollis",
    "Nasuella_meridensis",
    "Otaria_bryonia",
    "Pseudalopex_culpaeus",
    "Pseudalopex_fulvipes",
    "Pseudalopex_griseus",
    "Pseudalopex_gymnocercus",
    "Pseudalopex_sechurae",
    "Pseudalopex_vetulus",
    "Salanoia_durrelli"
  ),
  mdd_sciName = c(
    "Canis_lupaster",
    "Urva_edwardsii",
    "Hydrictis_maculicollis",
    "Nasua_olivacea",
    "Otaria_flavescens",
    "Lycalopex_culpaeus",
    "Lycalopex_fulvipes",
    "Lycalopex_grisea",
    "Lycalopex_gymnocerca",
    "Lycalopex_sechurae",
    "Lycalopex_vetula",
    "Salanoia_concolor"
  ),
  Rationale = c(
    "MDD v2.0 accepts Canis_lupaster; taxonomy notes indicate anthus is not the preferred epithet and may represent a nomen dubium.",
    "MDD v2.0 accepts Urva_edwardsii; MSW3_sciName links this species to the Herpestes-era name.",
    "MDD v2.0 accepts Hydrictis_maculicollis; originalNameCombination is Lutra maculicollis.",
    "MDD v2.0 treats meridensis as a nominal name under Nasua_olivacea.",
    "MDD v2.0 accepts Otaria_flavescens; nominalNames include byronia/byronii and iucnStatus notes use as Otaria byronia.",
    "MDD v2.0 accepts Lycalopex_culpaeus; tree/taxonomy retains the older Pseudalopex combination.",
    "MDD v2.0 accepts Lycalopex_fulvipes; tree/taxonomy retains the older Pseudalopex combination.",
    "MDD v2.0 accepts Lycalopex_griseus; tree/taxonomy retains the older Pseudalopex combination.",
    "MDD v2.0 accepts Lycalopex_gymnocerca; tree/taxonomy retains the older Pseudalopex combination.",
    "MDD v2.0 accepts Lycalopex_sechurae; tree/taxonomy retains the older Pseudalopex combination.",
    "MDD v2.0 accepts Lycalopex_vetula; tree/taxonomy retains the older Pseudalopex combination.",
    "MDD v2.0 treats durrelli as a nominal name under Salanoia_concolor."
  ),
  stringsAsFactors = FALSE
) %>%
  left_join(
    mdd_carn %>%
      select(
        mdd_sciName = sciName,
        mainCommonName,
        family
      ),
    by = "mdd_sciName"
  )

missing_manual_mdd <- carnivora_manual_synonyms %>%
  filter(is.na(mainCommonName) | is.na(family))

if (nrow(missing_manual_mdd) > 0) {
  stop(
    "Manual MDD synonym crosswalk failed to match accepted MDD names: ",
    paste(missing_manual_mdd$mdd_sciName, collapse = ", ")
  )
}

carnivora_manual_synonym_audit <- carnivora_manual_synonyms %>%
  transmute(
    `Tree/taxonomy name` = str_replace_all(matchName, "_", " "),
    `MDD matched name` = str_replace_all(mdd_sciName, "_", " "),
    `Common name from MDD` = mainCommonName,
    `Family from MDD` = toupper(family),
    `Treatment note` = Rationale
  )

write.csv(
  carnivora_manual_synonym_audit,
  carnivora_manual_synonym_output,
  row.names = FALSE
)

carnivora_manual_lookup <- carnivora_manual_synonyms %>%
  select(matchName, mainCommonName, family)

mdd_lookup <- bind_rows(
  carnivora_manual_lookup,
  mdd_lookup
) %>%
  distinct(matchName, .keep_all = TRUE)

make_carnivora_table <- function(lookup_tbl, use_tax_family_fallback = TRUE) {
  tree_species %>%
    left_join(selected_clean, by = "sciName") %>%
    mutate(
      `Subset membership` = ifelse(
        is.na(`Subset membership`),
        "Neither",
        `Subset membership`
      )
    ) %>%
    left_join(lookup_tbl, by = c("sciName" = "matchName")) %>%
    transmute(
      `Subset membership`,
      `Common name` = mainCommonName,
      `Latin name` = str_replace_all(sciName, "_", " "),
      Family = toupper(ifelse(
        is.na(family) & use_tax_family_fallback,
        tax_family,
        family
      ))
    ) %>%
    arrange(
      factor(`Subset membership`, levels = c("Clustered", "Dispersed", "Neither")),
      Family,
      `Latin name`
    )
}

final_tbl <- make_carnivora_table(
  lookup_tbl = mdd_lookup,
  use_tax_family_fallback = TRUE
)

final_tbl_unmatched_na <- make_carnivora_table(
  lookup_tbl = mdd_lookup %>%
    filter(!matchName %in% carnivora_manual_lookup$matchName),
  use_tax_family_fallback = FALSE
)

membership_counts <- table(final_tbl$`Subset membership`)
print(membership_counts)

unmatched <- final_tbl %>%
  filter(is.na(`Common name`) | is.na(Family))

print(unmatched)

write.csv(
  final_tbl,
  carnivora_output,
  row.names = FALSE
)

write.csv(
  final_tbl_unmatched_na,
  carnivora_unmatched_na_output,
  row.names = FALSE,
  na = ""
)

carnivora_formatted <- rbind(
  data.frame(
    `Subset membership` = "Supplementary Table S1 Species membership in the Carnivora empirical candidate pool and selected subsets.",
    `Common name` = "",
    `Latin name` = "",
    Family = "",
    check.names = FALSE
  ),
  data.frame(
    `Subset membership` = "",
    `Common name` = "",
    `Latin name` = "",
    Family = "",
    check.names = FALSE
  ),
  final_tbl
)

write.csv(
  carnivora_formatted,
  carnivora_formatted_output,
  row.names = FALSE,
  na = ""
)
