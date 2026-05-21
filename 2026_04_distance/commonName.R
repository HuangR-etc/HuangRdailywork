library(dplyr)
library(stringr)
setwd("D:/mammal_tre")
# 读取树文件
tree <- read.nexus("MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre")
tax <- read.csv("taxonomy_mamPhy_5911species_toPublish.csv")
sampled <- tax$tiplabel[tax$samp == "sampled"]

setdiff(tree$tip.label, sampled)
# "_Anolis_carolinensis"
# "Nyctinomops_kalinowskii_MOLOSSIDAE_CHIROPTERA"

drop.tips <- c(
  "_Anolis_carolinensis",
  "Nyctinomops_kalinowskii_MOLOSSIDAE_CHIROPTERA"
)

tr_4098 <- drop.tip(tree, drop.tips)

Ntip(tr_4098)

tax_carn <- tax[
  toupper(tax[["ord"]]) == "CARNIVORA" &
    tax[["extinct."]] == 0,
]

# 进一步保留 sampled 物种
if ("samp" %in% names(tax)) {
  tax_carn <- tax_carn[tax_carn$samp == "sampled", ]
}

# 6. 与树的 tip 做交集
candidate_tips <- intersect(tax_carn[["tiplabel"]], tree$tip.label)
carn_tree <- drop.tip(tree, setdiff(tree$tip.label, candidate_tips))

# 8. 候选池大小，即文中的 [N]
N <- Ntip(carn_tree)
#Carnivora食肉目：261物种
N


## -----------------------------
## 1. 你的已整理表
## -----------------------------
subset_tbl <- read.csv("subset_table.csv", stringsAsFactors = FALSE)
## -----------------------------
## 2. MDD 数据
## -----------------------------
mdd <- read.csv("MDD_v2.0_6759species.csv", stringsAsFactors = FALSE)
## -----------------------------
## 3. 从系统发育树 tip.label 提取物种名
## -----------------------------

tree_species <- data.frame(
  tip_label = carn_tree$tip.label,
  stringsAsFactors = FALSE
) %>%
  mutate(
    ## 去掉最后的 _FAMILY_ORDER
    sciName = str_remove(tip_label, "_[A-Z]+_CARNIVORA$")
  )


## -----------------------------
## 4. 整理已分类物种名称
## -----------------------------

subset_clean <- subset_tbl %>%
  mutate(
    sciName = str_replace_all(`Latin.name`, " ", "_")
  )


## -----------------------------
## 5. 找出树上剩下的 neither 物种
## -----------------------------

neither_tbl <- tree_species %>%
  anti_join(subset_clean, by = "sciName") %>%
  left_join(
    mdd %>%
      filter(order == "Carnivora") %>%
      select(
        sciName,
        mainCommonName,
        family
      ),
    by = "sciName"
  ) %>%
  transmute(
    `Subset.membership` = "Neither",
    `Common.name` = mainCommonName,
    `Latin.name` = str_replace_all(sciName, "_", " "),
    Family = toupper(family)
  )


## -----------------------------
## 6. 合并 Clustered / Dispersed / Neither
## -----------------------------

final_tbl <- subset_clean %>%
  select(
    `Subset.membership`,
    `Common.name`,
    `Latin.name`,
    Family
  ) %>%
  bind_rows(neither_tbl) %>%
  arrange(`Subset.membership`, Family, `Latin.name`)


## -----------------------------
## 7. 检查是否有没匹配到 MDD 的物种
## -----------------------------

unmatched <- final_tbl %>%
  filter(is.na(`Common.name`) | is.na(Family))

unmatched


## -----------------------------
## 8. 导出结果
## -----------------------------

write.csv(
  final_tbl,
  "Carnivora_clustered_dispersed_neither.csv",
  row.names = FALSE
)