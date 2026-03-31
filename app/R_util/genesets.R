library(jsonlite)

## Load gene sets for ssGSEA, GSEA and/or GESECA ##
app_dir = normalizePath(getwd())
pathways_dir = file.path(app_dir, "REF_DATA", "genesets_and_pathways")

to_list = function(df) split(df$gene_symbol, df$gs_name)

all_pathways = readRDS(file.path(pathways_dir, "human_pathways.rds"))
c2_raw = fromJSON(file.path(pathways_dir, "c2_all_genesets.json"))

boyault_data = c2_raw[grepl("BOYAULT", names(c2_raw))]
boyault_sets = lapply(boyault_data, function(entry) entry[["geneSymbols"]])

c2_gene_sets = lapply(c2_raw, function(entry) entry[["geneSymbols"]])


REACTOME_pathways = c2_gene_sets[grepl("REACTOME", names(c2_gene_sets), fixed = TRUE)]
KEGG_pathways = c2_gene_sets[grepl("KEGG", names(c2_gene_sets), fixed = TRUE)]

GOBP_pathways = all_pathways$c5[grepl("GOBP", names(all_pathways$c5), fixed = TRUE)]


#Freshly added but not working yet lol
# transcript_factor_target = all_pathways$c3[grepl("TARGET_GENES", names(all_pathways$c3), fixed = TRUE)]

stress_msig_names = c(
  # proximal adrenergic / second messenger
  "GOBP_ADRENERGIC_RECEPTOR_SIGNALING_PATHWAY", "BIOCARTA_CREB_PATHWAY", "GOBP_RESPONSE_TO_MONOAMINE",
  "REACTOME_NOREPINEPHRINE_NEUROTRANSMITTER_RELEASE_CYCLE",
  # inflammatory / CTRA-adjacent
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  # stress-linked tumor outputs
  "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_ANGIOGENESIS", "PID_LYMPH_ANGIOGENESIS_PATHWAY",
  # neural remodeling
  "GOBP_NEUROTROPHIN_TRK_RECEPTOR_SIGNALING_PATHWAY", "GOBP_AXON_GUIDANCE",
  # immune polarization / dysfunction
  # these DN sets are intentionally chosen because enrichment means "more M2" or "more exhausted"
  "COATES_MACROPHAGE_M1_VS_M2_DN", "GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN", "GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_DN", "GSE24026_PD1_LIGATION_VS_CTRL_IN_ACT_TCELL_LINE_UP"
)

paths_hallmark = intersect(stress_msig_names, names(all_pathways$h))
paths_c2 = intersect(stress_msig_names, names(c2_gene_sets))
paths_c5 = intersect(stress_msig_names, names(all_pathways$c5))

QoL_pathways = c(c2_gene_sets[paths_c2], all_pathways$c5[paths_c5], all_pathways$h[paths_hallmark])

