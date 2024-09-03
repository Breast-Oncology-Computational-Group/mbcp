get_mbcp_colors <- function(set) {

  receptor <- c("HR+/HER2-" = RColorBrewer::brewer.pal(6,"Set1")[3],
                "HER2+" = RColorBrewer::brewer.pal(6,"Set1")[5],
                "TNBC" = RColorBrewer::brewer.pal(6,"Set1")[4],
                "Unknown or indeterminate" = "#9D989E")

  pam50 <- c("LumA" = RColorBrewer::brewer.pal(6,"Set1")[3],
             "LumB" = RColorBrewer::brewer.pal(6,"Set1")[2],
             "Her2" = RColorBrewer::brewer.pal(6,"Set1")[1],
             "Basal"= RColorBrewer::brewer.pal(6,"Set1")[4],
             "Normal" = RColorBrewer::brewer.pal(6,"Set1")[6],
             "NC"= RColorBrewer::brewer.pal(8,"Dark2")[8],
             "No RNA-seq" = "#9D989E")

  dna_alts <- c("DeepDEL" = "#377eb8", "AMP" = "#ff6363", "HighAMP" = "#bc0000",
                "FocalHighAMP" = "#7f0000", "Gain_of_function" = "#d95f02",
                "Loss_of_function" = "#0173B2", "Putative_loss_of_function" = "#17c1e9ff",
                "Biallelic_inactivation"="darkblue", "Multi_hit"="black",
                "Hotspot" = "#008000","Missense_Mutation" = "#00e000",
                "Any Amp"="pink", "Other_Mutation"="#4dbd4d", "WT" = "#9D989E")

  mutation <- c("MUT"="#FC8D62","WT"="#8DA0CB")

  histology <- c("IDC" = RColorBrewer::brewer.pal(6,"Set1")[1],
                 "ILC" = RColorBrewer::brewer.pal(6,"Set1")[2],
                 "DCIS" = RColorBrewer::brewer.pal(6,"Set1")[3],
                 "MIXED_IDLC" = RColorBrewer::brewer.pal(6,"Set1")[4],
                 "ADENOCARCINOMA" = RColorBrewer::brewer.pal(6,"Set1")[5],
                 "CARCINOMA" = RColorBrewer::brewer.pal(8,"Set1")[8],
                 "OTHER_RARE_SUBTYPE" = "yellow",
                 "N/A_BLOOD_SAMPLE" = RColorBrewer::brewer.pal(7,"Set1")[7])

  treatment <- c("Primary Tx Naive" = "darkcyan", "MET Txed" = "orange")

  metastasis <- c("NO_METASTATIC_DISEASE_PRESENT" = "darkcyan", "METASTATIC_DISEASE_PRESENT" = "orange")

  expression <- c("Lower Decile" = "#5e3c99", "Lower Quartile" = "#b2abd2",
                       "Mid" = "#CECECE", "Upper Quartile" = "#fdb863",
                       "Upper Decile" = "#e66101")

  germline_alts <- c("DeepDEL" = "#9754CD", "LOF" = "#0173B2",
                     "LOH" = "#DE8F05", "VUS SNV"="#D55E00",
                     "None Detected" = "#029E73", "No Tumor WES" = "#666666")

  dna_signatures <- c("Clock-like" = "#DD310D", "APOBEC" = "#84C774", "HR" = "#FFB05C",
                      "SBS.12.unknown" = "#9379C8")

  # return vector with same name as set argument with switch
  return(switch(set,
                receptor = receptor,
                pam50 = pam50,
                dna_alts = dna_alts,
                mutation = mutation,
                histology = histology,
                treatment = treatment,
                metastasis = metastasis,
                expression = expression,
                germline_alts = germline_alts,
                dna_signatures = dna_signatures))

}
