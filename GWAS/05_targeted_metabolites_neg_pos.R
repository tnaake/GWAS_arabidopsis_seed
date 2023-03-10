################################################################################
################################# negative mode ################################
################################################################################

options(stringsAsFactors = FALSE)
setwd("~/GitHub/GWAS_arabidopsis_seed/")
rel_rep1 <- read.csv("rep1_negative_match.csv", sep = "\t")
rel_rep2 <- read.csv("rep2_negative_match.csv", sep = "\t")

## 4-methylsulfinylbutyl glucosinolate ## mz 436.041634
rel_rep1[rel_rep1[, "mz"] < 436.0421 & rel_rep1[, "mz"] > 436.0410, c("Name", "RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 436.0421 & rel_rep2[, "mz"] > 436.0411, c("RT", "mz", "mapping_rep1_rep2")]

## 5-methylsulfinylpentyl glucosinolate ## mz 450.0563348
rel_rep1[rel_rep1[, "mz"] < 450.0568 & rel_rep1[, "mz"] > 450.0557, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 450.0578 & rel_rep2[, "mz"] > 450.0557, c("RT", "mz", "mapping_rep1_rep2")]

## 6-methylsulfinylhexyl glucosinolate ## mz 464.0717755
rel_rep1[rel_rep1[, "mz"] < 464.0731 & rel_rep1[, "mz"] > 464.0713, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 464.0731 & rel_rep2[, "mz"] > 464.0713, c("RT", "mz", "mapping_rep1_rep2")]

## 7-methylsulfinylheptyl glucosinolate ## mz 478.088666
rel_rep1[rel_rep1[, "mz"] < 478.0888 & rel_rep1[, "mz"] > 478.0874, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 478.0888 & rel_rep2[, "mz"] > 478.0874, c("RT", "mz", "mapping_rep1_rep2")]

## 8-Methylsulfinyloctyl glucosinolate ## mz 492.1043394
rel_rep1[rel_rep1[, "mz"] < 492.1045 & rel_rep1[, "mz"] > 492.1031, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 492.1045 & rel_rep2[, "mz"] > 492.1031, c("RT", "mz", "mapping_rep1_rep2")]

## 3-hydroxypropyl-glucosinolate ## mz 376.0374431
rel_rep1[rel_rep1[, "mz"] < 376.0376 & rel_rep1[, "mz"] > 376.0362, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 376.0376 & rel_rep2[, "mz"] > 376.0372, c("RT", "mz", "mapping_rep1_rep2")]

## 4-hydroxybutylglucosinolate ## mz 390.0535754
rel_rep1[rel_rep1[, "mz"] < 390.0537 & rel_rep1[, "mz"] > 390.0533, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 390.0537 & rel_rep2[, "mz"] > 390.0533, c("RT", "mz", "mapping_rep1_rep2")]

## 3_methylthiopropylglucosinolate ## mz 406.03
rel_rep1[rel_rep1[, "mz"] < 406.0302 & rel_rep1[, "mz"] > 406.0298, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 406.0303 & rel_rep2[, "mz"] > 406.0290, c("RT", "mz", "mapping_rep1_rep2")]

## 4-methylthiobutyl glucosinolate ## mz 420.0461484
rel_rep1[rel_rep1[, "mz"] < 420.0463 & rel_rep1[, "mz"] > 420.0450, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 420.0463 & rel_rep2[, "mz"] > 420.0459, c("RT", "mz", "mapping_rep1_rep2")]

## 5-methylthiopentylglucosinolate ## mz 434.0618295
rel_rep1[rel_rep1[, "mz"] < 434.0620 & rel_rep1[, "mz"] > 434.0610, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 434.0620 & rel_rep2[, "mz"] > 434.0616, c("RT", "mz", "mapping_rep1_rep2")]

## 6-methylthiohexylglucosinolate ## mz 448.0777523
rel_rep1[rel_rep1[, "mz"] < 448.0779 & rel_rep1[, "mz"] > 448.0775, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 448.0779 & rel_rep2[, "mz"] > 448.0770, c("RT", "mz", "mapping_rep1_rep2")]

## 7-methylthioheptyl glucosinolate ## mz 462.0931615
rel_rep1[rel_rep1[, "mz"] < 462.0933 & rel_rep1[, "mz"] > 462.0929, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 462.0933 & rel_rep2[, "mz"] > 462.0929, c("RT", "mz", "mapping_rep1_rep2")]

## 8-Methylthiooctyl glucosinolate ## mz 476.1093415
rel_rep1[rel_rep1[, "mz"] < 476.1095 & rel_rep1[, "mz"] > 476.1031, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 476.1095 & rel_rep2[, "mz"] > 476.1031, c("RT", "mz", "mapping_rep1_rep2")]

## 3-indolymethyl glucosinolate ## mz 447.0537959 ## not found
rel_rep1[rel_rep1[, "mz"] < 447.0559 & rel_rep1[, "mz"] > 447.0515, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 447.0949 & rel_rep2[, "mz"] > 447.0225, c("RT", "mz", "mapping_rep1_rep2")]

## 1-methoxy-3-indolylmethyl-glucosinolate ## mz 477.0635223
## 4-methoxy-3-indolylmethyl-glucosinolate ## mz 477.0633195
rel_rep1[rel_rep1[, "mz"] < 477.0655 & rel_rep1[, "mz"] > 477.0611, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 477.0655 & rel_rep2[, "mz"] > 477.0611, c("RT", "mz", "mapping_rep1_rep2")]

## 3-benzoyloxypropyl-glucosinolate ## mz 480.0638973
rel_rep1[rel_rep1[, "mz"] < 480.0643 & rel_rep1[, "mz"] > 480.0631, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 480.0640 & rel_rep2[, "mz"] > 480.0636, c("RT", "mz", "mapping_rep1_rep2")]

## 4-benzoyloxybutylglucosinolate ## mz 494.0799077
rel_rep1[rel_rep1[, "mz"] < 494.0801 & rel_rep1[, "mz"] > 494.0788, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 494.0801 & rel_rep2[, "mz"] > 494.0788, c("RT", "mz", "mapping_rep1_rep2")]

## 5_Benzoyloxypentylglucosinolate ## mz 508.0947
rel_rep1[rel_rep1[, "mz"] < 508.0959 & rel_rep1[, "mz"] > 508.0935, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 508.0959 & rel_rep2[, "mz"] > 508.0935, c("RT", "mz", "mapping_rep1_rep2")]

## 6_Benzoyloxyhexylglucosinolate ## mz 522.1104
rel_rep1[rel_rep1[, "mz"] < 522.1118 & rel_rep1[, "mz"] > 522.1090, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 522.1118 & rel_rep2[, "mz"] > 522.1090, c("RT", "mz", "mapping_rep1_rep2")]

## 3-sinapoyloxypropylglucosinolate ## mz 582.095379
rel_rep1[rel_rep1[, "mz"] < 582.0965 & rel_rep1[, "mz"] > 582.0941, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 582.0965 & rel_rep2[, "mz"] > 582.0941, c("RT", "mz", "mapping_rep1_rep2")]

## 4-sinapoyloxybutylglucosinolate ## 596.1114532
rel_rep1[rel_rep1[, "mz"] < 596.1116 & rel_rep1[, "mz"] > 596.1112, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 596.1236 & rel_rep2[, "mz"] > 596.0982, c("RT", "mz", "mapping_rep1_rep2")]

## 3_methylsulfinylpropylglucosinolate ## mz 422.0249
rel_rep1[rel_rep1[, "mz"] < 422.0255 & rel_rep1[, "mz"] > 422.0245, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 422.0255 & rel_rep2[, "mz"] > 422.0245, c("RT", "mz", "mapping_rep1_rep2")] 

## 2_hydroxy_3_butenyl_glucosinolate ## mz 388.0372
rel_rep1[rel_rep1[, "mz"] < 388.0380 & rel_rep1[, "mz"] > 388.0370, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 388.0380 & rel_rep2[, "mz"] > 388.0370, c("RT", "mz", "mapping_rep1_rep2")]

## 3_butenyl_glucosinolate ## mz 372.0423
rel_rep1[rel_rep1[, "mz"] < 372.0429 & rel_rep1[, "mz"] > 372.0415, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 372.0429 & rel_rep2[, "mz"] > 372.0415, c("RT", "mz", "mapping_rep1_rep2")]

## 4_pentenyl_glucosinolate ## mz 386.0579
rel_rep1[rel_rep1[, "mz"] < 386.0683 & rel_rep1[, "mz"] > 386.0576, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 386.0683 & rel_rep2[, "mz"] > 386.0476, c("RT", "mz", "mapping_rep1_rep2")]

## (epi)Catechin_Hex_isomer1/2 ## mz 451.124
rel_rep1[rel_rep1[, "mz"] < 451.1250 & rel_rep1[, "mz"] > 451.1230, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 451.1250 & rel_rep2[, "mz"] > 451.1230, c("RT", "mz", "mapping_rep1_rep2")]

## sinapoylcholine_sinapoylcholine_dimer_isomer1/2 ## mz 294.1342, 354.1551
rel_rep1[rel_rep1[, "mz"] < 294.1344 & rel_rep1[, "mz"] > 294.1340, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 354.1575 & rel_rep1[, "mz"] > 354.1547, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 294.1364 & rel_rep2[, "mz"] > 294.1340, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 354.1575 & rel_rep2[, "mz"] > 354.1547, c("RT", "mz", "mapping_rep1_rep2")]

## 1_O_sinapoylglucose ## mz 385.1135
rel_rep1[rel_rep1[, "mz"] < 385.1137 & rel_rep1[, "mz"] > 385.1133, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 385.1137 & rel_rep2[, "mz"] > 385.1133, c("RT", "mz", "mapping_rep1_rep2")]

## N_((4'_O_glycosyl)_sinapoyl)_N'_sinapoylspermidine/N_Sinapoyl_N_(3,5_dimethoxy_4_hexosyloxycinnamoyl)spermidine_isomer2
## mz 718.3187
rel_rep1[rel_rep1[, "mz"] < 718.3193 & rel_rep1[, "mz"] > 718.3183, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 718.3193 & rel_rep2[, "mz"] > 718.3183, c("RT", "mz", "mapping_rep1_rep2")]

## isorhamnetin_DeoxyHex_Hex ## mz 623.1612
rel_rep1[rel_rep1[, "mz"] < 623.1614 & rel_rep1[, "mz"] > 623.1610, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 623.1614 & rel_rep2[, "mz"] > 623.1610, c("RT", "mz", "mapping_rep1_rep2")]

## N1_N8_di(sinapoyl)_spermidine/N_N_Bis(sinapoyl)spermidine ## mz 556.2659 ## not found
rel_rep1[rel_rep1[, "mz"] < 556.2663 & rel_rep1[, "mz"] > 556.2656, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 556.2663 & rel_rep2[, "mz"] > 556.2656, c("RT", "mz", "mapping_rep1_rep2")]

## 2_O_Sinapoylmalate ## mz 339.0716 223.0606
rel_rep1[rel_rep2[, "mz"] < 339.0719 & rel_rep1[, "mz"] > 339.0712, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep2[, "mz"] < 223.0607 & rel_rep1[, "mz"] > 223.0605, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 339.0719 & rel_rep2[, "mz"] > 339.0712, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 223.0607 & rel_rep2[, "mz"] > 223.0605, c("RT", "mz", "mapping_rep1_rep2")]

## isorhamnetin_(DeoxyHex)2 ## mz 607.1663 ## not found
rel_rep1[rel_rep1[, "mz"] < 607.1666 & rel_rep1[, "mz"] > 607.1660, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 607.1666 & rel_rep2[, "mz"] > 607.1660, c("RT", "mz", "mapping_rep1_rep2")]

## 3_[2_(4_Hydroxy_3_methoxyphenyl)_3_hydroxymethyl_7_methoxy_2_3_dihydrobenzo_furan_5_yl]acryloylcholine
## mz 456.2022 ## not found

rel_rep1[rel_rep1[, "mz"] < 456.2026 & rel_rep1[, "mz"] > 456.2018, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 456.2026 & rel_rep2[, "mz"] > 456.2018, c("RT", "mz", "mapping_rep1_rep2")]
## quercetin_deoxyhexose ## mz 895.197485 449.10784

rel_rep1[rel_rep1[, "mz"] < 895.1988 & rel_rep1[, "mz"] > 895.1965, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 449.1082 & rel_rep1[, "mz"] > 449.1074, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 895.1988 & rel_rep2[, "mz"] > 895.1965, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 449.1082 & rel_rep2[, "mz"] > 449.1074, c("RT", "mz", "mapping_rep1_rep2")]

## isorhamnetin_hexose ## mz 477.1033
rel_rep1[rel_rep1[, "mz"] < 477.1036 & rel_rep1[, "mz"] > 477.1030, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 477.1036 & rel_rep2[, "mz"] > 477.1030, c("RT", "mz", "mapping_rep1_rep2")]

## kaempferol_deoxyhex ## mz 431.0978
rel_rep1[rel_rep1[, "mz"] < 431.0982 & rel_rep1[, "mz"] > 431.0976, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 431.0982 & rel_rep2[, "mz"] > 431.0976, c("RT", "mz", "mapping_rep1_rep2")]


################################################################################
################################# positive mode ################################
################################################################################

options(stringsAsFactors = FALSE)
rel_rep1 <- read.csv("rep1_positive_match.csv", sep="\t")
rel_rep2 <- read.csv("rep2_positive_match.csv", sep="\t")

## 3_methylsulfinylpropylglucosinolate ## mz 461.9965 446.0225 366.0657 344.0838 182.0309 132.0483 114.0377
rel_rep1[rel_rep1[, "mz"] < 461.9969 & rel_rep1[, "mz"] > 461.9961, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 446.0227 & rel_rep1[, "mz"] > 446.0223, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 366.0659 & rel_rep1[, "mz"] > 366.0655, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 344.0840 & rel_rep1[, "mz"] > 344.0836, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 182.0311 & rel_rep1[, "mz"] > 182.0307, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 132.0485 & rel_rep1[, "mz"] > 132.0481, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 114.0379 & rel_rep1[, "mz"] > 114.0375, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 461.9969 & rel_rep2[, "mz"] > 461.9961, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 446.0227 & rel_rep2[, "mz"] > 446.0223, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 366.0659 & rel_rep2[, "mz"] > 366.0655, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 344.0840 & rel_rep2[, "mz"] > 344.0836, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 182.0311 & rel_rep2[, "mz"] > 182.0307, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 132.0485 & rel_rep2[, "mz"] > 132.0481, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 114.0379 & rel_rep2[, "mz"] > 114.0375, c("RT", "mz", "mapping_rep1_rep2")]

## 4-methylsulfinylbutyl glucosinolate ## mz 438.0556171
rel_rep1[rel_rep1[, "mz"] < 438.0558 & rel_rep1[, "mz"] > 438.0554, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 438.0558 & rel_rep2[, "mz"] > 438.0554, c("RT", "mz", "mapping_rep1_rep2")]

## 5-methylsulfinylpentyl glucosinolate ## mz 452.0714
rel_rep1[rel_rep1[, "mz"] < 452.0716 & rel_rep1[, "mz"] > 452.0712, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 452.0716 & rel_rep2[, "mz"] > 452.0712, c("RT", "mz", "mapping_rep1_rep2")]

## 6-methylsulfinylhexyl glucosinolate ## mz 466.0866
rel_rep1[rel_rep1[, "mz"] < 466.0868 & rel_rep1[, "mz"] > 466.0864, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 466.0868 & rel_rep2[, "mz"] > 466.0864, c("RT", "mz", "mapping_rep1_rep2")]

## 7-methylsulfinylheptyl glucosinolate ## mz 480.1026
rel_rep1[rel_rep1[, "mz"] < 480.1028 & rel_rep1[, "mz"] > 480.1024, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 480.1028 & rel_rep2[, "mz"] > 480.1024, c("RT", "mz", "mapping_rep1_rep2")]

## 8-Methylsulfinyloctyl glucosinolate ## mz 494.118
rel_rep1[rel_rep1[, "mz"] < 494.1182 & rel_rep1[, "mz"] > 494.1178, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 494.1182 & rel_rep2[, "mz"] > 494.1178, c("RT", "mz", "mapping_rep1_rep2")]

## 3-hydroxypropyl-glucosinolate ## mz 378.0518
rel_rep1[rel_rep1[, "mz"] < 378.0520 & rel_rep1[, "mz"] > 378.0516, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 378.0520 & rel_rep2[, "mz"] > 378.0516, c("RT", "mz", "mapping_rep1_rep2")]

## 4-hydroxybutylglucosinolate ## mz 392.0673
rel_rep1[rel_rep1[, "mz"] < 392.0675 & rel_rep1[, "mz"] > 392.0671, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 392.0675 & rel_rep2[, "mz"] > 392.0671, c("RT", "mz", "mapping_rep1_rep2")]

## 3_methylthiopropylglucosinolate ## mz 446.0016 430.0276 408.0457 350.0708 328.0889 166.036 148.0255 100.0221
rel_rep1[rel_rep1[, "mz"] < 446.0018 & rel_rep1[, "mz"] > 446.0014, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 430.0278 & rel_rep1[, "mz"] > 430.0274, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 408.0459 & rel_rep1[, "mz"] > 408.0455, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 350.0710 & rel_rep1[, "mz"] > 350.0706, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 328.0891 & rel_rep1[, "mz"] > 328.0887, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 166.0362 & rel_rep1[, "mz"] > 166.0358, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 148.0257 & rel_rep1[, "mz"] > 148.0253, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 100.0219 & rel_rep1[, "mz"] > 100.0223, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 446.0018 & rel_rep2[, "mz"] > 446.0014, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 430.0278 & rel_rep2[, "mz"] > 430.0274, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 408.0459 & rel_rep2[, "mz"] > 408.0455, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 350.0710 & rel_rep2[, "mz"] > 350.0706, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 328.0891 & rel_rep2[, "mz"] > 328.0887, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 166.0362 & rel_rep2[, "mz"] > 166.0358, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 148.0257 & rel_rep2[, "mz"] > 148.0253, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 100.0219 & rel_rep2[, "mz"] > 100.0223, c("RT", "mz", "mapping_rep1_rep2")]

## 4-methylthiobutyl glucosinolate ## mz 422.0598
rel_rep1[rel_rep1[, "mz"] < 422.0600 & rel_rep1[, "mz"] > 422.0596, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 422.0600 & rel_rep2[, "mz"] > 422.0596, c("RT", "mz", "mapping_rep1_rep2")]

## 5-methylthiopentylglucosinolate ## mz 436.0759
rel_rep1[rel_rep1[, "mz"] < 436.0761 & rel_rep1[, "mz"] > 436.0757, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 436.0761 & rel_rep2[, "mz"] > 436.0757, c("RT", "mz", "mapping_rep1_rep2")]

## 6-methylthiohexylglucosinolate ## mz 450.0916
rel_rep1[rel_rep1[, "mz"] < 450.0918 & rel_rep1[, "mz"] > 450.0914, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 450.0918 & rel_rep2[, "mz"] > 450.0914, c("RT", "mz", "mapping_rep1_rep2")]

## 7-methylthioheptyl glucosinolate ## mz 464.1082
rel_rep1[rel_rep1[, "mz"] < 464.1084 & rel_rep1[, "mz"] > 464.1080, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 464.1084 & rel_rep2[, "mz"] > 464.1080, c("RT", "mz", "mapping_rep1_rep2")]

## 8-Methylthiooctyl glucosinolate ## mz 478.1230
rel_rep1[rel_rep1[, "mz"] < 478.1232 & rel_rep1[, "mz"] > 478.1228, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 478.1232 & rel_rep2[, "mz"] > 478.1228, c("RT", "mz", "mapping_rep1_rep2")]

## 3-indolymethyl glucosinolate ## mz 449.0680
rel_rep1[rel_rep1[, "mz"] < 449.0682 & rel_rep1[, "mz"] > 449.0678, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 449.0682 & rel_rep2[, "mz"] > 449.0678, c("RT", "mz", "mapping_rep1_rep2")]

## 4-methylthiobutyldesulfoglucosinolate ## mz 342.1037
rel_rep1[rel_rep1[, "mz"] < 342.1039 & rel_rep1[, "mz"] > 342.1035, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 342.1039 & rel_rep2[, "mz"] > 342.1035, c("RT", "mz", "mapping_rep1_rep2")]

## 5-methylthiopentyldesulfoglucosinolate ## mz 356.1193
rel_rep1[rel_rep1[, "mz"] < 356.1195 & rel_rep1[, "mz"] > 356.1191, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 356.1195 & rel_rep2[, "mz"] > 356.1191, c("RT", "mz", "mapping_rep1_rep2")]

## 6-methylthiohexyldesulfoglucosinolate ## mz 370.1348
rel_rep1[rel_rep1[, "mz"] < 370.1350 & rel_rep1[, "mz"] > 370.1346, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 370.1350 & rel_rep2[, "mz"] > 370.1346, c("RT", "mz", "mapping_rep1_rep2")]

## 7-methylthioheptyldesulfoglucosinolate ## mz 384.1506
rel_rep1[rel_rep1[, "mz"] < 384.1508 & rel_rep1[, "mz"] > 384.1504, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 384.1508 & rel_rep2[, "mz"] > 384.1504, c("RT", "mz", "mapping_rep1_rep2")]

## 3-benzoyloxypropyl-glucosinolate ## mz 482.0788
rel_rep1[rel_rep1[, "mz"] < 482.0790 & rel_rep1[, "mz"] > 482.0786, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 482.0790 & rel_rep2[, "mz"] > 482.0786, c("RT", "mz", "mapping_rep1_rep2")]

## 4-benzoyloxybutylglucosinolate ## mz 496.0943
rel_rep1[rel_rep1[, "mz"] < 496.0945 & rel_rep1[, "mz"] > 496.0941, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 496.0945 & rel_rep2[, "mz"] > 496.0941, c("RT", "mz", "mapping_rep1_rep2")]

## 5_Benzoyloxypentylglucosinolate ## mz 548.0663 532.0923 452.1355 430.1536
rel_rep1[rel_rep1[, "mz"] < 548.0665 & rel_rep1[, "mz"] > 548.0661, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 532.0925 & rel_rep1[, "mz"] > 532.0921, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 452.1357 & rel_rep1[, "mz"] > 452.1353, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 430.1538 & rel_rep1[, "mz"] > 430.1534, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 548.0665 & rel_rep2[, "mz"] > 548.0661, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 532.0925 & rel_rep2[, "mz"] > 532.0921, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 452.1357 & rel_rep2[, "mz"] > 452.1353, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 430.1538 & rel_rep2[, "mz"] > 430.1534, c("RT", "mz", "mapping_rep1_rep2")]

## 6_Benzoyloxyhexylglucosinolate ## mz 562.0819 546.108 466.1512 444.1692
rel_rep1[rel_rep1[, "mz"] < 562.0821 & rel_rep1[, "mz"] > 562.0817, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 546.1082 & rel_rep1[, "mz"] > 546.1078, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 466.1514 & rel_rep1[, "mz"] > 466.1510, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 444.1694 & rel_rep1[, "mz"] > 444.1690, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 562.0821 & rel_rep2[, "mz"] > 562.0817, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 546.1082 & rel_rep2[, "mz"] > 546.1078, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 466.1514 & rel_rep2[, "mz"] > 466.1510, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 444.1694 & rel_rep2[, "mz"] > 444.1690, c("RT", "mz", "mapping_rep1_rep2")]

## 2_Hydroxy_3_butenyl_glucosinolate ## mz 428.0087 412.0348
rel_rep1[rel_rep1[, "mz"] < 428.0089 & rel_rep1[, "mz"] > 428.0085, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 412.0350 & rel_rep1[, "mz"] > 412.0346, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 428.0089 & rel_rep2[, "mz"] > 428.0085, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 412.0350 & rel_rep2[, "mz"] > 412.0346, c("RT", "mz", "mapping_rep1_rep2")]

## 4_pentenyl_glucosinolate ## mz 426.0295 410.0555
rel_rep1[rel_rep1[, "mz"] < 426.0297 & rel_rep1[, "mz"] > 426.0293, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 410.0557 & rel_rep1[, "mz"] > 410.0553, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 426.0297 & rel_rep2[, "mz"] > 426.0293, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 410.0557 & rel_rep2[, "mz"] > 410.0553, c("RT", "mz", "mapping_rep1_rep2")]

## Kaempferol 3-O-rhamnoside 7-O-rhamnoside ## mz 579.1696
rel_rep1[rel_rep1[, "mz"] < 579.1698 & rel_rep1[, "mz"] > 579.1694, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 579.1698 & rel_rep2[, "mz"] > 579.1694, c("RT", "mz", "mapping_rep1_rep2")]

## Kaempferol 3-O-glucoside 7-O-rhamnoside ## mz 595.1661
rel_rep1[rel_rep1[, "mz"] < 595.1663 & rel_rep1[, "mz"] > 595.1659, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 595.1663 & rel_rep2[, "mz"] > 595.1659, c("RT", "mz", "mapping_rep1_rep2")]

## Quercetin 3-O-rhamnoside 7-O-rhamnoside ## mz 595.1665
rel_rep1[rel_rep1[, "mz"] < 595.1667 & rel_rep1[, "mz"] > 595.1663, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 595.1667 & rel_rep2[, "mz"] > 595.1663, c("RT", "mz", "mapping_rep1_rep2")]

## Kaempferol 3-O-[2''-O-(rhamnosyl) glucoside] 7-O-rhamnoside ## mz 741.2225
rel_rep1[rel_rep1[, "mz"] < 741.2227 & rel_rep1[, "mz"] > 741.2223, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 741.2227 & rel_rep2[, "mz"] > 741.2223, c("RT", "mz", "mapping_rep1_rep2")]

## Kaempferol 3-O-glucosyl-glucoside 7-O-rhamnoside ## mz 757.2175
rel_rep1[rel_rep1[, "mz"] < 757.2177 & rel_rep1[, "mz"] > 757.2173, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 757.2177 & rel_rep2[, "mz"] > 757.2173, c("RT", "mz", "mapping_rep1_rep2")]

## Quercetin 3-O-[2''-O-(rhamnosyl) glucoside] 7-O-rhamnoside ## mz 757.2173
rel_rep1[rel_rep1[, "mz"] < 757.2175 & rel_rep1[, "mz"] > 757.2171, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 757.2175 & rel_rep2[, "mz"] > 757.2171, c("RT", "mz", "mapping_rep1_rep2")]

## Quercetin 3-O-glucoside 7-O-rhamnoside ## mz 611.1602
rel_rep1[rel_rep1[, "mz"] < 611.1604 & rel_rep1[, "mz"] > 611.1600, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 611.1604 & rel_rep2[, "mz"] > 611.1600, c("RT", "mz", "mapping_rep1_rep2")]

## glutamic acid ## mz 148.0604 130.0499 102.055
rel_rep1[rel_rep1[, "mz"] < 148.0606 & rel_rep1[, "mz"] > 148.0602, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 130.0500 & rel_rep1[, "mz"] > 130.0497, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 102.0552 & rel_rep1[, "mz"] > 102.0548, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 148.0606 & rel_rep2[, "mz"] > 148.0602, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 130.0500 & rel_rep2[, "mz"] > 130.0497, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 102.0552 & rel_rep2[, "mz"] > 102.0548, c("RT", "mz", "mapping_rep1_rep2")]

## valine ## mz 118.0863
rel_rep1[rel_rep1[, "mz"] < 118.0865 & rel_rep1[, "mz"] > 118.0861, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 118.0865 & rel_rep2[, "mz"] > 118.0861, c("RT", "mz", "mapping_rep1_rep2")]

## proline ## mz 116.0706
rel_rep1[rel_rep1[, "mz"] < 116.0708 & rel_rep1[, "mz"] > 116.0704, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 116.0708 & rel_rep2[, "mz"] > 116.0704, c("RT", "mz", "mapping_rep1_rep2")]

## 4_Hexosyloxybenzoylcholine ## mz 386.18094 165.0546
rel_rep1[rel_rep1[, "mz"] < 386.1811 & rel_rep1[, "mz"] > 386.1809, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 165.0548 & rel_rep1[, "mz"] > 165.0544, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 386.1811 & rel_rep2[, "mz"] > 386.1809, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 165.0548 & rel_rep2[, "mz"] > 165.0544, c("RT", "mz", "mapping_rep1_rep2")]

## 3_butenyl_glucosinolate ## mz 412.0138 396.0399
rel_rep1[rel_rep1[, "mz"] < 412.0140 & rel_rep1[, "mz"] > 412.0136, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 396.0401 & rel_rep1[, "mz"] > 396.0397, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 412.0140 & rel_rep2[, "mz"] > 412.0136, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 396.0401 & rel_rep2[, "mz"] > 396.0397, c("RT", "mz", "mapping_rep1_rep2")]

## 4-Hexosyloxy-3-methoxybenzoylcholine ## mz 416.19151 195.0652 151.039
rel_rep1[rel_rep1[, "mz"] < 416.1917 & rel_rep1[, "mz"] > 416.1913, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 195.0654 & rel_rep1[, "mz"] > 195.0650, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 151.0392 & rel_rep1[, "mz"] > 151.0388, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 416.1917 & rel_rep2[, "mz"] > 416.1913, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 195.0654 & rel_rep2[, "mz"] > 195.0650, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 151.0392 & rel_rep2[, "mz"] > 151.0388, c("RT", "mz", "mapping_rep1_rep2")]

## phenylalanine ## mz 166.0863 120.0808 103.0543
rel_rep1[rel_rep1[, "mz"] < 166.0865 & rel_rep1[, "mz"] > 166.0861, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 120.0810 & rel_rep1[, "mz"] > 120.0806, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 103.0545 & rel_rep1[, "mz"] > 103.0541, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 166.0865 & rel_rep2[, "mz"] > 166.0861, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 120.0810 & rel_rep2[, "mz"] > 120.0806, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 103.0545 & rel_rep2[, "mz"] > 103.0541, c("RT", "mz", "mapping_rep1_rep2")]

## 3,5-Dimethoxy-4-hexosyloxybenzoylcholine ## mz 446.2021 225.0758
rel_rep1[rel_rep1[, "mz"] < 446.2023 & rel_rep1[, "mz"] > 446.2019, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 225.0760 & rel_rep1[, "mz"] > 225.0756, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 446.2023 & rel_rep2[, "mz"] > 446.2019, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 225.0760 & rel_rep2[, "mz"] > 225.0756, c("RT", "mz", "mapping_rep1_rep2")]

## 4_Hexosyloxycinnamoylcholine ## mz 412.19659 147.0441
rel_rep1[rel_rep1[, "mz"] < 412.1967 & rel_rep1[, "mz"] > 412.1963, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 147.0443 & rel_rep1[, "mz"] > 147.0439, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 412.1967 & rel_rep2[, "mz"] > 412.1963, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 147.0443 & rel_rep2[, "mz"] > 147.0439, c("RT", "mz", "mapping_rep1_rep2")]

## 8_(Methylsulfinyl)octane_1_amine ## mz 192.1417 128.1434
rel_rep1[rel_rep1[, "mz"] < 192.1419 & rel_rep1[, "mz"] > 192.1415, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 128.1436 & rel_rep1[, "mz"] > 128.1432, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 192.1419 & rel_rep2[, "mz"] > 192.1415, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 128.1436 & rel_rep2[, "mz"] > 128.1432, c("RT", "mz", "mapping_rep1_rep2")]

## 4hexosyloxy_3_5_dimethoxy_cinnamoylcholine ## mz 472.21773 251.0914
rel_rep1[rel_rep1[, "mz"] < 472.2179 & rel_rep1[, "mz"] > 472.2175, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 251.0916 & rel_rep1[, "mz"] > 251.0912, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 472.2179 & rel_rep2[, "mz"] > 472.2175, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 251.0916 & rel_rep2[, "mz"] > 251.0912, c("RT", "mz", "mapping_rep1_rep2")]

## tryptophane ## mz 205.0972 188.0706 170.06 159.0917 146.06 144.0808 130.0651 118.065
rel_rep1[rel_rep1[, "mz"] < 205.0974 & rel_rep1[, "mz"] > 205.0970, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 188.0708 & rel_rep1[, "mz"] > 188.0704, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 170.0602 & rel_rep1[, "mz"] > 170.0598, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 159.0919 & rel_rep1[, "mz"] > 159.0915, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 146.0602 & rel_rep1[, "mz"] > 146.0598, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 144.0810 & rel_rep1[, "mz"] > 144.0806, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 130.0653 & rel_rep1[, "mz"] > 130.0649, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 118.0652 & rel_rep1[, "mz"] > 118.0648, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 205.0974 & rel_rep2[, "mz"] > 205.0970, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 188.0708 & rel_rep2[, "mz"] > 188.0704, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 170.0602 & rel_rep2[, "mz"] > 170.0598, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 159.0919 & rel_rep2[, "mz"] > 159.0915, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 146.0602 & rel_rep2[, "mz"] > 146.0598, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 144.0810 & rel_rep2[, "mz"] > 144.0806, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 130.0653 & rel_rep2[, "mz"] > 130.0649, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 118.0652 & rel_rep2[, "mz"] > 118.0648, c("RT", "mz", "mapping_rep1_rep2")]

## (epi)Catechin_Hex_isomer1/2 ## mz 453.1391 291.0863
rel_rep1[rel_rep1[, "mz"] < 453.1393 & rel_rep1[, "mz"] > 453.1389, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 291.0865 & rel_rep1[, "mz"] > 291.0861, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 453.1393 & rel_rep2[, "mz"] > 453.1389, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 291.0865 & rel_rep2[, "mz"] > 291.0861, c("RT", "mz", "mapping_rep1_rep2")]

## 4hexosyloxy_3_5_dimethoxy_cinnamoylcholine_isomer2 ## mz 472.21773 251.0914
rel_rep1[rel_rep1[, "mz"] < 472.2179 & rel_rep1[, "mz"] > 472.2175, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 251.0916 & rel_rep1[, "mz"] > 251.0912, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 472.2179 & rel_rep2[, "mz"] > 472.2175, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 251.0916 & rel_rep2[, "mz"] > 251.0912, c("RT", "mz", "mapping_rep1_rep2")]

## hydroxyferuloylcholine ## mz 296.14925 237.0758
rel_rep1[rel_rep1[, "mz"] < 296.1494 & rel_rep1[, "mz"] > 296.1490, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 237.0760 & rel_rep1[, "mz"] > 237.0756, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 296.1494 & rel_rep2[, "mz"] > 296.1490, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 237.0760 & rel_rep2[, "mz"] > 237.0756, c("RT", "mz", "mapping_rep1_rep2")]

## 4_[2_Hydroxy_2_(4_hydroxy_3_methoxyphenyl)_1_hydroxymethylethoxy]_3_methoxy_benzoylcholine_isomer1/2
## mz 450.21225 195.0652
rel_rep1[rel_rep1[, "mz"] < 450.2124 & rel_rep1[, "mz"] > 450.2120, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 195.0654 & rel_rep1[, "mz"] > 195.0650, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 450.2124 & rel_rep2[, "mz"] > 450.2120, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 195.0654 & rel_rep2[, "mz"] > 195.0650, c("RT", "mz", "mapping_rep1_rep2")]

## feruloylcholine ## mz 280.1543 221.0808
rel_rep1[rel_rep1[, "mz"] < 280.1545 & rel_rep1[, "mz"] > 280.1541, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 221.0810 & rel_rep1[, "mz"] > 221.0806, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 280.1545 & rel_rep2[, "mz"] > 280.1541, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 221.0810 & rel_rep2[, "mz"] > 221.0806, c("RT", "mz", "mapping_rep1_rep2")]

## sinapoylcholine_sinapoylcholine_dimer_isomer1/2 ## mz 310.16461 251.0914
rel_rep1[rel_rep1[, "mz"] < 310.1648 & rel_rep1[, "mz"] > 310.1644, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 251.0916 & rel_rep1[, "mz"] > 251.0912, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 310.1648 & rel_rep2[, "mz"] > 310.1644, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 251.0916 & rel_rep2[, "mz"] > 251.0912, c("RT", "mz", "mapping_rep1_rep2")]

## 2_phenylethylglucosinolate ## mz 462.0295 446.0555
rel_rep1[rel_rep1[, "mz"] < 462.0297 & rel_rep1[, "mz"] > 462.0293, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 446.0557 & rel_rep1[, "mz"] > 446.0553, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 462.0297 & rel_rep2[, "mz"] > 462.0293, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 446.0557 & rel_rep2[, "mz"] > 446.0553, c("RT", "mz", "mapping_rep1_rep2")]

## 1_O_Sinapoylglucose ## mz 409.1105 387.1286 207.0652
rel_rep1[rel_rep1[, "mz"] < 409.1107 & rel_rep1[, "mz"] > 409.1103, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 387.1288 & rel_rep1[, "mz"] > 387.1284, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 207.0654 & rel_rep1[, "mz"] > 207.0650, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 409.1107 & rel_rep2[, "mz"] > 409.1103, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 387.1288 & rel_rep2[, "mz"] > 387.1284, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 207.0654 & rel_rep2[, "mz"] > 207.0650, c("RT", "mz", "mapping_rep1_rep2")]

## N_((4'_O_glycosyl)_sinapoyl)_N'_sinapoylspermidine/N_Sinapoyl_N_(3,5_dimethoxy_4_hexosyloxycinnamoyl)spermidine_isomer2
## mz 720.3344 558.2815 264.123
rel_rep1[rel_rep1[, "mz"] < 720.3346 & rel_rep1[, "mz"] > 720.3342, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 558.2817 & rel_rep1[, "mz"] > 558.2813, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 264.1232 & rel_rep1[, "mz"] > 264.1228, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 720.3346 & rel_rep2[, "mz"] > 720.3342, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 558.2817 & rel_rep2[, "mz"] > 558.2813, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 264.1232 & rel_rep2[, "mz"] > 264.1228, c("RT", "mz", "mapping_rep1_rep2")]

## isorhamnetin_DeoxyHex_Hex ## mz 625.17631 463.1235 317.0656
rel_rep1[rel_rep1[, "mz"] < 625.1765 & rel_rep1[, "mz"] > 625.1763, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 463.1237 & rel_rep1[, "mz"] > 463.1233, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 317.0658 & rel_rep1[, "mz"] > 317.0654, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 625.1765 & rel_rep2[, "mz"] > 625.1763, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 463.1237 & rel_rep2[, "mz"] > 463.1233, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 317.0658 & rel_rep2[, "mz"] > 317.0654, c("RT", "mz", "mapping_rep1_rep2")]

## N1_N8_di(sinapoyl)_spermidine/N_N_Bis(sinapoyl)spermidine ## mz 558.2815 264.123
rel_rep1[rel_rep1[, "mz"] < 558.2817 & rel_rep1[, "mz"] > 558.2813, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 264.1232 & rel_rep1[, "mz"] > 264.1228, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 558.2817 & rel_rep2[, "mz"] > 558.2813, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 264.1232 & rel_rep2[, "mz"] > 264.1228, c("RT", "mz", "mapping_rep1_rep2")]

## 2_O_Sinapoylmalate ## mz 719.122 703.1481 379.0426 363.0687 175.039 147.0441
rel_rep1[rel_rep1[, "mz"] < 719.1222 & rel_rep1[, "mz"] > 719.1218, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 703.1483 & rel_rep1[, "mz"] > 703.1479, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 379.0428 & rel_rep1[, "mz"] > 379.0424, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 363.0689 & rel_rep1[, "mz"] > 363.0685, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 175.0392 & rel_rep1[, "mz"] > 175.0388, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 147.0443 & rel_rep1[, "mz"] > 147.0439, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 719.1222 & rel_rep2[, "mz"] > 719.1218, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 703.1483 & rel_rep2[, "mz"] > 703.1479, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 379.0428 & rel_rep2[, "mz"] > 379.0424, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 363.0689 & rel_rep2[, "mz"] > 363.0685, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 175.0392 & rel_rep2[, "mz"] > 175.0388, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 147.0443 & rel_rep2[, "mz"] > 147.0439, c("RT", "mz", "mapping_rep1_rep2")]

## isorhamnetin_(DeoxyHex)2 ## mz 609.1814 463.1235 317.068
rel_rep1[rel_rep1[, "mz"] < 609.1816 & rel_rep1[, "mz"] > 609.1812, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 463.1237 & rel_rep1[, "mz"] > 463.1233, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 317.0682 & rel_rep1[, "mz"] > 317.0678, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 609.1816 & rel_rep2[, "mz"] > 609.1812, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 463.1237 & rel_rep2[, "mz"] > 463.1233, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 317.0682 & rel_rep2[, "mz"] > 317.0678, c("RT", "mz", "mapping_rep1_rep2")]

## 3_[2_(4_Hydroxy_3_methoxyphenyl)_3_hydroxymethyl_7_methoxy_2_3_dihydrobenzo_furan_5_yl]acryloylcholine
## mz 458.2173 399.1438
rel_rep1[rel_rep1[, "mz"] < 458.2175 & rel_rep1[, "mz"] > 458.2171, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 399.1440 & rel_rep1[, "mz"] > 399.1436, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 458.2175 & rel_rep2[, "mz"] > 458.2171, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 399.1440 & rel_rep2[, "mz"] > 399.1436, c("RT", "mz", "mapping_rep1_rep2")]

## quercetin_deoxyhexose ## mz 449.10784 303.052
rel_rep1[rel_rep1[, "mz"] < 449.1090 & rel_rep1[, "mz"] > 449.1066, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 303.0522 & rel_rep1[, "mz"] > 303.0518, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 449.1090 & rel_rep2[, "mz"] > 449.1066, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 303.0532 & rel_rep2[, "mz"] > 303.0508, c("RT", "mz", "mapping_rep1_rep2")]

## isorhamnetin_hexose ## mz 501.1009 479.119 317.0656
rel_rep1[rel_rep1[, "mz"] < 501.1011 & rel_rep1[, "mz"] > 501.1007, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 479.1192 & rel_rep1[, "mz"] > 479.1188, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 317.0658 & rel_rep1[, "mz"] > 317.0654, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 501.1011 & rel_rep2[, "mz"] > 501.1007, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 479.1192 & rel_rep2[, "mz"] > 479.1188, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 317.0658 & rel_rep2[, "mz"] > 317.0654, c("RT", "mz", "mapping_rep1_rep2")]

## kaempferol_deoxyhex ## mz 433.1129 287.055
rel_rep1[rel_rep1[, "mz"] < 433.1131 & rel_rep1[, "mz"] > 433.1127, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 287.0552 & rel_rep1[, "mz"] > 287.0548, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 433.1131 & rel_rep2[, "mz"] > 433.1127, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 287.0552 & rel_rep2[, "mz"] > 287.0548, c("RT", "mz", "mapping_rep1_rep2")]

## isorhamnetin_deoxyhex ## mz 463.12349 317.0656
rel_rep1[rel_rep1[, "mz"] < 463.1236 & rel_rep1[, "mz"] > 463.1232, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep1[rel_rep1[, "mz"] < 317.0658 & rel_rep1[, "mz"] > 317.0654, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 463.1236 & rel_rep2[, "mz"] > 463.1232, c("RT", "mz", "mapping_rep1_rep2")]
rel_rep2[rel_rep2[, "mz"] < 317.0658 & rel_rep2[, "mz"] > 317.0654, c("RT", "mz", "mapping_rep1_rep2")]

