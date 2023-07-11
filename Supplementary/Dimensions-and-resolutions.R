#Thresholds,Dimensions and Resolutions

#QC cutoffs
datasetsUpper<-c("Mouse_T3_Mixed"= 6000,
                 "Mouse_A1_Mixed" = 5000,
                 "Mouse_T2_GFP" = 3000,
                 "Mouse_T1_dsRed" = 4000,
                 "Mouse_T1_GFP" = 5000,
                 "Mouse_U1_Unmarked" = 6500,
                 "Mouse_U2_Unmarked" = 6000,
                 "Mouse_U1_dsRed" = 7000,
                 "Mouse_U2_dsRed" = 6500,
                 "Mouse_A2_GFP" = 4500,
                 "Mouse_A2_dsRed" = 7000,
                 "Mouse_A3_GFP" = 5000,
                 "Mouse_A3_dsRed1" = 7500,
                 "Mouse_A3_dsRed2" = 7500)

datasetsHumanMito<-c("Mouse_T3_Mixed"= 15,
                     "Mouse_A1_Mixed" = 10,
                     "Mouse_T2_GFP" = 5,
                     "Mouse_T1_dsRed" = 5,
                     "Mouse_T1_GFP" = 15,
                     "Mouse_U1_Unmarked" = 15,
                     "Mouse_U2_Unmarked" = 15,
                     "Mouse_U1_dsRed" = 15,
                     "Mouse_U2_dsRed" = 15,
                     "Mouse_A2_GFP" = 5,
                     "Mouse_A2_dsRed" = 15,
                     "Mouse_A3_GFP" = 5,
                     "Mouse_A3_dsRed1" = 15,
                     "Mouse_A3_dsRed2" = 15)





#Pre-Subsetting using canonical cellmarkers
datasetsUMAP<-c("Mouse_T3_Mixed"= 30,
                "Mouse_A1_Mixed" = 30,
                "Mouse_T2_GFP" = 30,
                "Mouse_T1_dsRed" = 30,
                "Mouse_T1_GFP" = 30,
                "Mouse_U1_Unmarked" = 30,
                "Mouse_U2_Unmarked" = 30,
                "Mouse_U1_dsRed" = 30,
                "Mouse_U2_dsRed" = 30,
                "Mouse_A2_GFP" = 30,
                "Mouse_A2_dsRed" = 30,
                "Mouse_A3_GFP" = 30,
                "Mouse_A3_dsRed1" = 30,
                "Mouse_A3_dsRed2" = 30)

datasetsClusterRes<-c("Mouse_T3_Mixed"= 0.4,
                      "Mouse_A1_Mixed" = 0.7,
                      "Mouse_T2_GFP" = 0.6,
                      "Mouse_T1_dsRed" = 0.6,
                      "Mouse_T1_GFP" = 0.6,
                      "Mouse_U1_Unmarked" = 0.6,
                      "Mouse_U2_Unmarked" = 0.6,
                      "Mouse_U1_dsRed" = 0.4,
                      "Mouse_U2_dsRed" = 0.4,
                      "Mouse_A2_GFP" = 0.5,
                      "Mouse_A2_dsRed" = 0.5,
                      "Mouse_A3_GFP" = 0.5,
                      "Mouse_A3_dsRed1" = 0.5,
                      "Mouse_A3_dsRed2" = 0.5)




#Final
datasetsUMAP2<-c("Mouse_T3_Mixed"= 20,
                 "Mouse_A1_Mixed" = 20,
                 "Mouse_T2_GFP" = 20,
                 "Mouse_T1_dsRed" = 20,
                 "Mouse_T1_GFP" = 20,
                 "Mouse_U1_Unmarked" = 25,
                 "Mouse_U2_Unmarked" = 25,
                 "Mouse_U1_dsRed" = 30,
                 "Mouse_U2_dsRed" = 30,
                 "Mouse_A2_GFP" = 30,
                 "Mouse_A2_dsRed" = 30,
                 "Mouse_A3_GFP" = 30,
                 "Mouse_A3_dsRed1" = 30,
                 "Mouse_A3_dsRed2" = 30)

datasetsClusterRes2<-c("Mouse_T3_Mixed"= 0.5,
                       "Mouse_A1_Mixed" = 0.4,
                       "Mouse_T2_GFP" = 0.3,
                       "Mouse_T1_dsRed" = 0.3,
                       "Mouse_T1_GFP" = 0.3,
                       "Mouse_U1_Unmarked" = 0.5,
                       "Mouse_U2_Unmarked" = 0.5,
                       "Mouse_U1_dsRed" = 0.4,
                       "Mouse_U2_dsRed" = 0.4,
                       "Mouse_A2_GFP" = 0.4,
                       "Mouse_A2_dsRed" = 0.4,
                       "Mouse_A3_GFP" = 0.4,
                       "Mouse_A3_dsRed1" = 0.4,
                       "Mouse_A3_dsRed2" = 0.4)