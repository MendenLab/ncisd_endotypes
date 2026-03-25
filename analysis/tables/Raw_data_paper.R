library(DESeq2)
library(magrittr)
library(data.table)


# Select and rename columns
get_metadata <- function(df) {
  df_metadata <- df[, c("healthysamp_Pattern_v2", "healthysamp_diag_v2", 
                        "PatientID", "age", "Sex.x")]
  
  # Rename columns
  colnames(df_metadata) <- c("Pattern", "Disease", "Patient", "Age", "Sex")
  
  return(df_metadata)
}

get_infos <- function(df)
{
  # Sample: assuming your dataframe is named df and has columns:
  # Patient_ID, Pattern, diag, Age, Sex
  
  dt <- as.data.table(df)
  res <- dt[, .(
    Patients_n = uniqueN(Patient),
    Age_mean_sd = sprintf("%.1f ± %.1f", mean(Age, na.rm = TRUE), sd(Age, na.rm = TRUE)),
    Male_percent = round(mean(tolower(Sex) == "m", na.rm = TRUE) * 100, 1)
  ), by = .(Pattern, Disease)]
  
  # res <- df %>%
  #   group_by(Pattern, Disease) %>%
  #   summarise(
  #     `Patients (n)` = n_distinct(Patient),  # Count unique patients
  #     `Age (mean ± SD)` = sprintf("%.1f ± %.1f", mean(Age, na.rm = TRUE), sd(Age, na.rm = TRUE)),
  #     `Sex (% of male)` = round(mean(tolower(Sex) == "m", na.rm = TRUE) * 100, 1),
  #     .groups = "drop"
  #   )
  
  # Optional: Rename for final display
  setnames(res, c("Patients_n", "Age_mean_sd", "Male_percent"), 
           c("Patients (n)", "Age (mean ± SD)", "Sex (% of male)"))
  
  # Print the result
  print(res)
  
  return(res)
}

rename_nonlesion_samples_to_paired_lesion_name <- function(
    df, ref_obsname = "diag", obsname = "healthysamp_diag_v2", 
    group_name = "diag") 
{
  df[[ref_obsname]] <- as.character(df[[ref_obsname]])
  df[[obsname]] <- "Unknown"
  
  # Get unique diagnostic categories (prioritizing "non-lesional" first)
  categories <- unique(df[[ref_obsname]])
  categories <- c("non-lesional", setdiff(categories, "non-lesional"))
  
  for (diagnosis in categories) {
    # Subset patients with the current diagnosis
    patients_with_diag <- df[df[[ref_obsname]] == diagnosis, "PatientID"]
    
    # Match all samples from those patients where diag == diagnosis OR diag == non-lesional
    matched_rows <- df$PatientID %in% patients_with_diag & 
      (df[[group_name]] == diagnosis | df[[group_name]] == "non-lesional")
    
    # Assign diagnosis to these matched rows
    df[matched_rows, obsname] <- diagnosis
  }
  
  return(df)
}


# Load raw data
dds <- readRDS('/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/raw_data/count_matrix/data_freeze/raw/dds_v04.rds')

# Dimension before filtering
dim(dds) # genes: 18341   samples:727

length(unique(dds$sampleID))
length(unique(dds$PatientID))
length(unique(dds$diag))

df <- colData(dds)

index = c('MUC10222','MUC10223','MUC10224','MUC10225','MUC10226', 'MUC10227','MUC10228',
           'MUC10229','MUC10230','MUC10231','MUC10232','MUC10233','MUC10234', 'MUC10235',
           'MUC10236','MUC10237','MUC10238','MUC10239','MUC10240','MUC10242','MUC20403',
           'MUC20404','MUC20406','MUC20408',
           'MUC20409',
           'MUC20410',
           'MUC20411',
           'MUC20412',
           'MUC20413',
           'MUC20414',
           'MUC20415',
           'MUC20416',
           'MUC20417',
           'MUC20418',
           'MUC20419',
           'MUC20420',
           'MUC20421',
           'MUC20425',
           'MUC20426',
           'MUC20429',
           'MUC20430',
           'MUC20431',
           'MUC20432',
           'MUC20433',
           'MUC20434',
           'MUC20435',
           'MUC20436',
           'MUC20437',
           'MUC20438',
           'MUC20439',
           'MUC20440',
           'MUC20441',
           'MUC20442',
           'MUC20443',
           'MUC20444',
           'MUC20445',
           'MUC20446',
           'MUC20447',
           'MUC20448',
           'MUC20449',
           'MUC20452',
           'MUC20454',
           'MUC20455',
           'MUC20456',
           'MUC20457',
           'MUC20458',
           'MUC20459',
           'MUC20460',
           'MUC20461',
           'MUC20462',
           'MUC20463',
           'MUC20464',
           'MUC20465',
           'MUC20466',
           'MUC20467',
           'MUC20468',
           'MUC20469',
           'MUC20470',
           'MUC20471',
           'MUC20472',
           'MUC20473',
           'MUC20474',
           'MUC20475',
           'MUC20478',
           'MUC20479',
           'MUC20480',
           'MUC20481',
           'MUC20482',
           'MUC20483',
           'MUC20484',
           'MUC20485',
           'MUC20495',
           'MUC20499',
           'MUC20501',
           'MUC20503',
           'MUC20505',
           'MUC20508',
           'MUC20512',
           'MUC20516',
           'MUC20518',
           'MUC20520',
           'MUC20522',
           'MUC20524',
           'MUC20526',
           'MUC20528',
           'MUC20530',
           'MUC20709',
           'MUC20721',
           'MUC20722',
           'MUC20723',
           'MUC2674',
           'MUC2675',
           'MUC2676',
           'MUC2677',
           'MUC2678',
           'MUC2679',
           'MUC2680',
           'MUC2681',
           'MUC2682',
           'MUC2683',
           'MUC2684',
           'MUC2685',
           'MUC2686',
           'MUC2687',
           'MUC2688',
           'MUC2689',
           'MUC2690',
           'MUC2691',
           'MUC2692',
           'MUC2693',
           'MUC2694',
           'MUC2695',
           'MUC2696',
           'MUC2697',
           'MUC2698',
           'MUC2699',
           'MUC2700',
           'MUC2701',
           'MUC2702',
           'MUC2703',
           'MUC2704',
           'MUC2705',
           'MUC2706',
           'MUC2707',
           'MUC2708',
           'MUC2709',
           'MUC2710',
           'MUC2711',
           'MUC2712',
           'MUC2713',
           'MUC2714',
           'MUC2715',
           'MUC2716',
           'MUC2717',
           'MUC2718',
           'MUC2719',
           'MUC2720',
           'MUC2721',
           'MUC2722',
           'MUC2723',
           'MUC2724',
           'MUC2725',
           'MUC2726',
           'MUC2727',
           'MUC2728',
           'MUC2729',
           'MUC2730',
           'MUC2731',
           'MUC2732',
           'MUC2733',
           'MUC2734',
           'MUC2735',
           'MUC2736',
           'MUC2737',
           'MUC2738',
           'MUC2739',
           'MUC2740',
           'MUC2741',
           'MUC2742',
           'MUC2743',
           'MUC2744',
           'MUC2745',
           'MUC2746',
           'MUC2747',
           'MUC2748',
           'MUC2749',
           'MUC2750',
           'MUC2751',
           'MUC2752',
           'MUC2753',
           'MUC2754',
           'MUC2755',
           'MUC2756',
           'MUC2757',
           'MUC2758',
           'MUC2759',
           'MUC2760',
           'MUC2761',
           'MUC2762',
           'MUC2763',
           'MUC2764',
           'MUC2765',
           'MUC2766',
           'MUC2767',
           'MUC2768',
           'MUC2769',
           'MUC2770',
           'MUC2771',
           'MUC2772',
           'MUC2773',
           'MUC2774',
           'MUC2775',
           'MUC2776',
           'MUC2777',
           'MUC2778',
           'MUC2779',
           'MUC2780',
           'MUC2781',
           'MUC2782',
           'MUC2783',
           'MUC2784',
           'MUC2785',
           'MUC2786',
           'MUC2787',
           'MUC2788',
           'MUC2789',
           'MUC2790',
           'MUC2791',
           'MUC2792',
           'MUC2793',
           'MUC2794',
           'MUC2795',
           'MUC2796',
           'MUC2797',
           'MUC2798',
           'MUC2799',
           'MUC2800',
           'MUC2801',
           'MUC2802',
           'MUC2803',
           'MUC2804',
           'MUC2805',
           'MUC2806',
           'MUC2807',
           'MUC2808',
           'MUC2809',
           'MUC2810',
           'MUC2811',
           'MUC2812',
           'MUC2813',
           'MUC2814',
           'MUC2815',
           'MUC2816',
           'MUC2817',
           'MUC2818',
           'MUC2819',
           'MUC2820',
           'MUC2821',
           'MUC2822',
           'MUC2823',
           'MUC2824',
           'MUC2825',
           'MUC2826',
           'MUC2827',
           'MUC2828',
           'MUC2829',
           'MUC2830',
           'MUC2831',
           'MUC2832',
           'MUC2833',
           'MUC2834',
           'MUC2835',
           'MUC2836',
           'MUC2837',
           'MUC2838',
           'MUC2839',
           'MUC2840',
           'MUC2841',
           'MUC2842',
           'MUC2843',
           'MUC2844',
           'MUC2845',
           'MUC2846',
           'MUC2847',
           'MUC2848',
           'MUC2849',
           'MUC2850',
           'MUC2851',
           'MUC2852',
           'MUC2853',
           'MUC2854',
           'MUC2855',
           'MUC2856',
           'MUC2857',
           'MUC2858',
           'MUC2859',
           'MUC2860',
           'MUC2861',
           'MUC2862',
           'MUC2863',
           'MUC2864',
           'MUC2865',
           'MUC2866',
           'MUC2868',
           'MUC2869',
           'MUC2870',
           'MUC2871',
           'MUC2872',
           'MUC2873',
           'MUC2874',
           'MUC2875',
           'MUC2876',
           'MUC2877',
           'MUC2878',
           'MUC2879',
           'MUC2880',
           'MUC2881',
           'MUC2882',
           'MUC2883',
           'MUC2884',
           'MUC2885',
           'MUC2886',
           'MUC2887',
           'MUC2888',
           'MUC2889',
           'MUC2890',
           'MUC2891',
           'MUC2892',
           'MUC2893',
           'MUC2894',
           'MUC2895',
           'MUC2897',
           'MUC2898',
           'MUC2899',
           'MUC2900',
           'MUC2902',
           'MUC2903',
           'MUC2905',
           'MUC2906',
           'MUC2907',
           'MUC2908',
           'MUC2910',
           'MUC2911',
           'MUC2913',
           'MUC2914',
           'MUC2916',
           'MUC2917',
           'MUC2921',
           'MUC2922',
           'MUC2923',
           'MUC2924',
           'MUC2925',
           'MUC2926',
           'MUC2927',
           'MUC2928',
           'MUC2929',
           'MUC2930',
           'MUC2931',
           'MUC2932',
           'MUC2933',
           'MUC2934',
           'MUC2935',
           'MUC2936',
           'MUC2937',
           'MUC2938',
           'MUC2939',
           'MUC2940',
           'MUC2941',
           'MUC2942',
           'MUC2943',
           'MUC2944',
           'MUC2945',
           'MUC2946',
           'MUC2949',
           'MUC2950',
           'MUC2951',
           'MUC2952',
           'MUC2953',
           'MUC2954',
           'MUC2955',
           'MUC2956',
           'MUC2957',
           'MUC2958',
           'MUC2959',
           'MUC2960',
           'MUC2961',
           'MUC3667',
           'MUC3668',
           'MUC3669',
           'MUC3670',
           'MUC3671',
           'MUC3672',
           'MUC3673',
           'MUC3674',
           'MUC3675',
           'MUC3676',
           'MUC3677',
           'MUC3678',
           'MUC3679',
           'MUC3680',
           'MUC3681',
           'MUC3682',
           'MUC3683',
           'MUC3684',
           'MUC3685',
           'MUC3686',
           'MUC3687',
           'MUC3688',
           'MUC3689',
           'MUC3690',
           'MUC3691',
           'MUC3692',
           'MUC3693',
           'MUC3694',
           'MUC3695',
           'MUC3696',
           'MUC3697',
           'MUC3698',
           'MUC3699',
           'MUC3700',
           'MUC3701',
           'MUC3702',
           'MUC3703',
           'MUC3704',
           'MUC3705',
           'MUC3706',
           'MUC3707',
           'MUC3708',
           'MUC3709',
           'MUC3710',
           'MUC3711',
           'MUC3712',
           'MUC3713',
           'MUC3714',
           'MUC3715',
           'MUC3716',
           'MUC3717',
           'MUC3718',
           'MUC3719',
           'MUC3720',
           'MUC3721',
           'MUC3722',
           'MUC3723',
           'MUC3724',
           'MUC3725',
           'MUC3726',
           'MUC3727',
           'MUC3728',
           'MUC3729',
           'MUC3730',
           'MUC3731',
           'MUC3732',
           'MUC3733',
           'MUC3734',
           'MUC3735',
           'MUC3737',
           'MUC3738',
           'MUC3739',
           'MUC3740',
           'MUC3741',
           'MUC3742',
           'MUC3743',
           'MUC3744',
           'MUC3745',
           'MUC3746',
           'MUC3747',
           'MUC3748',
           'MUC3749',
           'MUC3750',
           'MUC3751',
           'MUC3752',
           'MUC3753',
           'MUC3754',
           'MUC3755',
           'MUC3756',
           'MUC3757',
           'MUC3758',
           'MUC3759',
           'MUC3760',
           'MUC3761',
           'MUC3762',
           'MUC4241',
           'MUC4242',
           'MUC4243',
           'MUC4244',
           'MUC4245',
           'MUC4246',
           'MUC4247',
           'MUC4248',
           'MUC4249',
           'MUC4250',
           'MUC4251',
           'MUC4252',
           'MUC4253',
           'MUC4254',
           'MUC4255',
           'MUC4256',
           'MUC4257',
           'MUC4258',
           'MUC4259',
           'MUC4260',
           'MUC4261',
           'MUC4262',
           'MUC4263',
           'MUC4264',
           'MUC4265',
           'MUC4266',
           'MUC4267',
           'MUC4268',
           'MUC4269',
           'MUC4270',
           'MUC4271',
           'MUC4272',
           'MUC4273',
           'MUC4274',
           'MUC4275',
           'MUC4276',
           'MUC4277',
           'MUC4278',
           'MUC4279',
           'MUC4280',
           'MUC4281',
           'MUC4282',
           'MUC4283',
           'MUC4284',
           'MUC4285',
           'MUC4286',
           'MUC4287',
           'MUC4288',
           'MUC4289',
           'MUC4290',
           'MUC4291',
           'MUC4292',
           'MUC4293',
           'MUC4294',
           'MUC4295',
           'MUC4296',
           'MUC4297',
           'MUC4298',
           'MUC4299',
           'MUC4300',
           'MUC4301',
           'MUC4302',
           'MUC4303',
           'MUC4304',
           'MUC4305',
           'MUC4306',
           'MUC4307',
           'MUC4308',
           'MUC4309',
           'MUC4310',
           'MUC4311',
           'MUC4312',
           'MUC4313',
           'MUC4314',
           'MUC4315',
           'MUC4316',
           'MUC4317',
           'MUC4318',
           'MUC4319',
           'MUC4320',
           'MUC4321',
           'MUC4322',
           'MUC4323',
           'MUC4324',
           'MUC4325',
           'MUC4326',
           'MUC4327',
           'MUC4328',
           'MUC4329',
           'MUC4330',
           'MUC4331',
           'MUC4332',
           'MUC4333',
           'MUC4334',
           'MUC4335',
           'MUC4336',
           'MUC7233',
           'MUC7234',
           'MUC7235',
           'MUC7236',
           'MUC7237',
           'MUC7238',
           'MUC7239',
           'MUC7240',
           'MUC7241',
           'MUC7242',
           'MUC7243',
           'MUC7244',
           'MUC7245',
           'MUC7246',
           'MUC7247',
           'MUC7248',
           'MUC7249',
           'MUC7250',
           'MUC7251',
           'MUC7252',
           'MUC7253',
           'MUC7254',
           'MUC7255',
           'MUC7256',
           'MUC7257',
           'MUC7258',
           'MUC7259',
           'MUC7260',
           'MUC7261',
           'MUC7262',
           'MUC7263',
           'MUC7264',
           'MUC7265',
           'MUC7266',
           'MUC7267',
           'MUC7268',
           'MUC7269',
           'MUC7270',
           'MUC7271',
           'MUC7272',
           'MUC7273',
           'MUC7274',
           'MUC7275',
           'MUC7276',
           'MUC7277',
           'MUC7278',
           'MUC7279',
           'MUC7280',
           'MUC7281',
           'MUC7282',
           'MUC7283',
           'MUC7284',
           'MUC7285',
           'MUC7286',
           'MUC7287',
           'MUC7288',
           'MUC7289',
           'MUC7290',
           'MUC7291',
           'MUC7292',
           'MUC7293',
           'MUC7294',
           'MUC7297',
           'MUC7298',
           'MUC7299',
           'MUC7300',
           'MUC7301',
           'MUC7302',
           'MUC7303',
           'MUC7304',
           'MUC7305',
           'MUC7306',
           'MUC7307',
           'MUC7308',
           'MUC7309',
           'MUC7310',
           'MUC7311',
           'MUC7312',
           'MUC7314',
           'MUC7315',
           'MUC7317',
           'MUC7318',
           'MUC7320',
           'MUC7321',
           'MUC7322',
           'MUC7324',
           'MUC7325',
           'MUC7326',
           'MUC7329',
           'MUC7330',
           'MUC7331',
           'MUC7332',
           'MUC7333',
           'MUC7334',
           'MUC7335',
           'MUC7336',
           'MUC7337',
           'MUC7338',
           'MUC7339',
           'MUC7340',
           'MUC7341',
           'MUC7342',
           'MUC7343',
           'MUC7344',
           'MUC7345',
           'MUC7346',
           'MUC7347',
           'MUC7348',
           'MUC7349',
           'MUC7350',
           'MUC7351','MUC7352','MUC7353','MUC7354','MUC7355','MUC7356','MUC7359',
           'MUC7360','MUC7361','MUC7362','MUC7363','MUC7364','MUC7365','MUC7366',
           'MUC7367','MUC7368','MUC7369',
           'MUC7370','MUC7371','MUC7372','MUC7375','MUC7376','MUC7377','MUC7378',
           'MUC7379','MUC7380','MUC7381','MUC7382',
           'MUC7383','MUC7384','MUC7385','MUC7386','MUC7387','MUC7388',
           'MUC7390','MUC7391','MUC7392','MUC7393','MUC7394','MUC7395')

# Reorder index by list
all(rownames(df[index, ]) == index)
df = df[index, ]

# Save MetaData of Raw Count data
# xlsx::write.xlsx(
#   file = '/Volumes/CH__data/Projects/data/Summary/Raw_MetaData_Overview.xlsx', 
#   x = df[c('PatientID', 'diag', 'Pattern', 'sampleType', 'age', 'batch', 'Sex.x')],
#   sheetName = 'BRAIN')
# PatientID sometimes contains strings like HG16.711_H and HG16.711_Lym 
# -> remove everything rubstring after "_" to make same patient
df$PatientID <- gsub("\\_.*","", df$PatientID)


crosstab <- table(df[c('PatientID', 'sampleType')])

# Total 390 Patients with 727 samples before preprocessing
# 1 Patient with 1 L and 2 NL samples (2, 1) = 3 samples
crosstab_LNL_more_NL_than_L <- crosstab[(crosstab[, 1] == 2 & crosstab[, 2] == 1), ]
# 4 Patients with 2 L and 1 NL samples (1, 2) = 12 samples
crosstab_LNL_more_L_than_NL <- dim(crosstab[(crosstab[, 1] == 1 & crosstab[, 2] == 2), ])
# 326 Patients with 1 L and 1 NL sample (1, 1) = 652 samples
crosstab_LNL_11 <- dim(crosstab[(crosstab[, 1] == 1 & crosstab[, 2] == 1), ])
# 58 Patients with L sample only (0, 1[57] | 2[1]) = 59 samples
crosstab_only_L <- dim(crosstab[crosstab[, 1] == 0 & crosstab[, 2] > 0, ])
# 1 Patient (1, 0) = 1 sample
crosstab_only_H <- crosstab[crosstab[, 1] == 1  & crosstab[, 2] == 0, ]


# Load latest version of metaData file
# '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/raw_data/clinical_data/20240730_patient_meta_data_final_PS.xlsx'
filename = '/Users/christina.hillig/Downloads/patient_meta_data_final (1).xlsx'
df_final_metadata <- xlsx::read.xlsx(filename, sheetIndex = 1)

# Change one patient's disease and pattern
mask_patient <- (df_final_metadata$Pseudo.ID == '10350') & (df_final_metadata$Helmholtz.identifyer == 'MUC4259')
# df_final_metadata[mask_patient, 'Pattern'] <- "1"
# df_final_metadata[mask_patient, 'diag'] <- "lichenoid drug reaction"

# PatientID sometimes contains strings like HG16.711_H and HG16.711_Lym 
# -> remove everything rubstring after "_" to make same patient
df_final_metadata$PatientID <- gsub("\\_.*","", df_final_metadata$Pseudo.ID)



# Apply function
df_final_metadata <- rename_nonlesion_samples_to_paired_lesion_name(
  df_final_metadata, ref_obsname = "diag", obsname = "healthysamp_diag_v2", 
  group_name = "diag")
df_final_metadata$Pattern[is.na(df_final_metadata$Pattern)] <- "non-lesional"
# df_final_metadata$Pattern <- gsub("0", "non-lesional", df_final_metadata$Pattern)
df_final_metadata <- rename_nonlesion_samples_to_paired_lesion_name(
  df_final_metadata, ref_obsname = "Pattern", obsname = "healthysamp_Pattern_v2", 
  group_name = "Pattern")

# Remove duplicate patients (keep last occurrence) 
df_unique <- df_final_metadata[!duplicated(df_final_metadata$PatientID, fromLast = TRUE), ]
# Get sub metadata
df_metadata <- get_metadata(df_unique)
# Read out Table S1 
df_res <- get_infos(df=df_metadata)
# Reorder by Pattern
df_res <- df_res[order(df_res$Pattern, decreasing = FALSE), ]

# Save to excel file
xlsx::write.xlsx(
  file = '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/Paper/20250526_Second_Draft/Table_S1_before_preprocessing.xlsx',
  x = df_res,
  sheetName = 'Table_S1_LSNL_samples_final_MetaData')


crosstab <- table(df_final_metadata[c('PatientID', 'sampleType')])

# Total 391 Patients with 729 samples before processing .. 2 more than were actually sequenced
# 1 Patient with 1 L and 2 NL samples (2, 1) = 3 samples
crosstab_LNL_more_NL_than_L <- crosstab[(crosstab[, 1] == 1 & crosstab[, 2] == 2), ]
# 4 Patients with 2 L and 1 NL samples (1, 2) = 12 samples
crosstab_LNL_more_L_than_NL <- crosstab[(crosstab[, 1] == 2 & crosstab[, 2] == 1), ]
# 327 Patients with 1 L and 1 NL sample (1, 1) = 654 samples
crosstab_LNL_11 <- dim(crosstab[(crosstab[, 1] == 1 & crosstab[, 2] == 1), ])
# 58 Patients with L sample only (0, 1) = 59 samples
crosstab_only_L <- dim(crosstab[crosstab[, 1] >= 1  & crosstab[, 2] == 0, ])
# 1 Patient with NL samples only (1, 0) = 1 sample
crosstab_only_H <- crosstab[crosstab[, 1] == 0 & crosstab[, 2] == 1, ]

# Welcher patient wurde nicht sequenziert? 
missing_patient_id <- setdiff(df_final_metadata$PatientID, df$PatientID)
# Welche sample fehlt?
missing_samples <- setdiff(df_final_metadata$Helmholtz.identifyer, df$sampleID)
# drop patient and read out Table S1 again
df_final_metadata_seq <- df_final_metadata[!(
  df_final_metadata$Helmholtz.identifyer %in% c("MUC20424", "MUC2867")), ]

# Remove duplicate patients (keep last occurrence) 
df_unique_seq <- df_final_metadata_seq[!duplicated(
  df_final_metadata_seq$PatientID, fromLast = TRUE), ]
# Get sub metadata
df_metadata_seq <- get_metadata(df_unique_seq)
# Read out Table S1 
df_res_seq <- get_infos(df=df_metadata_seq)
# Reorder by Pattern
df_res_seq <- df_res_seq[order(df_res_seq$Pattern, decreasing = FALSE), ]

library(openxlsx)
# Load the existing Excel file
wb <- loadWorkbook('/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/Paper/20250526_Second_Draft/Table_S1_before_preprocessing.xlsx')
# Add a new worksheet with your data frame
addWorksheet(wb, "Table_S1_LSNL_samples_sequenced")  # change to desired sheet name
writeData(wb, sheet = "Table_S1_LSNL_samples_sequenced", x = df_res_seq)
# Save back to the same file
saveWorkbook(wb, '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/Paper/20250526_Second_Draft/Table_S1_before_preprocessing.xlsx', 
             overwrite = TRUE)




