
process_artificial_data <- function(filename, sheetname) {
  
  df <- read_xlsx(filename,
                  sheet=sheetname) %>% 
    rename_all(tolower)
  names(df) <- gsub("-", "", names(df))
  names(df) <- gsub(" ", "_", names(df))
  names(df) <- gsub("__", "_", names(df))
  df <- df %>% 
    rename(
      hpb_normalised=normalized_data_hpb_negative_control,
      sherp_normalised=normalized_data_sherp_negative_control
    ) %>% 
    mutate(
      hpb_ct_mean=if_else(is.na(hpb_ct_mean), 40, hpb_ct_mean),
      sherp_ct_mean=if_else(is.na(sherp_ct_mean), 40, sherp_ct_mean)
    ) %>% 
    mutate(
      hpb_ct_mean=if_else(hpb_ct_mean > 40, 40, hpb_ct_mean),
      sherp_ct_mean=if_else(sherp_ct_mean > 40, 40, sherp_ct_mean)
    ) %>% 
    mutate(
      hpb_normalised=40-hpb_ct_mean,
      sherp_normalised=40-sherp_ct_mean
    ) %>% 
    mutate(bloodmeal_number=substr(blood_meal_event, 3, 4)) %>% 
    mutate(bloodmeal_number=as.numeric(bloodmeal_number)) %>% 
    mutate(timepoint=as.character(timepoint)) %>% 
    mutate(bloodmeal_type=if_else(bloodmeal_number>1, "multiple", "single"))
 df 
}
