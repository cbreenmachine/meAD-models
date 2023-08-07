


clean_master_samplesheet <- function(ifile){
    # 
    read.csv(args$master_file) %>%
        dplyr::arrange(sample_id) %>% 
        dplyr::mutate(sample_id = as.character(sample_id)) %>%
        dplyr::filter(diagnostic_group != "MCI") %>%
        dplyr::select(-starts_with("PC")) %>%
        dplyr::mutate(
            diagnostic_group_coded = 
            ifelse(diagnostic_group == "CONTROL", 0, 1))
}

get_valid_ids <- function(master_df, bed_dir){
    # Check that the samples in the directory have phenotype and vice versa
    idir_samples <- 
        unique(stringr::str_split_fixed(list.files(args$idir), "\\.", 3)[ ,2])

    # All samples that are in the directory and sampleshseet
    valid_samples <- intersect(master_df$sample_id, idir_samples)
    valid_samples
}


get_load_samples <- function(master_df, valid_samples){
    load_samples <- intersect(
        master_df$sample_id[master_df$diagnostic_group == "LOAD"], 
        valid_samples)
    return(load_samples)
}

get_control_samples <- function(master_df, valid_samples){
    ctrl_samples <- intersect(
        master_df$sample_id[master_df$diagnostic_group == "CONTROL"], 
        idir_samples
    )
    return(ctrl_samples)
}
