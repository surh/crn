
read_single_experiment <- function(od600_file, timepoints_file, syncoms_file = NULL,
                                   mpn_file = NULL, type = "syncom",
                                   batch_name =  basename(dirname(od600_file))){
  
  # od600_file <- "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/od600.tsv"
  # timepoints_file <- "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/timepoints.tsv"
  # syncoms_file <- "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/syncoms.tsv"
  # mpn_file <- "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/mpn.tsv"
  # type <- "syncom"
  # batch_name <- "NS1"
  
  
  # od600_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/single_strains/NS1/od600.tsv"
  # timepoints_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/single_strains/NS1/timepoints.tsv"
  # type <- "strain"
  # batch_name <- "SS_NS1"
  
  if(type == "syncom"){
    id_cols <- c("Community", "temp")
  }else if(type == "strain"){
    id_cols <- c("strain", "temp")
  }
  
  
  #' Read and proces OD
  OD <- read_tsv(od600_file)
  Timepoints <- read_tsv(timepoints_file)
  
  Dat <- OD %>%
    pivot_longer(-id_cols,
                 names_to = "timepoint",
                 values_to = "OD600") %>%
    left_join(Timepoints, by = c("timepoint", "temp")) %>%
    mutate(OD600 = (OD600 - blank_od) * od_factor) %>%
    mutate(OD600 = replace(OD600, timepoint == "t_0", 0.001)) 
  
  
  #' Add MPN if available
  if(!is.null(mpn_file) && file.exists(mpn_file)){
    MPN <- read_tsv(mpn_file, na = c("", "NA", "Error"))
    
    Dat <- MPN %>%
      pivot_longer(-id_cols,
                   names_to = "timepoint",
                   values_to = "MPN") %>%
      right_join(Dat, by = c(id_cols, "timepoint")) 
  }
  
  if(!is.null(batch_name)){
    Dat[["batch"]] <- batch_name
  }
  
  
  return(Dat)
  
}