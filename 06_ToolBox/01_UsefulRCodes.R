####
#### Useful R codes
#### 2021.10.17 Ushio
####

# ----------------------------------------- #
# Creat output folder
# ----------------------------------------- #
od_name <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od_name, end = -3), "Out")); rm(od_name)
dir.create(output_folder)
