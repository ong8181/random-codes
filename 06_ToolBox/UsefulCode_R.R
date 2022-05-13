####
#### Useful R codes
#### 2021.10.17 Ushio
####

# ----------------------------------------- #
# Creat output folder
# ----------------------------------------- #
## Option 1
od_name <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od_name, end = -3), "Out")); rm(od_name)
dir.create(output_folder)
## Option 2
ouput_folder <- rstudioapi::getSourceEditorContext()$path %>%
   basename %>% str_sub(end = -3) %>% paste0("Out")
dir.create(output_folder)


# ----------------------------------------- #
# tidyverse
# ----------------------------------------- #
# Change the number of rows/columns to be displayed 
options(tibble.print_min = 20)
options(tibble.width = Inf)


# ----------------------------------------- #
# ggplot2 tools
# ----------------------------------------- #
# Label exponents in ggplot
label_func <- function(x) {
  ifelse(x == 0, "0", 
         parse(text = gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x))))
         )
  }

# Create color palette
library(RColorBrewer)
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))
get_palette(20)

# Temperature
g1 + xlab(expression(paste("Temperature (", degree, "C)")))

# Add new line in expression()
g1 + xlab(expression(atop(paste("Water temperature (", degree, "C)"), "(Maximum value)")))

