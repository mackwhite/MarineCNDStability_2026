###project: LTER Marine Consumer Nutrient Dynamic Synthesis Working Group
###author(s): MW, AC, LK, WRJ
###goal(s): Wrangling raw CND data such that it is ready for model data joins
###date(s): February 2026
###note(s): 
## Here, we calculate CND stability and community structure (e.g., evenness, richness) metrics
## at the population, trophic, and community levels. Community level data is used in manuscript (starts on line 736)
## Population level calculations begin on line 22
## Trophic level calculations begin on line 378
## Community level calculations begin on line 736
## Please - direct any questions to mackwhiteecology@gmail.com

# Housekeeping ------------------------------------------------------------
### load necessary libraries
# install.packages("librarian")
librarian::shelf(tidyverse, vegan, readxl, dplyr, splitstackshape)

### set custom functions
nacheck <- function(df) {
      na_count_per_column <- sapply(df, function(x) sum(is.na(x)))
      print(na_count_per_column)
}

# Load and prepare data ----------------------------------------------------
strata_list <- readxl::read_excel(path = file.path("tier2", "strata_class.xlsx"),na=".") |> 
      ### remove decimals from numbered sites
      mutate(site = str_remove(site, "\\.0$"),
             subsite_level1 = str_remove(subsite_level1, "\\.0$"),
             subsite_level2 = str_remove(subsite_level2, "\\.0$"),
             subsite_level3 = str_remove(subsite_level3, "\\.0$"))
dt <- read.csv(file.path("tier2", "harmonized_consumer_excretion_CLEAN.csv"),stringsAsFactors = F,na.strings =".") |> 
      # tidy up column names
      janitor::clean_names() |> 
      # replace NAs with character data
      mutate(subsite_level1 = replace_na(subsite_level1, "Not Available"),
             subsite_level2 = replace_na(subsite_level2, "Not Available"),
             subsite_level3 = replace_na(subsite_level3, "Not Available")) |> 
      
      # filter out projects we are using in this manuscript
      filter(project %in% c('MCR', 'CoastalCA', 'SBC', 'VCR', 'FCE')) |> 
      
      # filter out beach habitat [talitrids at SBC]
      filter(habitat %in% c('estuary', 'ocean'))
glimpse(dt)
nacheck(dt)

### add phylum where necessary and update excretion values for those taxa
dta <- dt |> filter(!is.na(phylum))
dtb1 <- dt |> filter(is.na(phylum), density_num_m2 == 0) |> 
      mutate(phylum = 'Chordata')
dtb2 <- dt |> filter(is.na(phylum), density_num_m2 > 0) |> 
      mutate(phylum = 'Chordata') |> 
      mutate(n_vert_coef = if_else(phylum == "Chordata", 0.7804, 0),
             n_diet_coef = if_else(diet_cat == "algae_detritus", -0.0389,
                                   if_else(diet_cat == "invert", -0.2013,
                                           if_else(diet_cat == "fish", -0.0537,
                                                   if_else(diet_cat == "fish_invert", -0.1732, 
                                                           if_else(diet_cat == "algae_invert", 0,
                                                                   NA))))),
             nexc_log10  = ifelse(dmperind_g_ind > 0, 1.461 + 0.6840*(log10(dmperind_g_ind)) + 0.0246*temp_c + n_diet_coef + n_vert_coef,NA),
             nind_ug_hr  = 10^nexc_log10,
             nind_ug_hr  = ifelse(is.na(nind_ug_hr),0,nind_ug_hr)) |> 
      mutate(p_vert_coef = if_else(phylum == "Chordata", 0.7504, 0),
             p_diet_coef = if_else(diet_cat == "algae_detritus", 0.0173,
                                   if_else(diet_cat == "invert", -0.2480,
                                           if_else(diet_cat == "fish", -0.0337,
                                                   if_else(diet_cat == "fish_invert", -0.4525, 
                                                           if_else(diet_cat == "algae_invert",0,
                                                                   NA))))),
             pexc_log10  = ifelse(dmperind_g_ind >0, 0.6757 + 0.5656*(log10(dmperind_g_ind)) + 0.0194*temp_c + p_diet_coef + p_vert_coef, NA),
             pind_ug_hr  = 10^pexc_log10,
             pind_ug_hr  = ifelse(is.na(pind_ug_hr),0,pind_ug_hr)) |> 
      select(-n_vert_coef, -n_diet_coef, -nexc_log10, 
             -p_vert_coef, -p_diet_coef, -pexc_log10)
dtb <- rbind(dtb1, dtb2)
dt_ab <- rbind(dta, dtb)      

dt1 <- dt_ab |>
      
      # filter out organisms that are not fish
      mutate(
            order = case_when(
                  is.na(order) ~ 'missing',
                  TRUE ~ order
            )) |>
      filter(phylum == 'Chordata',
             order != 'Decapoda') |> 
      
      # set upper end for California Moray eel 
      mutate(
            dmperind_g_ind = case_when(
                  dmperind_g_ind > 9071 & scientific_name == "Gymnothorax mordax" ~ 9071,
                  TRUE ~ dmperind_g_ind
            )) |> 
      
      # pull out big schools of fishes 
      mutate(
            area = case_when(
                  project == 'CoastalCA' ~ 60,
                  project == 'SBC' ~ 80,
                  project == 'VCR' ~ 25,
                  project == 'MCR' & subsite_level3 == '1' ~ 50,
                  project == 'MCR' & subsite_level3 == '5' ~ 250,
                  TRUE ~ NA_real_
            ),
            count = area*density_num_m2
      ) |> 
      filter(project == "FCE" | count <10000) |> # FCE doesn't catch thousands of fish per transect
      select(-area, -count) |> 
      
      # pull out 'biomass buster' sharks and rays from CoastalCA, SBC, and MCR
      group_by(project, habitat) |> 
      mutate(
            mean_dmperind = mean(dmperind_g_ind, na.rm = TRUE),
            sd_dmperind   = sd(dmperind_g_ind, na.rm = TRUE),  
            lower_bound   = mean_dmperind - 5 * sd_dmperind,  
            upper_bound   = mean_dmperind + 5 * sd_dmperind,
            outlier       = dmperind_g_ind < lower_bound | dmperind_g_ind > upper_bound,
            sharkray      = grepl("\\bshark\\b|\\bray\\b", common_name, ignore.case = TRUE),
            elasmo        = class %in% c("Chondrichthyes", "Elasmobranchii")
      ) |> 
      ungroup() |> 
      filter(!(outlier & (sharkray | elasmo))) |> 
      select(-mean_dmperind, -sd_dmperind, -lower_bound, -upper_bound, -outlier, -sharkray, -elasmo) |> 
      
      # coalesce density columns 
      mutate(density = coalesce(density_num_m, density_num_m2)) |> 
      select(-density_num_m, -density_num_m2) |> 
      
      # calculate excretion per ind where phylum was missing in dataset
      mutate(
            
      )
glimpse(dt1)
nacheck(dt1)
head(dt1)
rm(dt, dta, dtb1, dtb2, dtb, dt_ab)

# Summarize data ----------------------------------------------------------
dt2 <- dt1 |> 
      group_by(project, habitat, year, month, 
               site, subsite_level1, subsite_level2, subsite_level3, 
               scientific_name) |> 
      ### sum across unique taxa at the transect level
      summarize(
            total_n = sum(nind_ug_hr*density, na.rm = TRUE),
            total_bm = sum(dmperind_g_ind*density, na.rm = TRUE),
            .groups = 'drop'
      )
glimpse(dt2)
head(dt2)
nacheck(dt2)
### notes - need to consider the fact that diet_category varies within some families... 
### for example, no longer in dt2 - but I don't know if we would want it in dt2...
### need to also calculate trophic richness and species richness
### also, need to remember that MCR needs 'fixed' and FCE needs to be brought to annual time_step