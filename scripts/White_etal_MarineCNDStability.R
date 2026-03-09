###project: LTER Marine Consumer Nutrient Dynamic Synthesis Working Group
###author(s): MW, AC, LK, WRJ
###goal(s): Wrangling and summarizing raw CND data such that it is ready for analysis
###date(s): Summer 2024, Revised Spring 2026
###note(s): 

# Housekeeping ------------------------------------------------------------
### load necessary libraries
# install.packages("librarian")
librarian::shelf(tidyverse, vegan, readxl, dplyr, splitstackshape, codyn, lavaan,
                 MuMIn, corrplot, performance, ggeffects, ggpubr, parameters, ggstats,
                 brms, mixedup, rstatix, sf, ggspatial)

### set custom functions
nacheck <- function(df) {
      na_count_per_column <- sapply(df, function(x) sum(is.na(x)))
      print(na_count_per_column)
}

# Load and prepare data ----------------------------------------------------
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
      filter(habitat %in% c('estuary', 'ocean')) |> 
      
      # filter out sites not sampled consistently at FCE
      filter(
            !(project == "FCE" & site == "TB" & subsite_level1 == "5"),
            !(project == "FCE" & site == "RB" & subsite_level1 %in% c("17", "19"))
      )
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
      select(-density_num_m, -density_num_m2, -density_num_m3) |> 
      
      # filter out high densities of large-bodied fishes that interrupt ts
      filter(!(density > 1 & nind_ug_hr > 20000))
glimpse(dt1)
nacheck(dt1)
head(dt1)
rm(dt, dta, dtb1, dtb2, dtb, dt_ab)


##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
# Summarize data ---------------------------------------------------------------------------------
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################

##################################################################################################
##################################################################################################
##################################################################################################
### Calculating CND, Biomass and Species Diversity metrics ---------------------------------------
##################################################################################################
##################################################################################################
##################################################################################################

### MCR subsite_level3 has two levels (1 and 5) that should be summed across, 
### as they are conducted within same space and time so split them out and go
### sum across subsite_level2 
dt2_mcr <- dt1 |> 
      filter(project == 'MCR') |> 
      group_by(project, habitat, year, month, 
               site, subsite_level1, subsite_level2, subsite_level3, 
               scientific_name) |> 
      summarize(
            n       = sum(nind_ug_hr*density, na.rm = TRUE),
            bm      = sum(dmperind_g_ind*density, na.rm = TRUE),
            dens    = sum(density, na.rm = TRUE),
            .groups = 'drop'
      ) |> 
      group_by(project, habitat, year, month, 
               site, subsite_level1, subsite_level2, 
               scientific_name) |> 
      summarize(
            species_n    = sum(n),
            species_bm   = sum(bm),
            species_dens = sum(dens),
            .groups = 'drop'
      ) |> 
      mutate(subsite_level3 = '15') |> 
      select(project, habitat, year, month, site, subsite_level1, subsite_level2, subsite_level3,
             scientific_name, species_n, species_bm, species_dens)

glimpse(dt2_mcr)
head(dt2_mcr)
nacheck(dt2_mcr)

### all other sites should be summed down subsite_level3 as these are
### considered the 'transect'
dt2_other <- dt1 |> 
      filter(project != 'MCR') |> 
      group_by(project, habitat, year, month, 
               site, subsite_level1, subsite_level2, subsite_level3, 
               scientific_name) |> 
      ### sum across unique taxa at the transect level
      ### coalesces multiple species observation to single n and bm value for each transect
      summarize(
            species_n    = sum(nind_ug_hr*density, na.rm = TRUE),
            species_bm   = sum(dmperind_g_ind*density, na.rm = TRUE),
            species_dens = sum(density, na.rm = TRUE),
            .groups = 'drop'
      ) |> 
      select(project, habitat, year, month, site, subsite_level1, subsite_level2, subsite_level3,
             scientific_name, species_n, species_bm, species_dens)

glimpse(dt2_other)
head(dt2_other)
nacheck(dt2_other)

### join mcr and other data back together now that they are at same 'transect' scale
dt2 <- rbind(dt2_other, dt2_mcr)
glimpse(dt2)
head(dt2)
nacheck(dt2)
rm(dt2_mcr, dt2_other)

### sum across species at the transect scale to get community-level excretion and biomass
### per unit area
dt3_cnd <- dt2 |> 
      group_by(project, habitat, year, month, 
               site, subsite_level1, subsite_level2, subsite_level3) |> 
      summarize(
            comm_n    = sum(species_n, na.rm = TRUE),
            comm_bm   = sum(species_bm, na.rm = TRUE),
            comm_dens = sum(species_dens, na.rm = TRUE),
            .groups   = 'drop'
      ) |> 
      mutate(project = case_when(
            project  == 'CoastalCA' & site == 'CENTRAL' ~ 'PCCC',
            project  == 'CoastalCA' & site == 'SOUTH' ~ 'PCCS',
            TRUE ~ project
      )) |> 
      ### take everything to the true 'site' level in which transects should be
      ### averaged across
      mutate(
            site_site = case_when(
                  project == 'SBC' ~ site,
                  project == 'FCE' ~ paste(site, subsite_level1, sep = ''),
                  project == 'VCR' ~ paste(site, subsite_level1, sep = ''),
                  project == 'MCR' ~ paste(subsite_level1, site, sep = ''),
                  project == 'PCCC' ~ subsite_level2,
                  project == 'PCCS' ~ subsite_level2,
            ) 
      ) |> 
      group_by(project, site_site, year) |> 
      summarize(
            comm_n    = mean(comm_n, na.rm = TRUE),
            comm_bm   = mean(comm_bm, na.rm = TRUE),
            comm_dens = mean(comm_dens, na.rm = TRUE),
            .groups   = 'drop'
      ) |> 
      rename(site = site_site)
      
glimpse(dt3_cnd)
head(dt3_cnd)
nacheck(dt3_cnd)

## Calculating Species Richness and Diversity ----
head(dt2)
glimpse(dt2)
nacheck(dt2)

dt3_sp <- dt2 |> 
      group_by(project, habitat, year, month, 
               site, subsite_level1, subsite_level2, subsite_level3) |> 
      summarize(
            species_richness = vegan::specnumber(species_dens),                 
            inv_simpson      = vegan::diversity(species_dens, index = "invsimpson"),
            total_dens       = sum(species_dens, na.rm = TRUE),
            .groups = "drop"
      ) |> 
      mutate(
            inv_simpson = case_when(
                  inv_simpson == Inf ~ 0,
                  TRUE ~ inv_simpson
            )) |> 
      mutate(project = case_when(
            project  == 'CoastalCA' & site == 'CENTRAL' ~ 'PCCC',
            project  == 'CoastalCA' & site == 'SOUTH' ~ 'PCCS',
            TRUE ~ project
      )) |> 
      ### take everything to the true 'site' level in which transects should be
      ### averaged across
      mutate(
            site_site = case_when(
                  project == 'SBC' ~ site,
                  project == 'FCE' ~ paste(site, subsite_level1, sep = ''),
                  project == 'VCR' ~ paste(site, subsite_level1, sep = ''),
                  project == 'MCR' ~ paste(subsite_level1, site, sep = ''),
                  project == 'PCCC' ~ subsite_level2,
                  project == 'PCCS' ~ subsite_level2,
            ) 
      ) |> 
      group_by(project, site_site, year) |> 
      summarize(
            s_rich    = mean(species_richness, na.rm = TRUE),
            s_div     = mean(inv_simpson, na.rm = TRUE),
            .groups   = 'drop'
      ) |> 
      rename(site = site_site)
glimpse(dt3_sp)
head(dt3_sp)
nacheck(dt3_sp)

##################################################################################################
##################################################################################################
##################################################################################################
### Calculating Trophic Diversity metrics --------------------------------------------------------
##################################################################################################
##################################################################################################
##################################################################################################

dt2_mcr_troph <- dt1 |> 
      filter(project == 'MCR') |> 
      group_by(project, habitat, year, month, 
               site, subsite_level1, subsite_level2, subsite_level3, 
               diet_cat) |> 
      summarize(
            n       = sum(nind_ug_hr*density, na.rm = TRUE),
            bm      = sum(dmperind_g_ind*density, na.rm = TRUE),
            dens    = sum(density, na.rm = TRUE),
            .groups = 'drop'
      ) |> 
      group_by(project, habitat, year, month, 
               site, subsite_level1, subsite_level2, 
               diet_cat) |> 
      summarize(
            troph_n    = sum(n),
            troph_bm   = sum(bm),
            troph_dens = sum(dens),
            .groups = 'drop'
      ) |> 
      mutate(subsite_level3 = '15') |> 
      select(project, habitat, year, month, site, subsite_level1, subsite_level2, subsite_level3,
             diet_cat, troph_n, troph_bm, troph_dens)

glimpse(dt2_mcr_troph)
head(dt2_mcr_troph)
nacheck(dt2_mcr_troph)

### all other sites should be summed down subsite_level3 as these are
### considered the 'transect'
dt2_other_troph <- dt1 |> 
      filter(project != 'MCR') |> 
      group_by(project, habitat, year, month, 
               site, subsite_level1, subsite_level2, subsite_level3, 
               diet_cat) |> 
      ### sum across unique taxa at the transect level
      ### coalesces multiple troph observation to single n and bm value for each transect
      summarize(
            troph_n    = sum(nind_ug_hr*density, na.rm = TRUE),
            troph_bm   = sum(dmperind_g_ind*density, na.rm = TRUE),
            troph_dens = sum(density, na.rm = TRUE),
            .groups = 'drop'
      ) |> 
      select(project, habitat, year, month, site, subsite_level1, subsite_level2, subsite_level3,
             diet_cat, troph_n, troph_bm, troph_dens)

glimpse(dt2_other_troph)
head(dt2_other_troph)
nacheck(dt2_other_troph)

### join mcr and other data back together now that they are at same 'transect' scale
dt2_troph <- rbind(dt2_other_troph, dt2_mcr_troph)
glimpse(dt2_troph)
head(dt2_troph)
nacheck(dt2_troph)
rm(dt2_mcr_troph, dt2_other_troph)

## Calculating Species Richness ----
head(dt2_troph)
glimpse(dt2_troph)
nacheck(dt2_troph)

dt3_troph <- dt2_troph |> 
      group_by(project, habitat, year, month, 
               site, subsite_level1, subsite_level2, subsite_level3) |> 
      summarize(
            trophic_richness       = vegan::specnumber(troph_dens),                 
            inv_simpson_troph      = vegan::diversity(troph_dens, index = "invsimpson"),
            total_dens_troph       = sum(troph_dens, na.rm = TRUE),
            .groups = "drop"
      ) |> 
      mutate(
            inv_simpson_troph = case_when(
                  inv_simpson_troph == Inf ~ 0,
                  TRUE ~ inv_simpson_troph
            )) |> 
      mutate(project = case_when(
            project  == 'CoastalCA' & site == 'CENTRAL' ~ 'PCCC',
            project  == 'CoastalCA' & site == 'SOUTH' ~ 'PCCS',
            TRUE ~ project
      )) |> 
      ### take everything to the true 'site' level in which transects should be
      ### averaged across
      mutate(
            site_site = case_when(
                  project == 'SBC' ~ site,
                  project == 'FCE' ~ paste(site, subsite_level1, sep = ''),
                  project == 'VCR' ~ paste(site, subsite_level1, sep = ''),
                  project == 'MCR' ~ paste(subsite_level1, site, sep = ''),
                  project == 'PCCC' ~ subsite_level2,
                  project == 'PCCS' ~ subsite_level2,
            ) 
      ) |> 
      group_by(project, site_site, year) |> 
      summarize(
            t_rich    = mean(trophic_richness, na.rm = TRUE),
            t_div     = mean(inv_simpson_troph, na.rm = TRUE),
            .groups   = 'drop'
      ) |> 
      rename(site = site_site)
glimpse(dt3_troph)
head(dt3_troph)
nacheck(dt3_troph)

### bring it all together ----
dt3 <- dt3_cnd |> 
      left_join(dt3_sp) |> 
      left_join(dt3_troph)
glimpse(dt3)
head(dt3)
nacheck(dt3)
cnd_ts_data <- dt3
# write_csv(cnd_ts_data, "local_data/annual-dt-for-summary.csv")

cnd_model_data <- cnd_ts_data |> 
      group_by(project, site) |> 
      summarize(
            comm_n_mean      = mean(comm_n, na.rm = TRUE),
            comm_n_sd        = sd(comm_n, na.rm = TRUE),
            comm_n_cv        = sd(comm_n, na.rm = TRUE)/mean(comm_n, na.rm = TRUE),
            comm_n_stability = 1/comm_n_cv,
            s_rich_mean      = mean(s_rich, na.rm = TRUE),
            s_div_mean       = mean(s_div, na.rm = TRUE),
            t_rich_mean      = mean(t_rich, na.rm = TRUE),
            t_div_mean       = mean(t_div, na.rm = TRUE),
            .groups = 'drop'
      )
glimpse(cnd_model_data)
head(cnd_model_data)
nacheck(cnd_model_data)
# write_csv(cnd_model_data, "local_data/community-level-nutrient-stability.csv")

##################################################################################################
##################################################################################################
##################################################################################################
### Calculating Species Dynamics -----------------------------------------------------------------
##################################################################################################
##################################################################################################
##################################################################################################
glimpse(dt2)
head(dt2)
nacheck(dt2)

dt3_sp_dyn <- dt2 |> 
      mutate(project = case_when(
            project  == 'CoastalCA' & site == 'CENTRAL' ~ 'PCCC',
            project  == 'CoastalCA' & site == 'SOUTH' ~ 'PCCS',
            TRUE ~ project
      )) |> 
      mutate(
            site_site = case_when(
                  project == 'SBC' ~ site,
                  project == 'FCE' ~ paste(site, subsite_level1, sep = ''),
                  project == 'VCR' ~ paste(site, subsite_level1, sep = ''),
                  project == 'MCR' ~ paste(subsite_level1, site, sep = ''),
                  project == 'PCCC' ~ subsite_level2,
                  project == 'PCCS' ~ subsite_level2,
            ) 
      ) |> 
      group_by(project, site_site, year, scientific_name) |> 
      summarize(
            spp_n    = mean(species_n, na.rm = TRUE),
            spp_bm   = mean(species_bm, na.rm = TRUE),
            spp_dens = mean(species_dens, na.rm = TRUE),
            .groups   = 'drop'
      ) |> 
      rename(site = site_site)
glimpse(dt3_sp_dyn)
head(dt3_sp_dyn)
nacheck(dt3_sp_dyn)

spp_dyn_model_data <- dt3_sp_dyn |> 
      group_by(project, site) |> 
      nest() |> 
      mutate(
            spp_turnover = map_dbl(data, ~{
                  beta_temp <- turnover(
                        df            = .x,
                        time.var      = "year",
                        abundance.var = "spp_dens",
                        species.var   = "scientific_name",
                        metric        = "total"
                  )
                  mean(beta_temp[, 1], na.rm = TRUE)
            }),
            spp_synchrony = map_dbl(data, ~{
                  synchrony(
                        df            = .x,
                        time.var      = "year",
                        species.var   = "scientific_name",
                        abundance.var = "spp_bm",
                        metric        = "Loreau",
                        replicate.var = NA
                  )
            })
      ) |> 
      ungroup() |> 
      select(project, site, spp_turnover, spp_synchrony)

##################################################################################################
##################################################################################################
##################################################################################################
### Calculating Trophic Dynamics -----------------------------------------------------------------
##################################################################################################
##################################################################################################
##################################################################################################
glimpse(dt2_troph)
head(dt2_troph)
nacheck(dt2_troph)

dt3_troph_dyn <- dt2_troph |> 
      mutate(project = case_when(
            project  == 'CoastalCA' & site == 'CENTRAL' ~ 'PCCC',
            project  == 'CoastalCA' & site == 'SOUTH' ~ 'PCCS',
            TRUE ~ project
      )) |> 
      mutate(
            site_site = case_when(
                  project == 'SBC' ~ site,
                  project == 'FCE' ~ paste(site, subsite_level1, sep = ''),
                  project == 'VCR' ~ paste(site, subsite_level1, sep = ''),
                  project == 'MCR' ~ paste(subsite_level1, site, sep = ''),
                  project == 'PCCC' ~ subsite_level2,
                  project == 'PCCS' ~ subsite_level2,
            ) 
      ) |> 
      group_by(project, site_site, year, diet_cat) |> 
      summarize(
            troph_n    = mean(troph_n, na.rm = TRUE),
            troph_bm   = mean(troph_bm, na.rm = TRUE),
            troph_dens = mean(troph_dens, na.rm = TRUE),
            .groups    = 'drop'
      ) |> 
      rename(site = site_site)
glimpse(dt3_troph_dyn)
head(dt3_troph_dyn)
nacheck(dt3_troph_dyn)

troph_dyn_model_data <- dt3_troph_dyn |> 
      group_by(project, site) |> 
      nest() |> 
      mutate(
            troph_turnover = map_dbl(data, ~{
                  beta_temp <- turnover(
                        df            = .x,
                        time.var      = "year",
                        abundance.var = "troph_dens",
                        species.var   = "diet_cat",
                        metric        = "total"
                  )
                  mean(beta_temp[, 1], na.rm = TRUE)
            }),
            troph_synchrony = map_dbl(data, ~{
                  synchrony(
                        df            = .x,
                        time.var      = "year",
                        species.var   = "diet_cat",
                        abundance.var = "troph_bm",
                        metric        = "Loreau",
                        replicate.var = NA
                  )
            })
      ) |> 
      ungroup() |> 
      select(project, site, troph_turnover, troph_synchrony)

### join everything together ----
glimpse(cnd_model_data)
glimpse(spp_dyn_model_data)
glimpse(troph_dyn_model_data)

model_data_all <- cnd_model_data |> 
      left_join(spp_dyn_model_data) |> 
      left_join(troph_dyn_model_data)
glimpse(model_data_all)
head(model_data_all)
nacheck(model_data_all)
# write_csv(model_data_all, "local_data/model-data-all.csv")

##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
# Run across system model ------------------------------------------------------------------------
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################

keep <- c("nacheck", "model_data_all", "cnd_ts_data")
rm(list = setdiff(ls(), keep))

glimpse(model_data_all)
dat_scaled <- model_data_all |> 
      rename(program = project) |> 
      select(program, site, comm_n_stability, everything()) |> 
      mutate(comm_n_stability = as.numeric(scale(comm_n_stability))) |>
      mutate(across(comm_n_mean:troph_synchrony, \(x) as.numeric(scale(x, center = TRUE))))
glimpse(dat_scaled)      

glimpse(dat_scaled)
dat_ready <- dat_scaled      
glimpse(dat_ready)

# path_model <- '
#   # Structural equations for the path model
#   comm_n_stability ~ s_rich_mean + troph_turnover + spp_turnover + spp_synchrony  # Stability regressed on Richness, Trophic Turnover, Pop. Turnover, Synchrony
#   spp_turnover     ~ s_rich_mean + t_rich_mean                                    # Population Turnover regressed on Species Richness and Trophic Richness
#   spp_synchrony    ~ s_div_mean + s_rich_mean + t_div_mean                        # Population Synchrony on Species Evenness, Species Richness, Trophic Evenness
# 
#   # Label specific paths to calculate indirect effects
#   spp_synchrony ~ a_se*s_div_mean                                                 # a_se: effect of Species Evenness on Synchrony
#   comm_n_stability ~ b_syn*spp_synchrony                                          # b_syn: effect of Synchrony on Stability
# 
#   # Define the indirect effect of Species Evenness on Stability via Synchrony
#   indirect_evenness := a_se * b_syn
# '
# 
# fit <- sem(path_model, data = dat_ready)
# summary(
#       fit,
#       standardized = TRUE,
#       fit.measures = TRUE,
#       rsquare = TRUE
# )

path_model <- '
   # Regressions
   comm_n_stability ~ cp*s_rich_mean + b1*spp_synchrony + b2*spp_turnover + b3*troph_turnover
   
   spp_synchrony    ~ a1*s_rich_mean
   
   spp_turnover     ~ a2*t_rich_mean + a4*s_rich_mean
   
   troph_turnover   ~ a3*t_rich_mean

   # Covariances
   s_rich_mean      ~~ t_rich_mean
   spp_synchrony    ~~ troph_turnover
   troph_turnover   ~~ s_rich_mean
   spp_turnover     ~~ troph_turnover

   # Indirect Effects 
   # ind_rich_sync := a1 * b1
   # ind_rich_turn := a4 * b2
   # ind_troph_turn := a3 * b2

   # total_rich_effect := cp + (a1 * b1) + (a4 * b2)
   # total_rich_effect := cp + ind_rich_sync + ind_rich_turn

'

fit <- sem(path_model, data = dat_ready, estimator = "MLR")
summary(fit, standardized = TRUE, fit.measures = TRUE)
# modificationIndices(fit, sort = TRUE, maximum.number = 3)

# run across system model -------------------------------------------------
keep <- c("nacheck", "model_data_all", "cnd_ts_data", "fit")
rm(list = setdiff(ls(), keep))

glimpse(model_data_all)
dat_scaled <- model_data_all |> 
      rename(program = project) |> 
      select(program, site, comm_n_stability, everything()) |> 
      mutate(comm_n_stability = as.numeric(scale(comm_n_stability))) |>
      group_by(program) |> 
      mutate(across(comm_n_mean:troph_synchrony, \(x) as.numeric(scale(x, center = TRUE)))) |> 
      ungroup()
glimpse(dat_scaled)      
glimpse(dat_scaled)
dat_ready <- dat_scaled      
glimpse(dat_ready)

### set priors following Lemoine (2019, Ecology)
pr = prior(normal(0, 1), class = 'b')

##################################################################################################
##################################################################################################
##################################################################################################
### Full Models ----------------------------------------------------------------------------------
##################################################################################################
##################################################################################################
##################################################################################################

test_corr <- dat_ready |> select(s_rich_mean, s_div_mean,
                                 t_rich_mean, t_div_mean,
                                 spp_turnover, spp_synchrony,
                                 troph_turnover, troph_synchrony)

matrix <- cor(test_corr, use = 'complete.obs')

corrplot(matrix, method = "number", type = "lower", tl.col = "black", tl.srt = 45)
glimpse(dat_ready)

### round one ------------------------------------------------------------------------------------
m1 <- brm(
      comm_n_stability ~ t_rich_mean + (t_rich_mean | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

m4 <- brm(
      comm_n_stability ~ spp_synchrony + (spp_synchrony | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

### best fit single term model
# saveRDS(m4, file = 'local_data/rds-single-synchrony.rds')

m5 <- brm(
      comm_n_stability ~ spp_turnover + (spp_turnover | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

m6 <- brm(
      comm_n_stability ~ troph_synchrony + (troph_synchrony | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

m7 <- brm(
      comm_n_stability ~ troph_turnover + (troph_turnover | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

m8 <- brm(
      comm_n_stability ~ s_rich_mean + (s_rich_mean | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

### model of interest, given across-ecosystem importance
# saveRDS(m8, file = 'local_data/rds-single-richness.rds')

m9 <- brm(
      comm_n_stability ~ s_div_mean + (s_div_mean | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

m10 <- brm(
      comm_n_stability ~ t_div_mean + (t_div_mean | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

model_table_all <- performance::compare_performance(m1,m4,m5,m6,m7,m8,m9,m10)

model_selection1 <- model_table_all |>
      mutate(dWAIC = WAIC - min(WAIC))

write_csv(model_selection1, "output/tables/brms-fullmodel-selection-table-roundone.csv")

keep <- c("nacheck", "model_data_all", "cnd_ts_data", "fit", "dat_ready", "pr", "palette", 'm4', 'model_selection1')
rm(list = setdiff(ls(), keep))

##################################################################################################
### round two ------------------------------------------------------------------------------------
##################################################################################################

m41 <- brm(
      comm_n_stability ~ t_rich_mean + spp_synchrony + (t_rich_mean + spp_synchrony | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

m45 <- brm(
      comm_n_stability ~ spp_turnover + spp_synchrony + (spp_turnover + spp_synchrony | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

m47 <- brm(
      comm_n_stability ~ troph_turnover + spp_synchrony + (troph_turnover + spp_synchrony | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

m48 <- brm(
      comm_n_stability ~ s_rich_mean + spp_synchrony + (s_rich_mean + spp_synchrony | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

m49 <- brm(
      comm_n_stability ~ t_div_mean + spp_synchrony + (t_div_mean + spp_synchrony | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

m40 <- brm(
      comm_n_stability ~ s_div_mean + spp_synchrony + (s_div_mean + spp_synchrony | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

model_table_all <- performance::compare_performance(m41,m4,m45,m47,m48,m49,m40)

model_selection2 <- model_table_all |>
      mutate(dWAIC = WAIC - min(WAIC))
write_csv(model_selection2, "output/tables/brms-fullmodel-selection-table-roundtwo.csv")

keep <- c("nacheck", "model_data_all", "cnd_ts_data", "fit", "dat_ready", "pr", "palette", 'm45', 'model_selection1', 'model_selection2')
rm(list = setdiff(ls(), keep))

##################################################################################################
### round three ----------------------------------------------------------------------------------
##################################################################################################

m451 <- brm(
      comm_n_stability ~ t_rich_mean + spp_turnover + spp_synchrony + (t_rich_mean + spp_turnover + spp_synchrony | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

m457 <- brm(
      comm_n_stability ~ troph_turnover + spp_turnover + spp_synchrony + (troph_turnover + spp_turnover + spp_synchrony | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

m458 <- brm(
      comm_n_stability ~ s_rich_mean + spp_turnover + spp_synchrony + (s_rich_mean + spp_turnover + spp_synchrony | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

m459 <- brm(
      comm_n_stability ~ t_div_mean + spp_turnover + spp_synchrony + (t_div_mean + spp_turnover + spp_synchrony | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

m450 <- brm(
      comm_n_stability ~ s_div_mean + spp_turnover + spp_synchrony + (s_div_mean + spp_turnover + spp_synchrony | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

model_table_all <- performance::compare_performance(m45,m451,m457,m458,m459,m450)

model_selection3 <- model_table_all |>
      mutate(dWAIC = WAIC - min(WAIC))
write_csv(model_selection3, "output/tables/brms-fullmodel-selection-table-roundthree.csv")

keep <- c("nacheck", "model_data_all", "cnd_ts_data", "fit", "dat_ready", "pr", "palette", 'm45', 'model_selection1', 'model_selection2', 'model_selection3')
rm(list = setdiff(ls(), keep))
full_model <- m45
saveRDS(full_model, file = 'local_data/rds-full-model.rds')

##################################################################################################
##################################################################################################
##################################################################################################
### Visualize Models -----------------------------------------------------------------------------
##################################################################################################
##################################################################################################
##################################################################################################

##################################################################################################
### Single Term Models ---------------------------------------------------------------------------
##################################################################################################

synch = readRDS("local_data/rds-single-synchrony.rds")
rich = readRDS("local_data/rds-single-richness.rds")
glimpse(dat_ready)

# stats ----
p_synch = posterior_samples(synch)
glimpse(p_synch)
mean(p_synch$`r_program[PCCC,spp_synchrony]` + p_synch$b_spp_synchrony  < 0)
mean(p_synch$`r_program[SBC,spp_synchrony]`  + p_synch$b_spp_synchrony  < 0)
mean(p_synch$`r_program[VCR,spp_synchrony]`  + p_synch$b_spp_synchrony  < 0)
mean(p_synch$`r_program[PCCS,spp_synchrony]` + p_synch$b_spp_synchrony  < 0)
mean(p_synch$`r_program[MCR,spp_synchrony]`  + p_synch$b_spp_synchrony  < 0)
mean(p_synch$`r_program[FCE,spp_synchrony]`  + p_synch$b_spp_synchrony  < 0)

p_sr = posterior_samples(rich)
mean(p_sr$`r_program[PCCC,s_rich_mean]` + p_sr$b_s_rich_mean > 0)
mean(p_sr$`r_program[SBC,s_rich_mean]`  + p_sr$b_s_rich_mean > 0)
mean(p_sr$`r_program[VCR,s_rich_mean]`  + p_sr$b_s_rich_mean > 0)
mean(p_sr$`r_program[PCCS,s_rich_mean]` + p_sr$b_s_rich_mean > 0)
mean(p_sr$`r_program[MCR,s_rich_mean]`  + p_sr$b_s_rich_mean > 0)
mean(p_sr$`r_program[FCE,s_rich_mean]`  + p_sr$b_s_rich_mean > 0)

### set color scheme ---
program_palette = c("Overall"="#000000",
                    "FCE"="#64a988",
                    "MCR"="#ff967d",
                    'PCCC'="#2A788EFF",
                    "PCCS"="#8b6b93",
                    'SBC'='#ff3f4c',
                    "VCR"="#9b9254")


prog = c('SBC', 'PCCC', 'VCR', 'PCCS', 'FCE', 'MCR')

# species richness 
# random effects
re95 = mixedup::extract_random_coefs(rich, ci_level = c(0.95))
re80 = mixedup::extract_random_coefs(rich, ci_level = c(0.8))

re_beta = left_join(re95, re80) |> 
      rename(term = effect,
             program = group)

# fixed effects
fe95 = mixedup::extract_fixed_effects(rich, ci_level = c(0.95))
fe80 = mixedup::extract_fixed_effects(rich, ci_level = c(0.8))

fe_beta = left_join(fe95, fe80) |> 
      mutate(program = 'Overall') 

# make data frame of all betas
df_beta = bind_rows(re_beta, fe_beta) |> 
      filter(term != 'Intercept') |> 
      mutate(program = factor(program, levels = prog))

# make equation data set
df_eq = bind_rows(re_beta, fe_beta) |> 
      select(term, program, value) |> 
      pivot_wider(names_from = term, values_from = value)  |> 
      rename(beta = s_rich_mean)

# load and get min and max values
### read in necessary data ---
ov <- model_data_all |>
      rename(program = project) |> 
      select(program, value = s_rich_mean) |>
      distinct() |> 
      group_by(program) |> 
      mutate(scaled = scale(value)) |> 
      ungroup() |> 
      slice(c(which.min(value), which.max(value))) |> 
      mutate(program = 'Overall')

dat <-model_data_all |>
      rename(program = project) |> 
      select(program, value = s_rich_mean) |>
      distinct() |> 
      group_by(program) |> 
      mutate(scaled = scale(value)) |>  
      slice(c(which.min(value), which.max(value))) |> 
      ungroup() |> 
      bind_rows(ov)

dat_scaled <- dat |>
              group_by(program, scaled) |> 
              mutate(i = row_number()) |> 
              select(program, scaled)

raw <-model_data_all |>
      rename(program = project) |> 
      select(program, 
             value = s_rich_mean,
             stab = comm_n_stability) 

df <- dat_scaled |> 
      left_join(df_eq) |> 
      mutate(pred = beta*scaled + Intercept,
             stab = pred*sd(raw$stab) + mean(raw$stab)) |> 
      left_join(dat) |> 
      mutate(program = factor(program, levels = prog))

a <- ggplot(df |> filter(program != 'Overall'), aes(value, stab, color = program)) +
      geom_point(data = raw, aes(value, stab, color = program), size = 2)+
      geom_line(linewidth = 1.75)+
      scale_color_manual(values = program_palette)+
      labs(y = 'CND Stability', x = 'Species Richness')+
      theme_classic()+
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 12),
            axis.text.y = element_text(face = "bold", color = "black", size = 12),
            axis.title.x = element_text(face = "bold", color = "black", size = 14),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            strip.text = element_blank(),
            strip.background = element_blank(),
            legend.position = "none",
            legend.text = element_text(face = "bold", color = "black", size = 12),
            legend.title = element_text(face = "bold", color = "black", size = 14))
a

b <- ggplot(df_beta|> filter(program != 'Overall'), aes(program, value, color = program)) +
      geom_hline(aes(yintercept = 0), linetype = "dashed", linewidth = 0.75) +
      geom_pointrange(aes(ymin = lower_10, ymax = upper_90), linewidth = 2) +
      geom_pointrange(aes(ymin = lower_2.5, ymax = upper_97.5), linewidth = 1, size = .9) +
      labs(y = 'Beta', x = NULL) +
      scale_color_manual(values = program_palette) +
      scale_y_continuous(limits = c(-0.35, 1.03)) +
      coord_flip() +
      theme_classic() +
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 12),
            axis.text.y = element_text(face = "bold", color = "black", size = 12),
            axis.title.x = element_text(face = "bold", color = "black", size = 14),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            strip.text = element_text(face = "bold", color = "black", size = 12),
            legend.position = "none",
            legend.background = element_blank(),
            legend.text = element_text(face = "bold", color = "black", size = 12),
            legend.title = element_text(face = "bold", color = "black", size = 14))
b

ggpubr::ggarrange(a,b, align = 'h')

sr = ggpubr::ggarrange(a,b, labels = c('c', 'd'), align = 'h', legend = 'none')

# synchrony 
# random effects
re95 = mixedup::extract_random_coefs(synch, ci_level = c(0.95))
re80 = mixedup::extract_random_coefs(synch, ci_level = c(0.8))

re_beta = left_join(re95, re80) |> 
      rename(term = effect,
             program = group)

# fixed effects
fe95 = mixedup::extract_fixed_effects(synch, ci_level = c(0.95))
fe80 = mixedup::extract_fixed_effects(synch, ci_level = c(0.8))

fe_beta = left_join(fe95, fe80) |> 
      mutate(program = 'Overall') 

# make data frame of all betas
df_beta = bind_rows(re_beta, fe_beta) |> 
      filter(term != 'Intercept') |> 
      mutate(program = factor(program, levels = prog))

# make equation data set
df_eq = bind_rows(re_beta, fe_beta) |> 
      select(term, program, value) |> 
      pivot_wider(names_from = term, values_from = value)  |> 
      rename(beta = spp_synchrony)

ov <- model_data_all |>
      rename(program = project) |> 
      select(program, value = spp_synchrony) |>
      distinct() |> 
      group_by(program) |> 
      mutate(scaled = scale(value)) |> 
      ungroup() |> 
      slice(c(which.min(value), which.max(value))) |> 
      mutate(program = 'Overall')

dat <-model_data_all |>
      rename(program = project) |> 
      select(program, value = spp_synchrony) |>
      distinct() |> 
      group_by(program) |> 
      mutate(scaled = scale(value)) |>  
      slice(c(which.min(value), which.max(value))) |> 
      ungroup() |> 
      bind_rows(ov)

dat_scaled <- dat |>
      group_by(program, scaled) |> 
      mutate(i = row_number()) |> 
      select(program, scaled)

raw <-model_data_all |>
      rename(program = project) |> 
      select(program, 
             value = spp_synchrony,
             stab = comm_n_stability) 

df <- dat_scaled |> 
      left_join(df_eq) |> 
      mutate(pred = beta*scaled + Intercept,
             stab = pred*sd(raw$stab) + mean(raw$stab)) |> 
      left_join(dat) |> 
      mutate(program = factor(program, levels = prog))

a <- ggplot(df |> filter(program != 'Overall'), aes(value, stab, color = program)) +
      geom_point(data = raw, aes(value, stab, color = program), size = 2)+
      geom_line(linewidth = 1.75)+
      scale_color_manual(values = program_palette)+
      labs(y = 'CND Stability', x = 'Species Synchrony')+
      theme_classic()+
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 12),
            axis.text.y = element_text(face = "bold", color = "black", size = 12),
            axis.title.x = element_text(face = "bold", color = "black", size = 14),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            strip.text = element_blank(),
            strip.background = element_blank(),
            legend.position = "none",
            legend.text = element_text(face = "bold", color = "black", size = 12),
            legend.title = element_text(face = "bold", color = "black", size = 14))
a

b <- ggplot(df_beta|> filter(program != 'Overall'), aes(program, value, color = program)) +
      geom_hline(aes(yintercept = 0), linetype = "dashed", linewidth = 0.75) +
      geom_pointrange(aes(ymin = lower_10, ymax = upper_90), linewidth = 2) +
      geom_pointrange(aes(ymin = lower_2.5, ymax = upper_97.5), linewidth = 1, size = .9) +
      labs(y = 'Beta', x = NULL) +
      scale_color_manual(values = program_palette) +
      # scale_y_continuous(limits = c(-0.35, 1.03)) +
      coord_flip() +
      theme_classic() +
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 12),
            axis.text.y = element_text(face = "bold", color = "black", size = 12),
            axis.title.x = element_text(face = "bold", color = "black", size = 14),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            strip.text = element_text(face = "bold", color = "black", size = 12),
            legend.position = "none",
            legend.background = element_blank(),
            legend.text = element_text(face = "bold", color = "black", size = 12),
            legend.title = element_text(face = "bold", color = "black", size = 14))
b

syn = ggpubr::ggarrange(a,b, labels = c('a', 'b', align = 'h', legend = 'none'))
ggarrange(syn, sr, align = 'v', nrow =2)

ggsave('output/fig3.png', dpi = 600, units= 'in', height = 6, width = 6)

##################################################################################################
### Full Model -----------------------------------------------------------------------------------
##################################################################################################

### read in necessary data ----
full_model = readRDS('local_data/rds-full-model.rds')

#summary stats -----
post = posterior_samples(full_model)

mean(post$`r_program[MCR,spp_synchrony]`+ post$b_spp_synchrony < 0)
mean(post$`r_program[MCR,spp_turnover]` + post$b_spp_turnover < 0)
mean(post$`r_program[PCCS,spp_synchrony]`+ post$b_spp_synchrony < 0)
mean(post$`r_program[PCCS,spp_turnover]` + post$b_spp_turnover < 0)
mean(post$`r_program[FCE,spp_synchrony]`+ post$b_spp_synchrony < 0)
mean(post$`r_program[FCE,spp_turnover]` + post$b_spp_turnover > 0)
mean(post$`r_program[PCCC,spp_synchrony]`+ post$b_spp_synchrony > 0)
mean(post$`r_program[PCCC,spp_turnover]` + post$b_spp_turnover < 0)
mean(post$`r_program[SBC,spp_synchrony]`+ post$b_spp_synchrony < 0)
mean(post$`r_program[SBC,spp_turnover]` + post$b_spp_turnover < 0)
mean(post$`r_program[VCR,spp_synchrony]`+ post$b_spp_synchrony < 0)
mean(post$`r_program[VCR,spp_turnover]` + post$b_spp_turnover < 0)

# random effects
re95 = mixedup::extract_random_coefs(full_model, ci_level = c(0.95))
re80 = mixedup::extract_random_coefs(full_model, ci_level = c(0.8))
re90 = mixedup::extract_random_coefs(full_model, ci_level = c(0.9))
re50 = mixedup::extract_random_coefs(full_model, ci_level = c(0.5))

re_beta = left_join(re95, re80) |> 
      left_join(re90) |> 
      left_join(re50)|> 
      rename(term = effect,
             program = group)

# fixed effects
fe95 = mixedup::extract_fixed_effects(full_model, ci_level = c(0.95))
fe80 = mixedup::extract_fixed_effects(full_model, ci_level = c(0.8))
fe90 = mixedup::extract_fixed_effects(full_model, ci_level = c(0.9))
fe50 = mixedup::extract_fixed_effects(full_model, ci_level = c(0.5))

fe_beta = left_join(fe95, fe80) |> 
      left_join(fe90) |> 
      left_join(fe50) |> 
      mutate(program = 'Overall') 

# make data frame of all betas
df_beta = bind_rows(re_beta, fe_beta) |> 
      filter(term != 'Intercept') |> 
      mutate(program = factor(program, levels = prog),
             term = factor(term, levels = c('spp_synchrony',
                                            'spp_turnover'),
                           labels = c('Species Synchrony',
                                      'Species Turnover')))

# make equation data set
df_eq = bind_rows(re_beta, fe_beta) |> 
      select(term, program, value) |> 
      pivot_wider(names_from = term, values_from = value) |> 
      pivot_longer(spp_turnover:spp_synchrony, names_to = 'term', values_to = 'beta')

ov =  model_data_all |>
      rename(program = project) |> 
      select(program,
             spp_turnover,
             spp_synchrony) |>
      distinct() |> 
      pivot_longer(spp_turnover:spp_synchrony, 
                   names_to = 'term', values_to = 'value') |> 
      group_by(program,term) |> 
      mutate(scaled = scale(value)) |> 
      group_by(term) |> 
      slice(c(which.min(value), which.max(value))) |> 
      mutate(program = 'Overall')

dat = model_data_all |>
      rename(program = project) |> 
      select(program,
             spp_turnover,
             spp_synchrony) |>
      distinct() |> 
      pivot_longer(spp_turnover:spp_synchrony, 
                   names_to = 'term', values_to = 'value') |> 
      group_by(program,term) |> 
      mutate(scaled = scale(value)) |>  
      slice(c(which.min(value), which.max(value))) |> 
      ungroup() |> 
      bind_rows(ov)

dat_scaled = dat |>
      group_by(program, term) |> 
      mutate(i = row_number()) |> 
      select(program, term, scaled)

df = dat_scaled |> 
      left_join(df_eq) |> 
      mutate(stab = beta*scaled + Intercept) |> 
      left_join(dat) |> 
      mutate(program = factor(program, levels = prog),
             term = factor(term, levels = c('spp_synchrony',
                                            'spp_turnover'),
                           labels = c('Species Synchrony',
                                      'Species Turnover')))

a = ggplot(df|> filter(program != 'Overall'), aes(value, stab, color = program))+
      geom_line(linewidth = 1.75)+
      scale_color_manual(values = program_palette)+
      facet_wrap(~term, ncol = 1, strip.position = 'right', scales = 'free')+
      labs(y = 'CND Stability', x = NULL)+
      theme_classic()+
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_text(face = "bold", color = "black", size = 12),
            axis.title.y = element_text(face = "bold", color = "black", size = 12),
            strip.text = element_blank(),
            strip.background = element_blank(),
            legend.position = "bottom",
            legend.text = element_text(face = "bold", color = "black", size = 12),
            legend.title = element_text(face = "bold", color = "black", size = 14))
a

# betas 
b = ggplot(df_beta|> filter(program != 'Overall'), aes(program, value, color = program))+
      geom_hline(aes(yintercept = 0), linetype = "dashed", linewidth = 0.75) +
      geom_pointrange(aes(ymin = lower_10, ymax = upper_90), linewidth = 2)+
      geom_pointrange(aes(ymin = lower_2.5, ymax = upper_97.5), linewidth = 1, size = .9)+
      labs(y = 'Beta', x = NULL)+
      scale_color_manual(values = program_palette)+
      coord_flip()+
      facet_wrap(~term, ncol = 1, strip.position = 'right')+
      theme_classic()+
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 10),
            axis.text.y = element_text(face = "bold", color = "black", size = 10),
            axis.title.x = element_text(face = "bold", color = "black", size = 12),
            axis.title.y = element_text(face = "bold", color = "black", size = 12),
            strip.text = element_text(face = "bold", color = "black", size = 10),
            legend.position = "none",
            legend.background = element_blank(),
            legend.text = element_text(face = "bold", color = "black", size = 12),
            legend.title = element_text(face = "bold", color = "black", size = 14))
b

plot = ggpubr::ggarrange(a,b,labels = c('a', 'b'), align = 'h', legend = 'none', label.x = -0.01)
plot

ggsave('output/fig4.png', plot = plot, dpi = 600, units= 'in', height = 5, width = 5)

##################################################################################################
##################################################################################################
##################################################################################################
### Supplemental Materials -----------------------------------------------------------------------
##################################################################################################
##################################################################################################
##################################################################################################

##################################################################################################
### Boxplot --------------------------------------------------------------------------------------
##################################################################################################
keep <- c("nacheck", "model_data_all", "cnd_ts_data", "program_palette")
rm(list = setdiff(ls(), keep))
dat <- model_data_all |> rename(program = project)
summ_test <- dat |>  
      group_by(program) |> 
      summarize(
            mean = mean(comm_n_stability, na.rm = TRUE),
            median = median(comm_n_stability, na.rm = TRUE)
      )
anova_mod <- aov(comm_n_stability ~ program, data = dat)
summary(anova_mod)
par(mfrow=c(2,2)); plot(anova_mod)
oneway.test(comm_n_stability ~ program, data = dat)
posthoc <- games_howell_test(dat, comm_n_stability ~ program)
pw <- posthoc |>
      transmute(group1, group2, p.adj)
pvec <- pw |>
      mutate(comparison = paste(group1, group2, sep = "-")) |>
      select(comparison, p.adj) |>
      deframe()
cld <- multcompLetters(pvec, threshold = 0.05)
letters_df <- tibble(
      program = names(cld$Letters),
      letters = cld$Letters) |>
      mutate(program = factor(program, levels = c("MCR","PCCS","FCE","PCCC","SBC","VCR"))) |>
      arrange(program) |> 
      mutate(
            letters = chartr(
                  old = "cabd",
                  new = "abcd",
                  letters
            )
      )
letters_df

letters <- letters_df |> 
      left_join(dat) |> 
      select(program, letters, comm_n_stability) |> 
      group_by(program, letters) |> 
      summarize(max = max(comm_n_stability, na.rm = TRUE),
                .groups = 'drop') |> 
      group_by(letters) |> 
      mutate(pos = max(max))

dat |> 
      mutate(program = factor(
            program,
            levels = c("MCR", "PCCS", "FCE", "PCCC", "SBC", "VCR")
      )) |>
      ggplot(aes(x=program, y=comm_n_stability, fill = program)) +
      geom_jitter(aes(color = program), shape = 16, size = 2, width = 0.2, alpha = 1.0)+
      geom_boxplot(outlier.shape = NA, alpha = 0.35) +
      scale_fill_manual(values = program_palette) + 
      scale_color_manual(values = program_palette) +
      labs(y = "CND Stability", 
           fill = "Program",
           x = NULL) +
      theme_classic() +
      geom_text(
            data = letters, 
            aes(x = program, y = pos + 0.2, label = letters), 
            inherit.aes = FALSE, 
            vjust = 0,
            fontface = "bold",
            size = 3.5,            
            color = "black"      
      ) +
      theme(axis.text.x = element_text(face = "bold", color = "black"),
            axis.text.y = element_text(face = "bold", color = "black"),
            axis.title.x = element_text(face = "bold", color = "black"),
            axis.title.y = element_text(face = "bold", color = "black"),
            legend.position = "none",
            legend.text = element_text(face = "bold", color = "black"),
            legend.title = element_text(face = "bold", color = "black"),
            strip.text = element_text(face = "bold", color = "black"))

ggsave("output/figure-one.png", units = "in", width = 4,
       height = 4, dpi =  600)

##################################################################################################
### Simple Regression ----------------------------------------------------------------------------
##################################################################################################
glimpse(dat)
model <- lm(log1p(comm_n_stability) ~ log1p(s_rich_mean), data = dat)
summary(model)$r.squared 
summary(model)
r2 <- summary(model)$r.squared
r2

summ <- dat |>
      group_by(program) |> 
      mutate(stability = mean(comm_n_stability),
             richness  = mean(s_rich_mean))

summ_model <- lm(log1p(stability) ~ log1p(richness), data = summ)
summary(summ_model)$r.squared 
summary(summ_model)
r2_summ <- summary(summ_model)$r.squared
r2_summ

dat |>
      ggplot(aes(x = log1p(s_rich_mean), y = log1p(comm_n_stability))) +
      # ggplot(aes(x = log(s_rich_mean), y = log(comm_n_stability))) +
      geom_smooth(method = "lm", size = 1.5, color = "black", linetype = "solid", se = FALSE) +
      geom_point(aes(color = program), size = 1.5, alpha = 0.30) +
      geom_point(aes(x = log1p(richness), y = log1p(stability), color = program), size = 5, dat = summ) +
      # geom_point(aes(x = log(richness), y = log(stability), color = program), size = 5, dat = summ) +
      labs(x = "log(Species Richness)",
           y = "log(CND Stability)",
           fill = 'Program',
           color = 'Program') +
      scale_y_continuous(breaks = seq(0.25, 1.75, by = 0.5)) +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 14),
            axis.text.y = element_text(face = "bold", color = "black", size = 14),
            axis.title.x = element_text(face = "bold", color = "black", size = 16),
            axis.title.y = element_text(face = "bold", color = "black", size = 16),
            legend.position = c(0.95, 0.05),
            legend.justification = c(1, 0),
            legend.text = element_text(face = "bold", color = "black"),
            legend.title = element_text(face = "bold", color = "black"))

ggsave("output/fig2-panelb.png", units = "in", width = 4.2,
       height = 4.2, dpi =  600)

##################################################################################################
### Map   ----------------------------------------------------------------------------------------
##################################################################################################

world <- st_read('../../qgis/continent/world-continents.shp')
glimpse(world)

coords <- read_csv("local_data/LTER_Site_Coords_w_information.csv") |> 
      filter(Site %in% c("FCE", "MCR", "SBC", "VCR")) |> 
      rename(program = Site,
             y = Latitude,
             x = Longitude) |> 
      dplyr::select(program, y, x)
glimpse(coords)

pisco <- tribble(
      ~program, ~y,       ~x,
      "PCCC",   38.47409, -123.2485,
      "PCCS",   35.5000,  -121.019020
)

dt <- rbind(coords, pisco) |> 
      mutate(program = as.factor(program)) |> 
      st_as_sf(coords = c("x", "y"), crs = 4326) |> 
      st_transform(crs = 26917)

all <- ggplot() +
      geom_sf(data = world, fill = "grey", color = NA) +
      geom_sf(data = dt, color = "white", size = 5) +
      geom_sf(data = dt, aes(color = program), size = 4) +  
      annotation_scale(location = "br", width_hint = 0.3, 
                       bar_cols = c("black", "aliceblue"), text_cex = 1.0) +
      annotation_north_arrow(location = "tl", which_north = "true", 
                             style = north_arrow_fancy_orienteering(fill = c("black", "aliceblue"),
                                                                    line_col = "black" ), 
                             height = unit(2.0, "cm"), width = unit(2.0, "cm")) +
      scale_color_manual(values = program_palette, breaks = levels(dt$program)) +
      coord_sf(xlim = c(-148.7, -70),
               ylim = c(-20.0, 45.0)) +
      theme_bw() +
      theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.ticks.length = unit(0, "pt"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "aliceblue", color = NA),
            legend.title = element_blank(),
            legend.background = element_rect(fill = 'aliceblue'),
            legend.position = c(0.11,0.5),
            legend.box.background = element_rect(fill = NA, color = 'black', linewidth = 2),
            legend.text = element_text(size = 12, face = 'bold')
      )

all

ggsave("output/sitemap.png", units = "in", width = 5,
       height = 5, dpi =  600)

##################################################################################################
### Map Time Series ------------------------------------------------------------------------------
##################################################################################################

ann_dt <- cnd_ts_data |> rename(program = project)
glimpse(ann_dt)

ann_dt |> 
      filter(program == "SBC") |> 
      ggplot(aes(x = year, y = comm_n, group = site, color = program)) +
      geom_line(alpha = 0.8, linewidth = 1) +
      labs(x = 'Year',
           y = expression(bold('Aggregate Nitrogen Supply (μg '*~m^-2~""~hr^-1*')'))) +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      scale_x_continuous(breaks = c(2000,2005,2010,2015,2020)) +
      theme(
            axis.text = element_text(face = "bold", size = 12, color = "black"),
            axis.title.y = element_text(face = "bold", size = 14, color = "black"),
            axis.title.x = element_blank(),
            axis.line = element_line("black"),
            legend.position = "none",
            legend.text = element_text(face = "bold", size = 14, color = "black"),
            legend.title = element_text(face = "bold", size = 14, color = "black"),
            panel.background = element_rect(fill = "white"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold", size = 12, color = "black"))

ggsave("output/map_insets/sbc-timeseries.png", units = "in", width = 5,
       height = 5, dpi =  600)

ann_dt |> 
      filter(program == "MCR") |> 
      ggplot(aes(x = year, y = comm_n, group = site, color = program)) +
      geom_line(alpha = 0.8, linewidth = 1) +
      labs(x = 'Year',
           y = expression(bold('Aggregate Nitrogen Supply (μg '*~m^-2~""~hr^-1*')'))) +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      scale_x_continuous(breaks = c(2006,2010,2014,2018,2022)) +
      theme(
            axis.text = element_text(face = "bold", size = 12, color = "black"),
            axis.title.y = element_text(face = "bold", size = 14, color = "black"),
            axis.title.x = element_blank(),
            axis.line = element_line("black"),
            legend.position = "none",
            legend.text = element_text(face = "bold", size = 14, color = "black"),
            legend.title = element_text(face = "bold", size = 14, color = "black"),
            panel.background = element_rect(fill = "white"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold", size = 12, color = "black"))

ggsave("output/map_insets/mcr-timeseries.png", units = "in", width = 5,
       height = 5, dpi =  600)

ann_dt |> 
      filter(program == "VCR") |> 
      ggplot(aes(x = year, y = comm_n, group = site, color = program)) +
      geom_line(alpha = 0.8, linewidth = 1) +
      labs(x = 'Year',
           y = expression(bold('Aggregate Nitrogen Supply (μg '*~m^-2~""~hr^-1*')'))) +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      scale_x_continuous(breaks = c(2012,2014,2016,2018)) +
      theme(
            axis.text = element_text(face = "bold", size = 12, color = "black"),
            axis.title.y = element_text(face = "bold", size = 14, color = "black"),
            axis.title.x = element_blank(),
            axis.line = element_line("black"),
            legend.position = "none",
            legend.text = element_text(face = "bold", size = 14, color = "black"),
            legend.title = element_text(face = "bold", size = 14, color = "black"),
            panel.background = element_rect(fill = "white"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold", size = 12, color = "black"))

ggsave("output/map_insets/vcr-timeseries.png", units = "in", width = 5,
       height = 5, dpi =  600)

ann_dt|> 
      filter(program == "FCE") |> 
      ggplot(aes(x = year, y = comm_n, group = site, color = program)) +
      geom_line(alpha = 0.8, linewidth = 1) +
      labs(x = 'Year',
           y = expression(bold('Aggregate Nitrogen Supply (μg '*~m^-1~""~hr^-1*')'))) +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      scale_x_continuous(breaks = c(2005,2008,2011,2014,2017,2020,2023)) +
      theme(
            axis.text = element_text(face = "bold", size = 12, color = "black"),
            axis.title.y = element_text(face = "bold", size = 14, color = "black"),
            axis.title.x = element_blank(),
            axis.line = element_line("black"),
            legend.position = "none",
            legend.text = element_text(face = "bold", size = 14, color = "black"),
            legend.title = element_text(face = "bold", size = 14, color = "black"),
            panel.background = element_rect(fill = "white"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold", size = 12, color = "black"))

ggsave("output/map_insets/fce-timeseries.png", units = "in", width = 5,
       height = 5, dpi =  600)

ann_dt |> 
      filter(program == "PCCC") |> 
      ggplot(aes(x = year, y = comm_n, group = site, color = program)) +
      geom_line(alpha = 0.8, linewidth = 1) +
      labs(x = 'Year',
           y = expression(bold('Aggregate Nitrogen Supply (μg '*~m^-2~""~hr^-1*')'))) +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      scale_x_continuous(breaks = c(2000,2005,2010,2015,2020)) +
      theme(
            axis.text = element_text(face = "bold", size = 12, color = "black"),
            axis.title.y = element_text(face = "bold", size = 14, color = "black"),
            axis.title.x = element_blank(),
            axis.line = element_line("black"),
            legend.position = "none",
            legend.text = element_text(face = "bold", size = 14, color = "black"),
            legend.title = element_text(face = "bold", size = 14, color = "black"),
            panel.background = element_rect(fill = "white"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold", size = 12, color = "black"))

ggsave("output/map_insets/pccc-timeseries.png", units = "in", width = 5,
       height = 5, dpi =  600)

ann_dt |> 
      filter(program == "PCCS") |> 
      ggplot(aes(x = year, y = comm_n, group = site, color = program)) +
      geom_line(alpha = 0.8, linewidth = 1) +
      labs(x = 'Year',
           y = expression(bold('Aggregate Nitrogen Supply (μg '*~m^-2~""~hr^-1*')'))) +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      scale_x_continuous(breaks = c(2000,2005,2010,2015,2020)) +
      theme(
            axis.text = element_text(face = "bold", size = 12, color = "black"),
            axis.title.y = element_text(face = "bold", size = 14, color = "black"),
            axis.title.x = element_blank(),
            axis.line = element_line("black"),
            legend.position = "none",
            legend.text = element_text(face = "bold", size = 14, color = "black"),
            legend.title = element_text(face = "bold", size = 14, color = "black"),
            panel.background = element_rect(fill = "white"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold", size = 12, color = "black"))

ggsave("output/map_insets/pccs-timeseries.png", units = "in", width = 5,
       height = 5, dpi =  600)