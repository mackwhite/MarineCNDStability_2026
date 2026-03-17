###project: LTER Marine Consumer Nutrient Dynamic Synthesis Working Group
###author(s): MW, AC, LK, WRJ
###goal(s): Wrangling and summarizing raw CND data such that it is ready for analysis
###date(s): Summer 2024, Revised Spring 2026
###note(s): 
# L759-1059 commented out - undo to run within-ecosystem models

# Housekeeping ------------------------------------------------------------
### load necessary libraries ----
# install.packages("librarian")
librarian::shelf(tidyverse, vegan, readxl, splitstackshape, codyn, lavaan,
                 MuMIn, corrplot, performance, ggeffects, ggpubr, parameters, ggstats,
                 brms, mixedup, rstatix, sf, ggspatial, waldo, multcompView, tidySEM,
                 lme4, glmmTMB)

### set custom functions ----
nacheck <- function(df) {
      na_count_per_column <- sapply(df, function(x) sum(is.na(x)))
      print(na_count_per_column)
}

# Load and clean data ----------------------------------------------------
dt <- read.csv(file.path("tier2", "harmonized_consumer_excretion_CLEAN.csv"),
               stringsAsFactors = F,
               ### all NAs were inititally transformed to '.'
               na.strings =".") |> 
      
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
      
      # set upper end for California Moray eel based on reported maximum weight
      mutate(
            dmperind_g_ind = case_when(
                  dmperind_g_ind > 9071 & scientific_name == "Gymnothorax mordax" ~ 9071,
                  TRUE ~ dmperind_g_ind
            )) |> 
      
      # remove extraordinary large schools of fish [lost < 0.0001 % of data]
      # FCE doesn't catch thousands of fish per transect, so fine for this given variable area measurements (i.e., electrofishing program)
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
      filter(project == "FCE" | count <10000) |> 
      select(-area, -count) |> 
      
      # remove 'biomass buster' sharks and rays from CoastalCA, SBC, and MCR [lost < 0.01% of data]
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
      
      # coalesce density columns and remove unnecessary '..m3' column
      mutate(density_num = coalesce(density_num_m, density_num_m2)) |> 
      select(-density_num_m, -density_num_m2, -density_num_m3) |> 
      
      # remove high density of large-bodied fish observations that skew entire time series [lost < 0.0001% of data]
      filter(!(density_num > 1 & nind_ug_hr > 20000))
glimpse(dt1)
nacheck(dt1)
head(dt1)
rm(dt, dta, dtb1, dtb2, dtb, dt_ab)


# Summarize data ---------------------------------------------------------------------------------
## Calculating CND and Biomass metrics -----------------------------------------------------------
### MCR subsite_level3 has two levels (1 and 5) that should be summed across, 
### as they are conducted within same space and time so split them out and go
### sum across subsite_level2 
dt2_mcr <- dt1 |> 
      filter(project == 'MCR') |> 
      group_by(project, habitat, year, month, 
               site, subsite_level1, subsite_level2, subsite_level3, 
               scientific_name) |> 
      summarize(
            n       = sum(nind_ug_hr*density_num, na.rm = TRUE),
            bm      = sum(dmperind_g_ind*density_num, na.rm = TRUE),
            dens    = sum(density_num, na.rm = TRUE),
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
      # subsite_level3 set to '15' to denote combined transects 1 and 5
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
            species_n    = sum(nind_ug_hr*density_num, na.rm = TRUE),
            species_bm   = sum(dmperind_g_ind*density_num, na.rm = TRUE),
            species_dens = sum(density_num, na.rm = TRUE),
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
      # Inf occurs at transects with zero observed species (density = 0 rows);
      # set to 0 to indicate absence of diversity
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


## Calculating Trophic Richness and Diversity  --------------------------------------
dt2_mcr_troph <- dt1 |> 
      filter(project == 'MCR') |> 
      group_by(project, habitat, year, month, 
               site, subsite_level1, subsite_level2, subsite_level3, 
               diet_cat) |> 
      summarize(
            n       = sum(nind_ug_hr*density_num, na.rm = TRUE),
            bm      = sum(dmperind_g_ind*density_num, na.rm = TRUE),
            dens    = sum(density_num, na.rm = TRUE),
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
      # subsite_level3 set to '15' to denote combined transects 1 and 5
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
            troph_n    = sum(nind_ug_hr*density_num, na.rm = TRUE),
            troph_bm   = sum(dmperind_g_ind*density_num, na.rm = TRUE),
            troph_dens = sum(density_num, na.rm = TRUE),
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

dt3_troph <- dt2_troph |> 
      group_by(project, habitat, year, month, 
               site, subsite_level1, subsite_level2, subsite_level3) |> 
      summarize(
            trophic_richness       = vegan::specnumber(troph_dens),                 
            inv_simpson_troph      = vegan::diversity(troph_dens, index = "invsimpson"),
            total_dens_troph       = sum(troph_dens, na.rm = TRUE),
            .groups = "drop"
      ) |> 
      # Inf occurs at transects with zero observed trophic groups (density = 0 rows);
      # set to 0 to indicate absence of diversity
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

## Join trophic and species level data together ---------------------------------------
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


## Calculate Species Dynamics -----------------------------------------------------------------
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
                  mean(beta_temp$total, na.rm = TRUE)
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


## Calculate Trophic Dynamics -----------------------------------------------------------------
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

## Join trophic and species level dynamics together -----------------------------------------------
glimpse(cnd_model_data)
glimpse(spp_dyn_model_data)
glimpse(troph_dyn_model_data)

## Join all data to create final model dataset ----------------------------------------------------
model_data_all <- cnd_model_data |> 
      left_join(spp_dyn_model_data) |> 
      left_join(troph_dyn_model_data)
glimpse(model_data_all)
head(model_data_all)
nacheck(model_data_all)
# write_csv(model_data_all, "local_data/model-data-all.csv")


# Run across system model ------------------------------------------------------------------------
keep <- c("nacheck", "model_data_all", "cnd_ts_data")
rm(list = setdiff(ls(), keep))

glimpse(model_data_all)
dat_scaled <- model_data_all |> 
      rename(program = project) |> 
      select(program, site, comm_n_stability, everything()) |> 
      
      # scale response metric
      mutate(comm_n_stability = as.numeric(scale(comm_n_stability))) |>
      
      # scale suite of predictor metrics 
      mutate(across(comm_n_mean:troph_synchrony, \(x) as.numeric(scale(x, center = TRUE))))
glimpse(dat_scaled)      
dat_ready <- dat_scaled      
glimpse(dat_ready)

## Check skewness and kurtosis of predictor metrics nested in path model -------------------------
# (Kline 2016; West, Finch & Curran 1995)
dat_ready |>
      select(comm_n_stability, spp_synchrony, 
             spp_turnover, troph_turnover) |>
      summarise(across(everything(), list(
            skew = skewness,
            kurt = kurtosis
      ))) |>
      pivot_longer(
            everything(),
            names_to  = c("variable", "metric"),
            names_sep = "_(?=[^_]+$)",
            values_to = "value"
      ) |>
      pivot_wider(
            names_from  = metric,
            values_from = value
      )

## Review correlation structure ------------------------------------------------------------------
num_vars <- dat_ready |> 
      dplyr::select(
            comm_n_stability,
            s_rich_mean, s_div_mean, spp_turnover, spp_synchrony,
            t_rich_mean, t_div_mean, troph_turnover, troph_synchrony
      )

m <- cor(num_vars, use = "complete.obs")
corrplot(m, 
         method = "circle",      
         type = "upper",         
         addCoef.col = "black",  
         tl.col = "black",       
         diag = FALSE)

## Specify model paths -------------------------------------------------------------------
path_model <- '
   # Regressions
   comm_n_stability ~ cp*s_rich_mean + b1*spp_synchrony + b2*spp_turnover + b3*troph_turnover
   
   spp_synchrony    ~ a1*s_rich_mean
   
   spp_turnover     ~ a3*s_rich_mean + a2*t_rich_mean
   
   troph_turnover   ~ a4*t_rich_mean

   # Covariances
   s_rich_mean      ~~ t_rich_mean
   spp_synchrony    ~~ troph_turnover
   troph_turnover   ~~ s_rich_mean
   spp_turnover     ~~ troph_turnover

   # Indirect Effects 
   ### indirect effect of species richness through species synchrony
   ind_srich_ssync := a1 * b1
   
   ### indirect effect of species richness through species turnover
   ind_srich_sturn := a3 * b2
   ind_trich_sturn := a2 * b2
   ### indirect effect of trophic richness through trophic turnover
   ind_trich_tturn := a4 * b3
   
   ### cumulative effect of species richness (direct and indirect paths)
   total_srich_effect := cp + ind_srich_ssync + ind_srich_sturn
'

## Examine model fit, summarize, and visualize -------------------------------------------------------
# MLR (Maximum Likelihood with Robust standard errors) estimator used for robustness to non-normality in scaled continuous outcomes
fit <- sem(path_model, data = dat_ready, estimator = "MLR")
summary(fit, standardized = TRUE, fit.measures = TRUE)
# no indices exceeded mi = 10, supporting retention of the - indicates no need to reform model structure
modindices(fit, sort = TRUE, maximum.number = 10)
# extract R-squared values for all endogenous variables
lavInspect(fit, "r2")

# visualize model here, but fit in .ppt for paper
lay <- get_layout(
      "s_rich_mean",    "",                  "t_rich_mean",               
      "spp_synchrony",  "spp_turnover",      "troph_turnover", 
      "",               "comm_n_stability",  "",             
      rows = 3
)

graph <- prepare_graph(
      model  = fit,
      layout = lay,
      label  = "est_sig_std"
)

edges <- graph$edges
edges <- edges |>
      filter(op != "~~" | lhs == rhs)
graph$edges <- edges

graph <- edit_edges(graph,
                    linetype  = "solid",
                    linewidth = ifelse(abs(est) > 0.3, 1.5, 0.8)
)

graph <- edit_nodes(graph,
                    fill  = "white",
                    color = "black"
)

p <- plot(graph) +
      theme_void() +
      theme(
            plot.background = element_rect(fill = "white", color = NA),
            plot.margin     = margin(20, 20, 20, 20)
      )
p


# Run within system model ------------------------------------------------------------------------
#
#
# keep <- c("nacheck", "model_data_all", "cnd_ts_data", "fit")
# rm(list = setdiff(ls(), keep))
# 
# glimpse(model_data_all)
# dat_scaled <- model_data_all |> 
#       rename(program = project) |> 
#       select(program, site, comm_n_stability, everything()) |> 
#       
#       # comm_n_stability scaled globally (outcome); predictors scaled within-program
#       mutate(comm_n_stability = as.numeric(scale(comm_n_stability))) |>
#       
#       # scaling predictors within program to control for between-program differences (i.e., across vs. within)
#       # in baseline community metrics; coefficients reflect within-program
#       # standardized effects
#       group_by(program) |> 
#       mutate(across(comm_n_mean:troph_synchrony, \(x) as.numeric(scale(x, center = TRUE)))) |> 
#       ungroup()
# glimpse(dat_scaled)      
# dat_ready <- dat_scaled      
# glimpse(dat_ready)
# 
## Check skewness and kurtosis of predictor metric in brms model (Kline 2016; West, Finch & Curran 1995) ----
# dat_ready |>
#       select(comm_n_stability) |>
#       summarise(across(everything(), list(
#             skew = skewness,
#             kurt = kurtosis
#       ))) |>
#       pivot_longer(
#             everything(),
#             names_to  = c("variable", "metric"),
#             names_sep = "_(?=[^_]+$)",  # splits on last underscore
#             values_to = "value"
#       ) |>
#       pivot_wider(
#             names_from  = metric,
#             values_from = value
#       )
# 
# # normal(0,1) weakly informative prior appropriate for standardized predictors - following Lemoine (2019, Ecology)
# pr = prior(normal(0, 1), class = 'b')
# 
#
## Review correlation structure --------------------------------------------------------------------
# test_corr <- dat_ready |> select(s_rich_mean, s_div_mean,
#                                  t_rich_mean, t_div_mean,
#                                  spp_turnover, spp_synchrony,
#                                  troph_turnover, troph_synchrony)
# 
# matrix <- cor(test_corr, use = 'complete.obs')
# 
# corrplot(matrix, method = "number", type = "lower", tl.col = "black", tl.srt = 45)
# glimpse(dat_ready)
# 
## Fit models using forward selection --------------------------------------------------------------
### round one ------------------------------------------------------------------------------------
#
# 
# ### round one: single-term models to identify best individual predictor
# m1 <- brm(
#       comm_n_stability ~ t_rich_mean + (t_rich_mean | program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4,
#       seed = 20
# )
# 
# m4 <- brm(
#       comm_n_stability ~ spp_synchrony + (spp_synchrony | program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4,
#       seed = 20
# )
# 
# ### best fit single term model
# # saveRDS(m4, file = 'local_data/rds-single-synchrony.rds')
# 
# m5 <- brm(
#       comm_n_stability ~ spp_turnover + (spp_turnover | program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4,
#       seed = 20
# )
# 
# m6 <- brm(
#       comm_n_stability ~ troph_synchrony + (troph_synchrony | program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4,
#       seed = 20
# )
# 
# m7 <- brm(
#       comm_n_stability ~ troph_turnover + (troph_turnover | program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4,
#       seed = 20
# )
# 
# m8 <- brm(
#       comm_n_stability ~ s_rich_mean + (s_rich_mean | program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4,
#       seed = 20
# )
# 
# ### model of interest, given across-ecosystem importance
# # saveRDS(m8, file = 'local_data/rds-single-richness.rds')
# 
# m9 <- brm(
#       comm_n_stability ~ s_div_mean + (s_div_mean | program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4,
#       seed = 20
# )
# 
# m10 <- brm(
#       comm_n_stability ~ t_div_mean + (t_div_mean | program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4,
#       seed = 20
# )
# 
# model_table_all <- performance::compare_performance(m1,m4,m5,m6,m7,m8,m9,m10)
# 
# model_selection1 <- model_table_all |>
#       mutate(dWAIC = WAIC - min(WAIC))
# 
# # write_csv(model_selection1, "output/tables/brms-fullmodel-selection-table-roundone.csv")
# 
# keep <- c("nacheck", "model_data_all", "cnd_ts_data", "fit", "dat_ready", "pr", 'm4', 'model_selection1')
# rm(list = setdiff(ls(), keep))
# 
#
### round two ------------------------------------------------------------------------------------
#
# ### round two: add spp_synchrony (best round one model) to all other predictors
# m41 <- brm(
#       comm_n_stability ~ t_rich_mean + spp_synchrony + (t_rich_mean + spp_synchrony | program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4,
#       seed = 20
# )
# 
# m45 <- brm(
#       comm_n_stability ~ spp_turnover + spp_synchrony + (spp_turnover + spp_synchrony | program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4,
#       seed = 20
# )
# 
# m47 <- brm(
#       comm_n_stability ~ troph_turnover + spp_synchrony + (troph_turnover + spp_synchrony | program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4,
#       seed = 20
# )
# 
# m48 <- brm(
#       comm_n_stability ~ s_rich_mean + spp_synchrony + (s_rich_mean + spp_synchrony | program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4,
#       seed = 20
# )
# 
# m49 <- brm(
#       comm_n_stability ~ t_div_mean + spp_synchrony + (t_div_mean + spp_synchrony | program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4,
#       seed = 20
# )
# 
# m40 <- brm(
#       comm_n_stability ~ s_div_mean + spp_synchrony + (s_div_mean + spp_synchrony | program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4,
#       seed = 20
# )
# 
# model_table_all <- performance::compare_performance(m41,m4,m45,m47,m48,m49,m40)
# 
# model_selection2 <- model_table_all |>
#       mutate(dWAIC = WAIC - min(WAIC))
# write_csv(model_selection2, "output/tables/brms-fullmodel-selection-table-roundtwo.csv")
# 
# keep <- c("nacheck", "model_data_all", "cnd_ts_data", "fit", "dat_ready", "pr", 'm45', 'model_selection1', 'model_selection2')
# rm(list = setdiff(ls(), keep))
# 
#
### round three ----------------------------------------------------------------------------------
#
#
# ### round three: add spp_turnover (best round two model) to spp_synchrony + one additional predictor
# 
# m451 <- brm(
#       comm_n_stability ~ t_rich_mean + spp_turnover + spp_synchrony + (t_rich_mean + spp_turnover + spp_synchrony | program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4,
#       seed = 20
# )
# 
# m457 <- brm(
#       comm_n_stability ~ troph_turnover + spp_turnover + spp_synchrony + (troph_turnover + spp_turnover + spp_synchrony | program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4,
#       seed = 20
# )
# 
# m458 <- brm(
#       comm_n_stability ~ s_rich_mean + spp_turnover + spp_synchrony + (s_rich_mean + spp_turnover + spp_synchrony | program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4,
#       seed = 20
# )
# 
# m459 <- brm(
#       comm_n_stability ~ t_div_mean + spp_turnover + spp_synchrony + (t_div_mean + spp_turnover + spp_synchrony | program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4,
#       seed = 20
# )
# 
# m450 <- brm(
#       comm_n_stability ~ s_div_mean + spp_turnover + spp_synchrony + (s_div_mean + spp_turnover + spp_synchrony | program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4,
#       seed = 20
# )
# 
# model_table_all <- performance::compare_performance(m45,m451,m457,m458,m459,m450)
# 
# model_selection3 <- model_table_all |>
#       mutate(dWAIC = WAIC - min(WAIC))
# write_csv(model_selection3, "output/tables/brms-fullmodel-selection-table-roundthree.csv")
# 
# keep <- c("nacheck", "model_data_all", "cnd_ts_data", "fit", "dat_ready", "pr", 'm45', 'model_selection1', 'model_selection2', 'model_selection3')
# rm(list = setdiff(ls(), keep))
# full_model <- m45
# 
# saveRDS(full_model, file = 'local_data/rds-full-model.rds')


## Examine model fit, summarize, and visualize --------------------------------------------------
### single term models --------------------------------------------------------------------------

#### read in necessary models ----
synch = readRDS("local_data/rds-single-synchrony.rds")
rich = readRDS("local_data/rds-single-richness.rds")
glimpse(dat_ready)
performance(synch)
performance(rich)

#### summary stats ----
p_synch = as_draws_df(synch)
glimpse(p_synch)
mean(p_synch$`r_program[PCCC,spp_synchrony]` + p_synch$b_spp_synchrony  < 0)
mean(p_synch$`r_program[SBC,spp_synchrony]`  + p_synch$b_spp_synchrony  < 0)
mean(p_synch$`r_program[VCR,spp_synchrony]`  + p_synch$b_spp_synchrony  < 0)
mean(p_synch$`r_program[PCCS,spp_synchrony]` + p_synch$b_spp_synchrony  < 0)
mean(p_synch$`r_program[MCR,spp_synchrony]`  + p_synch$b_spp_synchrony  < 0)
mean(p_synch$`r_program[FCE,spp_synchrony]`  + p_synch$b_spp_synchrony  < 0)

p_sr = as_draws_df(rich)
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
      # extract min and max predictor values per program
      # sufficient for linear prediction lines
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
             # back-transform predictions to raw stability scale
             # reverses z-scoring: pred * SD + mean
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


### forward selected model -----------------------------------------------------
#### read in necessary models --------------------------------------------------
full_model = readRDS('local_data/rds-full-model.rds')
performance(full_model)

#### summary stats -----
post = as_draws_df(full_model)
mean(post$`r_program[MCR,spp_synchrony]`+ post$b_spp_synchrony < 0)
mean(post$`r_program[MCR,spp_turnover]` + post$b_spp_turnover < 0)

mean(post$`r_program[PCCS,spp_synchrony]`+ post$b_spp_synchrony < 0)
mean(post$`r_program[PCCS,spp_turnover]` + post$b_spp_turnover < 0)

mean(post$`r_program[FCE,spp_synchrony]`+ post$b_spp_synchrony < 0)
mean(post$`r_program[FCE,spp_turnover]` + post$b_spp_turnover < 0)

mean(post$`r_program[PCCC,spp_synchrony]`+ post$b_spp_synchrony < 0)
mean(post$`r_program[PCCC,spp_turnover]` + post$b_spp_turnover < 0)

mean(post$`r_program[SBC,spp_synchrony]`+ post$b_spp_synchrony < 0)
mean(post$`r_program[SBC,spp_turnover]` + post$b_spp_turnover < 0)

mean(post$`r_program[VCR,spp_synchrony]`+ post$b_spp_synchrony < 0)
mean(post$`r_program[VCR,spp_turnover]` + post$b_spp_turnover < 0)

# random effects
re95 = mixedup::extract_random_coefs(full_model, ci_level = c(0.95))
re80 = mixedup::extract_random_coefs(full_model, ci_level = c(0.8))

re_beta = left_join(re95, re80) |> 
      rename(term = effect,
             program = group)

# fixed effects
fe95 = mixedup::extract_fixed_effects(full_model, ci_level = c(0.95))
fe80 = mixedup::extract_fixed_effects(full_model, ci_level = c(0.8))

fe_beta = left_join(fe95, fe80) |> 
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


# Additional Analyses and Summary ----------------------------------------------
## Figure One Map --------------------------------------------------------------
keep <- c("nacheck", "model_data_all", "cnd_ts_data", "program_palette")

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


## Figure Two Boxplot and ANOVA -----------------------------------------------------------------
keep <- c("nacheck", "model_data_all", "cnd_ts_data", "program_palette")
rm(list = setdiff(ls(), keep))
dat <- model_data_all |> rename(program = project)
summ_test <- dat |>  
      group_by(program) |> 
      summarize(
            mean = mean(comm_n_stability, na.rm = TRUE),
            sd = sd(comm_n_stability, na.rm = TRUE),
            median = median(comm_n_stability, na.rm = TRUE),
            .groups = 'drop'
      )

# test for homogeneity of variance across programs
# significant result justifies Welch's ANOVA and Games-Howell post-hoc
car::leveneTest(comm_n_stability ~ factor(program), data = dat)
oneway.test(comm_n_stability ~ program, data = dat)

# Games-Howell post-hoc — appropriate for unequal variances and sample sizes
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
      # remap CLD letters to alphabetical order for readability
      # WARNING: if data changes, verify cld$Letters output before updating chartr mapping
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


## Figure Three Part B Simple Regression [CND Stability ~ Richness] ------------
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
      geom_smooth(method = "lm", size = 1.5, color = "black", linetype = "solid", se = FALSE) +
      geom_point(aes(color = program), size = 1.5, alpha = 0.30) +
      geom_point(aes(x = log1p(richness), y = log1p(stability), color = program), size = 5, dat = summ) +
      labs(x = "log(Species Richness + 1)",
           y = "log(CND Stability + 1)",
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


# Supporting Materials ----------------------------------------------------
## Supplemental Figure One Model Validation Regression -------------------------
keep <- c("nacheck", "model_data_all", "cnd_ts_data", "program_palette")
rm(list = setdiff(ls(), keep))

emp <- read_csv('local_data/empirical_excretion_kelp.csv') |> 
      dplyr::rename(order = TAXON_ORDER,
                    family = TAXON_FAMILY,
                    genus = TAXON_GENUS,
                    species = TAXON_SPECIES,
                    wetmass_g = WM_g,
                    diet_cat = `Functional group`,
                    nind_umol_hr = `exc_rate_NH4_umol_hour-1`,
                    region = Region) |> 
      mutate(scientific_name = paste(genus, species, sep = " "),
             nind_ug_hr = 18.042*nind_umol_hr) |> 
      dplyr::select(scientific_name, wetmass_g, nind_ug_hr, nind_umol_hr, order, family, genus, species, diet_cat, region)
glimpse(emp)

kelp_diet <- read_csv('local_data/kelp_empirical_species_list_updated.csv')
glimpse(kelp_diet)      
kelp <- emp |> dplyr::select(-diet_cat)
glimpse(kelp)
kelp_all <- kelp |> left_join(kelp_diet) |> distinct()

### clean up environment 
keep <- c("nacheck", "model_data_all", "cnd_ts_data", "program_palette", "kelp_all")
rm(list = setdiff(ls(), keep))

kelp1 <- kelp_all |> 
      rename(nind_ug_hr_emp = nind_ug_hr,
             nind_umol_hr_emp = nind_umol_hr) |> 
      mutate(drymass_g = wetmass_g*dm_conv,
             phylum = "Chordata") |> 
      mutate(temp = case_when(
            region == "Central California" ~ 13,
            region == "Southern California" ~ 17
      ))
glimpse(kelp1)

kelp2 <- kelp1 |> 
      ### vertebrate coefficient classification
      mutate(vert_coef = if_else(phylum == "Chordata", 0.7804, 0),
             vert_coef_upper = if_else(phylum == "Chordata", 0.7804 + 0.0655, 0),
             vert_coef_lower = if_else(phylum == "Chordata", 0.7804 - 0.0655, 0)) |> 
      ### diet coefficient classification
      mutate(diet_coef = case_when(
            diet_cat == "algae_detritus" ~ -0.0389,
            diet_cat == "invert" ~ -0.2013,
            diet_cat == "fish" ~ -0.0537,
            diet_cat == "fish_invert" ~ -0.1732,
            diet_cat == "algae_invert" ~ 0,
            TRUE ~ NA)) |> 
      mutate(diet_coef_upper = case_when(
            diet_cat == "algae_detritus" ~ -0.0389 + 0.0765,
            diet_cat == "invert" ~ -0.2013 + 0.0771,
            diet_cat == "fish" ~ -0.0537 + 0.2786,
            diet_cat == "fish_invert" ~ -0.1732 + 0.1384,
            diet_cat == "algae_invert" ~ 0,
            TRUE ~ NA)) |> 
      mutate(diet_coef_lower = case_when(
            diet_cat == "algae_detritus" ~ -0.0389 - 0.0765,
            diet_cat == "invert" ~ -0.2013 - 0.0771,
            diet_cat == "fish" ~ -0.0537 - 0.2786,
            diet_cat == "fish_invert" ~ -0.1732 - 0.1384,
            diet_cat == "algae_invert" ~ 0,
            TRUE ~ NA)) |> 
      ### temperature coefficient classification
      mutate(temp_coef = 0.0246,
             temp_coef_upper = 0.0246 + 0.0014,
             temp_coef_lower = 0.0246 - 0.0014) |> 
      ### dry mass coefficient classification
      mutate(dm_coef = 0.6840,
             dm_coef_upper = 0.6840 + 0.0177,
             dm_coef_lower = 0.6840 - 0.0177) |> 
      ### intercept coefficient classification
      mutate(int_coef = 1.4610,
             int_coef_upper = 1.4610 + 0.0897,
             int_coef_lower = 1.4610 - 0.0897)
glimpse(kelp2)

kelp3 <- kelp2 |> 
      mutate(n10 = int_coef + dm_coef*(log10(drymass_g)) + temp_coef*temp + diet_coef + vert_coef,
             n10_lower = int_coef_lower + dm_coef_lower*(log10(drymass_g)) + temp_coef_lower*temp + diet_coef_lower + vert_coef_lower,
             n10_upper = int_coef_upper + dm_coef_upper*(log10(drymass_g)) + temp_coef_upper*temp + diet_coef_upper + vert_coef_upper) |> 
      mutate(nind_ug_hr = 10^n10,
             nind_ug_hr_lower = 10^n10_lower,
             nind_ug_hr_upper = 10^n10_upper)

kelp4 <- kelp3 |> 
      dplyr::select(scientific_name, wetmass_g, drymass_g, diet_cat,
                    nind_ug_hr_emp, nind_ug_hr, nind_ug_hr_lower, nind_ug_hr_upper,
                    phylum, order, family, genus, species) |> 
      mutate(
            n10_emp = log10(nind_ug_hr_emp),
            n10_mod = log10(nind_ug_hr),
            n10_mod_low = log10(nind_ug_hr_lower),
            n10_mod_upp = log10(nind_ug_hr_lower),
            n10_dm = log10(drymass_g)
      ) |> 
      group_by(scientific_name) |> mutate(n = n()) |> ungroup() |> 
      filter(n > 2) |> 
      dplyr::select(-n) |> 
      rename(diet = diet_cat) |> 
      group_by(scientific_name) |> mutate(species_n = n()) |> ungroup() |> 
      group_by(genus) |> mutate(genus_n = n()) |> ungroup() |> 
      group_by(family) |> mutate(family_n = n()) |> ungroup()

mod <- glmmTMB(
      n10_emp ~ n10_mod, family = gaussian(), data = kelp4
)
summary(mod)
performance::performance(mod)

kelp4 |> 
      ggplot(aes(x = n10_mod, y = n10_emp)) +
      geom_point(size = 2, alpha = 0.2, color = '#2A788EFF') +
      geom_smooth(method = lm, color = "black", fill = "black") +
      theme_classic() +
      scale_y_continuous(breaks = c(1,2,3,4,5), limits = c(1.13,5.233)) +
      scale_x_continuous(breaks = c(1,2,3,4,5), limits = c(1.13,5.233)) +
      ylab(expression(bold("Empirical Log" [10] * " N Excretion (" * mu * "g" %.% ind^-1 %.% hr^-1 * ")"))) +
      xlab(expression(bold("Modeled Log" [10] * " N Excretion (" * mu * "g" %.% ind^-1 %.% hr^-1 * ")"))) +
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 14),
            axis.text.y = element_text(face = "bold", color = "black", size = 14),
            axis.title.x = element_text(face = "bold", color = "black", size = 14),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            legend.position = "none",
            legend.text = element_text(face = "bold", color = "black"),
            legend.title = element_text(face = "bold", color = "black"))

ggsave("output/smf1-model-validation.png", units = "in", width = 8,
       height = 4, dpi = 600)

## Supplemental Figure Two POR Effect on CND Stability Regression --------------
keep <- c("nacheck", "model_data_all", "cnd_ts_data", "program_palette")
rm(list = setdiff(ls(), keep))

summary <- cnd_ts_data |> 
      group_by(project, site) |> 
      summarize(years = n_distinct(year)) |> 
      rename(program = project,
             site = site)

dat_summary <- model_data_all |> 
      dplyr::select(site, comm_n_stability)

all <- left_join(dat_summary, summary, by = "site")  
model <- lm(comm_n_stability ~ years, data = all)
summary(model)$r.squared 
r2 <- summary(model)$r.squared
summary(model)

a <- all |>
      # filter(Program != "VCR") |> 
      ggplot(aes(x = years, y = comm_n_stability)) +
      geom_point(aes(color = program), size = 2) +  # Adds the scatter plot points
      geom_smooth(method = "lm", size = 2, color = "black", linetype = "solid", se = FALSE) +
      labs(x = "Period of Record (years)",
           y = "CND Stability",
           color = 'Program') +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      theme(axis.text.x = element_text(face = "bold", color = "black"),
            axis.text.y = element_text(face = "bold", color = "black"),
            axis.title.x = element_text(face = "bold", color = "black"),
            axis.title.y = element_text(face = "bold", color = "black"),
            legend.position = c(0.15,0.75),
            legend.text = element_text(face = "bold", color = "black"),
            legend.title = element_text(face = "bold", color = "black"))
a

all_short <- all |> filter(program != "VCR")
model_short <- lm(comm_n_stability ~ years, data = all_short)
summary(model_short)$r.squared 
r2_short <- summary(model_short)$r.squared
summary(model_short)

b <- all |>
      filter(program != "VCR") |>
      ggplot(aes(x = years, y = comm_n_stability)) +
      geom_point(aes(color = program), size = 2) +
      geom_smooth(method = "lm", size = 2, color = "black", linetype = "dashed", se = FALSE) +
      labs(x = "Period of Record (years)",
           y = "CND Stability") +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      theme(axis.text.x = element_text(face = "bold", color = "black"),
            axis.text.y = element_text(face = "bold", color = "black"),
            axis.title.x = element_text(face = "bold", color = "black"),
            axis.title.y = element_text(face = "bold", color = "black"),
            legend.position = "none",
            legend.text = element_text(face = "bold", color = "black"),
            legend.title = element_text(face = "bold", color = "black"))
b
ggarrange(a, b,
          labels = c('a)','b)'),
          ncol = 2, vjust = 1.3, align = "h")

ggsave("output/smf2-por-effect.png", units = "in", width = 8,
       height = 4, dpi = 600)