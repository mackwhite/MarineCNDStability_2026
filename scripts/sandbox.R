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

# # Load and clean data ----------------------------------------------------
# dt <- read.csv(file.path("tier2", "harmonized_consumer_excretion_CLEAN.csv"),
#                stringsAsFactors = F,
#                ### all NAs were inititally transformed to '.'
#                na.strings =".") |> 
#       
#       # tidy up column names
#       janitor::clean_names() |> 
#       
#       # replace NAs with character data
#       mutate(subsite_level1 = replace_na(subsite_level1, "Not Available"),
#              subsite_level2 = replace_na(subsite_level2, "Not Available"),
#              subsite_level3 = replace_na(subsite_level3, "Not Available")) |> 
#       
#       # filter out projects we are using in this manuscript
#       filter(project %in% c('MCR', 'CoastalCA', 'SBC', 'VCR', 'FCE')) |> 
#       
#       # filter out beach habitat [talitrids at SBC]
#       filter(habitat %in% c('estuary', 'ocean')) |> 
#       
#       # filter out sites not sampled consistently at FCE
#       filter(
#             !(project == "FCE" & site == "TB" & subsite_level1 == "5"),
#             !(project == "FCE" & site == "RB" & subsite_level1 %in% c("17", "19"))
#       )
# glimpse(dt)
# nacheck(dt)
# 
# ### add phylum where necessary and update excretion values for those taxa
# dta <- dt |> filter(!is.na(phylum))
# dtb1 <- dt |> filter(is.na(phylum), density_num_m2 == 0) |> 
#       mutate(phylum = 'Chordata')
# dtb2 <- dt |> filter(is.na(phylum), density_num_m2 > 0) |> 
#       mutate(phylum = 'Chordata') |> 
#       mutate(n_vert_coef = if_else(phylum == "Chordata", 0.7804, 0),
#              n_diet_coef = if_else(diet_cat == "algae_detritus", -0.0389,
#                                    if_else(diet_cat == "invert", -0.2013,
#                                            if_else(diet_cat == "fish", -0.0537,
#                                                    if_else(diet_cat == "fish_invert", -0.1732, 
#                                                            if_else(diet_cat == "algae_invert", 0,
#                                                                    NA))))),
#              nexc_log10  = ifelse(dmperind_g_ind > 0, 1.461 + 0.6840*(log10(dmperind_g_ind)) + 0.0246*temp_c + n_diet_coef + n_vert_coef,NA),
#              nind_ug_hr  = 10^nexc_log10,
#              nind_ug_hr  = ifelse(is.na(nind_ug_hr),0,nind_ug_hr)) |> 
#       mutate(p_vert_coef = if_else(phylum == "Chordata", 0.7504, 0),
#              p_diet_coef = if_else(diet_cat == "algae_detritus", 0.0173,
#                                    if_else(diet_cat == "invert", -0.2480,
#                                            if_else(diet_cat == "fish", -0.0337,
#                                                    if_else(diet_cat == "fish_invert", -0.4525, 
#                                                            if_else(diet_cat == "algae_invert",0,
#                                                                    NA))))),
#              pexc_log10  = ifelse(dmperind_g_ind >0, 0.6757 + 0.5656*(log10(dmperind_g_ind)) + 0.0194*temp_c + p_diet_coef + p_vert_coef, NA),
#              pind_ug_hr  = 10^pexc_log10,
#              pind_ug_hr  = ifelse(is.na(pind_ug_hr),0,pind_ug_hr)) |> 
#       select(-n_vert_coef, -n_diet_coef, -nexc_log10, 
#              -p_vert_coef, -p_diet_coef, -pexc_log10)
# dtb <- rbind(dtb1, dtb2)
# dt_ab <- rbind(dta, dtb)      
# 
# dt1 <- dt_ab |>
#       
#       # filter out organisms that are not fish
#       mutate(
#             order = case_when(
#                   is.na(order) ~ 'missing',
#                   TRUE ~ order
#             )) |>
#       filter(phylum == 'Chordata',
#              order != 'Decapoda') |> 
#       
#       # set upper end for California Moray eel based on reported maximum weight
#       mutate(
#             dmperind_g_ind = case_when(
#                   dmperind_g_ind > 9071 & scientific_name == "Gymnothorax mordax" ~ 9071,
#                   TRUE ~ dmperind_g_ind
#             )) |> 
#       
#       # remove extraordinary large schools of fish [lost < 0.0001 % of data]
#       # FCE doesn't catch thousands of fish per transect, so fine for this given variable area measurements (i.e., electrofishing program)
#       mutate(
#             area = case_when(
#                   project == 'CoastalCA' ~ 60,
#                   project == 'SBC' ~ 80,
#                   project == 'VCR' ~ 25,
#                   project == 'MCR' & subsite_level3 == '1' ~ 50,
#                   project == 'MCR' & subsite_level3 == '5' ~ 250,
#                   TRUE ~ NA_real_
#             ),
#             count = area*density_num_m2
#       ) |> 
#       filter(project == "FCE" | count <10000) |> 
#       select(-area, -count) |> 
#       
#       # remove 'biomass buster' sharks and rays from CoastalCA, SBC, and MCR [lost < 0.01% of data]
#       group_by(project, habitat) |> 
#       mutate(
#             mean_dmperind = mean(dmperind_g_ind, na.rm = TRUE),
#             sd_dmperind   = sd(dmperind_g_ind, na.rm = TRUE),  
#             lower_bound   = mean_dmperind - 5 * sd_dmperind,  
#             upper_bound   = mean_dmperind + 5 * sd_dmperind,
#             outlier       = dmperind_g_ind < lower_bound | dmperind_g_ind > upper_bound,
#             sharkray      = grepl("\\bshark\\b|\\bray\\b", common_name, ignore.case = TRUE),
#             elasmo        = class %in% c("Chondrichthyes", "Elasmobranchii")
#       ) |> 
#       ungroup() |> 
#       filter(!(outlier & (sharkray | elasmo))) |> 
#       select(-mean_dmperind, -sd_dmperind, -lower_bound, -upper_bound, -outlier, -sharkray, -elasmo) |> 
#       
#       # coalesce density columns and remove unnecessary '..m3' column
#       mutate(density_num = coalesce(density_num_m, density_num_m2)) |> 
#       select(-density_num_m, -density_num_m2, -density_num_m3) |> 
#       
#       # remove high density of large-bodied fish observations that skew entire time series [lost < 0.0001% of data]
#       filter(!(density_num > 1 & nind_ug_hr > 20000))
# glimpse(dt1)
# nacheck(dt1)
# head(dt1)
# rm(dt, dta, dtb1, dtb2, dtb, dt_ab)
# 
# 
# # Summarize data ---------------------------------------------------------------------------------
# ## Calculating CND and Biomass metrics -----------------------------------------------------------
# ### MCR subsite_level3 has two levels (1 and 5) that should be summed across, 
# ### as they are conducted within same space and time so split them out and go
# ### sum across subsite_level2 
# dt2_mcr <- dt1 |> 
#       filter(project == 'MCR') |> 
#       group_by(project, habitat, year, month, 
#                site, subsite_level1, subsite_level2, subsite_level3, 
#                scientific_name) |> 
#       summarize(
#             n       = sum(nind_ug_hr*density_num, na.rm = TRUE),
#             bm      = sum(dmperind_g_ind*density_num, na.rm = TRUE),
#             dens    = sum(density_num, na.rm = TRUE),
#             .groups = 'drop'
#       ) |> 
#       group_by(project, habitat, year, month, 
#                site, subsite_level1, subsite_level2, 
#                scientific_name) |> 
#       summarize(
#             species_n    = sum(n),
#             species_bm   = sum(bm),
#             species_dens = sum(dens),
#             .groups = 'drop'
#       ) |> 
#       # subsite_level3 set to '15' to denote combined transects 1 and 5
#       mutate(subsite_level3 = '15') |> 
#       select(project, habitat, year, month, site, subsite_level1, subsite_level2, subsite_level3,
#              scientific_name, species_n, species_bm, species_dens)
# 
# glimpse(dt2_mcr)
# head(dt2_mcr)
# nacheck(dt2_mcr)
# 
# ### all other sites should be summed down subsite_level3 as these are
# ### considered the 'transect'
# dt2_other <- dt1 |> 
#       filter(project != 'MCR') |> 
#       group_by(project, habitat, year, month, 
#                site, subsite_level1, subsite_level2, subsite_level3, 
#                scientific_name) |> 
#       ### sum across unique taxa at the transect level
#       ### coalesces multiple species observation to single n and bm value for each transect
#       summarize(
#             species_n    = sum(nind_ug_hr*density_num, na.rm = TRUE),
#             species_bm   = sum(dmperind_g_ind*density_num, na.rm = TRUE),
#             species_dens = sum(density_num, na.rm = TRUE),
#             .groups = 'drop'
#       ) |> 
#       select(project, habitat, year, month, site, subsite_level1, subsite_level2, subsite_level3,
#              scientific_name, species_n, species_bm, species_dens)
# 
# glimpse(dt2_other)
# head(dt2_other)
# nacheck(dt2_other)
# 
# ### join mcr and other data back together now that they are at same 'transect' scale
# dt2 <- rbind(dt2_other, dt2_mcr)
# glimpse(dt2)
# head(dt2)
# nacheck(dt2)
# rm(dt2_mcr, dt2_other)
# 
# ### sum across species at the transect scale to get community-level excretion and biomass
# ### per unit area
# dt3_cnd <- dt2 |> 
#       group_by(project, habitat, year, month, 
#                site, subsite_level1, subsite_level2, subsite_level3) |> 
#       summarize(
#             comm_n    = sum(species_n, na.rm = TRUE),
#             comm_bm   = sum(species_bm, na.rm = TRUE),
#             comm_dens = sum(species_dens, na.rm = TRUE),
#             .groups   = 'drop'
#       ) |> 
#       mutate(project = case_when(
#             project  == 'CoastalCA' & site == 'CENTRAL' ~ 'PCCC',
#             project  == 'CoastalCA' & site == 'SOUTH' ~ 'PCCS',
#             TRUE ~ project
#       )) |> 
#       ### take everything to the true 'site' level in which transects should be
#       ### averaged across
#       mutate(
#             site_site = case_when(
#                   project == 'SBC' ~ site,
#                   project == 'FCE' ~ paste(site, subsite_level1, sep = ''),
#                   project == 'VCR' ~ paste(site, subsite_level1, sep = ''),
#                   project == 'MCR' ~ paste(subsite_level1, site, sep = ''),
#                   project == 'PCCC' ~ subsite_level2,
#                   project == 'PCCS' ~ subsite_level2,
#             ) 
#       ) |> 
#       group_by(project, site_site, year) |> 
#       summarize(
#             comm_n    = mean(comm_n, na.rm = TRUE),
#             comm_bm   = mean(comm_bm, na.rm = TRUE),
#             comm_dens = mean(comm_dens, na.rm = TRUE),
#             .groups   = 'drop'
#       ) |> 
#       rename(site = site_site)
# 
# glimpse(dt3_cnd)
# head(dt3_cnd)
# nacheck(dt3_cnd)
# 
# ## Calculating Species Richness and Diversity ----
# head(dt2)
# glimpse(dt2)
# nacheck(dt2)
# 
# dt3_sp <- dt2 |> 
#       group_by(project, habitat, year, month, 
#                site, subsite_level1, subsite_level2, subsite_level3) |> 
#       summarize(
#             species_richness = vegan::specnumber(species_dens),                 
#             inv_simpson      = vegan::diversity(species_dens, index = "invsimpson"),
#             total_dens       = sum(species_dens, na.rm = TRUE),
#             .groups = "drop"
#       ) |> 
#       # Inf occurs at transects with zero observed species (density = 0 rows);
#       # set to 0 to indicate absence of diversity
#       mutate(
#             inv_simpson = case_when(
#                   inv_simpson == Inf ~ 0,
#                   TRUE ~ inv_simpson
#             )) |> 
#       mutate(project = case_when(
#             project  == 'CoastalCA' & site == 'CENTRAL' ~ 'PCCC',
#             project  == 'CoastalCA' & site == 'SOUTH' ~ 'PCCS',
#             TRUE ~ project
#       )) |> 
#       ### take everything to the true 'site' level in which transects should be
#       ### averaged across
#       mutate(
#             site_site = case_when(
#                   project == 'SBC' ~ site,
#                   project == 'FCE' ~ paste(site, subsite_level1, sep = ''),
#                   project == 'VCR' ~ paste(site, subsite_level1, sep = ''),
#                   project == 'MCR' ~ paste(subsite_level1, site, sep = ''),
#                   project == 'PCCC' ~ subsite_level2,
#                   project == 'PCCS' ~ subsite_level2,
#             ) 
#       ) |> 
#       group_by(project, site_site, year) |> 
#       summarize(
#             s_rich    = mean(species_richness, na.rm = TRUE),
#             s_div     = mean(inv_simpson, na.rm = TRUE),
#             .groups   = 'drop'
#       ) |> 
#       rename(site = site_site)
# glimpse(dt3_sp)
# head(dt3_sp)
# nacheck(dt3_sp)
# 
# 
# ## Calculating Trophic Richness and Diversity  --------------------------------------
# dt2_mcr_troph <- dt1 |> 
#       filter(project == 'MCR') |> 
#       group_by(project, habitat, year, month, 
#                site, subsite_level1, subsite_level2, subsite_level3, 
#                diet_cat) |> 
#       summarize(
#             n       = sum(nind_ug_hr*density_num, na.rm = TRUE),
#             bm      = sum(dmperind_g_ind*density_num, na.rm = TRUE),
#             dens    = sum(density_num, na.rm = TRUE),
#             .groups = 'drop'
#       ) |> 
#       group_by(project, habitat, year, month, 
#                site, subsite_level1, subsite_level2, 
#                diet_cat) |> 
#       summarize(
#             troph_n    = sum(n),
#             troph_bm   = sum(bm),
#             troph_dens = sum(dens),
#             .groups = 'drop'
#       ) |> 
#       # subsite_level3 set to '15' to denote combined transects 1 and 5
#       mutate(subsite_level3 = '15') |> 
#       select(project, habitat, year, month, site, subsite_level1, subsite_level2, subsite_level3,
#              diet_cat, troph_n, troph_bm, troph_dens)
# 
# glimpse(dt2_mcr_troph)
# head(dt2_mcr_troph)
# nacheck(dt2_mcr_troph)
# 
# ### all other sites should be summed down subsite_level3 as these are
# ### considered the 'transect'
# dt2_other_troph <- dt1 |> 
#       filter(project != 'MCR') |> 
#       group_by(project, habitat, year, month, 
#                site, subsite_level1, subsite_level2, subsite_level3, 
#                diet_cat) |> 
#       ### sum across unique taxa at the transect level
#       ### coalesces multiple troph observation to single n and bm value for each transect
#       summarize(
#             troph_n    = sum(nind_ug_hr*density_num, na.rm = TRUE),
#             troph_bm   = sum(dmperind_g_ind*density_num, na.rm = TRUE),
#             troph_dens = sum(density_num, na.rm = TRUE),
#             .groups = 'drop'
#       ) |> 
#       select(project, habitat, year, month, site, subsite_level1, subsite_level2, subsite_level3,
#              diet_cat, troph_n, troph_bm, troph_dens)
# 
# glimpse(dt2_other_troph)
# head(dt2_other_troph)
# nacheck(dt2_other_troph)
# 
# ### join mcr and other data back together now that they are at same 'transect' scale
# dt2_troph <- rbind(dt2_other_troph, dt2_mcr_troph)
# glimpse(dt2_troph)
# head(dt2_troph)
# nacheck(dt2_troph)
# rm(dt2_mcr_troph, dt2_other_troph)
# 
# dt3_troph <- dt2_troph |> 
#       group_by(project, habitat, year, month, 
#                site, subsite_level1, subsite_level2, subsite_level3) |> 
#       summarize(
#             trophic_richness       = vegan::specnumber(troph_dens),                 
#             inv_simpson_troph      = vegan::diversity(troph_dens, index = "invsimpson"),
#             total_dens_troph       = sum(troph_dens, na.rm = TRUE),
#             .groups = "drop"
#       ) |> 
#       # Inf occurs at transects with zero observed trophic groups (density = 0 rows);
#       # set to 0 to indicate absence of diversity
#       mutate(
#             inv_simpson_troph = case_when(
#                   inv_simpson_troph == Inf ~ 0,
#                   TRUE ~ inv_simpson_troph
#             )) |> 
#       mutate(project = case_when(
#             project  == 'CoastalCA' & site == 'CENTRAL' ~ 'PCCC',
#             project  == 'CoastalCA' & site == 'SOUTH' ~ 'PCCS',
#             TRUE ~ project
#       )) |> 
#       ### take everything to the true 'site' level in which transects should be
#       ### averaged across
#       mutate(
#             site_site = case_when(
#                   project == 'SBC' ~ site,
#                   project == 'FCE' ~ paste(site, subsite_level1, sep = ''),
#                   project == 'VCR' ~ paste(site, subsite_level1, sep = ''),
#                   project == 'MCR' ~ paste(subsite_level1, site, sep = ''),
#                   project == 'PCCC' ~ subsite_level2,
#                   project == 'PCCS' ~ subsite_level2,
#             ) 
#       ) |> 
#       group_by(project, site_site, year) |> 
#       summarize(
#             t_rich    = mean(trophic_richness, na.rm = TRUE),
#             t_div     = mean(inv_simpson_troph, na.rm = TRUE),
#             .groups   = 'drop'
#       ) |> 
#       rename(site = site_site)
# glimpse(dt3_troph)
# head(dt3_troph)
# nacheck(dt3_troph)
# 
# ## Join trophic and species level data together ---------------------------------------
# dt3 <- dt3_cnd |> 
#       left_join(dt3_sp) |> 
#       left_join(dt3_troph)
# glimpse(dt3)
# head(dt3)
# nacheck(dt3)
# cnd_ts_data <- dt3
# # write_csv(cnd_ts_data, "local_data/annual-dt-for-summary.csv")
# 
# cnd_model_data <- cnd_ts_data |> 
#       group_by(project, site) |> 
#       summarize(
#             
#             ### nutrients
#             comm_n_mean      = mean(comm_n, na.rm = TRUE),
#             comm_n_sd        = sd(comm_n, na.rm = TRUE),
#             comm_n_cv        = sd(comm_n, na.rm = TRUE)/mean(comm_n, na.rm = TRUE),
#             comm_n_stability = 1/comm_n_cv,
#             
#             ### biomass
#             comm_bm_mean      = mean(comm_bm, na.rm = TRUE),
#             comm_bm_sd        = sd(comm_bm, na.rm = TRUE),
#             comm_bm_cv        = sd(comm_bm, na.rm = TRUE)/mean(comm_bm, na.rm = TRUE),
#             comm_bm_stability = 1/comm_bm_cv,
#             
#             ### diversity
#             s_rich_mean      = mean(s_rich, na.rm = TRUE),
#             s_div_mean       = mean(s_div, na.rm = TRUE),
#             t_rich_mean      = mean(t_rich, na.rm = TRUE),
#             t_div_mean       = mean(t_div, na.rm = TRUE),
#             .groups = 'drop'
#       )
# glimpse(cnd_model_data)
# head(cnd_model_data)
# nacheck(cnd_model_data)
# # write_csv(cnd_model_data, "local_data/community-level-nutrient-stability.csv")
# 
# 
# ## Calculate Species Dynamics -----------------------------------------------------------------
# glimpse(dt2)
# head(dt2)
# nacheck(dt2)
# 
# dt3_sp_dyn <- dt2 |> 
#       mutate(project = case_when(
#             project  == 'CoastalCA' & site == 'CENTRAL' ~ 'PCCC',
#             project  == 'CoastalCA' & site == 'SOUTH' ~ 'PCCS',
#             TRUE ~ project
#       )) |> 
#       mutate(
#             site_site = case_when(
#                   project == 'SBC' ~ site,
#                   project == 'FCE' ~ paste(site, subsite_level1, sep = ''),
#                   project == 'VCR' ~ paste(site, subsite_level1, sep = ''),
#                   project == 'MCR' ~ paste(subsite_level1, site, sep = ''),
#                   project == 'PCCC' ~ subsite_level2,
#                   project == 'PCCS' ~ subsite_level2,
#             ) 
#       ) |> 
#       group_by(project, site_site, year, scientific_name) |> 
#       summarize(
#             spp_n    = mean(species_n, na.rm = TRUE),
#             spp_bm   = mean(species_bm, na.rm = TRUE),
#             spp_dens = mean(species_dens, na.rm = TRUE),
#             .groups   = 'drop'
#       ) |> 
#       rename(site = site_site)
# glimpse(dt3_sp_dyn)
# head(dt3_sp_dyn)
# nacheck(dt3_sp_dyn)
# 
# spp_dyn_model_data <- dt3_sp_dyn |> 
#       group_by(project, site) |> 
#       nest() |> 
#       mutate(
#             spp_turnover = map_dbl(data, ~{
#                   beta_temp <- turnover(
#                         df            = .x,
#                         time.var      = "year",
#                         abundance.var = "spp_dens",
#                         species.var   = "scientific_name",
#                         metric        = "total"
#                   )
#                   mean(beta_temp$total, na.rm = TRUE)
#             }),
#             spp_synchrony = map_dbl(data, ~{
#                   synchrony(
#                         df            = .x,
#                         time.var      = "year",
#                         species.var   = "scientific_name",
#                         abundance.var = "spp_bm",
#                         metric        = "Loreau",
#                         replicate.var = NA
#                   )
#             })
#       ) |> 
#       ungroup() |> 
#       select(project, site, spp_turnover, spp_synchrony)
# 
# 
# ## Calculate Trophic Dynamics -----------------------------------------------------------------
# glimpse(dt2_troph)
# head(dt2_troph)
# nacheck(dt2_troph)
# 
# dt3_troph_dyn <- dt2_troph |> 
#       mutate(project = case_when(
#             project  == 'CoastalCA' & site == 'CENTRAL' ~ 'PCCC',
#             project  == 'CoastalCA' & site == 'SOUTH' ~ 'PCCS',
#             TRUE ~ project
#       )) |> 
#       mutate(
#             site_site = case_when(
#                   project == 'SBC' ~ site,
#                   project == 'FCE' ~ paste(site, subsite_level1, sep = ''),
#                   project == 'VCR' ~ paste(site, subsite_level1, sep = ''),
#                   project == 'MCR' ~ paste(subsite_level1, site, sep = ''),
#                   project == 'PCCC' ~ subsite_level2,
#                   project == 'PCCS' ~ subsite_level2,
#             ) 
#       ) |> 
#       group_by(project, site_site, year, diet_cat) |> 
#       summarize(
#             troph_n    = mean(troph_n, na.rm = TRUE),
#             troph_bm   = mean(troph_bm, na.rm = TRUE),
#             troph_dens = mean(troph_dens, na.rm = TRUE),
#             .groups    = 'drop'
#       ) |> 
#       rename(site = site_site)
# glimpse(dt3_troph_dyn)
# head(dt3_troph_dyn)
# nacheck(dt3_troph_dyn)
# 
# troph_dyn_model_data <- dt3_troph_dyn |> 
#       group_by(project, site) |> 
#       nest() |> 
#       mutate(
#             troph_turnover = map_dbl(data, ~{
#                   beta_temp <- turnover(
#                         df            = .x,
#                         time.var      = "year",
#                         abundance.var = "troph_dens",
#                         species.var   = "diet_cat",
#                         metric        = "total"
#                   )
#                   mean(beta_temp[, 1], na.rm = TRUE)
#             }),
#             troph_synchrony = map_dbl(data, ~{
#                   synchrony(
#                         df            = .x,
#                         time.var      = "year",
#                         species.var   = "diet_cat",
#                         abundance.var = "troph_bm",
#                         metric        = "Loreau",
#                         replicate.var = NA
#                   )
#             })
#       ) |> 
#       ungroup() |> 
#       select(project, site, troph_turnover, troph_synchrony)
# 
# ## Join trophic and species level dynamics together -----------------------------------------------
# glimpse(cnd_model_data)
# glimpse(spp_dyn_model_data)
# glimpse(troph_dyn_model_data)
# 
# ## Join all data to create final model dataset ----------------------------------------------------
# model_data_all <- cnd_model_data |> 
#       left_join(spp_dyn_model_data) |> 
#       left_join(troph_dyn_model_data)
# glimpse(model_data_all)
# head(model_data_all)
# nacheck(model_data_all)
# write_csv(model_data_all, "local_data/model-data-all-sandbox.csv")

# Run within system model ------------------------------------------------------------------------
model_data_all <- read_csv('local_data/model-data-all-sandbox.csv')
glimpse(model_data_all)
dat_scaled <- model_data_all |> 
      rename(program = project) |> 
      select(program, site, comm_n_stability, everything()) |> 
      
      # scale response metric
      mutate(comm_n_stability = as.numeric(scale(comm_n_stability))) |>
      
      # scale suite of predictor metrics 
      mutate(across(comm_n_mean:troph_synchrony, \(x) as.numeric(scale(x, center = TRUE))))
glimpse(dat_scaled)      
dat_ready_across <- dat_scaled      
glimpse(dat_ready_across)

num_vars <- dat_ready_across |> 
      dplyr::select(
            comm_n_stability, comm_bm_mean,
            s_rich_mean, s_div_mean, spp_turnover, spp_synchrony,
            t_rich_mean, t_div_mean, troph_turnover, troph_synchrony
      )

m_across <- cor(num_vars, use = "complete.obs")
corrplot(m_across, 
         method = "circle",      
         type = "upper",         
         addCoef.col = "black",  
         tl.col = "black",       
         diag = FALSE)

dat_scaled <- model_data_all |>
      rename(program = project) |>
      select(program, site, comm_n_stability, everything()) |>

      # comm_n_stability scaled globally (outcome); predictors scaled within-program
      mutate(comm_n_stability = as.numeric(scale(comm_n_stability))) |>

      # scaling predictors within program to control for between-program differences (i.e., across vs. within)
      # in baseline community metrics; coefficients reflect within-program
      # standardized effects
      group_by(program) |>
      mutate(across(comm_n_mean:troph_synchrony, \(x) as.numeric(scale(x, center = TRUE)))) |>
      ungroup()
glimpse(dat_scaled)
dat_ready_within <- dat_scaled
glimpse(dat_ready_within)

num_vars <- dat_ready_within |> 
      dplyr::select(
            comm_n_stability, comm_bm_mean,
            s_rich_mean, s_div_mean, spp_turnover, spp_synchrony,
            t_rich_mean, t_div_mean, troph_turnover, troph_synchrony
      )

m_within <- cor(num_vars, use = "complete.obs")
corrplot(m_within, 
         method = "circle",      
         type = "upper",         
         addCoef.col = "black",  
         tl.col = "black",       
         diag = FALSE)

# check models with biomass -----------------------------------------------

# # normal(0,1) weakly informative prior appropriate for standardized predictors - following Lemoine (2019, Ecology)
pr = prior(normal(0, 1), class = 'b')

### round one: single-term models to identify best individual predictor
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

m11 <- brm(
      comm_n_stability ~ comm_bm_mean + (comm_bm_mean | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

m12 <- brm(
      comm_n_stability ~ comm_n_mean + (comm_n_mean | program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

model_table_all <- performance::compare_performance(m1,m4,m5,m6,m7,m8,m9,m10,m11,m12)

model_selection1 <- model_table_all |>
      mutate(dWAIC = WAIC - min(WAIC))

write_csv(model_selection1, "output/tables/brms-fullmodel-selection-table-roundone-sandbox.csv")

keep <- c("nacheck", "model_data_all", "cnd_ts_data", "fit", "dat_ready", "pr", 'm4', 'model_selection1')
rm(list = setdiff(ls(), keep))

# Run across system model ------------------------------------------------------------------------
# Review correlation structure ------------------------------------------------------------------
num_vars <- dat_ready |> 
      dplyr::select(
            comm_n_stability, comm_bm_stability, comm_n_mean, comm_bm_mean,
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

m1 <- brm(
      comm_n_stability ~ spp_synchrony + comm_n_mean + (spp_synchrony + comm_n_mean | program),
      data = dat_ready_within,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

m2 <- brm(
      comm_n_stability ~ spp_synchrony + comm_bm_mean + (spp_synchrony + comm_bm_mean | program),
      data = dat_ready_within,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4,
      seed = 20
)

model_table_all <- performance::compare_performance(m0,m1,m2)

model_selection1 <- model_table_all |>
      mutate(dWAIC = WAIC - min(WAIC))

# trophic comp supplemental table -----------------------------------------

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
            troph_n  = mean(troph_n, na.rm = TRUE),
            .groups  = 'drop'
      ) |> 
      rename(site = site_site) |> 
      left_join(dt3_cnd) |> 
      dplyr::select(project, site, year, diet_cat, troph_n, comm_n) |> 
      
      complete(
            nesting(project, site, year, comm_n),
            diet_cat,
            fill = list(troph_n = 0)
      ) |> 
      
      mutate(prop = troph_n / comm_n) |> 
      
      group_by(project, diet_cat) |> 
      summarize(prop_m  = mean(prop, na.rm = TRUE),
                prop_sd = sd(prop, na.rm   = TRUE),
                prop_cv = prop_sd/prop_m,
                prop_stability = 1/prop_cv,
                .groups = 'drop') |> 
      dplyr::select(project, diet_cat, prop_m, prop_sd)
glimpse(dt3_troph)

dt3_troph_wide <- dt3_troph |>
      mutate(prop = paste0(
            formatC(prop_m, format = "f", digits = 2),
            " \u00b1 ",
            formatC(prop_sd, format = "f", digits = 2)
      )) |>
      select(project, diet_cat, prop) |>
      pivot_wider(names_from = diet_cat, values_from = prop) |> 
      dplyr::select(project, fish, fish_invert, invert, algae_invert, algae_detritus)
glimpse(dt3_troph_wide)
# write_csv(dt3_troph_wide, 'output/supplemental-table-trophic-proportions.csv')
