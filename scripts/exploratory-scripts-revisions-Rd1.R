# summary stats -----------------------------------------------------------

# supplemental analyses ---------------------------------------------------
glimpse(model_data_all)
sync_turn_cors <- model_data_all |> 
      group_by(project) |> 
      summarise(
            r = cor(spp_synchrony, spp_turnover, method = "pearson"),
            p = cor.test(spp_synchrony, spp_turnover)$p.value,
            n = n()
      ) |> 
      arrange(r)
print(sync_turn_cors)

project_order <- sync_turn_cors |> 
      arrange(r) |> 
      pull(project)

model_data_all |> 
      mutate(project = factor(project, levels = project_order)) |> 
      ggplot(aes(x = spp_synchrony, y = spp_turnover, color = project)) +
      geom_point(alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE) +
      facet_wrap(~ project, scales = "free") +
      theme_bw() +
      labs(
            x = "Species Synchrony",
            y = "Species Turnover",
            title = "Within-site covariance: Synchrony vs. Turnover"
      ) +
      theme(legend.position = "none")

##################################################################################################
### Stability ~ Magnitude ------------------------------------------------------------------------
##################################################################################################

summ <- model_data_all |>
      group_by(project) |> 
      mutate(stability = mean(comm_n_stability),
             magnitude  = mean(comm_n_mean))

summ_model <- lm(log1p(stability) ~ log1p(magnitude), data = summ)
summary(summ_model)$r.squared 
summary(summ_model)
r2_summ <- summary(summ_model)$r.squared
r2_summ

model_data_all |>
      ggplot(aes(x = log1p(comm_n_mean), y = log1p(comm_n_stability))) +
      geom_smooth(method = "lm", size = 1.5, color = "black", linetype = "solid", se = FALSE) +
      geom_point(aes(color = project), size = 1.5, alpha = 0.30) +
      geom_point(aes(x = log1p(magnitude), y = log1p(stability), color = project), size = 5, dat = summ) +
      labs(x = "log(CND Supply + 1)",
           y = "log(CND Stability + 1)",
           color = 'Program') +
      scale_y_continuous(breaks = seq(0.25, 1.75, by = 0.5)) +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 14),
            axis.text.y = element_text(face = "bold", color = "black", size = 14),
            axis.title.x = element_text(face = "bold", color = "black", size = 16),
            axis.title.y = element_text(face = "bold", color = "black", size = 16),
            legend.position = c(0.20, 0.65),
            legend.justification = c(1, 0),
            legend.text = element_text(face = "bold", color = "black"),
            legend.title = element_text(face = "bold", color = "black"))

##################################################################################################
### Magnitude ~ Richness -------------------------------------------------------------------------
##################################################################################################

summ <- dat |>
      group_by(program) |> 
      mutate(richness = mean(s_rich_mean),
             magnitude  = mean(comm_n_mean))

summ_model <- lm(log1p(magnitude) ~ log1p(richness), data = summ)
summary(summ_model)$r.squared 
summary(summ_model)
r2_summ <- summary(summ_model)$r.squared
r2_summ

dat |>
      ggplot(aes(x = log1p(s_rich_mean), y = log1p(comm_n_mean))) +
      geom_smooth(aes(color = program), method = "lm", se = FALSE) + 
      geom_point(aes(color = program), size = 1.5, alpha = 0.30) +
      # geom_point(aes(x = log1p(richness), y = log1p(magnitude), color = program), size = 5, dat = summ) +
      labs(x = "log(Species Richness + 1)",
           y = "log(CND Supply + 1)",
           color = 'Program') +
      # scale_y_continuous(breaks = seq(0.25, 1.75, by = 0.5)) +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 14),
            axis.text.y = element_text(face = "bold", color = "black", size = 14),
            axis.title.x = element_text(face = "bold", color = "black", size = 16),
            axis.title.y = element_text(face = "bold", color = "black", size = 16),
            legend.position = c(0.20, 0.65),
            legend.justification = c(1, 0),
            legend.text = element_text(face = "bold", color = "black"),
            legend.title = element_text(face = "bold", color = "black"))

dat |>
      ggplot(aes(x = log1p(s_rich_mean), 
                 y = log1p(comm_n_mean),
                 color = program)) +
      geom_smooth(method = "lm", se = FALSE) +
      geom_point(alpha = 0.30) +
      facet_wrap(~program, scales = "free") +
      scale_color_manual(values = program_palette) +
      theme_classic()

lmm_slopes <- lmer(
      log1p(comm_n_mean) ~ log1p(s_rich_mean) + 
            (log1p(s_rich_mean) | program), 
      data = dat
)
summary(lmm_slopes)

# extract R2 - marginal (fixed effects only) and conditional (full model)
r.squaredGLMM(lmm_slopes)
