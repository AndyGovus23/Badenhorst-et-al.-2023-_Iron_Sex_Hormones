#Paper Title: Does Chronic Oral Contraceptive Use Detrimentally Affect C-Reactive Protein 
#or Iron Status for Aerobically Trained Women?
#Date: 21/01/2023
#Purpose: To compare of iron and inflammatory parameters
#across the menstrual cycle for naturally menstruating women and OCP users.
#Outcome Variables: Iron, Ferritin, Transferrin, E2 (Estrogen), P4 (Progesterone), CRP, IL6
#Explanatory Variables: Time (4 levels: 0-3), Condition (2 levels: OCP Users, OCP non-users)

####LOAD LIBRARIES####

#Manipulate Excel Files
library(readxl)#Read in MS Excel .xlsx files

#Linear & Linear Mixed Modelling#
library(nlme) #Linear mixed modelling
library(lme4) #Linear mixed modelling
library(lmerTest) #Allow Kenward Rogers Approximation for DF
library(sjPlot) #Visualisation of mixed models
library(car) #Companion to Applied Regression
library(rmcorr) #Repeated measures correlation
library(rstatix)
library(mgcv) #Generalised additive modelling
library(itsadug) #Interpreting Time Series and Autocorrelated Data Using GAMMs

#Robust Modelling
library(robustlmm) #Robust linear mixed modelling
library(WRS2) #Robust linear modelling

#Visualisation & Descriptives#
library(jmv) #Jamovi package
library(dplyr)
library(tidyverse)
library(lattice)
library(ggResidpanel) #Check model residuals
library(performance) #Check model performance
library(cowplot) #For panelled plots

#Multiple Comparisons#
library(effects)
library(multcomp)
library(emmeans)
library(afex)

#Easy Stats Functions#
library(easystats)

####SECTION 0 - DATA WRANGLING####

####Open Data Source####
data.long.all <-read_excel("C:/Users/agovus/Desktop/Massey_Iron/Massey_Iron Parameters.xlsx", sheet = "Data")

#Extract required data
data.long.all_factors <- data.long.all[,1:3]
data.long.all_var <-  data.long.all[,13:20]
data.long.all_var <- dplyr::select(data.long.all_var, -TEMP)
data.long.all <- as_tibble(cbind(data.long.all_factors, data.long.all_var))

####Creation of Fixed Factors - Long Format####
data.long.all$ID <- factor(data.long.all$ID)
data.long.all$DAY <- factor(data.long.all$DAY)
data.long.all$GROUP <- factor(data.long.all$GROUP,
                          labels=c("Natural", "OCP"))

####Split data by OCP users vs. Non-OCP Users
data.long.ocp <- data.long.all %>% 
  group_by(GROUP) %>%
  filter(GROUP=="OCP")

#Drop unused level of GROUP#
data.long.ocp$GROUP <- droplevels(data.long.ocp$GROUP)
data.long.ocp$ID <- droplevels(data.long.ocp$ID)

####Data Adjustments####

#E2: 11 - Zero. Added 0.01 to create a non-zero value to allow regression estimation
#P4: 10, 11, 15 - Zero. Added 0.01 to create a non-zero value to allow regression estimation

####DATA TRANSFORMATIONS###

#Log transformation for later modelling - All Variables
data.long.all_2 <- log(data.long.all[,4:10])
colnames(data.long.all_2) <- paste("LN", colnames(data.long.all_2), sep="_")
data.long.all <- as_tibble(cbind(data.long.all,data.long.all_2))

####SECTION 1 - EXPLORATORY DATA ANALYSIS####

#Part 1 - Descriptive statistics#
all.desc <- data.long.all %>%
  group_by(GROUP, DAY) %>%
  summarize(across(P4:FE, list(mean=mean, 
                               sd=sd,
                               min=min,
                               max=max)))
                   
#Part 2 - Outlier Detection
out.p4 <- data.long.all %>%
  group_by(GROUP, DAY) %>%
  identify_outliers(P4)

out.e2 <- data.long.all %>%
  group_by(GROUP, DAY) %>%
  identify_outliers(E2)

out.fer <- data.long.all %>%
  group_by(GROUP, DAY) %>%
  identify_outliers(FER)

out.fe <- data.long.all %>%
  group_by(GROUP, DAY) %>%
  identify_outliers(FE)

out.trans <- data.long.all %>%
  group_by(GROUP, DAY) %>%
  identify_outliers(TRANS)

out.il6 <- data.long.all %>%
  group_by(GROUP, DAY) %>%
  identify_outliers(IL6)

out.crp <- data.long.all %>%
  group_by(GROUP, DAY) %>%
  identify_outliers(CRP)

####SECTION 2 - EXPLORATORY DATA VISUALISATION####

####PROGESTERONE (P4)####

####Plot - Individual Trajectories####
p4_plot_indiv <- ggplot(data.long.all,aes(DAY, P4, group=1,colour=GROUP)) +
  geom_point()+
  geom_line()+
  facet_wrap(~ID)+
  labs(x ="Day", y="P4 (pg/mL)")

####Plot - Individual Trajectories - Split by Group####
p4_plot_id_facet <- ggplot(data.long.all,aes(DAY, P4, group=ID, colour=ID)) +
  geom_point()+
  geom_line()+
  facet_wrap(~GROUP)+
  labs(x ="Day", y="P4 (pg/mL)")

###Plot - Time Series Boxplots - Faceted by Group###
p4_plot_id_boxplot_facet <- ggplot(data.long.all,aes(DAY, P4)) +
  geom_boxplot(outlier.alpha = 0)+
  geom_point(position=position_dodge(0.2), shape=18, size=3, aes(group=ID, colour=ID))+
  stat_summary(fun=median, colour="black", geom="line", linewidth=1, linetype=2, aes(group=1))+
  stat_summary(fun=median, colour="red", size=4, shape=18, geom="point",aes(group = 1))+
  facet_wrap(~GROUP)+
  labs(x ="Day", y="P4 (pg/mL)")+
  theme_bw()

####Plot - Group Average Trajectories - Faceted by Group####
p4_plot_group_av_facet <- ggplot(p4.desc,aes(DAY, mean_p4, group=1)) +
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=mean_p4-sd_p4,
                    ymax=mean_p4+sd_p4), width=.2,
                position=position_dodge(0.05))+
  facet_wrap(~GROUP)+
  labs(x ="Day", y="P4 (pg/mL)")

####Plot - Group Average Trajectory - No Facet####
p4_plot_group_av_nofacet <- ggplot(p4.desc, aes(x=DAY, y=mean_p4, group=GROUP, color=GROUP))+
  geom_line()+
  geom_point(position=position_dodge(0.05))+
  geom_errorbar(aes(ymin=mean_p4-sd_p4,
                    ymax=mean_p4+sd_p4), width=.2,
                position=position_dodge(0.05))+
  labs(x ="Day", y="P4 (pg/mL)")

####ESTROGEN (E2)####
e2_plot_indiv <- ggplot(data.long.all,aes(DAY, E2, group=1,colour=GROUP)) +
  geom_point()+
  geom_line()+
  facet_wrap(~ID)+
  labs(x ="Day", y="E2 (pg/mL)")

####Plot - Individual Trajectories####
e2_plot_indiv <- ggplot(data.long.all,aes(DAY, E2, group=1,colour=GROUP)) +
  geom_point()+
  geom_line()+
  facet_wrap(~ID)+
  labs(x ="Day", y="E2 (pg/mL)")

####Plot - Individual Trajectories - Split by Group####
e2_plot_id_facet <- ggplot(data.long.all,aes(DAY, E2, group=ID, colour=ID)) +
  geom_point()+
  geom_line()+
  facet_wrap(~GROUP)+
  labs(x ="Day", y="E2 (pg/mL)")

###Plot - Time Series Boxplots - Faceted by Group###
e2_plot_id_boxplot_facet <- ggplot(data.long.all,aes(DAY, E2)) +
  geom_boxplot(outlier.alpha = 0)+
  geom_point(position=position_dodge(0.2), shape=18, size=3, aes(group=ID, colour=ID))+
  stat_summary(fun=median, colour="black", geom="line", linewidth=1, linetype=2, aes(group=1))+
  stat_summary(fun=median, colour="red", size=4, shape=18, geom="point",aes(group = 1))+
  facet_wrap(~GROUP)+
  labs(x ="Day", y="E2 (pg/mL)")+
  theme_bw()

####Plot - Group Average Trajectories - Faceted by Group####
e2_plot_group_av_facet <- ggplot(e2.desc,aes(DAY, mean_e2, group=1)) +
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=mean_e2-sd_e2,
                    ymax=mean_e2+sd_e2), width=.2,
                position=position_dodge(0.05))+
  facet_wrap(~GROUP)+
  labs(x ="Day", y="E2 (pg/mL)")

####Plot - Group Average Trajectory - No Facet####
e2_plot_group_av_nofacet <- ggplot(e2.desc, aes(x=DAY, y=mean_e2, group=GROUP, color=GROUP))+
  geom_line()+
  geom_point(position=position_dodge(0.05))+
  geom_errorbar(aes(ymin=mean_e2-sd_e2,
                    ymax=mean_e2+sd_e2), width=.2,
                position=position_dodge(0.05))+
  labs(x ="Day", y="E2 (pg/mL)")

####FERRITIN####

####Plot - Individual Trajectories####
fer_plot_indiv <- ggplot(data.long.all,aes(DAY, FER, group=1,colour=GROUP)) +
  geom_point()+
  geom_line()+
  facet_wrap(~ID)+
  labs(x ="Day", y="Serum Ferritin (ng/L)")

####Plot - Individual Trajectories - Split by Group####
fer_plot_id_facet <- ggplot(data.long.all,aes(DAY, FER, group=ID, colour=ID)) +
  geom_point()+
  geom_line()+
  facet_wrap(~GROUP)+
  labs(x ="Day", y="Serum Ferritin (ng/L)")

###Plot - Time Series Boxplots - Faceted by Group###
fer_plot_id_boxplot_facet <- ggplot(data.long.all,aes(DAY, FER)) +
  geom_boxplot(outlier.alpha = 0)+
  geom_point(position=position_dodge(0.2), shape=18, size=3, aes(group=ID, colour=ID))+
  stat_summary(fun=median, colour="black", geom="line", linewidth=1, linetype=2, aes(group=1))+
  stat_summary(fun=median, colour="red", size=4, shape=18, geom="point",aes(group = 1))+
  facet_wrap(~GROUP)+
  labs(x ="Day", y="Serum Ferritin (ng/L)")+
  theme_bw()

####Plot - Group Average Trajectories - Faceted by Group####
fer_plot_group_av_facet <- ggplot(fer.desc,aes(DAY, mean_fer, group=1)) +
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=mean_fer-sd_fer,
                    ymax=mean_fer+sd_fer), width=.2,
                position=position_dodge(0.05))+
  facet_wrap(~GROUP)
labs(x ="Day", y="Serum Ferritin (ng/L)")

####IRON####

####Plot - Individual Trajectories####
fe_plot_indiv <- ggplot(data.long.all,aes(DAY, FE, group=1,colour=GROUP)) +
  geom_point()+
  geom_line()+
  facet_wrap(~ID)+
  labs(x ="Day", y="FE (mg/L)")

####Plot - Individual Trajectories - Split by Group####
fe_plot_id_facet <- ggplot(data.long.all,aes(DAY, FE, group=ID, colour=ID)) +
  geom_point()+
  geom_line()+
  facet_wrap(~GROUP)+
  labs(x ="Day", y="FE (mg/L)")

###Plot - Time Series Boxplots - Faceted by Group###
fe_plot_id_boxplot_facet <- ggplot(data.long.all,aes(DAY, FE)) +
  geom_boxplot(outlier.alpha = 0)+
  geom_point(position=position_dodge(0.2), shape=18, size=3, aes(group=ID, colour=ID))+
  stat_summary(fun=median, colour="black", geom="line", linewidth=1, linetype=2, aes(group=1))+
  stat_summary(fun=median, colour="red", size=4, shape=18, geom="point",aes(group = 1))+
  facet_wrap(~GROUP)+
  labs(x ="Day", y="FE (mg/L)")+
  theme_bw()

####Plot - Group Average Trajectories - Faceted by Group####
fe_plot_group_av_facet <- ggplot(fe.desc,aes(DAY, mean_fe, group=1)) +
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=mean_fe-sd_fe,
                    ymax=mean_fe+sd_fe), width=.2,
                position=position_dodge(0.05))+
  facet_wrap(~GROUP)+
  labs(x ="Day", y="FE (mg/L)")

####TRANSFERRIN####

####Plot - Individual Trajectories####
  fer_plot_indiv <- ggplot(data.long.all,aes(DAY, FER, group=1,colour=GROUP)) +
    geom_point()+
    geom_line()+
    facet_wrap(~ID)+
    labs(x ="Day", y="Serum Ferritin (ng/L)")
  
####Plot - Individual Trajectories - Split by Group####
  fer_plot_id_facet <- ggplot(data.long.all,aes(DAY, FER, group=ID, colour=ID)) +
    geom_point()+
    geom_line()+
    facet_wrap(~GROUP)+
    labs(x ="Day", y="Serum Ferritin (ng/L)")
  
###Plot - Time Series Boxplots - Faceted by Group###
  fer_plot_id_boxplot_facet <- ggplot(data.long.all,aes(DAY, FER)) +
    geom_boxplot(outlier.alpha = 0)+
    geom_point(position=position_dodge(0.2), shape=18, size=3, aes(group=ID, colour=ID))+
    stat_summary(fun=median, colour="black", geom="line", linewidth=1, linetype=2, aes(group=1))+
    stat_summary(fun=median, colour="red", size=4, shape=18, geom="point",aes(group = 1))+
    facet_wrap(~GROUP)+
    labs(x ="Day", y="Serum Ferritin (ng/L)")+
    theme_bw()
  
####Plot - Group Average Trajectories - Faceted by Group####
  fer_plot_group_av_facet <- ggplot(fer.desc,aes(DAY, mean_fer, group=1)) +
    geom_point()+
    geom_line()+
    geom_errorbar(aes(ymin=mean_fer-sd_fer,
                      ymax=mean_fer+sd_fer), width=.2,
                  position=position_dodge(0.05))+
    facet_wrap(~GROUP)
  labs(x ="Day", y="Serum Transferrin (mg/L)")
  
####IL6####
  
####Plot - Individual Trajectories####
  il6_plot_indiv <- ggplot(data.long.all,aes(DAY, IL6, group=1,colour=GROUP)) +
    geom_point()+
    geom_line()+
    facet_wrap(~ID)+
    labs(x ="Day", y="IL6 (pg/mL)")
  
####Plot - Individual Trajectories - Split by Group####
  il6_plot_id_facet <- ggplot(data.long.all,aes(DAY, IL6, group=ID, colour=ID)) +
    geom_point()+
    geom_line()+
    facet_wrap(~GROUP)+
    labs(x ="Day", y="IL6 (pg/mL)")
  
###Plot - Time Series Boxplots - Faceted by Group###
  il6_plot_id_boxplot_facet <- ggplot(data.long.all,aes(DAY, IL6)) +
    geom_boxplot(outlier.alpha = 0)+
    geom_point(position=position_dodge(0.2), shape=18, size=3, aes(group=ID, colour=ID))+
    stat_summary(fun=median, colour="black", geom="line", linewidth=1, linetype=2, aes(group=1))+
    stat_summary(fun=median, colour="red", size=4, shape=18, geom="point",aes(group = 1))+
    facet_wrap(~GROUP)+
    labs(x ="Day", y="IL6 (pg/mL)")+
    theme_bw()
  
####Plot - Group Average Trajectories - Faceted by Group####
  il6_plot_group_av_facet <- ggplot(il6.desc,aes(DAY, mean_il6, group=1)) +
    geom_point()+
    geom_line()+
    geom_errorbar(aes(ymin=mean_il6-sd_il6,
                      ymax=mean_il6+sd_il6), width=.2,
                  position=position_dodge(0.05))+
    facet_wrap(~GROUP)+
    labs(x ="Day", y="IL6 (pg/mL)")

####CRP####
  
####Plot - Individual Trajectories####
  crp_plot_indiv <- ggplot(data.long.all,aes(DAY, CRP, group=1,colour=GROUP)) +
    geom_point()+
    geom_line()+
    facet_wrap(~ID)+
    labs(x ="Day", y="CRP (pg/mL)")
  
####Plot - Individual Trajectories - Split by Group####
  crp_plot_id_facet <- ggplot(data.long.all,aes(DAY, CRP, group=ID, colour=ID)) +
    geom_point()+
    geom_line()+
    facet_wrap(~GROUP)+
    labs(x ="Day", y="CRP (pg/mL)")
  
###Plot - Time Series Boxplots - Faceted by Group###
  crp_plot_id_boxplot_facet <- ggplot(data.long.all,aes(DAY, CRP)) +
    geom_boxplot(outlier.alpha = 0)+
    geom_point(position=position_dodge(0.2), shape=18, size=3, aes(group=ID, colour=ID))+
    stat_summary(fun=median, colour="black", geom="line", linewidth=1, linetype=2, aes(group=1))+
    stat_summary(fun=median, colour="red", size=4, shape=18, geom="point",aes(group = 1))+
    facet_wrap(~GROUP)+
    labs(x ="Day", y="CRP (pg/mL)")+
    theme_bw()
  
####Plot - Group Average Trajectories - Faceted by Group####
  crp_plot_group_av_facet <- ggplot(crp.desc,aes(DAY, mean_crp, group=1)) +
    geom_point()+
    geom_line()+
    geom_errorbar(aes(ymin=mean_crp-sd_crp,
                      ymax=mean_crp+sd_crp), width=.2,
                  position=position_dodge(0.05))+
    facet_wrap(~GROUP)+
    labs(x ="Day", y="CRP (pg/mL)")

####SECTION 3 - LINEAR MIXED MODELLING####

###Model strategy###
#1. Fit all models with a random intercept for ID
#2. Examine the need for 1) heterogeneous variance for Day x Time
#3. Examine the need to autocorrelated errors
#4. Refit model with additional covariates
#5. Simplify model to remove non-influential covariates

#Specify Optimisation Algorithm & Increase Iteration Limit#
ctrl <- lmeControl(opt="optim", msMaxIter=100)

#Change to effects contrast coding to sum contrasts due to categorical interaction terms#
#NB - Contrast should sum to zero, rather than 1
contrasts(data.long.all$GROUP) <- contr.sum(2)
contrasts(data.long.all$DAY) <- contr.sum(4)

####PROGESTERONE (P4)####

####Part 1 - General Linear Mixed Models####

####Model 0####
p4_int_0 <- lme(P4~DAY*GROUP,
                random=~1|ID,
                method = "ML",
                control=ctrl,
                na.action=na.omit,
                data=data.long.all)

####Model 1####
p4_int_1 <- lme(P4~DAY*GROUP,
                random=~1|ID,
                corr=corExp(form=~DAY),
                weights = varIdent(form=~1|DAY*GROUP),
                method = "ML",
                control=ctrl,
                na.action=na.omit,
                data=data.long.all)

#Does not converge - Singularity at backsolve.

####Model 2####
p4_int_2 <- lme(P4~DAY*GROUP,
                random=~1|ID,
                corr=corExp(form=~DAY),
                weights = varIdent(form=~1|GROUP),
                method = "ML",
                control=ctrl,
                na.action=na.omit,
                data=data.long.all)

####Model 3####
p4_int_3 <- lme(P4~DAY*GROUP,
                random=~1|ID,
                corr=corExp(form=~DAY),
                method = "ML",
                control=ctrl,
                na.action=na.omit,
                data=data.long.all)

####Model 4####
p4_int_4 <- lme(P4~DAY*GROUP,
                random=~1|ID,
                weights = varIdent(form=~1|DAY*GROUP),
                method = "ML",
                control=ctrl,
                na.action=na.omit,
                data=data.long.all)

#Check - ANOVA Table#
anova(p4_int_0,
      p4_int_2,
      p4_int_3,
      p4_int_4)
#Autocorrelation - Does not improve model fit
#Heterogeneous variances - Improves model fit
#Model 4: Parsimonious

####Model 5 - Refit Parsimonious Model with REML####
p4_int_4_reml <- lme(P4~DAY*GROUP,
                     random=~1|ID,
                     weights = varIdent(form=~1|DAY*GROUP),
                     method = "REML",
                     control=ctrl,
                     na.action=na.omit,
                     data=data.long.all)

#Check - Model ANOVA Table
Anova(p4_int_4_reml, type = 3)
#Sig: Day, Group, Day x Group

#Plot - Model Residuals
resid_panel(p4_int_4_reml,qqbands=TRUE)

#Plot - Model Effects
plot(allEffects(p4_int_4_reml))
#Model Fit: Poor - Long tailed residuals and high and low end
#Consider log transformation or Gamma distribution instead

####Model 6 - Log Transformation####
ln_p4_int_6_reml <- lme(log(P4)~DAY*GROUP,
                        random=~1|ID,
                        weights = varIdent(form=~1|DAY*GROUP),
                        method = "REML",
                        control=ctrl,
                        na.action=na.omit,
                        data=data.long.all)

#Check - Model ANOVA Table
anova(ln_p4_int_6_reml, type="marginal")
#Sig: Day, Group, Day x Group

#Plot - Model Residuals
resid_panel(ln_p4_int_6_reml,qqbands=TRUE)
#Fit: Acceptable model fit

#Plot _ Model Diagnostics
check_model(ln_p4_int_6_reml)

#Plot - Model Effects
plot(allEffects(ln_p4_int_6_reml))

#Check - Model Effects
summary(allEffects(ln_p4_int_6_reml))

#Check - Post Tests - Pairwise Contrasts
post_hoc_p4_int <- confint(pairs(emmeans(ln_p4_int_6_reml, "GROUP", by = "DAY")))
print(post_hoc_p4_int)

post_hoc_p4_int_r <- confint(pairs(emmeans(ln_p4_int_6_reml, "GROUP", by = "DAY", type="response")))

#Check - Estimated Marginal Means#
emm_post_hoc_p4_int <-  emmeans(ln_p4_int_6_reml, "GROUP", by = "DAY")
tidy(emm_post_hoc_p4_int)

####Model 7 - Log Transformation####
ln_p4_crp_int_7_reml <- lme(LN_P4~DAY*GROUP+CRP,
                            random=~1|ID,
                            weights = varIdent(form=~1|DAY*GROUP),
                            method = "REML",
                            control=ctrl,
                            na.action=na.omit,
                            data=data.long.all)

#Check - Model ANOVA Table
anova(ln_p4_crp_int_7_reml, type = "marginal")
#Sig: Day, Group, Day x Group
#NS: CRP

#Plot - Model Residuals
resid_panel(ln_p4_crp_int_7_reml,qqbands=TRUE)
#Fit: Acceptable model fit

#Check - Model Tables
tab_model(ln_p4_int_6_reml,
          ln_p4_crp_int_7_reml,
          show.aicc=TRUE)

BIC(ln_p4_int_6_reml,
    ln_p4_crp_int_7_reml)

####Part 2 - Generalised Linear Mixed Model####

####Model 1 - Gamma Model - Unadjusted####
gamma_p4_int_1_reml <- glmer(P4~DAY*GROUP+(1|ID),
                             family=Gamma(link="log"),
                             na.action=na.omit,
                             data=data.long.all)

#Check - ANOVA Table
Anova(gamma_p4_int_1_reml, type=3)
#Sig: Group, Day, Group x Day

#Plot - Model Residuals
resid_panel(gamma_p4_int_1_reml,qqbands=TRUE)
#Fit: Residuals approximately normal

#Plot _ Model Diagnostics
check_model(gamma_p4_int_1_reml)

#Check - Model Performance
model_performance(gamma_p4_int_1_reml)

#Check - Model Effects
summary(allEffects(gamma_p4_int_1_reml))

#Plot - Model Effects
plot(allEffects(gamma_p4_int_1_reml))

#Check - Post Tests - Pairwise Contrasts
post_hoc_p4_int_gamma <- confint(pairs(emmeans(gamma_p4_int_1_reml, "GROUP", by = "DAY")))
print(post_hoc_p4_int_gamma)

#Check - Estimated Marginal Means#
emm_post_hoc_p4_int_gamma <-  emmeans(gamma_p4_int_1_reml, "GROUP", by = "DAY")
print(emm_post_hoc_p4_int_gamma)

####Model 2 - Gamma Model - Adjusted for CRP####
gamma_p4_int_2_reml <- glmer(P4~DAY*GROUP+CRP+(1|ID),
                             family=Gamma(link="log"),
                             na.action=na.omit,
                             data=data.long.all)

#Check - ANOVA Table
Anova(gamma_p4_int_2_reml, type=3)
#Sig: Day, Group, Day x Group
#NS: CRP

#Plot - Model Residuals
resid_panel(gamma_p4_int_2_reml,qqbands=TRUE)
#Fit: Residuals approximately normal

#Check - Model Tables
tab_model(gamma_p4_int_1_reml,
          gamma_p4_int_2_reml)

BIC(gamma_p4_int_1_reml,
    gamma_p4_int_2_reml)

####ESTROGEN (E2)####

####Model 0 - Null Model####
e2_int_0 <- lme(E2~DAY*GROUP,
                random=~1|ID,
                method = "ML",
                control=ctrl,
                na.action=na.omit,
                data=data.long.all)

####Model 1 - Exponential Autocorrelation + Heterogeneous variances for Day & Group####
e2_int_1 <- lme(E2~DAY*GROUP,
                random=~1|ID,
                corr=corExp(form=~DAY),
                weights = varIdent(form=~1|DAY*GROUP),
                method = "ML",
                control=ctrl,
                na.action=na.omit,
                data=data.long.all)

####Model 2 - Exponential Autocorrelation + Heterogeneous variances for Group####
e2_int_2 <- lme(E2~DAY*GROUP,
                random=~1|ID,
                corr=corExp(form=~DAY),
                weights = varIdent(form=~1|GROUP),
                method = "ML",
                control=ctrl,
                na.action=na.omit,
                data=data.long.all)

####Model 3 - Exponential Autocorrelation####
e2_int_3 <- lme(E2~DAY*GROUP,
                random=~1|ID,
                corr=corExp(form=~DAY),
                method = "ML",
                control=ctrl,
                na.action=na.omit,
                data=data.long.all)

####Model 4 - Heterogeneous variances for Day & Group####
e2_int_4 <- lme(E2~DAY*GROUP,
                random=~1|ID,
                weights = varIdent(form=~1|DAY*GROUP),
                method = "ML",
                control=ctrl,
                na.action=na.omit,
                data=data.long.all)

#Check - ANOVA Table#
anova(e2_int_0,
      e2_int_1,
      e2_int_2,
      e2_int_3,
      e2_int_4)

#Model 1: Parsimonious

####Model 5 - Parsimonious Model - REML####
ln_e2_int_1_reml <- lme(LN_E2~DAY*GROUP,
                     random=~1|ID,
                     corr=corExp(form=~DAY),
                     weights = varIdent(form=~1|DAY*GROUP),
                     method = "REML",
                     control=ctrl,
                     na.action=na.omit,
                     data=data.long.all)

#Check - Model ANOVA Table
anova(ln_e2_int_1_reml, type="marginal")
#Sig: Day, Group
#NS: Day x Group

#Plot - Model Residuals
resid_panel(ln_e2_int_1_reml,qqbands=TRUE) 
#Long-tailed on the left due to low E2 value (near zero)

###Model 6 - Day x Time adjusted for CRP - REML####
ln_e2_int_2_crp_reml <- lme(LN_E2~DAY*GROUP+CRP,
                         random=~1|ID,
                         corr=corExp(form=~DAY),
                         weights = varIdent(form=~1|DAY*GROUP),
                         method = "REML",
                         control=ctrl,
                         na.action=na.omit,
                         data=data.long.all)

#Check - Model ANOVA Table
anova(ln_e2_int_2_crp_reml, type="marginal")
#Sig: Group
#NS: CRP, Day x Group

#Plot - Model Residuals
resid_panel(ln_e2_int_2_crp_reml,qqbands=TRUE) 
#Long-tailed on the left due to low E2 value (near zero)
#Switch to a Gamma distribution instead

#Plot - Model Effects
plot(allEffects(ln_e2_int_2_crp_reml))

#Check - Model  Effects
summary(allEffects(ln_e2_int_2_crp_reml))

#Check - Model Tables
tab_model(e2_int_1_reml,
          e2_int_2_crp_reml,
          show.aicc = TRUE)

#Check - BIC
BIC(e2_int_1_reml,
    e2_int_2_crp_reml)

####Part 2 - Generalised Linear Mixed Model####

####Model 1 - Gamma Model - Unadjusted####
gamma_e2_int_1_reml <- glmer(E2~DAY*GROUP+(1|ID),
                             family=Gamma(link="log"),
                             na.action=na.omit,
                             data=data.long.all)

#Check - Model ANOVA
Anova(gamma_e2_int_1_reml, type=3)
#Sig: Group
#NS: Day, Day x Group

#Plot - Model Residuals
resid_panel(gamma_e2_int_1_reml,qqbands=TRUE)
#Fit: Residuals are approximately normal, but long-tailed at low end

####Model 2 - Gamma Model - Adjusted for CRP####
gamma_e2_int_2_reml <- glmer(E2~DAY*GROUP+CRP+(1|ID),
                             family=Gamma(link="log"),
                             na.action=na.omit,
                             data=data.long.all)

#NB - Convergence issues

#Check - Model ANOVA
Anova(gamma_e2_int_2_reml, type = 3)
#Sig: Group
#NS: Day, Day x Group, CRP

#Plot - Model Residuals
resid_panel(gamma_e2_int_2_reml,qqbands=TRUE)
#Fit: Residuals are approximately normal, but long-tailed at low end

#Check - Model Tables
tab_model(gamma_e2_int_1_reml,
          gamma_e2_int_2_reml)

BIC(gamma_e2_int_1_reml,
          gamma_e2_int_2_reml)

####SERUM FERRITIN####

####Part 1 - Linear Mixed Model####

####Model 1 - Random Intercept - Heterogeneous variances: Day*Group, Exponential Autocorrelation####
ln_fer_int_1 <- lme(log(FER)~DAY*GROUP,
                                random=~1|ID,
                                corr=corExp(form=~DAY),
                                weights = varIdent(form=~1|DAY*GROUP),
                                method = "REML",
                                control=ctrl,
                                na.action=na.omit,
                                data=data.long.all)

#Check - Model ANOVA
anova(ln_fer_int_1, type="marginal")
#NS: Everything

#Check - Model Residuals
resid_panel(ln_fer_int_1, qqbands=T)
#Fit: Residuals not approximately normal

####Model 2 - Add Covariates####
ln_fer_int_2 <- lme(LN_FER~DAY*GROUP+E2+P4+IL6+CRP,
                 random=~1|ID,
                 corr=corExp(form=~DAY),
                 weights = varIdent(form=~1|DAY*GROUP),
                 method = "REML",
                 control=ctrl,
                 na.action=na.omit,
                 data=data.long.all)

#Check - Model ANOVA Table
anova(ln_fer_int_2, type="marginal") 
#Sig: CRP
#NS: Everything else

#Check - Model Fit
resid_panel(ln_fer_int_2, qqbands=T)
#Fit: Residuals are not apprixmately normal, despite log transformation
#Fit: Consider a Gamma distribution instead

####Part 2 - Generalised Linear Mixed Model - Gamma Distribution####

####Model 1 - Random Intercepts Model - Log link####
fer_int_gamma <- glmer(FER~DAY*GROUP+
                            (1|ID),
                          family=Gamma(link="log"),
                          na.action=na.omit,
                          data=data.long.all)

#Check - Model ANOVA
Anova(fer_int_gamma, type=3)
#NS: Everything

#Check - Model Diagnostics
resid_panel(fer_int_gamma, qqbands=TRUE)

#Check - Model Table
tab_model(fer_int_gamma)

####Model 2 - Random Intercepts Model - Covariate Adjustment####
fer_int_gamma_cov <- glmer(FER~DAY*GROUP+E2+P4+IL6+CRP+
                          (1|ID),
                          family=Gamma(link="log"),
                          na.action=na.omit,
                          data=data.long.all)

#NB - Identifiability issues - Large eigenvalues

#Check - Model ANOVA
Anova(fer_int_gamma_cov)
#Sig: CRP
#NS: All other variables
#NB: Simplify model or scale excoriates

#Check - Model Diagnostics
resid_panel(fer_int_gamma_cov, qqbands=TRUE)
#NB: Acceptable residual distribution

####Model 3 - Random Intercepts Model - Reduced Covariate Model####
fer_int_gamma_crp <- glmer(FER~DAY*GROUP+CRP+
                          (1|ID),
                          family=Gamma(link="log"),
                          na.action=na.omit,
                          data=data.long.all)

#Check - Model ANOVA
Anova(fer_int_gamma_crp, type=3)
#NS: All variables, CRP (P = 0.054)

#Check - Model Diagnostics
resid_panel(fer_int_gamma_crp, qqbands=TRUE)
#Fit: Acceptable

#Check - Model Comparison
tab_model(fer_int_gamma,
          fer_int_gamma_cov,
          fer_int_gamma_crp)

BIC(fer_int_gamma,
          fer_int_gamma_cov,
          fer_int_gamma_crp)

#Check - Model Effects
summary(allEffects(fer_int_gamma_crp))

#Plot - Model Effects
plot(allEffects(fer_int_gamma_crp))

####SERUM IRON####

####Part 1 - Linear Mixed Model####

####Model 1 - Random Intercepts, Heterogeneous variance by Group####
fe_int_1 <- lme(FE~DAY*GROUP,
                 random=~1|ID,
                 method = "ML",
                 na.action=na.omit,
                 control=ctrl,
                 data=data.long.all)

####Model 2 - Random Intercept - Heterogeneous Variances by day split by group####
fe_int_2 <- lme(FE~DAY*GROUP,
                               random=~1|ID,
                               weights = varIdent(form=~1|DAY*GROUP),
                               method = "ML",
                               na.action=na.omit,
                               control=ctrl,
                               data=data.long.all)

####Model 3 - Random Intercept - Exponential Autocorrelation####
fe_int_3 <- lme(FE~DAY*GROUP,
                    random=~1|ID,
                    corr=corExp(form=~DAY),
                    method = "ML",
                    na.action=na.omit,
                    control=ctrl,
                    data=data.long.all)

####Model 4 - Random Intercept - Heterogeneous variances: Day*Group, Exponential Autocorrelation####
fe_int_4 <- lme(FE~DAY*GROUP,
                                  random=~1|ID,
                                  corr=corExp(form=~DAY),
                                  weights = varIdent(form=~1|DAY*GROUP),
                                  method = "ML",
                                  control=ctrl,
                                  na.action=na.omit,
                                  data=data.long.all)

#Check - Model Comparison
anova(fe_int_1,
      fe_int_2,
      fe_int_3,
      fe_int_4)
#Heterogeneous variances and autocorrelation - Not required

####Model 5 - GROUP x DAY - Log-Transformed - REML####
ln_int_5_reml <- lmer(LN_FE~DAY*GROUP+
                            (1|ID),
                          na.action=na.omit,
                          data=data.long.all)

#Check - Model Effects
anova(ln_int_5_reml, type="marginal", test="F")
#NS: Everything

#Plot - Model Residuals
resid_panel(ln_int_5_reml, qqbands=TRUE)
#Fit: Acceptable

####Model 6 - ANCOVA - All Covariates - Log Transformed - REML####
ln_int_6_reml <- lmer(log(FE)~DAY*GROUP+CRP+IL6+P4+E2+
                                (1|ID),
                              na.action=na.omit,
                              data=data.long.all)

#Check - Model Effects
anova(ln_int_6_reml, type="marginal", test="F")
#NS: Everything, CRP (P = 0.08)

#Plot - Model Residuals
resid_panel(ln_fe_int_cov_ri_reml, qqbands=TRUE)
#Fit: Acceptable.

####Model 7 - ANCOVA - Group*Day Only####
ln_int_7_reml <- lmer(LN_FE~DAY*GROUP+CRP+
                                (1|ID),
                              na.action=na.omit,
                              data=data.long.all)

#Check - Model Effects
anova(ln_int_7_reml, type="marginal")
#NS: All factors

#Plot - Model Residuals
resid_panel(ln_int_7_reml, qqbands=TRUE)
#Fit: Acceptable

#Check- Model Fit
check_model(ln_int_7_reml)
model_performance(ln_int_7_reml)

####Part 2 - Generalised Linear Model####

####Model 1 - GROUP x DAY - Gamma log link - REML####
gamma_fe_int_ri_reml <- glmer(FE~DAY*GROUP+
                                (1|ID),
                              family=Gamma(link="log"),
                              na.action=na.omit,
                              data=data.long.all)

#Check - Model Effects
Anova(gamma_fe_int_ri_reml, type=3)
#NS: Everything

#Plot - Model Residuals
resid_panel(gamma_fe_int_ri_reml, qqbands=TRUE)
#Fit: Acceptable

####Model 2 - ANCOVA - All Covariates - Gamma log link - REML####
gamma_fe_int_cov_ri_reml <- glmer(FE~DAY*GROUP+CRP+IL6+P4+E2+
                                    (1|ID),
                                  family=Gamma(link="log"),
                                  na.action=na.omit,
                                  data=data.long.all)

#NB: Identifiability issues - large eigenvalues

#Check - Model Effects
Anova(gamma_fe_int_cov_ri_reml, type=3)
#NS: Everything, CRP (P = 0.04), Day x Group (p = 0.085)

#Plot - Model Residuals
resid_panel(gamma_fe_int_cov_ri_reml, qqbands=TRUE)
#Fit: Acceptable.

####Model 3 - ANCOVA - Reduced Model - Gamma log link####
gamma_fe_int_crp_ri_reml <- glmer(FE~DAY*GROUP+CRP+
                                    (1|ID),
                                  family=Gamma(link="log"),
                                  na.action=na.omit,
                                  data=data.long.all)

#Check - Model Effects
Anova(gamma_fe_int_crp_ri_reml, type=3)
#NS: Everything, CRP (P = 0.08)

#Plot - Model Residuals
resid_panel(gamma_fe_int_crp_ri_reml, qqbands=TRUE)
#Fit: Acceptable.

#Check- Model Table
tab_model(gamma_fe_int_ri_reml,
          gamma_fe_int_cov_ri_reml,
          gamma_fe_int_crp_ri_reml)

BIC(gamma_fe_int_ri_reml,
    gamma_fe_int_cov_ri_reml,
    gamma_fe_int_crp_ri_reml)

AIC(gamma_fe_int_ri_reml,
    gamma_fe_int_cov_ri_reml,
    gamma_fe_int_crp_ri_reml)

####TRANSFERRIN####

####Model 0 - Basic Interaction Model####
trans_int_0_ri <- lme(TRANS~DAY*GROUP,
                           random=~1|ID,
                           method = "ML",
                           control=ctrl,
                           na.action=na.omit,
                           data=data.long.all)

####Model 1 - Heterogeneous Variances per Day and Group over Time####
trans_int_1_het_exp <- lme(TRANS~DAY*GROUP,
                   random=~1|ID,
                   corr=corExp(form=~DAY),
                   weights = varIdent(form=~1|DAY*GROUP),
                   method = "ML",
                   control=ctrl,
                   na.action=na.omit,
                   data=data.long.all)

####Model 2 - Heterogeneous Variances - Day x Group####
trans_int_2_het <- lme(TRANS~DAY*GROUP,
                   random=~1|ID,
                   weights = varIdent(form=~1|DAY*GROUP),
                   method = "ML",
                   control=ctrl,
                   na.action=na.omit,
                   data=data.long.all)

####Model 3 - Heterogeneous Variance by Day####
trans_int_2_het_day <- lme(TRANS~DAY*GROUP,
                       random=~1|ID,
                       weights = varIdent(form=~1|DAY),
                       method = "ML",
                       control=ctrl,
                       na.action=na.omit,
                       data=data.long.all)

####Model 4 - Exponential Autocorrelation####
trans_int_3_exp <- lme(TRANS~DAY*GROUP,
                   random=~1|ID,
                   corr=corExp(form=~DAY),
                   method = "ML",
                   control=ctrl,
                   na.action=na.omit,
                   data=data.long.all)

#Check - Model Comparison
anova(trans_int_0_ri,
      trans_int_1_het_exp,
      trans_int_2_het,
      trans_int_2_het_day,
      trans_int_3_exp)
#Heterogeneous variances = Parsimonious
#Random effects become uninterpretable with covariance structures, hence removed.

####Model 5 - All Covariates####
lmer_trans_int_1 <- lmer(TRANS~DAY*GROUP+CRP+IL6+E2+P4+
                      (1|ID),
                    na.action=na.omit,
                    data=data.long.all)

#Check - ANOVA Table
Anova(lmer_trans_int_1, type=3)
#Sig: Day, Group 

#Check - Model Fit
resid_panel(lmer_trans_int_1, qqbands=TRUE)
#Fit: Approximately normal

####Model 6 - Reduced Model####
lmer_trans_int_2 <- lmer(TRANS~DAY*GROUP+E2+
                      (1|ID),
                    na.action=na.omit,
                    data=data.long.all)

#Check - ANOVA Table
Anova(lmer_trans_int_2,type=3)
#Sig: Day, Group, E2

#Check - Model Fit
resid_panel(lmer_trans_int_2, qqbands=TRUE)
#Fit: Approximately normal

####Model 7 - Reduced Model####
lmer_trans_int_3 <- lmer(TRANS~DAY*GROUP+
                           (1|ID),
                         na.action=na.omit,
                         data=data.long.all)

#Check - ANOVA Table
Anova(lmer_trans_int_3, type=3)
#Sig: Day, Group

#Check - Model Fit
resid_panel(lmer_trans_int_3, qqbands=TRUE)
#Fit: Approximately normal

#Check - Model Tables
tab_model(lmer_trans_int_3,
          lmer_trans_int_1,
          lmer_trans_int_2)

BIC(lmer_trans_int_3,
          lmer_trans_int_1,
          lmer_trans_int_2)

####IL6####

####Part 1 - Linear Mixed Modelling####

####Model 1 - Full Model####
il6_int_cov_ri_het_ac <- lme(IL6~DAY*GROUP+E2+P4,
                             random=~1|ID,
                             corr=corExp(form=~DAY),
                             weights = varIdent(form=~1|DAY*GROUP),
                             method = "ML",
                             control=ctrl,
                             na.action=na.omit,
                             data=data.long.all)

####Model 2 - Remove Heterogeneous variances####
il6_int_cov_ri_ac <- lme(IL6~DAY*GROUP+E2+P4,
                         random=~1|ID,
                         corr=corExp(form=~DAY),
                         method = "ML",
                         control=ctrl,
                         na.action=na.omit,
                         data=data.long.all)

####Model 3 - Remove Autocorrelation####
il6_int_cov_ri_het <- lme(IL6~DAY*GROUP+E2+P4,
                          random=~1|ID,
                          weights = varIdent(form=~1|DAY*GROUP),
                          method = "ML",
                          control=ctrl,
                          na.action=na.omit,
                          data=data.long.all)

###Model 4 - Remove AC and Heterogenous Autocorrelation####
il6_int_cov_ri <- lme(IL6~DAY*GROUP+E2+P4,
                      random=~1|ID,
                      method = "ML",
                      control=ctrl,
                      na.action=na.omit,
                      data=data.long.all)

#Check - Model Comparison
anova(il6_int_cov_ri,
      il6_int_cov_ri_het_ac,
      il6_int_cov_ri_ac,
      il6_int_cov_ri_het)
#Heterogeneous variances - Improves model

#Check - Model Comparison
anova(il6_int_cov_ri_het_ac,
      il6_int_cov_ri_het)
#Autocorrelation - Does not improve model

####Model 5 - Interaction Only Model####
il6_int_ri_het_reml <- lme(IL6~DAY*GROUP,
                           random=~1|ID,
                           weights = varIdent(form=~1|DAY*GROUP),
                           method = "REML",
                           control=ctrl,
                           na.action=na.omit,
                           data=data.long.all)

#Check - Model ANOVA
Anova(il6_int_ri_het_reml)
#Sig: Group & Day

#Check - Model Residual
resid_panel(il6_int_ri_het_reml, qqband=TRUE)
#Fit: Poor, long right tail

####Model 6 - Full Covariate Model####
il6_int_ri_het_e2_p4_reml <- lme(IL6~DAY*GROUP+E2+P4,
                                 random=~1|ID,
                                 weights = varIdent(form=~1|DAY*GROUP),
                                 method = "REML",
                                 control=ctrl,
                                 na.action=na.omit,
                                 data=data.long.all)

#Check - Model ANOVA
Anova(il6_int_ri_het_e2_p4_reml)
#Sig: Group & P4

#Check - Model Residual
resid_panel(il6_int_ri_het_e2_p4_reml, qqband=TRUE)
#Fit: Poor, long right tail

####Model 7 - Reduced Model fit with REML####
il6_int_ri_het_p4_reml <- lme(IL6~DAY*GROUP+P4,
                              random=~1|ID,
                              weights = varIdent(form=~1|DAY*GROUP),
                              method = "REML",
                              control=ctrl,
                              na.action=na.omit,
                              data=data.long.all)

#Check - Model ANOVA
Anova(il6_int_ri_het_p4_reml)
#Sig: Group, P4

#Check - Model Residual
resid_panel(il6_int_ri_het_p4_reml, qqband=TRUE)
#Fit: Poor, long right tail

####Part 2 - Linear Mixed Model - Log Transformed####

####Model 1 - Interaction Only Model####
ln_il6_int_ri_het_reml <- lme(log(IL6)~DAY*GROUP,
                              random=~1|ID,
                              method = "REML",
                              control=ctrl,
                              na.action=na.omit,
                              data=data.long.all)

#Check - Model ANOVA
anova(ln_il6_int_ri_het_reml, type="marginal")
#Sig: Group

#Check - Model Residual
resid_panel(ln_il6_int_ri_het_reml, qqband=TRUE)
#Fit: Acceptable

####Model 2 - Full Covariate Model####
ln_il6_int_ri_het_e2_p4_reml <- lme(log(IL6)~DAY*GROUP+E2+P4,
                                    random=~1|ID,
                                    method = "REML",
                                    control=ctrl,
                                    na.action=na.omit,
                                    data=data.long.all)

#Check - Model ANOVA
anova(ln_il6_int_ri_het_e2_p4_reml, type="marginal")
#Sig: None

#Check - Model Residual
resid_panel(ln_il6_int_ri_het_e2_p4_reml, qqband=TRUE)
#Fit: Acceptable

####Model 3 - Reduced Model fit with REML####
ln_il6_int_ri_het_p4_reml <- lme(log(IL6)~DAY*GROUP+P4,
                                 random=~1|ID,
                                 method = "REML",
                                 control=ctrl,
                                 na.action=na.omit,
                                 data=data.long.all)

#Check - Model ANOVA
anova(ln_il6_int_ri_het_p4_reml, type="marginal")
#Sig: None.

#Check - Model Residual
resid_panel(ln_il6_int_ri_het_p4_reml, qqband=TRUE)
#Fit: Acceptable.

#Check - Model Tables
tab_model(ln_il6_int_ri_het_reml,
          ln_il6_int_ri_het_e2_p4_reml,
          ln_il6_int_ri_het_p4_reml,
          show.aicc=TRUE)

BIC(ln_il6_int_ri_het_reml,
    ln_il6_int_ri_het_e2_p4_reml,
    ln_il6_int_ri_het_p4_reml)

####Part 3 - Generalised Linear Mixed Model - Gamma Distribution####

####Model 1 - Gamma Model - ANCOVA Model####
il6_int_ri_gamma <- glmer(IL6~DAY*GROUP+
                            (1|ID),
                          family=Gamma(link="log"),
                          na.action=na.omit,
                          data=data.long.all)

#Check - ANOVA Table
Anova(il6_int_ri_gamma, type=3)
#Sig: Group, Day x Group (P = 0.07)

#Check - Model Diagnostics
resid_panel(il6_int_ri_gamma)
#Fit: Approximately normal

#Plot - Model Effects
plot(allEffects(il6_int_ri_gamma))

#Check - Model Table
tab_model(il6_int_ri_gamma)

#Plot - Model Fixed Effects
plot_model (il6_int_ri_gamma, type="pred",terms=c("DAY", "GROUP"))

####Model 2 - Gamma Model - ANCOVA, adjusted for P4 & E2####
il6_int_cov_ri_gamma_e2_p4 <- glmer(IL6~DAY*GROUP+P4+E2+
                                      (1|ID),
                                    family=Gamma(link="log"),
                                    na.action=na.omit,
                                    data=data.long.all)

#NB - Model identifiability issues

#Check - Model ANOVA
Anova(il6_int_cov_ri_gamma_e2_p4, type=3)
#NS: Everything, P4 (P = 0.09)

#Check - Model Diagnostics
resid_panel(il6_int_cov_ri_gamma_e2_p4, qqbands=T)
#Fit: Approximately normal

####Model 3 - Gamma Model - ANCOVA, adjusted for P4 & E2####
il6_int_cov_ri_gamma_p4 <- glmer(IL6~DAY*GROUP+P4+
                                   (1|ID),
                                 family=Gamma(link="log"),
                                 na.action=na.omit,
                                 data=data.long.all)


#Check - Model ANOVA
Anova(il6_int_cov_ri_gamma_p4, type=3)
#NS: Everything, P4 (P = 0.09)

#Check - Model Diagnostics
resid_panel(il6_int_cov_ri_gamma_p4, qqbands=T)
#Fit: Approximately normal

#Check - Model Table
tab_model(il6_int_ri_gamma,
          il6_int_cov_ri_gamma_e2_p4,
          il6_int_cov_ri_gamma_p4)

BIC(il6_int_ri_gamma,
    il6_int_cov_ri_gamma_e2_p4,
    il6_int_cov_ri_gamma_p4)

####CRP####

####Part 1 - Linear Mixed Model - Raw Data####

####Model 1 - Heterogeneous variances (Day x Time) and Exponential Autocorrelation####
crp_int_1_het_exp <- lme(CRP~DAY*GROUP,
                         random=~1|ID,
                         corr=corExp(form=~DAY),
                         weights = varIdent(form=~1|DAY*GROUP),
                         method = "ML",
                         control=ctrl,
                         na.action=na.omit,
                         data=data.long.all)

####Model 2 - Heterogeneous Variances (Day x Time)####
crp_int_2_het <- lme(CRP~DAY*GROUP,
                     random=~1|ID,
                     weights = varIdent(form=~1|DAY*GROUP),
                     method = "ML",
                     control=ctrl,
                     na.action=na.omit,
                     data=data.long.all)

####Model 3 - Exponential Autocorrelation####
crp_int_3_exp <- lme(CRP~DAY*GROUP,
                     random=~1|ID,
                     corr=corExp(form=~DAY),
                     method = "ML",
                     control=ctrl,
                     na.action=na.omit,
                     data=data.long.all)

####Model 4 - Random Intercepts Model####
crp_int_4_ml <- lme(CRP~DAY*GROUP,
                    random=~1|ID,
                    method = "ML",
                    control=ctrl,
                    na.action=na.omit,
                    data=data.long.all)

#Check - Model Comaprison
anova(crp_int_1_het_exp,
      crp_int_2_het,
      crp_int_3_exp,
      crp_int_4_ml)
#Parsimonious: Heterogeneous variances and exponential autocorrelation

####Model 5 - Time x Day - REML####
crp_int_2_reml <- lme(CRP~DAY*GROUP,
                      random=~1|ID,
                      corr=corExp(form=~DAY),
                      weights = varIdent(form=~1|DAY*GROUP),
                      method = "REML",
                      control=ctrl,
                      na.action=na.omit,
                      data=data.long.all)

####Model 6 - Time x Day + E2 + P4 - REML####
crp_int_3_reml <- lme(CRP~DAY*GROUP+E2+P4,
                      random=~1|ID,
                      weights = varIdent(form=~1|DAY*GROUP),
                      method = "REML",
                      control=ctrl,
                      na.action=na.omit,
                      data=data.long.all)

####Model 7 - Time x Day + P4 - REML####
crp_int_4_reml <- lme(CRP~DAY*GROUP+P4,
                      random=~1|ID,
                      weights = varIdent(form=~1|DAY*GROUP),
                      method = "REML",
                      control=ctrl,
                      na.action=na.omit,
                      data=data.long.all)

#Check - Model Table
tab_model(crp_int_2_reml,
          crp_int_3_reml,
          crp_int_4_reml,
          show.aicc=TRUE)

#No covariate model - Parsimonious

#Check - Model ANOVA
Anova(crp_int_2_reml)  #Sig: Group
Anova(crp_int_3_reml)  #Sig: Group
Anova(crp_int_4_reml)  #Sig: Group

#Check - Model Summary Effects
summary(allEffects(crp_int_2_reml))

#Check - Model Diagnostics
resid_panel(crp_int_2_reml)
#NB - Long tailed at the extremes.
#NB - Outliers may need treating/removing
#NB - Attempt log or square root transformation

#Box Cox Transformation#
lm.crp <- lm(CRP~DAY*GROUP, data=data.long)
boxcox(lm.crp)
lm.crp$x[which.max(lm.crp$y)]
#Power of 0.4-0.6 seems optimal
#Consider SQRT distribution or Gamma GLM.

####Part 2 - Linear Mixed Model - Log Transformation

####Model 1 - Unadjusted Model - Log Transformed####
ln_crp_int_1_reml <- lme(LN_CRP~DAY*GROUP,
                         random=~1|ID,
                         corr=corExp(form=~DAY),
                         weights = varIdent(form=~1|DAY*GROUP),
                         method = "REML",
                         control=ctrl,
                         na.action=na.omit,
                         data=data.long.all)

#Check - Model ANOVA
anova(ln_crp_int_1_reml, type="marginal")
#Sig: Group

#Check - Model Residuals
resid_panel(ln_crp_int_1_reml, qqbands=T)
#Fit: Approximately normal residuals

####Model 2 - ANCOVA Model - Log Transformed####
ln_crp_int_2_reml <- lme(log(CRP)~DAY*GROUP+E2+P4,
                         random=~1|ID,
                         corr=corExp(form=~DAY),
                         weights = varIdent(form=~1|DAY*GROUP),
                         method = "REML",
                         control=ctrl,
                         na.action=na.omit,
                         data=data.long.all)

#Check - Model ANOVA
anova(ln_crp_int_2_reml, type="marginal")
#Sig: Group

#Check - Model Residuals
resid_panel(ln_crp_int_2_reml, qqbands=T)
#Fit: Reasonable

####Model 3 - ANCOVA Model - Log Transformed####
ln_crp_int_3_reml <- lme(log(CRP)~DAY*GROUP+P4,
                         random=~1|ID,
                         corr=corExp(form=~DAY),
                         weights = varIdent(form=~1|DAY*GROUP),
                         method = "REML",
                         control=ctrl,
                         na.action=na.omit,
                         data=data.long.all)

#Check - Model ANOVA
anova(ln_crp_int_3_reml, type="marginal")
#Sig: Group

#Check - Model Residuals
resid_panel(ln_crp_int_3_reml, qqbands=T)
#Fit: Reasonable

#Check - Model Table
tab_model(ln_crp_int_1_reml,
          ln_crp_int_2_reml,
          ln_crp_int_3_reml)


BIC(ln_crp_int_1_reml,
    ln_crp_int_2_reml,
    ln_crp_int_3_reml)

####Part 3 - Generalised Linear Mixed Model - Gamma Model####

####Model 1 - Gamma with log-link####
crp_int_1_gamma <- glmer(CRP~DAY*GROUP+
                           (1|ID),
                         family=Gamma(link="log"),
                         na.action=na.omit,
                         data=data.long.all)

#Check - Model ANOVA
Anova(crp_int_1_gamma,type=3)
#Sig: Group

#Check - Model Residuals
resid_panel(crp_int_1_gamma, qqbands=T)
#Fit: Reasonable

####Model 2 - Gamma with log-link - Covariate Adjusted####
crp_int_2_gamma <- glmer(CRP~DAY*GROUP+E2+P4+
                           (1|ID),
                         family=Gamma(link="log"),
                         na.action=na.omit,
                         data=data.long.all)

#NB - Model convergence issues

#Check - Model ANOVA
Anova(crp_int_2_gamma, type=3)

#Check - Model Residuals
resid_panel(crp_int_2_gamma, qqbands=T)
#Fit: Reasonable

#Check - Model Table
tab_model(crp_int_1_gamma,
          crp_int_2_gamma)

BIC(crp_int_1_gamma,
    crp_int_2_gamma)

####All Final Models####

#Model table of the final models#
#Model 1 - Day x Group + All Covariates
#Model 2 - Reduced Model/Parsimonious Model

####P4####
Anova(gamma_p4_int_2_reml)
tab_model(gamma_p4_int_2_reml)

####E2####
Anova(e2_int_1_reml)
tab_model(e2_int_1_reml)

####Serum Iron####
Anova(fe_int_ri_reml)
Anova(fe_int_cov_ri_reml)
Anova(fe_int_ri_crp_reml)

tab_model(fe_int_ri_reml,
          fe_int_cov_ri_reml,
          fe_int_ri_crp_reml,
          show.aicc=TRUE)

####Serum Ferritin####
Anova(fer_int_gamma)
Anova(fer_int_gamma_cov)
Anova(fer_int_gamma_red)

tab_model(fer_int_gamma,
          fer_int_gamma_cov,
          fer_int_gamma_red,
          show.aicc=TRUE)

####Transferrin####
Anova(lmer_trans_int_1)
Anova(lmer_trans_int_2)

tab_model(lmer_trans_int_1,
          lmer_trans_int_2,
          show.aicc=TRUE)

####CRP####
Anova(crp_int_1_gamma)
Anova(crp_int_2_gamma)
Anova(crp_int_3_gamma)

tab_model(crp_int_1_gamma,
          crp_int_2_gamma,
          crp_int_3_gamma,
          show.aicc=TRUE)

####IL6####
Anova(il6_int_ri_gamma)
Anova(il6_int_cov_ri_gamma_e2_p4)
Anova(il6_int_cov_ri_gamma_p4)

tab_model(il6_int_ri_gamma,
          il6_int_cov_ri_gamma_e2_p4,
          il6_int_cov_ri_gamma_p4,
          show.aicc=TRUE)