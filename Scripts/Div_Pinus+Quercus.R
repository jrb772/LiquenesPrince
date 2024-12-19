rm(list=ls())
set.seed(10)
pacman::p_load(readxl, dplyr, vegan,
               glmmTMB, fitdistrplus, performance, DHARMa, car, emmeans, ggeffects,
               ggplot2, ggpubr)

setwd("D:/CloudDrive/OneDrive - Universitat de les Illes Balears/Investigación/Liquenes Prince")

# DATA ####
Data <- read_excel("Fichas líquenes2.xlsx", sheet="COMPLETA")
Data$Tampon <- as.numeric(Data$Tampon)
Data <- Data %>% mutate_if(is.character, as.factor)
Data2 <- subset(Data, Species == "Q. robur" | Species == "P. pinaster")

alpha <- cbind(
  Data2[, 1:11],
  # Shannon = diversity(Data2[, 12:length(Data2)]),
  # Simpson = diversity(Data2[, 12:length(Data2)], "simpson"),
  # InvSimpson = diversity(Data2[, 12:length(Data2)], "inv"),
  Rich = specnumber(Data2[, 12:length(Data2)])
)

NatAr <- read_excel("Fichas líquenes2.xlsx", sheet="SIOSE")
NatAr$Type <- as.factor(NatAr$Type)

SIOSE_H <- read_excel("Fichas líquenes2.xlsx", sheet="SIOSE_H")

Water_dist <- read_excel("Fichas líquenes2.xlsx", sheet="Water_dist")

dataset <- merge(alpha, subset(NatAr, NatAr$Type == "Natural"), by="ID")
dataset <- merge(dataset, SIOSE_H, by="ID")
dataset <- merge(dataset, Water_dist, by="ID")

# Plots alpha ####
plot.rich <-ggplot(alpha, aes(x = Species, y = Rich, colour = Species)) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  ylab("Species Richness") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.rich

# Models alpha ####
fitdistrplus::descdist(alpha$Rich, discrete = T)

m.Rich <- glmmTMB(Rich ~ Species * (H_1Km + Water_min) + (1|ID), 
              data=dataset, family=poisson())
m.Rich <- glmmTMB(Rich ~ Species * (H_500m + Water_min) + (1|ID), 
                  data=dataset, family=poisson())
m.Rich <- glmmTMB(Rich ~ Species * (H_250m + Water_min + Water_Ret + pH) + (1|ID), 
                  data=dataset, family=poisson())
drop1(m.Rich)

check_overdispersion(m.Rich)
check_collinearity(m.Rich)
check_model(m.Rich)
simulateResiduals(m.Rich, plot=T, quantreg=T)
model_performance(m.Rich)

summary(m.Rich)
Anova(m.Rich, type=3)

ggpredict(m.Rich, c("H_250m", "Species")) |> plot()
ggpredict(m.Rich, c("pH")) |> plot()
ggpredict(m.Rich, c("Water_min")) |> plot()
ggpredict(m.Rich, c("Water_Ret")) |> plot()
ggpredict(m.Rich, c("Water_Ret", "Species")) |> plot()

ph1 <- emmeans(m.Rich, ~ Species, type="response", adjust="tukey")
CLD1 <- multcomp::cld(ph1, alpha = 0.05, Letters = letters, decreasing = T)
CLD1$.group <- gsub("\\s+","",CLD1$.group,perl=T)
CLD1
