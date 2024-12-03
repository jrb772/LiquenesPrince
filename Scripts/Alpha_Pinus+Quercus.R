rm(list=ls())
set.seed(10)
pacman::p_load(readxl, dplyr, vegan,
               glmmTMB, fitdistrplus, performance, DHARMa, car, emmeans, ggeffects,
               ggplot2, ggpubr)

# setwd("D:/CloudDrive/OneDrive - Universitat de les Illes Balears/Investigación/Liquenes Prince")

# DATA ####
Data <- read_excel("Fichas líquenes2.xlsx", sheet="COMPLETA")
Data$Tampon <- as.numeric(Data$Tampon)
Data <- Data %>% mutate_if(is.character, as.factor)
Data2 <- subset(Data, Species == "Q. robur" | Species == "P. pinaster")

alpha <- cbind(
  Data2[, 1:11],
  Shannon = diversity(Data2[, 12:length(Data2)]),
  Simpson = diversity(Data2[, 12:length(Data2)], "simpson"),
  InvSimpson = diversity(Data2[, 12:length(Data2)], "inv"),
  Rich = specnumber(Data2[, 12:length(Data2)])
)

# Plots alpha ####
plot.shan <- ggplot(alpha, aes(x = Species, y = Shannon, colour = Species)) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  ylab("Shannon's H'") + 
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))

plot.simp <- ggplot(alpha, aes(x = Species, y = Simpson, colour = Species)) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  ylab("Simpson's diversity") + 
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))

plot.inv <- ggplot(alpha, aes(x = Species, y = InvSimpson, colour = Species)) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  ylab("Inverse Simpson's") + 
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))

plot.rich <-ggplot(alpha, aes(x = Species, y = Rich, colour = Species)) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  ylab("Species Richness") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))

legend <- get_legend(plot.even)

ggarrange(plot.shan + theme(legend.position = "none"),
          plot.simp + theme(legend.position = "none"),
          plot.inv + theme(legend.position = "none"),
          plot.rich + theme(legend.position = "none"),
          ncol = 4)


# Models alpha ####
fitdistrplus::descdist(alpha$Shannon, discrete = FALSE)
fitdistrplus::descdist(alpha$Simpson, discrete = FALSE)
fitdistrplus::descdist(alpha$InvSimpson, discrete = FALSE)
fitdistrplus::descdist(alpha$Rich, discrete = FALSE)

H1 <- glmmTMB(Shannon ~ Species + (pH + Antrop + Water) + (1|ID), 
              data=alpha, family=gaussian())

check_overdispersion(H1)
check_collinearity(H1)
check_model(H1)
simulateResiduals(H1, plot=T, quantreg=T)
model_performance(H1)

summary(H1)
Anova(H1, type=3)
