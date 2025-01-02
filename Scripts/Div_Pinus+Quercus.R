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

# Models alpha ####
fitdistrplus::descdist(alpha$Rich, discrete = T)

m.Rich_1Km <- glmmTMB(Rich ~ Species * (H_1Km + Water_min) + (1|ID), 
              data=dataset, family=poisson())
m.Rich_500m <- glmmTMB(Rich ~ Species * (H_500m + Water_min) + (1|ID), 
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
Anova(m.Rich_500m, type=3)
Anova(m.Rich_1Km, type=3)

ggpredict(m.Rich, c("H_250m", "Species"), bias_correction = T) |> plot()
ggpredict(m.Rich, c("pH")) |> plot()
ggpredict(m.Rich, c("Water_min")) |> plot()
ggpredict(m.Rich, c("Water_Ret")) |> plot()
ggpredict(m.Rich, c("Water_Ret", "Species")) |> plot()

ph1 <- emmeans(m.Rich, ~ Species, type="response", adjust="tukey")
CLD1 <- multcomp::cld(ph1, alpha = 0.05, Letters = letters, decreasing = T)
CLD1$.group <- gsub("\\s+","",CLD1$.group,perl=T)
CLD1

plot.rich <- ggplot() +
  geom_violin(data = dataset, aes(x = Species, y = Rich, fill = Species),
                                  alpha = .4, trim = FALSE) +
  geom_point(data = CLD1, aes(x = Species, y = rate),
             size = 3, color="black") +
  geom_errorbar(data = CLD1, aes(x = Species, ymin=rate-SE, ymax=rate+SE),
                width=0.2, color="black") +
  geom_text(data = CLD1, aes(x = Species, y=rate+SE, label=.group),
            size=5, vjust=-1.5, color="black") +
  ylab("Species Richness") + xlab("") +
  scale_y_continuous(limits=c(0,40), breaks=seq(0,50, by=5)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() + theme(legend.position = "none")
plot.rich

# NMDS ####
com <- as.matrix(Data2[, 12:length(Data2)])
env <- merge(Data2[,1:11], subset(NatAr, NatAr$Type == "Natural"), by="ID")
env <- merge(env, SIOSE_H, by="ID")
env <- merge(env, Water_dist, by="ID")

nmds <- metaMDS(com, distance="jaccard", k=2, trymax=100)
plot(nmds, display="sites")
fit <- envfit(nmds, env, permutations = 9999, na.rm = TRUE)
fit
plot(fit)

data.scores <- as.data.frame(scores(nmds)$sites)
data.scores <- cbind(env,data.scores)

en_coord_cont <- as.data.frame(scores(fit, "vectors"))[c(2:4,11),] * ordiArrowMul(fit)
en_coord_cat <- as.data.frame(scores(fit, "factors")) * ordiArrowMul(fit)

plot_beta <- ggplot(data.scores, aes(NMDS1, NMDS2)) +
  geom_point(size=4, aes(shape=Species, colour=Species)) +
  geom_segment(aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
               data=en_coord_cont, size=1, alpha=0.5, colour="black") +
  geom_text(data=en_coord_cont, aes(x=NMDS1, y=NMDS2), colour="black",
            fontface = "bold", label = row.names(en_coord_cont)) +
  theme_minimal()+
  labs(x="NMDS1", y="NMDS2", shape="Species") +   theme(legend.position = "top")
plot_beta

# CCA ####
spe.rda <- rda(com ~ Species * (Area_1km + pH + Crust_estab + Crust_text + Tampon +
                 Water_Ret + Water_min), data=env)
summary(spe.rda)

# Forward selection of variables:
fwd.sel <- ordiR2step(rda(com ~ 1, data = env), # lower model limit (simple!)
                      scope = formula(spe.rda), # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = TRUE, # can't surpass the "full" model's R2
                      pstep = 10000,
                      trace = FALSE) # change to TRUE to see the selection process!
# Check the new model with forward-selected variables
fwd.sel$call
# Write our new model
spe.rda.signif <- cca(com ~  Species + (Area_1km + pH + Tampon +
                                          Water_Ret + Water_min),
                      data = env)
# check the adjusted R2 (corrected for the number of
# explanatory variables)
RsquareAdj(spe.rda.signif)
anova.cca(spe.rda.signif, step = 1000)
anova.cca(spe.rda.signif, step = 1000, by = "term")

# Type 1 scaling
ordiplot(spe.rda.signif, scaling = 1)
# Type 2 scaling
ordiplot(spe.rda.signif, scaling = 2)

site.scores <- scores(spe.rda.signif, display = "sites", choices = c(1,2))
species.scores <- scores(spe.rda.signif, display = "species", choices = c(1,2))
env.scores <- as.data.frame(scores(spe.rda.signif, display = "bp", choices = c(1,2))) * ordiArrowMul(spe.rda.signif, display="bp")

# Create a data frame for plotting
plot.data <- data.frame(
  CCA1 = site.scores[,1],
  CCA2 = site.scores[,2],
  Species = env$Species
)

ggplot() +
  geom_point(data = plot.data, aes(x = CCA1, y = CCA2, color = Species), size = 3) +
  geom_segment(data = env.scores, aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.2, "cm"))) + # Add vector arrows for environmental variables
  geom_text(data = env.scores, aes(x = CCA1, y = CCA2, label = rownames(env.scores)),
            hjust = 1, vjust = -0.5, size = 3) + # Add labels for environmental variables
  xlab("CCA1") +
  ylab("CCA2") +
  theme_minimal() +
  scale_x_continuous(limits=c(-5,5), breaks=seq(-10,10, by=2)) +
  scale_y_continuous(limits=c(-5,5), breaks=seq(-10,10, by=2.5)) +
  theme(legend.position = "top")

