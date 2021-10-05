library(ggplot2)
library(drc)
library(dplyr)
library(extrafont)
library(gridExtra)
library(RGraphics)
library(data.table)
library(ggpubr)
library(viridis)
library(scales)
font_import()

theme_science <- function (base_size = 12, base_family = "Arial Black") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.border = element_blank(), axis.line = element_line(colour = "black", size=2), 
          panel.grid.major = element_line(), panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour= "black", size=1),  axis.line.y = element_line(colour= "black", size=1),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_line(), 
          panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), 
          strip.background = element_rect(colour = "black", 
                                          size = 0.5), legend.key = element_blank())
}



scaleFUN <- function(x) sprintf("%.2f", x)

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbviridis <- c("#440154FF", "#31688EFF", "#35B779FF", "#E69F00")


q_colors =  4 # for no particular reason
v_colors =  viridis(q_colors, option = "D")

setwd("/Users/cottr/Box Sync/ADAR_inhibitors/")

#read in viability data for 8-azaadenosine
combined_aza_overexpression <- read.delim("combined_aza_overexpression.txt", header = TRUE)

#fit ll.4 model to data
aza <- drm(Viability ~ Concentration, Cell_line, data = combined_aza_overexpression, fct = LL.4(),
           pmodels = data.frame(Cell_line, Cell_line, 1, Cell_line))
summary(aza)
#find EC50 from model
ED(aza, 50)

#group viability data by cell line and concentration
combined_aza_overexpression <- dplyr::group_by(combined_aza_overexpression, Cell_line, Concentration)

#summarise data
combined_aza_overexpression_sum <- dplyr::summarise(combined_aza_overexpression, Viability_mean = mean(Viability), SD = sd(Viability))

#plot dose response curves and save
ggplot(combined_aza_overexpression,aes(Concentration, Viability, colour = Cell_line)) + geom_point() + 
  geom_smooth(method = drm, method.args = list(fct = L.4(fixed = c(NA,NA,1,NA))), se = FALSE) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 50), labels=c(0.001, 0.01, 0.1, 1, 10, 50)) + theme_science() + 
  geom_pointrange(data = combined_aza_overexpression_sum, aes(x = Concentration, y = Viability_mean, ymin = Viability_mean - SD, ymax = Viability_mean + SD), size = 0.7) +
  scale_colour_manual(values = cbviridis) + 
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab(expression(paste("8-azaadenosine (", mu, "M)"))) + ylab("Relative Viability") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1))

ggplot(combined_aza_overexpression,aes(Concentration, Viability, colour = Cell_line)) + geom_point() + 
  geom_smooth(method = drm, method.args = list(fct = L.4(fixed = c(NA,NA,1,NA))), se = FALSE) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 50), labels=c(0.001, 0.01, 0.1, 1, 10, 50)) + theme_science() + 
  geom_pointrange(data = combined_aza_overexpression_sum, aes(x = Concentration, y = Viability_mean, ymin = Viability_mean - SD, ymax = Viability_mean + SD), size = 0.7) +
  scale_colour_manual(values = cbviridis, labels = c("HCC1806\nEC50 = 0.91", "MDA-MB-468\nEC50 = 1.64", "MCF-7\nEC50 = 1.27", "SKBR-3\nEC50 = 0.92")) + 
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab(expression(paste("8-azaadenosine (", mu, "M)"))) + ylab("Relative Viability") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) 
ggsave("aza_combined.tiff", height = 4, width = 5.9)

#read in viability data for 8-azaadenosine
combined_chloro <- read.delim("combined_chloro.txt", header = TRUE)

#fit ll.4 model to data
chloro <- drm(Viability ~ Concentration, Cell_line, data = combined_chloro, fct = LL.4(), 
           pmodels=data.frame(Cell_line, Cell_line, 1, Cell_line))
summary(chloro)

#find EC50 from model
ED(chloro, 50)

#group viability data by cell line and concentration
combined_chloro <- dplyr::group_by(combined_chloro, Cell_line, Concentration)

#summarise data
combined_chloro_sum <- dplyr::summarise(combined_chloro, Viability_mean = mean(Viability), SD = sd(Viability))

#plot dose response curves and save
ggplot(combined_chloro,aes(Concentration, Viability, colour = Cell_line)) + geom_point() + 
  geom_smooth(method = drm, method.args = list(fct = L.4(fixed = c(NA, NA, 1,NA))), se = FALSE) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 50), labels=c(0.001, 0.01, 0.1, 1, 10, 50)) + theme_science() + 
  geom_pointrange(data = combined_chloro_sum, aes(x = Concentration, y = Viability_mean, ymin = Viability_mean - SD, ymax = Viability_mean + SD), size = 0.7) +
  scale_colour_manual(values = cbviridis) + 
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab(expression(paste("8-chloroadenosine (", mu, "M)"))) + ylab("Relative Viability") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) 

ggplot(combined_chloro,aes(Concentration, Viability, colour = Cell_line)) + geom_point() + 
  geom_smooth(method = drm, method.args = list(fct = L.4(fixed = c(NA, NA, 1,NA))), se = FALSE) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 50), labels=c(0.001, 0.01, 0.1, 1, 10, 50)) + theme_science() + 
  geom_pointrange(data = combined_chloro_sum, aes(x = Concentration, y = Viability_mean, ymin = Viability_mean - SD, ymax = Viability_mean + SD), size = 0.7) +
  scale_colour_manual(values = cbviridis, labels = c("HCC1806\nEC50 = 0.60", "MDA-MB-468\nEC50 = 0.55", "MCF-7\nEC50 = 0.35", "SK-BR-3\nEC50 = 0.46")) + 
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab(expression(paste("8-chloroadenosine (", mu, "M)"))) + ylab("Relative Viability") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) 
ggsave("chloro_combined.tiff", height = 4, width = 5.9)


#read in viability data for SKBR3
skbr3_combined <- read.delim("combined_skbr3.txt", header = TRUE)

#subset based on Drug
skbr3_combined_aza_overexpression <- subset(skbr3_combined, Drug == "aza")
skbr3_combined_chloro <- subset(skbr3_combined, Drug == "chloro")

#fit LL.4 model to data
chloro <- drm(Viability ~ Concentration, Cell_line, data = skbr3_combined_chloro, fct = LL.4(), 
              pmodels=data.frame(Cell_line, Cell_line, 1, Cell_line, Cell_line))
summary(chloro)

#find EC50 from model
ED(chloro, 50)

#group and summarise data
skbr3_combined_chloro_sum <- dplyr::group_by(skbr3_combined_chloro, Cell_line, Concentration)

skbr3_combined_chloro_sum <- dplyr::summarise(skbr3_combined_chloro_sum, Viability_mean = mean(Viability), SD = sd(Viability))

#plot dose response curves and save
ggplot(skbr3_combined_chloro, aes(Concentration, Viability, colour = Cell_line)) + geom_point(alpha = 0.5) + 
  geom_smooth(method = drm, method.args = list(fct = L.4(fixed = c(NA,NA,1,NA))), se = FALSE) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 50), labels=c(0.001, 0.01, 0.1, 1, 10, 50)) + theme_science() +
  geom_pointrange(data = skbr3_combined_chloro_sum, aes(x = Concentration, y = Viability_mean, ymin = Viability_mean - SD, ymax = Viability_mean + SD), size = 0.7) +
  scale_colour_manual(values = c("black", "grey"), labels = c("Scr" = "shSCR\nEC50 = 0.81", "ADAR" = "shADAR\nEC50 = 0.38")) + 
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab(expression(paste("SK-BR-3 8-chloroadenosine (", mu, "M)"))) + ylab("Relative Viability") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) 
ggsave("skbr3_chloro.tiff", height = 4, width = 5.7)

#fit LL.4 model to data
aza <- drm(Viability ~ Concentration, Cell_line, data = skbr3_combined_aza_overexpression, fct = LL.4(), 
              pmodels=data.frame(Cell_line, Cell_line, 1, Cell_line, Cell_line))
summary(aza)

#find EC50 from model
ED(aza, 50)

#group and summarise data
skbr3_combined_aza_overexpression_sum <- dplyr::group_by(skbr3_combined_aza_overexpression, Cell_line, Concentration)

skbr3_combined_aza_overexpression_sum <- dplyr::summarise(skbr3_combined_aza_overexpression_sum, Viability_mean = mean(Viability), SD = sd(Viability))

#plot dose response curves and save
ggplot(skbr3_combined_aza_overexpression, aes(Concentration, Viability, colour = Cell_line)) + geom_point(alpha = 0.5) + 
  geom_smooth(method = drm, method.args = list(fct = L.4(fixed = c(NA,NA,1,NA))), se = FALSE) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 50), labels=c(0.001, 0.01, 0.1, 1, 10, 50)) + theme_science() +
  geom_pointrange(data = skbr3_combined_aza_overexpression_sum, aes(x = Concentration, y = Viability_mean, ymin = Viability_mean - SD, ymax = Viability_mean + SD), size = 0.7) +
  scale_colour_manual(values = c("black", "grey"), labels = c("Scr" = "shSCR\nEC50 = 0.74", "ADAR" = "shADAR\nEC50 = 0.71")) + 
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab(expression(paste("SK-BR-3 8-azaadenosine (", mu, "M)"))) + ylab("Relative Viability") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) 
ggsave("skbr3_aza.tiff", height = 4, width = 5.7)


#read in viability data for SKBR3
mcf7_combined <- read.delim("combined_mcf7.txt", header = TRUE)

#subset based on drug
mcf7_combined_aza_overexpression <- subset(mcf7_combined, Drug == "aza")
mcf7_combined_chloro <- subset(mcf7_combined, Drug == "chloro")

#fit LL.4 model to data
chloro <- drm(Viability ~ Concentration, Cell_line, data = mcf7_combined_chloro, fct = LL.4(), 
              pmodels=data.frame(Cell_line, Cell_line, 1, Cell_line, Cell_line))
summary(chloro)

#find EC50 from model
ED(chloro, 50)

#group and summarise data
mcf7_combined_chloro_sum <- dplyr::group_by(mcf7_combined_chloro, Cell_line, Concentration)

mcf7_combined_chloro_sum <- dplyr::summarise(mcf7_combined_chloro_sum, Viability_mean = mean(Viability), SD = sd(Viability))

#plot dose response curves and save
ggplot(mcf7_combined_chloro, aes(Concentration, Viability, colour = Cell_line)) + geom_point(alpha = 0.5) + 
  geom_smooth(method = drm, method.args = list(fct = L.4(fixed = c(NA,NA,1,NA))), se = FALSE) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 50), labels=c(0.001, 0.01, 0.1, 1, 10, 50)) + theme_science() +
  geom_pointrange(data = mcf7_combined_chloro_sum, aes(x = Concentration, y = Viability_mean, ymin = Viability_mean - SD, ymax = Viability_mean + SD), size = 0.7) +
  scale_colour_manual(values = c("black", "grey"), labels = c("scr" = "shSCR\nEC50 = 0.46", "ADAR" = "shADAR\nEC50 = 0.41")) + 
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab(expression(paste("MCF-7 8-chloroadenosine (", mu, "M)"))) + ylab("Relative Viability") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) 
ggsave("mcf7_chloro.tiff", height = 4, width = 5.7)


#fit LL.4 model to data
aza <- drm(Viability ~ Concentration, Cell_line, data = mcf7_combined_aza_overexpression, fct = LL.4(), 
           pmodels=data.frame(Cell_line, Cell_line, 1, Cell_line, Cell_line))
summary(aza)

#find EC50 from model
ED(aza, 50)

#group and summarise the data
mcf7_combined_aza_overexpression_sum <- dplyr::group_by(mcf7_combined_aza_overexpression, Cell_line, Concentration)

mcf7_combined_aza_overexpression_sum <- dplyr::summarise(mcf7_combined_aza_overexpression_sum, Viability_mean = mean(Viability), SD = sd(Viability))

#plot dose response curves and save
ggplot(mcf7_combined_aza_overexpression, aes(Concentration, Viability, colour = Cell_line)) + geom_point(alpha = 0.5) + 
  geom_smooth(method = drm, method.args = list(fct = L.4(fixed = c(NA,NA,1,NA))), se = FALSE) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 50), labels=c(0.001, 0.01, 0.1, 1, 10, 50)) + theme_science() +
  geom_pointrange(data = mcf7_combined_aza_overexpression_sum, aes(x = Concentration, y = Viability_mean, ymin = Viability_mean - SD, ymax = Viability_mean + SD), size = 0.7) +
  scale_colour_manual(values = c("black", "grey"), labels = c("scr" = "shSCR\nEC50 = 1.81", "ADAR" = "shADAR\nEC50 = 2.04")) + 
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab(expression(paste("MCF-7 8-azaadenosine (", mu, "M)"))) + ylab("Relative Viability") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) 
ggsave("mcf7_aza.tiff", height = 4, width = 5.7)




#read in viability data for 8-azaadenosine
combined_aza_overexpression <- read.delim("combined_aza_overexpression.txt", header = TRUE)

#fit ll.4 model to data
aza <- drm(Viability ~ Concentration, Cell_line, data = combined_aza_overexpression, fct = LL.4(),
           pmodels = data.frame(Cell_line, Cell_line, 1, Cell_line))
summary(aza)
#find EC50 from model
ED(aza, 50)

#group viability data by cell line and concentration
combined_aza_overexpression <- dplyr::group_by(combined_aza_overexpression, Cell_line, Concentration)

#summarise data
combined_aza_overexpression_sum <- dplyr::summarise(combined_aza_overexpression, Viability_mean = mean(Viability), SD = sd(Viability))

#plot dose response curves and save
ggplot(combined_aza_overexpression,aes(Concentration, Viability, colour = Cell_line)) + geom_point() + 
  geom_smooth(method = drm, method.args = list(fct = L.4(fixed = c(NA,NA,1,NA))), se = FALSE) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 50), labels=c(0.001, 0.01, 0.1, 1, 10, 50)) + theme_science() + 
  geom_pointrange(data = combined_aza_overexpression_sum, aes(x = Concentration, y = Viability_mean, ymin = Viability_mean - SD, ymax = Viability_mean + SD), size = 0.7) +
  scale_colour_manual(values = cbviridis) + 
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab(expression(paste("8-azaadenosine (", mu, "M)"))) + ylab("Relative Viability") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1))

ggplot(combined_aza_overexpression,aes(Concentration, Viability, colour = Cell_line)) + geom_point() + 
  geom_smooth(method = drm, method.args = list(fct = L.4(fixed = c(NA,NA,1,NA))), se = FALSE) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 50), labels=c(0.001, 0.01, 0.1, 1, 10, 50)) + theme_science() + 
  geom_pointrange(data = combined_aza_overexpression_sum, aes(x = Concentration, y = Viability_mean, ymin = Viability_mean - SD, ymax = Viability_mean + SD), size = 0.7) +
  scale_colour_manual(values = cbviridis, labels = c("EV = 0.52", "p110 = 0.61", "p150 = 0.58")) + 
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab(expression(paste("8-azaadenosine (", mu, "M)"))) + ylab("Relative Viability") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) 
ggsave("aza_combined_overexpression.tiff", height = 4, width = 5.9)




#read in viability data for 8-chloroadenosine
combined_chloro_overexpression <- read.delim("combined_chloro_overexpression.txt", header = TRUE)

#fit ll.4 model to data
chloro <- drm(Viability ~ Concentration, Cell_line, data = combined_chloro_overexpression, fct = LL.4(),
           pmodels = data.frame(Cell_line, Cell_line, 1, Cell_line))
summary(chloro)
#find EC50 from model
ED(chloro, 50)

#group viability data by cell line and concentration
combined_chloro_overexpression <- dplyr::group_by(combined_chloro_overexpression, Cell_line, Concentration)

#summarise data
combined_chloro_overexpression_sum <- dplyr::summarise(combined_chloro_overexpression, Viability_mean = mean(Viability), SD = sd(Viability))

#plot dose response curves and save
ggplot(combined_chloro_overexpression,aes(Concentration, Viability, colour = Cell_line)) + geom_point() + 
  geom_smooth(method = drm, method.args = list(fct = L.4(fixed = c(NA,NA,1,NA))), se = FALSE) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 50), labels=c(0.001, 0.01, 0.1, 1, 10, 50)) + theme_science() + 
  geom_pointrange(data = combined_chloro_overexpression_sum, aes(x = Concentration, y = Viability_mean, ymin = Viability_mean - SD, ymax = Viability_mean + SD), size = 0.7) +
  scale_colour_manual(values = cbviridis) + 
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab(expression(paste("8-chloroadenosine (", mu, "M)"))) + ylab("Relative Viability") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1))

ggplot(combined_chloro_overexpression,aes(Concentration, Viability, colour = Cell_line)) + geom_point() + 
  geom_smooth(method = drm, method.args = list(fct = L.4(fixed = c(NA,NA,1,NA))), se = FALSE) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 50), labels=c(0.001, 0.01, 0.1, 1, 10, 50)) + theme_science() + 
  geom_pointrange(data = combined_chloro_overexpression_sum, aes(x = Concentration, y = Viability_mean, ymin = Viability_mean - SD, ymax = Viability_mean + SD), size = 0.7) +
  scale_colour_manual(values = cbviridis, labels = c("EV = 0.77", "p110 = 0.88", "p150 = 0.91")) + 
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab(expression(paste("8-chloroadenosine (", mu, "M)"))) + ylab("Relative Viability") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) 
ggsave("chloro_combined_overexpression.tiff", height = 4, width = 5.9)



################################################################################################
#read in foci formation data
FF <- read.delim("FF_analysis.txt")

#make an empty data.frame for DMSO_mean
DMSO_mean <- NULL

#determine the mean foci area (percent.area) for DMSO for each cell line
for (i in unique(FF$Cell_line)) {
  m <- mean(FF$Percent.Area[FF$Drug == "DMSO" & FF$Cell_line == i])
  DMSO_mean <- rbind(DMSO_mean, data.frame(i, m))
}

#subset FF data by cell line
FF_SKBR3 <- subset(FF, Cell_line == "SKBR3")
FF_MCF7 <- subset(FF, Cell_line == "MCF7")
FF_MB468 <- subset(FF, Cell_line == "MB468")
FF_HCC1806 <- subset(FF, Cell_line == "HCC1806")

#set percent.area relative to mean percent area for DMSO for each cell line
FF_SKBR3$relative <- FF_SKBR3$Percent.Area/DMSO_mean$m[DMSO_mean$i == "SKBR3"]
FF_MCF7$relative <- FF_MCF7$Percent.Area/DMSO_mean$m[DMSO_mean$i == "MCF7"]
FF_MB468$relative <- FF_MB468$Percent.Area/DMSO_mean$m[DMSO_mean$i == "MB468"]
FF_HCC1806$relative <- FF_HCC1806$Percent.Area/DMSO_mean$m[DMSO_mean$i == "HCC1806"]

#bind together data from all cell lines
FF_relative <- rbind(FF_SKBR3, FF_MCF7, FF_HCC1806, FF_MB468)

#group and summarise data
FF_relative_sum <- dplyr::group_by(FF_relative, Cell_line, Drug)

FF_relative_sum <- dplyr::summarise(FF_relative_sum, Relative_mean = mean(relative), Relative_sd = sd(relative))


#make plots and save
ggplot(FF_relative, aes(Drug, relative, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("DMSO", "AZA", "Chloro"), labels = c("DMSO", expression(paste("1 ",mu,M," 8-aza")), expression(paste("0.5 ",mu,M, " 8-chloro")))) + labs(x = "", y = "Relative Foci Area", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,1.5), expand = c(0, 0)) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  scale_colour_manual(values = cbviridis) +
  geom_col(data = FF_relative_sum, aes(Drug, Relative_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = FF_relative_sum, aes(x = Drug, y = Relative_mean, ymin = Relative_mean - Relative_sd, ymax = Relative_mean + Relative_sd), width=0.2, position = position_dodge(width = 0.7))

ggplot(FF_relative, aes(Drug, relative, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("DMSO", "AZA", "Chloro"), labels = c("DMSO", expression(paste("1 ",mu,M," 8-aza")), expression(paste("0.5 ",mu,M, " 8-chloro")))) + labs(x = "", y = "Relative Foci Area", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,1.5), expand = c(0, 0)) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  scale_colour_manual(values = cbviridis, labels = c("HCC1806", "MDA-MB-468", "MCF-7", "SK-BR-3")) +
  geom_col(data = FF_relative_sum, aes(Drug, Relative_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = FF_relative_sum, aes(x = Drug, y = Relative_mean, ymin = Relative_mean - Relative_sd, ymax = Relative_mean + Relative_sd), width=0.2, position = position_dodge(width = 0.7))
ggsave("FF.tiff", height = 4, width = 6)


#read western (immunoblot) data for knockdowns
western_knockdown <- read.delim("western_knockdowns.txt")

western_knockdown <- na.omit(western_knockdown)

#group and summarise data
western_knockdown_sum <- dplyr::group_by(western_knockdown, Cell_line, Sample)

western_knockdown_sum <- dplyr::summarise(western_knockdown_sum, ADAR_mean = mean(ADAR), 
                                          ADAR_sd = sd(ADAR), pPKR_mean = mean(pPKR.PKR), 
                                          pPKR_sd = sd(pPKR.PKR),
                                          PKR_mean = mean(PKR), PKR_sd = sd(PKR))

#make plots and save
ggplot(western_knockdown, aes(Sample, ADAR, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shADAR")) + labs(x = "", y = "Relative ADAR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,1.5), expand = c(0, 0)) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = western_knockdown_sum, aes(Sample, ADAR_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = western_knockdown_sum, aes(x = Sample, y = ADAR_mean, ymin = ADAR_mean - ADAR_sd, ymax = ADAR_mean + ADAR_sd), width=0.2, position = position_dodge(width = 0.7))
ggsave("western_knockdowns_adar.tiff", height = 4, width = 4)




ggplot(western_knockdown, aes(Sample, pPKR.PKR, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shADAR")) + labs(x = "", y = "Relative pPKR/PKR", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,25), expand = c(0, 0)) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = western_knockdown_sum, aes(Sample, pPKR_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = western_knockdown_sum, aes(x = Sample, y = pPKR_mean, ymin = pPKR_mean - pPKR_sd, ymax = pPKR_mean + pPKR_sd), width=0.2, position = position_dodge(width = 0.7))
ggsave("western_knockdowns_pPKR.tiff", height = 4, width = 4)


ggplot(western_knockdown, aes(Sample, PKR, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shADAR")) + labs(x = "", y = "Relative PKR Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,1.5), expand = c(0, 0)) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = western_knockdown_sum, aes(Sample, PKR_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = western_knockdown_sum, aes(x = Sample, y = PKR_mean, ymin = PKR_mean - PKR_sd, ymax = PKR_mean + PKR_sd), width=0.2, position = position_dodge(width = 0.7))
ggsave("western_knockdowns_PKR.tiff", height = 4, width = 4)


#read western (immunoblot) data for 8-azaadenosine
mb468_aza_western <- read.delim("western_aza_mb468.txt")
hcc1806_aza_western <- read.delim("western_aza_hcc1806.txt")

#bind together data from both cell lines
aza_western <- rbind(mb468_aza_western, hcc1806_aza_western)

aza_western <- na.omit(aza_western)

#group and summarise data
western_aza_sum <- dplyr::group_by(aza_western, Cell_line, Sample)

western_aza_sum <- dplyr::summarise(western_aza_sum, ADAR_mean = mean(ADAR), 
                                    ADAR_sd = sd(ADAR), pPKR_mean = mean(pPKR.PKR), 
                                    pPKR_sd = sd(pPKR.PKR),
                                    PKR_mean = mean(PKR), PKR_sd = sd(PKR))
#make plots and save
ggplot(aza_western, aes(Sample, ADAR, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("DMSO", "1_aza", "10_aza"), labels = c("DMSO", expression(paste("1 ", mu, "M")), expression(paste("10 ", mu, "M")))) + labs(x = "8-azaadenosine", y = "Relative ADAR Abundance", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,2), expand = c(0, 0)) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = western_aza_sum, aes(Sample, ADAR_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = western_aza_sum, aes(x = Sample, y = ADAR_mean, ymin = ADAR_mean - ADAR_sd, ymax = ADAR_mean + ADAR_sd), width=0.2, position = position_dodge(width = 0.7))
ggsave("western_aza_adar.tiff", height = 4, width = 4)

ggplot(aza_western, aes(Sample, pPKR.PKR, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("DMSO", "1_aza", "10_aza"), labels = c("DMSO", expression(paste("1 ", mu, "M")), expression(paste("10 ", mu, "M")))) + labs(x = "8-azaadenosine", y = "Relative pPKR/PKR", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,2), expand = c(0, 0)) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = western_aza_sum, aes(Sample, pPKR_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = western_aza_sum, aes(x = Sample, y = pPKR_mean, ymin = pPKR_mean - pPKR_sd, ymax = pPKR_mean + pPKR_sd), width=0.2, position = position_dodge(width = 0.7))
ggsave("western_aza_pPKR.tiff", height = 4, width = 4)

ggplot(aza_western, aes(Sample, PKR, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("DMSO", "1_aza", "10_aza"), labels = c("DMSO", expression(paste("1 ", mu, "M")), expression(paste("10 ", mu, "M")))) + labs(x = "8-azaadenosine", y = "Relative PKR Abundance", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,2), expand = c(0, 0)) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = western_aza_sum, aes(Sample, PKR_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = western_aza_sum, aes(x = Sample, y = PKR_mean, ymin = PKR_mean - PKR_sd, ymax = PKR_mean + PKR_sd), width=0.2, position = position_dodge(width = 0.7))
ggsave("western_aza_PKR.tiff", height = 4, width = 4)



#read western (immunoblot) data for 8-chloroadenosine
mb468_chloro_western <- read.delim("western_chloro_mb468.txt")
hcc1806_chloro_western <- read.delim("western_chloro_hcc1806.txt")

#bind data together for both cell lines
chloro_western <- rbind(mb468_chloro_western, hcc1806_chloro_western)

#group and summarise data
western_chloro_sum <- dplyr::group_by(chloro_western, Cell_line, Sample)

western_chloro_sum <- dplyr::summarise(western_chloro_sum, ADAR_mean = mean(ADAR), 
                                       ADAR_sd = sd(ADAR), pPKR_mean = mean(pPKR.PKR), 
                                       pPKR_sd = sd(pPKR.PKR),
                                       PKR_mean = mean(PKR), PKR_sd = sd(PKR))
#make plots and save
ggplot(chloro_western, aes(Sample, ADAR, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("DMSO", "1_chloro", "10_chloro"), labels = c("DMSO", expression(paste("1 ", mu, "M")), expression(paste("10 ", mu, "M")))) + labs(x = "8-chloroadenosine", y = "Relative ADAR Abundance", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,2), expand = c(0, 0)) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = western_chloro_sum, aes(Sample, ADAR_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = western_chloro_sum, aes(x = Sample, y = ADAR_mean, ymin = ADAR_mean - ADAR_sd, ymax = ADAR_mean + ADAR_sd), width=0.2, position = position_dodge(width = 0.7))
ggsave("western_chloro_adar.tiff", height = 4, width = 4)

ggplot(chloro_western, aes(Sample, pPKR.PKR, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("DMSO", "1_chloro", "10_chloro"), labels = c("DMSO", expression(paste("1 ", mu, "M")), expression(paste("10 ", mu, "M")))) + labs(x = "8-chloroadenosine", y = "Relative pPKR/PKR", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,2), expand = c(0, 0)) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = western_chloro_sum, aes(Sample, pPKR_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = western_chloro_sum, aes(x = Sample, y = pPKR_mean, ymin = pPKR_mean - pPKR_sd, ymax = pPKR_mean + pPKR_sd), width=0.2, position = position_dodge(width = 0.7))
ggsave("western_chloro_pPKR.tiff", height = 4, width = 4)

t.test(chloro_western$pPKR.PKR[chloro_western$Sample == "DMSO" & chloro_western$Cell_line == "HCC1806"], 
       chloro_western$pPKR.PKR[chloro_western$Sample == "10_chloro" & chloro_western$Cell_line == "HCC1806"])

ggplot(chloro_western, aes(Sample, PKR, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("DMSO", "1_chloro", "10_chloro"), labels = c("DMSO", expression(paste("1 ", mu, "M")), expression(paste("10 ", mu, "M")))) + labs(x = "8-chloroadenosine", y = "Relative PKR Abundance", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,2), expand = c(0, 0)) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = western_chloro_sum, aes(Sample, PKR_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = western_chloro_sum, aes(x = Sample, y = PKR_mean, ymin = PKR_mean - PKR_sd, ymax = PKR_mean + PKR_sd), width=0.2, position = position_dodge(width = 0.7))
ggsave("western_chloro_PKR.tiff", height = 4, width = 4)




#read editing data for knockdowns
editing_knockdown_qsv <- fread("editing_knockdown_qsv.txt")

#group and summarise
editing_knockdown_qsv_sum <- dplyr::group_by(editing_knockdown_qsv, Cell_line, Sample)

editing_knockdown_qsv_sum <- dplyr::summarise(editing_knockdown_qsv_sum, Editing_mean = mean(Percent_edit), Editing_sd = sd(Percent_edit))

#make plots and save
ggplot(editing_knockdown_qsv, aes(Sample, Percent_edit, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shADAR")) + labs(x = "", y = "Percent Editing BPNT1", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), legend.position = "bottom", axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = editing_knockdown_qsv_sum, aes(Sample, Editing_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = editing_knockdown_qsv_sum, aes(x = Sample, y = Editing_mean, ymin = Editing_mean - Editing_sd, ymax = Editing_mean + Editing_sd), width=0.2, position = position_dodge(width = 0.7))
ggsave("editing_knockdown_qsv.tiff", height = 4, width = 3)


#read editing data for knockdowns
editing_knockdown_chr13 <- fread("editing_knockdown_chr13.txt")

#group and summarise
editing_knockdown_chr13_sum <- dplyr::group_by(editing_knockdown_chr13, Cell_line, Sample)

editing_knockdown_chr13_sum <- dplyr::summarise(editing_knockdown_chr13_sum, Editing_mean = mean(Percent_edit), Editing_sd = sd(Percent_edit))

#make plots and save
ggplot(editing_knockdown_chr13, aes(Sample, Percent_edit, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shADAR")) + labs(x = "", y = "Percent Editing ZDHHC20", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), legend.position = "bottom", axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = editing_knockdown_chr13_sum, aes(Sample, Editing_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = editing_knockdown_chr13_sum, aes(x = Sample, y = Editing_mean, ymin = Editing_mean - Editing_sd, ymax = Editing_mean + Editing_sd), width=0.2, position = position_dodge(width = 0.7))
ggsave("editing_knockdown_chr13.tiff", height = 4, width = 3)

t.test(editing_knockdown_chr13$Percent_edit[editing_knockdown_chr13$Sample == "shSCR" & editing_knockdown_chr13$Cell_line == "HCC1806"], 
       editing_knockdown_chr13$Percent_edit[editing_knockdown_chr13$Sample == "shADAR" & editing_knockdown_chr13$Cell_line == "HCC1806"]) 


#read editing data for knockdowns
editing_knockdown_chr10 <- fread("editing_knockdown_chr10.txt")

#group and summarise
editing_knockdown_chr10_sum <- dplyr::group_by(editing_knockdown_chr10, Cell_line, Sample)

editing_knockdown_chr10_sum <- dplyr::summarise(editing_knockdown_chr10_sum, Editing_mean = mean(Percent_edit), Editing_sd = sd(Percent_edit))

#make plots and save
ggplot(editing_knockdown_chr10, aes(Sample, Percent_edit, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shADAR")) + labs(x = "", y = "Percent Editing MRPS16", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), legend.position = "bottom", axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(0,32), expand = c(0, 0)) +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = editing_knockdown_chr10_sum, aes(Sample, Editing_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = editing_knockdown_chr10_sum, aes(x = Sample, y = Editing_mean, ymin = Editing_mean - Editing_sd, ymax = Editing_mean + Editing_sd), width=0.2, position = position_dodge(width = 0.7))
ggsave("editing_knockdown_chr10.tiff", height = 4, width = 3)

t.test(editing_knockdown_chr10$Percent_edit[editing_knockdown_chr10$Sample == "shSCR" & editing_knockdown_chr10$Cell_line == "HCC1806"], 
       editing_knockdown_chr10$Percent_edit[editing_knockdown_chr10$Sample == "shADAR" & editing_knockdown_chr10$Cell_line == "HCC1806"]) 




#read editing data for 8-azaadenosine
editing_aza_qsv <- fread("editing_aza_qsv.txt")

#group and summarise
editing_aza_qsv_sum <- dplyr::group_by(editing_aza_qsv, Cell_line, Sample)

editing_aza_qsv_sum <- dplyr::summarise(editing_aza_qsv_sum, Editing_mean = mean(Percent_edit), Editing_sd = sd(Percent_edit))

#make plots and save
ggplot(editing_aza_qsv, aes(Sample, Percent_edit, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("DMSO", "1A", "10A"), labels = c("DMSO", expression(paste("1 ", mu, "M")), expression(paste("10 ", mu, "M")))) + labs(x = "8-azaadenosine", y = "Percent Editing BPNT1", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = editing_aza_qsv_sum, aes(Sample, Editing_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = editing_aza_qsv_sum, aes(x = Sample, y = Editing_mean, ymin = Editing_mean - Editing_sd, ymax = Editing_mean + Editing_sd), width=0.2, position = position_dodge(width = 0.7))
ggsave("editing_aza_qsv.tiff", height = 4, width = 4)


#read editing data for azas
editing_aza_chr13 <- fread("aza_chr13.txt")

#group and summarise
editing_aza_chr13_sum <- dplyr::group_by(editing_aza_chr13, Cell_line, Sample)

editing_aza_chr13_sum <- dplyr::summarise(editing_aza_chr13_sum, Editing_mean = mean(Percent_edit), Editing_sd = sd(Percent_edit))

#make plots and save
ggplot(editing_aza_chr13, aes(Sample, Percent_edit, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("DMSO", "1A", "10A"), labels = c("DMSO", expression(paste("1 ", mu, "M")), expression(paste("10 ", mu, "M")))) + labs(x = "8-azaadenosine", y = "Percent Editing ZDHHC20", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = editing_aza_chr13_sum, aes(Sample, Editing_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = editing_aza_chr13_sum, aes(x = Sample, y = Editing_mean, ymin = Editing_mean - Editing_sd, ymax = Editing_mean + Editing_sd), width=0.2, position = position_dodge(width = 0.7))
ggsave("editing_aza_chr13.tiff", height = 4, width = 4)

t.test(editing_aza_chr13$Percent_edit[editing_aza_chr13$Sample == "DMSO" & editing_aza_chr13$Cell_line == "HCC1806"], 
       editing_aza_chr13$Percent_edit[editing_aza_chr13$Sample == "1A" & editing_aza_chr13$Cell_line == "HCC1806"]) 
t.test(editing_aza_chr13$Percent_edit[editing_aza_chr13$Sample == "DMSO" & editing_aza_chr13$Cell_line == "HCC1806"], 
       editing_aza_chr13$Percent_edit[editing_aza_chr13$Sample == "10A" & editing_aza_chr13$Cell_line == "HCC1806"]) 


#read editing data for azas
editing_aza_chr10 <- fread("aza_chr10.txt")

#group and summarise
editing_aza_chr10_sum <- dplyr::group_by(editing_aza_chr10, Cell_line, Sample)

editing_aza_chr10_sum <- dplyr::summarise(editing_aza_chr10_sum, Editing_mean = mean(Percent_edit), Editing_sd = sd(Percent_edit))

#make plots and save
ggplot(editing_aza_chr10, aes(Sample, Percent_edit, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("DMSO", "1A", "10A"), labels = c("DMSO", expression(paste("1 ", mu, "M")), expression(paste("10 ", mu, "M")))) + labs(x = "8-azaadenosine", y = "Percent Editing MRPS16", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,40), expand = c(0, 0)) +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = editing_aza_chr10_sum, aes(Sample, Editing_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = editing_aza_chr10_sum, aes(x = Sample, y = Editing_mean, ymin = Editing_mean - Editing_sd, ymax = Editing_mean + Editing_sd), width=0.2, position = position_dodge(width = 0.7))
ggsave("editing_aza_chr10.tiff", height = 4, width = 4)



#read editing data for 8-chloroadenosine
editing_chloro_qsv <- fread("editing_chloro_qsv.txt")

#group and summarise
editing_chloro_qsv_sum <- dplyr::group_by(editing_chloro_qsv, Cell_line, Sample)

editing_chloro_qsv_sum <- dplyr::summarise(editing_chloro_qsv_sum, Editing_mean = mean(Percent_edit), Editing_sd = sd(Percent_edit))

#make plots and save
ggplot(editing_chloro_qsv, aes(Sample, Percent_edit, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("DMSO", "1C", "10C"), labels = c("DMSO", expression(paste("1 ", mu, "M")), expression(paste("10 ", mu, "M")))) + labs(x = "8-chloroadenosine", y = "Percent Editing BPNT1", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = editing_chloro_qsv_sum, aes(Sample, Editing_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = editing_chloro_qsv_sum, aes(x = Sample, y = Editing_mean, ymin = Editing_mean - Editing_sd, ymax = Editing_mean + Editing_sd), width=0.2, position = position_dodge(width = 0.7))
ggsave("editing_chloro_qsv.tiff", height = 4, width = 4)


#read editing data for chloros
editing_chloro_chr13 <- fread("editing_chloro_chr13.txt")

#group and summarise
editing_chloro_chr13_sum <- dplyr::group_by(editing_chloro_chr13, Cell_line, Sample)

editing_chloro_chr13_sum <- dplyr::summarise(editing_chloro_chr13_sum, Editing_mean = mean(Percent_edit), Editing_sd = sd(Percent_edit))

#make plots and save
ggplot(editing_chloro_chr13, aes(Sample, Percent_edit, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("DMSO", "1C", "10C"), labels = c("DMSO", expression(paste("1 ", mu, "M")), expression(paste("10 ", mu, "M")))) + labs(x = "8-chloroadenosine", y = "Percent Editing ZDHHC20", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = editing_chloro_chr13_sum, aes(Sample, Editing_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = editing_chloro_chr13_sum, aes(x = Sample, y = Editing_mean, ymin = Editing_mean - Editing_sd, ymax = Editing_mean + Editing_sd), width=0.2, position = position_dodge(width = 0.7))
ggsave("editing_chloro_chr13.tiff", height = 4, width = 4)

t.test(editing_chloro_chr13$Percent_edit[editing_chloro_chr13$Sample == "shSCR" & editing_chloro_chr13$Cell_line == "HCC1806"], 
       editing_chloro_chr13$Percent_edit[editing_chloro_chr13$Sample == "shADAR" & editing_chloro_chr13$Cell_line == "HCC1806"]) 


#read editing data for chloros
editing_chloro_chr10 <- fread("editing_chloro_chr10.txt")

#group and summarise
editing_chloro_chr10_sum <- dplyr::group_by(editing_chloro_chr10, Cell_line, Sample)

editing_chloro_chr10_sum <- dplyr::summarise(editing_chloro_chr10_sum, Editing_mean = mean(Percent_edit), Editing_sd = sd(Percent_edit))

#make plots and save
ggplot(editing_chloro_chr10, aes(Sample, Percent_edit, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("DMSO", "1C", "10C"), labels = c("DMSO", expression(paste("1 ", mu, "M")), expression(paste("10 ", mu, "M")))) + labs(x = "8-chloroadenosine", y = "Percent Editing MRPS16", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,40), expand = c(0, 0)) +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = editing_chloro_chr10_sum, aes(Sample, Editing_mean, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = editing_chloro_chr10_sum, aes(x = Sample, y = Editing_mean, ymin = Editing_mean - Editing_sd, ymax = Editing_mean + Editing_sd), width=0.2, position = position_dodge(width = 0.7))
ggsave("editing_chloro_chr10.tiff", height = 4, width = 4)






#read qPCR data for knockdowns
qPCR_knockdown <- read.delim("knockdown_IFN_qPCR.txt")

#group and summarise data
qPCR_knockdown_sum <- dplyr::group_by(qPCR_knockdown, Cell_line, Sample, Gene)

qPCR_knockdown_sum <- dplyr::summarise(qPCR_knockdown_sum, mean_FC = mean(Fold_Change), sd_FC = sd(Fold_Change))

qPCR_knockdown_CMPK2 <- subset(qPCR_knockdown, Gene == "CMPK2")

qPCR_knockdown_sum_CMPK2 <- subset(qPCR_knockdown_sum, Gene == "CMPK2")

qPCR_knockdown_CXCL10 <- subset(qPCR_knockdown, Gene == "CXCL10")

qPCR_knockdown_sum_CXCL10 <- subset(qPCR_knockdown_sum, Gene == "CXCL10")


#make plots and save
ggplot(qPCR_knockdown_CMPK2, aes(Sample, Fold_Change, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shADAR")) + labs(x = "", y = "Relative CMPK2 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,4), expand = c(0, 0)) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = qPCR_knockdown_sum_CMPK2, aes(Sample, mean_FC, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = qPCR_knockdown_sum_CMPK2, aes(x = Sample, y = mean_FC, ymin = mean_FC - sd_FC, ymax = mean_FC + sd_FC), width=0.2, position = position_dodge(width = 0.7))
ggsave("qPCR_knockdowns_CMPK2.tiff", height = 4, width = 4)

t.test(qPCR_knockdown_CMPK2$Fold_Change[qPCR_knockdown_CMPK2$Sample == "shSCR" & qPCR_knockdown_CMPK2$Cell_line == "HCC1806"], 
       qPCR_knockdown_CMPK2$Fold_Change[qPCR_knockdown_CMPK2$Sample == "shADAR" & qPCR_knockdown_CMPK2$Cell_line == "HCC1806"])
t.test(qPCR_knockdown_CMPK2$Fold_Change[qPCR_knockdown_CMPK2$Sample == "shSCR" & qPCR_knockdown_CMPK2$Cell_line == "MB468"], 
       qPCR_knockdown_CMPK2$Fold_Change[qPCR_knockdown_CMPK2$Sample == "shADAR" & qPCR_knockdown_CMPK2$Cell_line == "MB468"])


ggplot(qPCR_knockdown_CXCL10, aes(Sample, Fold_Change, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("shSCR", "shADAR")) + labs(x = "", y = "Relative CXCL10 Abundance", colour = "", fill = "") + 
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,6.5), expand = c(0, 0)) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = qPCR_knockdown_sum_CXCL10, aes(Sample, mean_FC, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = qPCR_knockdown_sum_CXCL10, aes(x = Sample, y = mean_FC, ymin = mean_FC - sd_FC, ymax = mean_FC + sd_FC), width=0.2, position = position_dodge(width = 0.7))
ggsave("qPCR_knockdowns_CXCL10.tiff", height = 4, width = 4)

t.test(qPCR_knockdown_CXCL10$Fold_Change[qPCR_knockdown_CXCL10$Sample == "shSCR" & qPCR_knockdown_CXCL10$Cell_line == "HCC1806"], 
       qPCR_knockdown_CXCL10$Fold_Change[qPCR_knockdown_CXCL10$Sample == "shADAR" & qPCR_knockdown_CXCL10$Cell_line == "HCC1806"])
t.test(qPCR_knockdown_CXCL10$Fold_Change[qPCR_knockdown_CXCL10$Sample == "shSCR" & qPCR_knockdown_CXCL10$Cell_line == "MB468"], 
       qPCR_knockdown_CXCL10$Fold_Change[qPCR_knockdown_CXCL10$Sample == "shADAR" & qPCR_knockdown_CXCL10$Cell_line == "MB468"])




#read qPCR data for chloros
qPCR_chloro <- read.delim("chloro_ISG_qPCR.txt")

#group and summarise data
qPCR_chloro_sum <- dplyr::group_by(qPCR_chloro, Cell_line, Sample, Gene)

qPCR_chloro_sum <- dplyr::summarise(qPCR_chloro_sum, mean_FC = mean(Fold_Change), sd_FC = sd(Fold_Change))

qPCR_chloro_CXCL10 <- subset(qPCR_chloro, Gene == "CXCL10")

qPCR_chloro_sum_CXCL10 <- subset(qPCR_chloro_sum, Gene == "CXCL10")

qPCR_chloro_CMPK2 <- subset(qPCR_chloro, Gene == "CMPK2")

qPCR_chloro_sum_CMPK2 <- subset(qPCR_chloro_sum, Gene == "CMPK2")

#make plots and save
ggplot(qPCR_chloro_CXCL10, aes(Sample, Fold_Change, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("DMSO", "1C", "10C"), labels = c("DMSO", expression(paste("1 ", mu, "M")), expression(paste("10 ", mu, "M")))) +
  labs(x = "8-chloroadenosine", y = "Relative CXCL10 Abundance", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,6.5), expand = c(0, 0)) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = qPCR_chloro_sum_CXCL10, aes(Sample, mean_FC, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = qPCR_chloro_sum_CXCL10, aes(x = Sample, y = mean_FC, ymin = mean_FC - sd_FC, ymax = mean_FC + sd_FC), width=0.2, position = position_dodge(width = 0.7))
ggsave("qPCR_chloros_CXCL10.tiff", height = 4, width = 4)

t.test(qPCR_chloro_CXCL10$Fold_Change[qPCR_chloro_CXCL10$Sample == "DMSO" & qPCR_chloro_CXCL10$Cell_line == "HCC1806"], 
       qPCR_chloro_CXCL10$Fold_Change[qPCR_chloro_CXCL10$Sample == "10C" & qPCR_chloro_CXCL10$Cell_line == "HCC1806"])
t.test(qPCR_chloro_CXCL10$Fold_Change[qPCR_chloro_CXCL10$Sample == "DMSO" & qPCR_chloro_CXCL10$Cell_line == "MB468"], 
       qPCR_chloro_CXCL10$Fold_Change[qPCR_chloro_CXCL10$Sample == "10C" & qPCR_chloro_CXCL10$Cell_line == "MB468"])



#make plots and save
ggplot(qPCR_chloro_CMPK2, aes(Sample, Fold_Change, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("DMSO", "1C", "10C"), labels = c("DMSO", expression(paste("1 ", mu, "M")), expression(paste("10 ", mu, "M")))) +
  labs(x = "8-chloroadenosine", y = "Relative CMPK2 Abundance", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,4), expand = c(0, 0)) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = qPCR_chloro_sum_CMPK2, aes(Sample, mean_FC, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = qPCR_chloro_sum_CMPK2, aes(x = Sample, y = mean_FC, ymin = mean_FC - sd_FC, ymax = mean_FC + sd_FC), width=0.2, position = position_dodge(width = 0.7))
ggsave("qPCR_chloros_CMPK2.tiff", height = 4, width = 4)


t.test(qPCR_chloro_CMPK2$Fold_Change[qPCR_chloro_CMPK2$Sample == "DMSO" & qPCR_chloro_CMPK2$Cell_line == "HCC1806"], 
       qPCR_chloro_CMPK2$Fold_Change[qPCR_chloro_CMPK2$Sample == "10C" & qPCR_chloro_CMPK2$Cell_line == "HCC1806"])
t.test(qPCR_chloro_CMPK2$Fold_Change[qPCR_chloro_CMPK2$Sample == "DMSO" & qPCR_chloro_CMPK2$Cell_line == "HCC1806"], 
       qPCR_chloro_CMPK2$Fold_Change[qPCR_chloro_CMPK2$Sample == "1C" & qPCR_chloro_CMPK2$Cell_line == "HCC1806"])
t.test(qPCR_chloro_CMPK2$Fold_Change[qPCR_chloro_CMPK2$Sample == "DMSO" & qPCR_chloro_CMPK2$Cell_line == "MB468"], 
       qPCR_chloro_CMPK2$Fold_Change[qPCR_chloro_CMPK2$Sample == "10C" & qPCR_chloro_CMPK2$Cell_line == "MB468"])


#read qPCR data for azas
qPCR_aza <- read.delim("aza_qpcr.txt")

#group and summarise data
qPCR_aza_sum <- dplyr::group_by(qPCR_aza, Cell_line, Sample, Gene)

qPCR_aza_sum <- dplyr::summarise(qPCR_aza_sum, mean_FC = mean(Fold_Change), sd_FC = sd(Fold_Change))

qPCR_aza_CXCL10 <- subset(qPCR_aza, Gene == "CXCL10")

qPCR_aza_sum_CXCL10 <- subset(qPCR_aza_sum, Gene == "CXCL10")

qPCR_aza_CMPK2 <- subset(qPCR_aza, Gene == "CMPK2")

qPCR_aza_sum_CMPK2 <- subset(qPCR_aza_sum, Gene == "CMPK2")

#make plots and save
ggplot(qPCR_aza_CXCL10, aes(Sample, Fold_Change, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("DMSO", "1A", "10A"), labels = c("DMSO", expression(paste("1 ", mu, "M")), expression(paste("10 ", mu, "M")))) +
  labs(x = "8-azaadenosine", y = "Relative CXCL10 Abundance", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,6.5), expand = c(0, 0)) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = qPCR_aza_sum_CXCL10, aes(Sample, mean_FC, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = qPCR_aza_sum_CXCL10, aes(x = Sample, y = mean_FC, ymin = mean_FC - sd_FC, ymax = mean_FC + sd_FC), width=0.2, position = position_dodge(width = 0.7))
ggsave("qPCR_azas_CXCL10.tiff", height = 4, width = 4)

t.test(qPCR_aza_CXCL10$Fold_Change[qPCR_aza_CXCL10$Sample == "DMSO" & qPCR_aza_CXCL10$Cell_line == "HCC1806"], 
       qPCR_aza_CXCL10$Fold_Change[qPCR_aza_CXCL10$Sample == "10A" & qPCR_aza_CXCL10$Cell_line == "HCC1806"])
t.test(qPCR_aza_CXCL10$Fold_Change[qPCR_aza_CXCL10$Sample == "DMSO" & qPCR_aza_CXCL10$Cell_line == "HCC1806"], 
       qPCR_aza_CXCL10$Fold_Change[qPCR_aza_CXCL10$Sample == "1A" & qPCR_aza_CXCL10$Cell_line == "HCC1806"])
t.test(qPCR_aza_CXCL10$Fold_Change[qPCR_aza_CXCL10$Sample == "DMSO" & qPCR_aza_CXCL10$Cell_line == "MB468"], 
       qPCR_aza_CXCL10$Fold_Change[qPCR_aza_CXCL10$Sample == "10A" & qPCR_aza_CXCL10$Cell_line == "MB468"])



#make plots and save
ggplot(qPCR_aza_CMPK2, aes(Sample, Fold_Change, colour = Cell_line)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05)) + theme_science(base_size = 13) + 
  scale_x_discrete(limits = c("DMSO", "1A", "10A"), labels = c("DMSO", expression(paste("1 ", mu, "M")), expression(paste("10 ", mu, "M")))) +
  labs(x = "8-azaadenosine", y = "Relative CMPK2 Abundance", colour = "", fill = "") + 
  theme(legend.key.size = unit(0.3, "cm"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0,4), expand = c(0, 0)) + geom_hline(yintercept = 1, size = 1, colour = "grey", linetype = "dashed") +
  scale_colour_manual(values = cbviridis, labels= c("HCC1806", "MDA-MB-468")) +
  geom_col(data = qPCR_aza_sum_CMPK2, aes(Sample, mean_FC, colour = Cell_line), width = 0.5, position_dodge(width = 0.7), fill = "white", alpha = 0) +
  geom_errorbar(data = qPCR_aza_sum_CMPK2, aes(x = Sample, y = mean_FC, ymin = mean_FC - sd_FC, ymax = mean_FC + sd_FC), width=0.2, position = position_dodge(width = 0.7))
ggsave("qPCR_azas_CMPK2.tiff", height = 4, width = 4)


t.test(qPCR_aza_CMPK2$Fold_Change[qPCR_aza_CMPK2$Sample == "DMSO" & qPCR_aza_CMPK2$Cell_line == "HCC1806"], 
       qPCR_aza_CMPK2$Fold_Change[qPCR_aza_CMPK2$Sample == "10A" & qPCR_aza_CMPK2$Cell_line == "HCC1806"])
t.test(qPCR_aza_CMPK2$Fold_Change[qPCR_aza_CMPK2$Sample == "DMSO" & qPCR_aza_CMPK2$Cell_line == "MB468"], 
       qPCR_aza_CMPK2$Fold_Change[qPCR_aza_CMPK2$Sample == "10A" & qPCR_aza_CMPK2$Cell_line == "MB468"])
