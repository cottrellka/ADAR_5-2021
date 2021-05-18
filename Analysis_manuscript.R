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

setwd("/Users/Kyle/Box Sync/ADAR_inhibitors/")

#read in viability data for 8-azaadenosine
combined_aza <- read.delim("combined_aza.txt", header = TRUE)

#fit ll.4 model to data
aza <- drm(Viability ~ Concentration, Cell_line, data = combined_aza, fct = LL.4(),
           pmodels = data.frame(Cell_line, Cell_line, 1, Cell_line))
summary(aza)
#find EC50 from model
ED(aza, 50)

#group viability data by cell line and concentration
combined_aza <- dplyr::group_by(combined_aza, Cell_line, Concentration)

#summarise data
combined_aza_sum <- dplyr::summarise(combined_aza, Viability_mean = mean(Viability), SD = sd(Viability))

#plot dose response curves and save
ggplot(combined_aza,aes(Concentration, Viability, colour = Cell_line)) + geom_point() + 
  geom_smooth(method = drm, method.args = list(fct = L.4(fixed = c(NA,NA,1,NA))), se = FALSE) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 50), labels=c(0.001, 0.01, 0.1, 1, 10, 50)) + theme_science() + 
  geom_pointrange(data = combined_aza_sum, aes(x = Concentration, y = Viability_mean, ymin = Viability_mean - SD, ymax = Viability_mean + SD), size = 0.7) +
  scale_colour_manual(values = cbviridis) + 
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab(expression(paste("8-azaadenosine (", mu, "M)"))) + ylab("Relative Viability") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1))

ggplot(combined_aza,aes(Concentration, Viability, colour = Cell_line)) + geom_point() + 
  geom_smooth(method = drm, method.args = list(fct = L.4(fixed = c(NA,NA,1,NA))), se = FALSE) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 50), labels=c(0.001, 0.01, 0.1, 1, 10, 50)) + theme_science() + 
  geom_pointrange(data = combined_aza_sum, aes(x = Concentration, y = Viability_mean, ymin = Viability_mean - SD, ymax = Viability_mean + SD), size = 0.7) +
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
skbr3_combined_aza <- subset(skbr3_combined, Drug == "aza")
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
aza <- drm(Viability ~ Concentration, Cell_line, data = skbr3_combined_aza, fct = LL.4(), 
              pmodels=data.frame(Cell_line, Cell_line, 1, Cell_line, Cell_line))
summary(aza)

#find EC50 from model
ED(aza, 50)

#group and summarise data
skbr3_combined_aza_sum <- dplyr::group_by(skbr3_combined_aza, Cell_line, Concentration)

skbr3_combined_aza_sum <- dplyr::summarise(skbr3_combined_aza_sum, Viability_mean = mean(Viability), SD = sd(Viability))

#plot dose response curves and save
ggplot(skbr3_combined_aza, aes(Concentration, Viability, colour = Cell_line)) + geom_point(alpha = 0.5) + 
  geom_smooth(method = drm, method.args = list(fct = L.4(fixed = c(NA,NA,1,NA))), se = FALSE) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 50), labels=c(0.001, 0.01, 0.1, 1, 10, 50)) + theme_science() +
  geom_pointrange(data = skbr3_combined_aza_sum, aes(x = Concentration, y = Viability_mean, ymin = Viability_mean - SD, ymax = Viability_mean + SD), size = 0.7) +
  scale_colour_manual(values = c("black", "grey"), labels = c("Scr" = "shSCR\nEC50 = 0.74", "ADAR" = "shADAR\nEC50 = 0.71")) + 
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab(expression(paste("SK-BR-3 8-azaadenosine (", mu, "M)"))) + ylab("Relative Viability") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) 
ggsave("skbr3_aza.tiff", height = 4, width = 5.7)


#read in viability data for SKBR3
mcf7_combined <- read.delim("combined_mcf7.txt", header = TRUE)

#subset based on drug
mcf7_combined_aza <- subset(mcf7_combined, Drug == "aza")
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
aza <- drm(Viability ~ Concentration, Cell_line, data = mcf7_combined_aza, fct = LL.4(), 
           pmodels=data.frame(Cell_line, Cell_line, 1, Cell_line, Cell_line))
summary(aza)

#find EC50 from model
ED(aza, 50)

#group and summarise the data
mcf7_combined_aza_sum <- dplyr::group_by(mcf7_combined_aza, Cell_line, Concentration)

mcf7_combined_aza_sum <- dplyr::summarise(mcf7_combined_aza_sum, Viability_mean = mean(Viability), SD = sd(Viability))

#plot dose response curves and save
ggplot(mcf7_combined_aza, aes(Concentration, Viability, colour = Cell_line)) + geom_point(alpha = 0.5) + 
  geom_smooth(method = drm, method.args = list(fct = L.4(fixed = c(NA,NA,1,NA))), se = FALSE) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 50), labels=c(0.001, 0.01, 0.1, 1, 10, 50)) + theme_science() +
  geom_pointrange(data = mcf7_combined_aza_sum, aes(x = Concentration, y = Viability_mean, ymin = Viability_mean - SD, ymax = Viability_mean + SD), size = 0.7) +
  scale_colour_manual(values = c("black", "grey"), labels = c("scr" = "shSCR\nEC50 = 1.81", "ADAR" = "shADAR\nEC50 = 2.04")) + 
  theme(legend.position = "bottom", legend.title = element_blank()) + xlab(expression(paste("MCF-7 8-azaadenosine (", mu, "M)"))) + ylab("Relative Viability") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) 
ggsave("mcf7_aza.tiff", height = 4, width = 5.7)


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

#group and summarise data
western_knockdown_sum <- dplyr::group_by(western_knockdown, Cell_line, Sample)

western_knockdown_sum <- dplyr::summarise(western_knockdown_sum, ADAR_mean = mean(ADAR), ADAR_sd = sd(ADAR), pPKR_mean = mean(pPKR.PKR), pPKR_sd = sd(pPKR.PKR))

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


#read western (immunoblot) data for 8-azaadenosine
mb468_aza_western <- read.delim("western_aza_mb468.txt")
hcc1806_aza_western <- read.delim("western_aza_hcc1806.txt")

#bind together data from both cell lines
aza_western <- rbind(mb468_aza_western, hcc1806_aza_western)

#group and summarise data
western_aza_sum <- dplyr::group_by(aza_western, Cell_line, Sample)

western_aza_sum <- dplyr::summarise(western_aza_sum, ADAR_mean = mean(ADAR), ADAR_sd = sd(ADAR), pPKR_mean = mean(pPKR.PKR), pPKR_sd = sd(pPKR.PKR))

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


#read western (immunoblot) data for 8-chloroadenosine
mb468_chloro_western <- read.delim("western_chloro_mb468.txt")
hcc1806_chloro_western <- read.delim("western_chloro_hcc1806.txt")

#bind data together for both cell lines
chloro_western <- rbind(mb468_chloro_western, hcc1806_chloro_western)

#group and summarise data
western_chloro_sum <- dplyr::group_by(chloro_western, Cell_line, Sample)

western_chloro_sum <- dplyr::summarise(western_chloro_sum, ADAR_mean = mean(ADAR), ADAR_sd = sd(ADAR), pPKR_mean = mean(pPKR.PKR), pPKR_sd = sd(pPKR.PKR))

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

