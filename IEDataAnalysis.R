#==============================Libraries
library(dplyr)
library(ggplot2)
library(foreach)
library(doParallel)
library(deSolve)
library(reticulate)
library(FME)
library(ABCoptim)
library(caRamel)
library(locfit)
library(nsga2R)
library(CEoptim)
#=============================ggplot theme
theme_min <- readRDS("/mnt/580CBE2464C5F83D/pracovni/helpfull_R_Matlab_script/ggtheme.rds")
#=============================DATA
IE <- read.csv("DataIncubation.csv")
names(IE)
str(IE)
summary(IE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Plots of variables over time~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#==================Chemical parameters
#Glucose concentration
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Time < 20) %>% 
  summarise(y = mean(Gl, na.rm = T), ySD = sd(Gl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Glucose (", mu, "mol C g ", DW^{-1}, ")"))) 
#pH
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(pH, na.rm = T), ySD = sd(pH, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("pH")))

#Respiration rate
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(R, na.rm = T), ySD = sd(R, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Respiration rate (", mu, "mol C g ", DW^{-1}~d^{-1}, ")")))
#Oxygen consumption rate
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(O2, na.rm = T), ySD = sd(O2, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(O[2], " cosumption rate (", mu, "mol ", O[2], " g ", DW^{-1}~d^{-1}, ")")))
#Respiration quotient
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & RQ <= 1) %>% 
  summarise(y = mean(RQ, na.rm = T), ySD = sd(RQ, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(Respiration~quotient)))
#Cumulative respiration
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(CumulativeR, na.rm = T), ySD = sd(CumulativeR, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Cumulative respiration (", mu, "mol C g ", DW^{-1}, ")")))
#DOCw
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(DOCw, na.rm = T), ySD = sd(DOCw, na.rm = T),
            y2 = mean(Gl, na.rm = T), y2SD = sd(Gl, na.rm = T),
            y3 = mean(Cpotass, na.rm = T), y3SD = sd(Cpotass, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  #geom_point(cex = 6, pch = 21, fill = "red", aes(x = Time, y = y2), alpha = 0.5) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD)) + 
  #geom_point(cex = 6, pch = 21, fill = "green", aes(x = Time, y = y3), alpha = 0.5) + geom_errorbar(aes(ymin = y3 - y3SD, ymax = y3 + y3SD)) + 
  xlab("Time (days)") + ylab(expression(paste("DOC in water (", mu, "mol C g ", DW^{-1}, ")")))
#DONw
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(DONw, na.rm = T), ySD = sd(DONw, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("DON in water (", mu, "mol N g ", DW^{-1}, ")")))
#NH4w
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(NH4w, na.rm = T), ySD = sd(NH4w, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(NH[4], " in water (", mu, "mol N g ", DW^{-1}, ")")))
#NO3w
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(NO3w, na.rm = T), ySD = sd(NO3w, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(NO[3], " in water (", mu, "mol N g ", DW^{-1}, ")")))
#DOPw
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(DOPw, na.rm = T), ySD = sd(DOPw, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(DOP, " in water (", mu, "mol P g ", DW^{-1}, ")")))
#SRPw
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(PO4w, na.rm = T), ySD = sd(PO4w, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(SRP, " in water (", mu, "mol P g ", DW^{-1}, ")")))
#Al
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(Al, na.rm = T), ySD = sd(Al, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(Al, " in water (", mu, "mol Al g ", DW^{-1}, ")")))
#Fe
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(Fe, na.rm = T), ySD = sd(Fe, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(Fe, " in water (", mu, "mol Fe g ", DW^{-1}, ")")))
#S
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(S, na.rm = T), ySD = sd(S, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(S, " in water (", mu, "mol S g ", DW^{-1}, ")")))

#==================Biochemical parameters
#Beta-glucosidase Vmax
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(VmaxBG, na.rm = T), ySD = sd(VmaxBG, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(V[MAX], " of BG (", mu, "mol MUB-C g ", DW^{-1}~hour^{-1}, ")")))
#Beta-glucosidase Km
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(KmBG, na.rm = T), ySD = sd(KmBG, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(K[M], " of BG (", mu, "mol MUB-C g ", DW^{-1}, ")")))
#Beta-glucosidase Ki
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(KiBG, na.rm = T), ySD = sd(KiBG, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(K[i], " of BG (", mu, "mol MUB-C g ", DW^{-1}, ")")))
#Cellobiosidase Vmax
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(VmaxCL, na.rm = T), ySD = sd(VmaxCL, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(V[MAX], " of CL (", mu, "mol MUB-C g ", DW^{-1}~hour^{-1}, ")")))
#Cellobiosidase Km
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(KmCL, na.rm = T), ySD = sd(KmCL, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(K[M], " of CL (", mu, "mol MUB-C g ", DW^{-1}, ")")))
#Cellobiosidase Ki
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(KiCL, na.rm = T), ySD = sd(KiCL, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(K[i], " of CL (", mu, "mol MUB-C g ", DW^{-1}, ")")))
#Chitinase Vmax
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(VmaxCh, na.rm = T), ySD = sd(VmaxCh, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(V[MAX], " of Ch (", mu, "mol MUB-C g ", DW^{-1}~hour^{-1}, ")")))
#Chitinase Km
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(KmCh, na.rm = T), ySD = sd(KmCh, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(K[M], " of Ch (", mu, "mol MUB-C g ", DW^{-1}, ")")))
#Chitinase Ki
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(KiCh, na.rm = T), ySD = sd(KiCh, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(K[i], " of Ch (", mu, "mol MUB-C g ", DW^{-1}, ")")))
#Leucine-aminopeptidase Vmax
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(VmaxL, na.rm = T), ySD = sd(VmaxL, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(V[MAX], " of L (", mu, "mol MUB-C g ", DW^{-1}~hour^{-1}, ")")))
#Leucine-aminopeptidase Km
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(KmL, na.rm = T), ySD = sd(KmL, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(K[M], " of L (", mu, "mol MUB-C g ", DW^{-1}, ")")))
#Phosphatase Vmax
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(VmaxP, na.rm = T), ySD = sd(VmaxP, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(V[MAX], " of P (", mu, "mol MUB-C g ", DW^{-1}~hour^{-1}, ")")))
#Phosphatase Km
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(KmP, na.rm = T), ySD = sd(KmP, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(K[M], " of P (", mu, "mol MUB-C g ", DW^{-1}, ")")))
#Phosphatase Ki
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(KiP, na.rm = T), ySD = sd(KiP, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(K[i], " of P (", mu, "mol MUB-C g ", DW^{-1}, ")")))

#==================Biological parameters
#Cflush
IE$CflushOut <- c("False")
IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic" & IE$Time == 2), "CflushOut"] <- c("True")

IE %>% group_by(Time, Treatment, Soil) %>% #filter(CflushOut == "False") %>% 
  summarise(y = mean(Cflush, na.rm = T), ySD = sd(Cflush, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(CHCl[3], " flush (", mu, "mol C g ", DW^{-1}, ")")))
#Nflush
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(Nflush, na.rm = T), ySD = sd(Nflush, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(CHCl[3], " flush (", mu, "mol N g ", DW^{-1}, ")")))
#Pflush
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(Pflush, na.rm = T), ySD = sd(Pflush, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(CHCl[3], " flush (", mu, "mol P g ", DW^{-1}, ")")))
#C to P in flush
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(Cflush/Pflush, na.rm = T), ySD = sd(Cflush/Pflush, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(frac(C, P)," in ", CHCl[3], " flush (mol/mol)")))
#C to N in flush
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(Cflush/Nflush, na.rm = T), ySD = sd(Cflush/Nflush, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(frac(C, N)," in ", CHCl[3], " flush (mol/mol)")))
#N to P in flush
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(Cflush/Nflush, na.rm = T), ySD = sd(Cflush/Nflush, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(frac(N, P)," in ", CHCl[3], " flush (mol/mol)")))
#DNA
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(DNA, na.rm = T), ySD = sd(DNA, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("DNA (", mu, "mol C-DNA g ", DW^{-1}, ")")))
#RNA
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(RNA, na.rm = T), ySD = sd(RNA, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free_y") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("RNA (", mu, "mol C-RNA g ", DW^{-1}, ")")))
#ATP
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(ATP, na.rm = T), ySD = sd(ATP, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free_y") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("ATP (", mu, "mol C-ATP g ", DW^{-1}, ")")))
#PLFA
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(PLFA, na.rm = T), ySD = sd(PLFA, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free_y") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("PLFA (", mu, "mol C-PLFA g ", DW^{-1}, ")")))
#CLC to DNA
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(Cflush/DNA, na.rm = T), ySD = sd(Cflush/DNA, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(" ")))
#CLC to ATP
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(ATP/Cflush, na.rm = T), ySD = sd(Cflush/DNA, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(" ")))
#====================================Functional groups
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(C20_3w6/PLFA*100, na.rm = T), ySD = sd(C20_3w6/PLFA*100, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(" ")))
#==================================We will only define bacteria and fungi
#Bacteria
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(DNAb, na.rm = T), ySD = sd(DNAb, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Bacterial genes (", 10^{9}, " copies ", mu, "mol C-DNA)")))
#Fungi
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(DNAf, na.rm = T), ySD = sd(DNAf, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Fungal genes (", 10^{9}, " copies ", mu, "mol C-DNA)")))

#Fungi to bacteria ratio
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(DNAf/DNAb, na.rm = T), ySD = sd(DNAf/DNAb, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Fungi to bacteria ratio (unitless)")))
#=============================PLFA
#certainly fungal PLFA is 18_2w6
IE$PLFAf <- IE$C18_2w6
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(PLFAf, na.rm = T), ySD = sd(PLFAf, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Fungal PLFA (", mu, "mol(C)", g(DW)^{-1}, ")")))

#Bacterial PLFA are selected according to Frostegard and Baath (1996) so the derived conversion factors can be used
##those are C15, iC15, aC15, iC16, iC17, aC17, C16_1w7, C18_1w7, and cyc19 (C17, cyc17 and C16_1w9 are missing)
IE$PLFAb <- with(IE, C15+iC15+aC15+iC17+aC17+C16_1w7+C18_1w7+cycC19)
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(PLFAb, na.rm = T), ySD = sd(PLFAb, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Bacterial PLFA (", mu, "mol(C)", g(DW)^{-1}, ")")))

#Fungi to bacteria ratio - PLFA based
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(PLFAf/PLFAb, na.rm = T), ySD = sd(PLFAf/PLFAb, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Fungi to bacteria ratio (unitless)")))

#Comparing DNA a PLFA based fungi to bacteria ratio
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(PLFAf/PLFAb, na.rm = T), ySD = sd(PLFAf/PLFAb, na.rm = T),
            y2 = mean(DNAf/DNAb, na.rm = T), ySD2 = sd(PLFAf/PLFAb, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  geom_point(cex = 6, pch = 21, col = "red", aes(Time, y2)) + 
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  geom_errorbar(aes(ymin = y2 - ySD2, ymax = y2 + ySD2), col = "red") + 
  xlab("Time (days)") + ylab(expression(paste("Fungi to bacteria ratio (unitless)")))

#======================Isotopes
#delta 13C of CO2
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(delta13CO2, na.rm = T), ySD = sd(delta13CO2, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste({}^{13},CO[2])))

#Respiration from glucose
##Fill missing data in PL and CT aerobic data
###=========================Plesne
which(is.na(IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic"), "delta13CO2"]))
IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic"), "delta13CO2"][45] <- mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic"), "delta13CO2"][46:48])
IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic"), "delta13CO2"][56] <- mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic"), "delta13CO2"][53:55])
IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic"), "delta13CO2"][67:68] <- mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic"), "delta13CO2"][65:66])
###==============================
###========================Certovo
which(is.na(IE[(IE$Soil == "Certovo" & IE$Treatment == "Aerobic"), "delta13CO2"]))
IE[(IE$Soil == "Certovo" & IE$Treatment == "Aerobic"), "delta13CO2"][54] <- mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Aerobic"), "delta13CO2"][c(53, 55, 56)])
###==============================
#Delta values to atm%
IE$CO2atm <- with(IE, ((delta13CO2/1000 + 1)*0.0111802)/(1 + ((delta13CO2/1000 + 1)*0.0111802))) 
#CO2 from glucose in umol C-CO2/g DW/h
IE$CO2gl <- numeric(length = nrow(IE))

for(i in 1:nrow(IE)){
  if(IE$Soil[i] == "Plesne" & IE$Treatment[i] == "Aerobic"){
    IE$CO2gl[i] <- (IE$CO2atm[i] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic" & IE$Time == 0), "CO2atm"]))/
      (IE$at13Gl[i] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic" & IE$Time == 0), "CO2atm"]))*IE$R[i]
  }else{
    if(IE$Soil[i] == "Plesne" & IE$Treatment[i] == "Anaerobic"){
      IE$CO2gl[i] <- (IE$CO2atm[i] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Anaerobic" & IE$Time == 0), "CO2atm"]))/
        (IE$at13Gl[i] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Anaerobic" & IE$Time == 0), "CO2atm"]))*IE$R[i]
    }else{
      if(IE$Soil[i] == "Certovo" & IE$Treatment[i] == "Anaerobic"){
        IE$CO2gl[i] <- (IE$CO2atm[i] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Anaerobic" & IE$Time == 0), "CO2atm"]))/
          (IE$at13Gl[i] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Anaerobic" & IE$Time == 0), "CO2atm"]))*IE$R[i]
      }else{
        IE$CO2gl[i] <- (IE$CO2atm[i] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Aerobic" & IE$Time == 0), "CO2atm"]))/
          (IE$at13Gl[i] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Aerobic" & IE$Time == 0), "CO2atm"]))*IE$R[i]
      }
    }
  }
}
IE[IE$Time == 0, "CO2gl"] <- 0
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(CO2gl, na.rm = T), ySD = sd(CO2gl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste({}^{13},CO[2], " (", mu, "mol ", g^{-1}~h^{-1}, ")"))) + scale_y_log10()
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(R - CO2gl, na.rm = T), ySD = sd(R - CO2gl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste({}^{13},CO[2]))) 

#Cumulative respiration from glucose and soil
##=============================Total cumulative respiration recalculated for better accuracy
CumulativeR <- function(Soil, Treatment, variable){
  if(Treatment == "Aerobic"){
    CO2t <- rep(0, 4)
    for(i in 2:length(unique(IE[(IE$Soil == Soil & IE$Treatment == Treatment), "Time"]))){
      xdt <- IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time == unique(IE[(IE$Soil == Soil & IE$Treatment == Treatment), "Time"])[i]), "Time"] -
        IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time == unique(IE[(IE$Soil == Soil & IE$Treatment == Treatment), "Time"])[i-1]), "Time"]
      rt0 <- IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time == unique(IE[(IE$Soil == Soil & IE$Treatment == Treatment), "Time"])[i-1]), variable]
      drdt <- (IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time == unique(IE[(IE$Soil == Soil & IE$Treatment == Treatment), "Time"])[i]), variable] -
                 IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time == unique(IE[(IE$Soil == Soil & IE$Treatment == Treatment), "Time"])[i-1]), variable])
      CO2t <- append(CO2t, xdt*24*(rt0 + drdt/2))
    }
    xRC <- numeric(0)
    for(i in 1:length(seq(1, 68, by = 4))){
      xRC <- append(xRC, c(cumsum(CO2t[seq(1, 68, by = 4)])[i],
                           cumsum(CO2t[seq(2, 68, by = 4)])[i],
                           cumsum(CO2t[seq(3, 68, by = 4)])[i],
                           cumsum(CO2t[seq(4, 68, by = 4)])[i]))
    }
  }else{
    CO2t <- rep(0, 4)
    for(i in 2:length(unique(IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time != 0.17 & IE$Time != 0.19), "Time"]))){
      xdt <- IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time == unique(IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time != 0.17 & IE$Time != 0.19), "Time"])[i]), "Time"] -
        IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time == unique(IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time != 0.17 & IE$Time != 0.19), "Time"])[i-1]), "Time"]
      rt0 <- IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time == unique(IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time != 0.17 & IE$Time != 0.19), "Time"])[i-1]), variable]
      drdt <- (IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time == unique(IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time != 0.17 & IE$Time != 0.19), "Time"])[i]), variable] -
                 IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time == unique(IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time != 0.17 & IE$Time != 0.19), "Time"])[i-1]), variable])
      CO2t <- append(CO2t, xdt*24*(rt0 + drdt/2))
    }
    
    xRC <- numeric(0)
    for(i in 1:length(seq(1, 64, by = 4))){
      xRC <- append(xRC, c(cumsum(CO2t[seq(1, 64, by = 4)])[i],
                           cumsum(CO2t[seq(2, 64, by = 4)])[i],
                           cumsum(CO2t[seq(3, 64, by = 4)])[i],
                           cumsum(CO2t[seq(4, 64, by = 4)])[i]))
    }
    
    xRC <- c(xRC[1:4], NA, NA, NA, NA, xRC[5:64])
  }
  return(xRC)
}

IE$CumulativeR <- numeric(length = nrow(IE))

for(i in unique(IE$Soil)){
  for(n in unique(IE$Treatment)){
    IE[(IE$Soil == i & IE$Treatment == n), "CumulativeR"] <- CumulativeR(Soil = i, Treatment = n, variable = "R")
  }
}

##========================================================================
##=============================Cumulative respiration from glucose
IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic"), c("CO2gl", "Time")]
IE$CO2s <- IE$R - IE$CO2gl
CumulativeRg <- function(Soil, Treatment, variable){
  CO2t <- rep(0, 4)
  for(i in 2:length(unique(IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time != 0.17 & IE$Time != 0.19 & IE$Time != 4), "Time"]))){
    xdt <- IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time == unique(IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time != 0.17 & IE$Time != 0.19 & IE$Time != 4), "Time"])[i]), "Time"] -
      IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time == unique(IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time != 0.17 & IE$Time != 0.19 & IE$Time != 4), "Time"])[i-1]), "Time"]
    rt0 <- IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time == unique(IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time != 0.17 & IE$Time != 0.19 & IE$Time != 4), "Time"])[i-1]), variable]
    drdt <- (IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time == unique(IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time != 0.17 & IE$Time != 0.19 & IE$Time != 4), "Time"])[i]), variable] -
               IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time == unique(IE[(IE$Soil == Soil & IE$Treatment == Treatment & IE$Time != 0.17 & IE$Time != 0.19 & IE$Time != 4), "Time"])[i-1]), variable])
    CO2t <- append(CO2t, xdt*24*(rt0 + drdt/2))
  }
  
  xRC <- numeric(0)
  for(i in 1:length(seq(1, 60, by = 4))){
    xRC <- append(xRC, c(cumsum(CO2t[seq(1, 64, by = 4)])[i],
                         cumsum(CO2t[seq(2, 64, by = 4)])[i],
                         cumsum(CO2t[seq(3, 64, by = 4)])[i],
                         cumsum(CO2t[seq(4, 64, by = 4)])[i]))
  }
  
  xRC <- c(xRC[1:4], NA, NA, NA, NA, xRC[5:16], NA, NA, NA, NA, xRC[17:60])
  return(xRC)
}
#Cumulative respiration from glucose in umol C/g DW
IE$CumulativeRg <- numeric(length = nrow(IE))
#Cumulative respiration from soil in umol C/g DW
IE$CumulativeRs <- numeric(length = nrow(IE))
for(i in unique(IE$Soil)){
  for(n in unique(IE$Treatment)){
    IE[(IE$Soil == i & IE$Treatment == n), "CumulativeRg"] <- CumulativeRg(Soil = i, Treatment = n, variable = "CO2gl")
    IE[(IE$Soil == i & IE$Treatment == n), "CumulativeRs"] <- CumulativeRg(Soil = i, Treatment = n, variable = "CO2s")
  }
}

IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(CumulativeRg, na.rm = T), ySD = sd(CumulativeRg, na.rm = T),
            y2 = mean(CumulativeRs, na.rm = T), y2SD = sd(CumulativeRs, na.rm = T),
            y3 = mean(CumulativeR, na.rm = T), y3SD = sd(CumulativeR, na.rm = T),
            y4 = mean(CumulativeRg + CumulativeRs, na.rm = T), y4SD = sd(CumulativeRg + CumulativeRs, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  geom_point(cex = 6, pch = 21, fill = "red", aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "red") + 
  xlab("Time (days)") + ylab(expression(paste("Cumulative ", CO[2]^{"Glucose"}, " production (", mu, "mol C ", g~(DW)^{-1}, ")"))) 

#Residual glucose in soil
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean((Gl0 - CumulativeRg)/Gl0*100, na.rm = T), ySD = sd((Gl0 - CumulativeRg)/Gl0*100, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab("Residual glucose (%)") 
#CUE from glucose and cumulative respiration
IE$CUE_RgGl <- 1 - IE$CumulativeRg/(500 - IE$Cpotassgl)
IE$CUE_RgGl <- ifelse(IE$CUE_RgGl < 0, NA, IE$CUE_RgGl)
IE$CUE_RgGl <- ifelse(IE$CUE_RgGl > 1, NA, IE$CUE_RgGl)
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(CUE_RgGl, na.rm = T), ySD = sd(CUE_RgGl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab("CUE") 

#delta 13C of K2SO4 extract
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(CpotassAtm, na.rm = T), ySD = sd(CpotassAtm, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste({}^{13},EC-K[2],SO[4], " (at%)")))

#Glucose C in K2SO4 extract in umol C/g DW
IE$Cpotassgl <- numeric(length = nrow(IE))

for(i in 1:nrow(IE)){
  if(IE$Soil[i] == "Plesne" & IE$Treatment[i] == "Aerobic"){
    IE$Cpotassgl[i] <- (IE$CpotassAtm[i] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic" & IE$Time == 0), "CpotassAtm"]))/
      (IE$at13Gl[i] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic" & IE$Time == 0), "CpotassAtm"]))*IE$Cpotass[i]
  }else{
    if(IE$Soil[i] == "Plesne" & IE$Treatment[i] == "Anaerobic"){
      IE$Cpotassgl[i] <- (IE$CpotassAtm[i] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Anaerobic" & IE$Time == 0), "CpotassAtm"]))/
        (IE$at13Gl[i] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Anaerobic" & IE$Time == 0), "CpotassAtm"]))*IE$Cpotass[i]
    }else{
      if(IE$Soil[i] == "Certovo" & IE$Treatment[i] == "Anaerobic"){
        IE$Cpotassgl[i] <- (IE$CpotassAtm[i] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Anaerobic" & IE$Time == 0), "CpotassAtm"]))/
          (IE$at13Gl[i] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Anaerobic" & IE$Time == 0), "CpotassAtm"]))*IE$Cpotass[i]
      }else{
        IE$Cpotassgl[i] <- (IE$CpotassAtm[i] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Aerobic" & IE$Time == 0), "CpotassAtm"]))/
          (IE$at13Gl[i] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Aerobic" & IE$Time == 0), "CpotassAtm"]))*IE$Cpotass[i]
      }
    }
  }
}
IE[IE$Time == 0, "Cpotassgl"] <- 0

IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(Cpotassgl, na.rm = T), ySD = sd(Cpotassgl, na.rm = T),
            y2 = mean(Cpotass, na.rm = T), y2SD = sd(Cpotass, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  geom_point(cex = 6, pch = 21, fill = "red", aes(x = Time, y = y2), alpha = 0.5) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), col = "red") +
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in ",EC-K[2],SO[4], " (", mu, "mol C ",g~(DW)^{-1}, ")" )))

#Glucose concentration vs glucose content of K2SO4 extract
ggplot(IE, aes(Cpotassgl, Gl)) + geom_point(cex = 6, pch = 21) + theme_min +
  xlab(expression(paste("Glucose-C in ",EC-K[2],SO[4], " (", mu, "mol C ",g~(DW)^{-1}, ")" ))) +
  ylab(expression(paste("Glucose-C (", mu, "mol C ",g~(DW)^{-1}, ")" ))) + geom_abline(intercept = 0, slope = 1) +
  facet_grid(Treatment~Soil)

#delta 13C of water extract
IE$DOCwAtm <- (IE$delta13DOCw/1000 + 1)*0.0111802
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(CpotassAtm, na.rm = T), ySD = sd(DOCwAtm, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste({}^{13},C, " in  DOC (at%)")))

#Glucose C in water extract in umol C/g DW
IE$DOCwgl <- numeric(length = nrow(IE))

for(i in 1:nrow(IE)){
  if(IE$Soil[i] == "Plesne" & IE$Treatment[i] == "Aerobic"){
    IE$DOCwgl[i] <- (IE$DOCwAtm[i] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic" & IE$Time == 0), "DOCwAtm"]))/
      (IE$at13Gl[i] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic" & IE$Time == 0), "DOCwAtm"]))*IE$DOCw[i]
  }else{
    if(IE$Soil[i] == "Plesne" & IE$Treatment[i] == "Anaerobic"){
      IE$DOCwgl[i] <- (IE$DOCwAtm[i] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Anaerobic" & IE$Time == 0), "DOCwAtm"]))/
        (IE$at13Gl[i] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Anaerobic" & IE$Time == 0), "DOCwAtm"]))*IE$DOCw[i]
    }else{
      if(IE$Soil[i] == "Certovo" & IE$Treatment[i] == "Anaerobic"){
        IE$DOCwgl[i] <- (IE$DOCwAtm[i] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Anaerobic" & IE$Time == 0), "DOCwAtm"]))/
          (IE$at13Gl[i] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Anaerobic" & IE$Time == 0), "DOCwAtm"]))*IE$DOCw[i]
      }else{
        IE$DOCwgl[i] <- (IE$DOCwAtm[i] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Aerobic" & IE$Time == 0), "DOCwAtm"]))/
          (IE$at13Gl[i] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Aerobic" & IE$Time == 0), "DOCwAtm"]))*IE$DOCw[i]
      }
    }
  }
}
IE[IE$Time == 0, "DOCwgl"] <- 0

IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(DOCwgl, na.rm = T), ySD = sd(DOCwgl, na.rm = T),
            y2 = mean(DOCw, na.rm = T), y2SD = sd(DOCw, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  geom_point(cex = 6, pch = 21, fill = "red", aes(x = Time, y = y2), alpha = 0.5) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), col = "red") +
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in DOC (", mu, "mol C ",g~(DW)^{-1}, ")" )))

#Glucose concentration vs glucose content of K2SO4 extract
ggplot(IE, aes(Cpotassgl, Gl)) + geom_point(cex = 6, pch = 21) + theme_min +
  xlab(expression(paste("Glucose-C in ",EC-K[2],SO[4], " (", mu, "mol C ",g~(DW)^{-1}, ")" ))) +
  ylab(expression(paste("Glucose-C (", mu, "mol C ",g~(DW)^{-1}, ")" ))) + geom_abline(intercept = 0, slope = 1) +
  facet_grid(Treatment~Soil)

#Glucose concentration vs glucose content of DOC extract
ggplot(IE, aes(DOCwgl, Glraw)) + geom_point(cex = 6, pch = 21) + theme_min +
  xlab(expression(paste("Glucose-C in DOC (", mu, "mol C ",g~(DW)^{-1}, ")" ))) +
  ylab(expression(paste("Glucose-C (", mu, "mol C ",g~(DW)^{-1}, ")" ))) + geom_abline(intercept = 0, slope = 1) +
  facet_grid(Treatment~Soil)

#Glucose C in Cflush extract in umol C/g DW
IE$CFlushgl <- numeric(length = nrow(IE))
for(i in 1:nrow(IE)){
  if(IE$Soil[i] == "Plesne" & IE$Treatment[i] == "Aerobic"){
    IE$CFlushgl[i] <- (IE$CflushAtm[i] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic" & IE$Time == 0), "CflushAtm"]))/
      (IE$at13Gl[i] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic" & IE$Time == 0), "CflushAtm"]))*IE$Cflush[i]
  }else{
    if(IE$Soil[i] == "Plesne" & IE$Treatment[i] == "Anaerobic"){
      IE$CFlushgl[i] <- (IE$CflushAtm[i] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Anaerobic" & IE$Time == 0), "CflushAtm"]))/
        (IE$at13Gl[i] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Anaerobic" & IE$Time == 0), "CflushAtm"]))*IE$Cflush[i]
    }else{
      if(IE$Soil[i] == "Certovo" & IE$Treatment[i] == "Anaerobic"){
        IE$CFlushgl[i] <- (IE$CflushAtm[i] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Anaerobic" & IE$Time == 0), "CflushAtm"]))/
          (IE$at13Gl[i] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Anaerobic" & IE$Time == 0), "CflushAtm"]))*IE$Cflush[i]
      }else{
        IE$CFlushgl[i] <- (IE$CflushAtm[i] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Aerobic" & IE$Time == 0), "CflushAtm"]))/
          (IE$at13Gl[i] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Aerobic" & IE$Time == 0), "CflushAtm"]))*IE$Cflush[i]
      }
    }
  }
}

IE[IE$Time == 0, "CFlushgl"] <- 0
IE$CflushOut[21] <- "True"
IE$CflushOut[17] <- "True"

IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(CflushAtm, na.rm = T), ySD = sd(CflushAtm, na.rm = T),
            y2 = mean(Gl, na.rm = T), y2SD = sd(Gl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in ",CHCl[3]~flush, " (", mu, "mol C ",g~(DW)^{-1}, ")" )))

IE %>% group_by(Time, Treatment, Soil) %>% #filter(CflushOut == "False") %>% 
  summarise(y = mean(CFlushgl, na.rm = T), ySD = sd(CFlushgl, na.rm = T),
            y2 = mean(Cflush, na.rm = T), y2SD = sd(Cflush, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  geom_point(cex = 6, pch = 21, fill = "green", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "green") +
  #geom_point(cex = 6, pch = 21, fill = "blue", alpha = 0.5, aes(x = Time, y = y3)) + geom_errorbar(aes(ymin = y3 - y3SD, ymax = y3 + y3SD), color = "blue") +
  #geom_point(cex = 6, pch = 21, fill = "red", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "red") +
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in ",CHCl[3]~flush, " (", mu, "mol C ",g~(DW)^{-1}, ")" )))

#Comparing glucose in Cflush with residual glucose C (everything that was not respired)
IE %>% group_by(Time, Treatment, Soil) %>% #filter(CflushOut == "False") %>% 
  summarise(y = mean(CFlushgl, na.rm = T), ySD = sd(CFlushgl, na.rm = T),
            y2 = mean((Gl0 - Gl - CumulativeRg), na.rm = T), y2SD = sd((Gl0 - Gl - CumulativeRg), na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  #geom_point(cex = 6, pch = 21, fill = "red", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "red") +
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in ",CHCl[3]~flush, " (", mu, "mol C ",g~(DW)^{-1}, ")" )))

#Glucose C in PLFAs in umol C/g DW
NC <- names(IE)[59:82]
for(i in 1:length(NC)){NC[i] <- paste0(NC[i], "gl")}
NCframe <- as.data.frame(matrix(nrow = nrow(IE), ncol = length(NC)))
colnames(NCframe) <- NC
IE <- cbind(IE, NCframe)

NC2 <- names(IE)[59:82]
NC3 <- names(IE)[84:107]

NCall <- data.frame(NC, NC2, NC3)

for(i in 1:nrow(IE)){
  for(l in 1:nrow(NCall)){
    if(IE$Soil[i] == "Plesne" & IE$Treatment[i] == "Aerobic"){
      IE[i, NCall[l, 1]] <- (IE[i, NCall[l, 3]] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic" & IE$Time == 0), NCall[l, 3]]))/
        (IE$at13Gl[i] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic" & IE$Time == 0), NCall[l, 3]]))*IE[i, NCall[l, 2]]
    }else{
      if(IE$Soil[i] == "Plesne" & IE$Treatment[i] == "Anaerobic"){
        IE[i, NCall[l, 1]] <- (IE[i, NCall[l, 3]] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Anaerobic" & IE$Time == 0), NCall[l, 3]]))/
          (IE$at13Gl[i] - mean(IE[(IE$Soil == "Plesne" & IE$Treatment == "Anaerobic" & IE$Time == 0), NCall[l, 3]]))*IE[i, NCall[l, 2]]
      }else{
        if(IE$Soil[i] == "Certovo" & IE$Treatment[i] == "Anaerobic"){
          IE[i, NCall[l, 1]] <- (IE[i, NCall[l, 3]] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Anaerobic" & IE$Time == 0), NCall[l, 3]]))/
            (IE$at13Gl[i] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Anaerobic" & IE$Time == 0), NCall[l, 3]]))*IE[i, NCall[l, 2]]
        }else{
          IE[i, NCall[l, 1]] <- (IE[i, NCall[l, 3]] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Aerobic" & IE$Time == 0), NCall[l, 3]]))/
            (IE$at13Gl[i] - mean(IE[(IE$Soil == "Certovo" & IE$Treatment == "Aerobic" & IE$Time == 0), NCall[l, 3]]))*IE[i, NCall[l, 2]]
        }
      }
    }
  }
}

IE[IE$Time == 0, NC] <- 0
for(i in NC){
  IE[(!is.na(IE[, i]) & IE[, i] < 0), i] <- 0 
}

#Total PLFA
IE$PLFAgl <- numeric(length = nrow(IE))
for(i in 1:nrow(IE)){
  IE$PLFAgl[i] <- sum(IE[i, c(121:144)], na.rm = T)
}

IE[(IE$PLFAgl == 0 & IE$Time > 0), "PLFAgl"] <- NA
summary(IE$PLFAgl)
#===================Total PLFA
IE %>% group_by(Time, Treatment, Soil) %>% #filter(CflushOut == "False") %>% 
  summarise(y2 = mean(CFlushgl, na.rm = T), y2SD = sd(CFlushgl, na.rm = T),
            y = mean(PLFAgl, na.rm = T), ySD = sd(PLFAgl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  #geom_point(cex = 6, pch = 21, fill = "green", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "green") +
  #geom_point(cex = 6, pch = 21, fill = "blue", alpha = 0.5, aes(x = Time, y = y3)) + geom_errorbar(aes(ymin = y3 - y3SD, ymax = y3 + y3SD), color = "blue") +
  #geom_point(cex = 6, pch = 21, fill = "red", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "red") +
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in PLFA (", mu, "mol C ",g~(DW)^{-1}, ")" )))

#Fungi
IE$Fungl <- IE$C18_2w6gl
#Bacteria
IE$Bacgl <- with(IE, C15gl+iC15gl+aC15gl+iC17gl+aC17gl+C16_1w7gl+C18_1w7gl+cycC19gl)
#======================Fungi
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(Fungl, na.rm = T), ySD = sd(Fungl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in fungal PLFA (", mu, "mol C ",g~(DW)^{-1}, ")" )))
#removing outliers
IE$FunglAll <- IE$Fungl
IE[(IE$FunglAll > 0.45 & IE$Treatment == "Aerobic" & IE$Soil == "Certovo" & !is.na(IE$FunglAll)), "Fungl"] <- NA
#======================Bacteria
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(Bacgl, na.rm = T), ySD = sd(Bacgl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in bacterial PLFA (", mu, "mol C ",g~(DW)^{-1}, ")" )))
#removing outliers
IE$BacglAll <- IE$Bacgl
IE[(IE$BacglAll > 0.18 & IE$Treatment == "Aerobic" & IE$Soil == "Certovo" & !is.na(IE$BacglAll)), "Bacgl"] <- NA
IE[(IE$BacglAll == max(IE[(IE$Treatment == "Aerobic" & IE$Soil == "Plesne" & IE$Time == 14), "Bacgl"]) & !is.na(IE$BacglAll)), "Bacgl"] <- NA
#======================Relative abundance of fungal and bacterial PLFA
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic") %>% 
  summarise(y = mean((Fungl+Bacgl)/PLFAgl, na.rm = T), ySD = sd((Fungl+Bacgl)/PLFAgl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in PLFA (", mu, "mol C ",g~(DW)^{-1}, ")" )))
#======================Actinobacteria
IE$Actgl <- with(IE, C18_1w9gl)
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(Actgl/PLFAgl, na.rm = T), ySD = sd(PLFAgl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in actinobacterial PLFA (", mu, "mol C ",g~(DW)^{-1}, ")" )))

#================pH sensitivity
pHbac <- read.csv2("pHbacteria.csv", header = F, sep = ";", col.names = c("pH", "growth"))
pHbac$group <- "Bacteria"
pHbac$growth <- pHbac$growth*24/1000/6

summary(lm(growth~pH, pHbac))
coef(lm(growth~pH, pHbac))[2]/1.6

phfun<- read.csv2("pHfungi.csv", header = F, sep = ";", col.names = c("pH", "growth"))
phfun$group <- "Fungi"
phfun$growth <- phfun$growth*24/1000/2

summary(lm(growth~pH, phfun))
coef(lm(growth~pH, phfun))[2]/1.6

phall <- rbind(pHbac, phfun)

ggplot(phall, aes(pH, growth)) + geom_point(cex = 6, pch = 21, aes(fill = group)) +
  theme_min + stat_smooth(method = lm, aes(color = group))

#Anderson and Domsch 1986, Anderson and Grey 1990 - substrate inhibition
AD <- data.frame(U = c(23.92/968, 30.63/968, 33.16/968, 32.53/968, 38.73/968, 40.51/968, 28.48/968),
                 S0 = c(90.8/12.01, 181.7/12.01, 254.3/12.01, 363.3/12.01, 726.6/12.01, 1453.3/12.01, 2179.9/12.01),
                 Legend = rep("Phaozem", 7))

AG <- read.csv("AndersonGrey1990.csv")

AD <- rbind(AD, AG)

ggplot(AD[AD$U < 0.15, ], aes(S0, U)) + geom_point(cex = 6, pch = 21, aes(fill = Legend)) + theme_min #+
  stat_function(fun = function(x){coef(SI_AD)[1]*x/(coef(SI_AD)[2] + x + x^2/coef(SI_AD)[3])})

SI_AD1 <- nls(U ~ Vmax*S0/(Km + S0 + S0^2/Ki), data = AD, start = list(Vmax = 0.1, Km = 50, Ki = 70),
              subset = Legend == "Phaozem")
summary(SI_AD1)
SI_AD2 <- nls(U ~ Vmax*S0/(Km + S0 + S0^2/Ki), data = AD, start = list(Vmax = 0.1, Km = 50, Ki = 70),
              subset = Legend == "Monoculture")
summary(SI_AD2)
SI_AD3 <- nls(U ~ Vmax*S0/(Km + S0 + S0^2/Ki), data = AD, start = list(Vmax = 0.1, Km = 50, Ki = 70),
              subset = Legend == "Rotation")
summary(SI_AD3)

SI_AD <- nls(U ~ Vmax*S0/(Km + S0 + S0^2/Ki), data = AD, start = list(Vmax = 0.1, Km = 50, Ki = 70))
summary(SI_AD)
#==============================================SubMicrobial model==============================================
#================================================Basic model==================================================#
#===================Basic model and respective solver
source_python("Python/subM0.py")
source_python("Python/subMSolver0.py")
#===================================================
#Objective function and initial parameter guess
source("Python/subMObjective0.R")
p0 <- c(6.5, 221, 0.95, 3, 2, 4e-4, 0.41, 0.24)#, 0.05
pl <- c(0.01, 0.1, 0.1, 0.5, 1e-7, 1e-7, 0.3, 0.04)#, 0.003
pu <- c(20, 400, 1, 25, 5, 1e-2, 1, 0.6)#, 0.1
#==================================================Certovo - aerobic conditions
#==================================================#
#~~~~~~~~~~~~~Parameters estimation~~~~~~~~~~~~~~~~#
#==================================================#
##First guess by MCMC 
CAguess <- modMCMC(subMObjective0, p = p0, lower = pl, upper = pu, niter = 30000, 
                   IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"))
summary(CAguess)

#Estimate
CAOpt <- abc_optim(fn = subMObjective0,
                   par = as.numeric(apply(CAguess$pars, 2, quantile, 0.5)), 
                   lb = as.numeric(apply(CAguess$pars, 2, quantile, 0.05)), 
                   ub = as.numeric(apply(CAguess$pars, 2, quantile, 0.95)), 
                   IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"))
round(CAOpt$par, 6)
#=============================================Model fit
source("Python/subMFit0.R")
CAout <- subMFit0(CAOpt$par, "Certovo", "Aerobic", 
                  IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"))
CAout$errors
CAout$R2all

#Glucose concentration
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Certovo") %>% 
  summarise(y = mean(Gl, na.rm = T), ySD = sd(Gl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Glucose (", mu, "mol C g ", DW^{-1}, ")"))) +
  geom_line(data = CAout$Simulation, aes(x = Time, y = Gl))
#Cumulative respiration from glucose
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Certovo") %>% 
  summarise(y = mean(CumulativeRg, na.rm = T), ySD = sd(CumulativeRg, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Cumulative ", CO[2]^{"Glucose"}, " production (", mu, "mol C ", g~(DW)^{-1}, ")"))) +
  geom_line(data = CAout$Simulation, aes(x = Time, y = CO2))
#Chloroform flush
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Certovo" & CFlushgl < 100) %>% 
  summarise(y = mean(CFlushgl, na.rm = T), ySD = sd(CFlushgl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in ",CHCl[3]~flush, " (", mu, "mol C ",g~(DW)^{-1}, ")" ))) +
  geom_line(data = CAout$Simulation, aes(x = Time, y = Cflush))

#===============================Fermentation products
##===================Certovo 
#1. Using CLC, CO2, glucose and pH (i.e. function of acetic acid production)
#pH change due to CO2 dissolution
IE$ExtraH <- with(IE, 10^-6.4187*(0.04554*96*pCO2/1e6*9.87e-3)/(10^-pH + 10^-6.4187))
for(i in c("Plesne", "Certovo")){
  for(n in c("Aerobic", "Anaerobic")){
    IE[(IE$Soil == i & IE$Treatment == n), "ExtraH"] <- IE[(IE$Soil == i & IE$Treatment == n), "ExtraH"] - 
      mean(as.numeric(IE[(IE$Soil == i & IE$Treatment == n & IE$Time == 0), "ExtraH"]), na.rm = T)
  }
}

IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Certovo") %>% 
  summarise(y = mean(pH, na.rm = T), ySD = sd(pH, na.rm = T),
            ExtraH = mean(ExtraH, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab("pH") +
  geom_point(aes(x = Time, y = -log10(10^-4.3+ExtraH)), color = "green")

#Spline model to fill out extra H+ from CO2 dissolution into simulation
lc1 <- locfit(ExtraH ~ lp(Time, nn = 0.1), data = subset(IE, Soil == "Certovo" & Treatment == "Aerobic"))
summary(lc1)
lc2 <- locfit(pH ~ lp(Time, nn = 0.1), data = subset(IE, Soil == "Certovo" & Treatment == "Aerobic"))
summary(lc2)
#===================Basic model and respective solver
source_python("Python/subMFerm.py")
source_python("Python/subMSolverFerm.py")
#===================================================
#Objective function and initial parameter guess
source("Python/subMObjectiveFerm.R")
##Optimization
###Model parameters (Initial guess, lower and upper bound) 
Im = c(1, 1e-2, 20)
Km = c(25, 0.1, 500)
yA = c(0.9, 0, 1)
Gm = c(1, 1e-3, 1e2)
m = c(1e-3, 1e-8, 1)
edemand = c(0.9, 0.1, 10)
emax = c(0.19, 0.01, 1)
etaf = c(2.55, 0.1, 20)
etar = c(1.33, 0.1, 20)
pb = c(0.71, 0, 1)
pg = c(0.61, 0, 1)
ng = c(0.8, 0, 1)
nb = c(0.3, 0, 1)
k = m

ParmsFerm = rbind(Im, Km, yA, Gm, m, emax, ng, nb)

##First guess by MCMC 
CTAGuess0 <- modMCMC(subMObjectiveFerm, p = ParmsFerm[,1], lower = ParmsFerm[, 2], upper = ParmsFerm[, 3], niter = 30000,
                     IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"))
summary(CTAGuess0)
##Estimate
CTAP0 <- abc_optim(fn = subMObjectiveFerm, IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"),
                   par = as.numeric(summary(CTAGuess0)[c("mean"), ]), 
                   lb = as.numeric(summary(CTAGuess0)[c("min"), ]), 
                   ub = as.numeric(summary(CTAGuess0)[c("max"), ]))
round(CTAP0$par, 6)
# #Uncertainty
# FermPU <- modMCMC(FermObjective, p = FermP$par, 
#                   lower = ParmsFerm[1:8, 2], upper = ParmsFerm[1:8, 3], niter = 5000)
# summary(FermPU)
#Goodness of fit and simulations
source("Python/subMFitFerm.R")
SimCTA0 <- subMFitFerm(CTAP0$par, "Certovo", "Aerobic", IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"))
SimCTA0$errors
SimCTA0$R2all

#Glucose concentration
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Certovo") %>% 
  summarise(y = mean(Gl, na.rm = T), ySD = sd(Gl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Glucose (", mu, "mol C g ", DW^{-1}, ")"))) +
  geom_line(data = SimCTA0$Simulation, aes(x = Time, y = Gl)) +
  geom_line(data = CAout$Simulation, aes(x = Time, y = Gl), color = "red")
#Cumulative respiration from glucose
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Certovo") %>% 
  summarise(y = mean(CumulativeRg, na.rm = T), ySD = sd(CumulativeRg, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Cumulative ", CO[2]^{"Glucose"}, " production (", mu, "mol C ", g~(DW)^{-1}, ")"))) +
  geom_line(data = SimCTA0$Simulation, aes(x = Time, y = CO2))+
  geom_line(data = CAout$Simulation, aes(x = Time, y = CO2), color = "red")
#Chloroform flush
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Certovo") %>% 
  summarise(y = mean(CFlushgl, na.rm = T), ySD = sd(CFlushgl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in ",CHCl[3]~flush, " (", mu, "mol C ",g~(DW)^{-1}, ")" ))) +
  geom_line(data = SimCTA0$Simulation, aes(x = Time, y = Cflush)) +
  geom_line(data = CAout$Simulation, aes(x = Time, y = Cflush), color = "red")
#pH
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Certovo") %>% 
  summarise(y = mean(pH, na.rm = T), ySD = sd(pH, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab("pH") +
  geom_line(data = SimCTA0$Simulation, aes(x = Time, y = pH))
#=======================fermentation products and PLFA data
#=======================Adding PLFA data according to microbial groups
#===================Basic model and respective solver
source_python("Python/subMFermGroups.py")
#source_python("Python/subMSolverFermGroups.py")
#===================================================
#Objective function and initial parameter guess
source("Python/subMObjectiveFermGroups.R")
##Optimization
###Model parameters (Initial guess, lower and upper bound) 
Im = c(6, 0.1, 15)
Km = c(221, 1, 400)
Ki = c(100, 10, 1000)
yA = c(0.9, 0.24, 1)
Gm = c(3, 0.5, 10)
m = c(1e-3, 1e-8, 1)
edemand = c(1.6, 0.8, 5)
emax = c(0.19, 0.01, 1)
etaf = c(2.55, 0.1, 20)
etar = c(1.33, 0.1, 20)
pb = c(0.71, 0, 1)
pg = c(0.61, 0, 1)
ng = c(0.4, 1e-3, 1)
nb = c(0.2, 1e-3, 1)
k = c(1e-3, 1e-8, 1)
kpf = c(0.05, 0, 1)
kpb = kpf
kp = c(0.14, 1e-6, 1)
s = c(0.04, 0, 1)
phm = c(0.03, 0, 0.2)
eta = c(0.5, 0.1, 3)
a = c(0.5, 0.001, 0.999)

ParmsFermGroups = rbind(Im, yA, k, eta, a)
#=============================Testing physiological differences between bacteria and fungi
source_python("Python/subMFermGroups.py")
#source_python("Python/subMSolverFermGroups.py")
source("Python/subMObjectiveFermGroups.R")
source("Python/subMFitFermGroups.R")

#Null model - all parameters fixed for bacteria and fungi
##First guess by MCMC 
nullGuess <- modMCMC(subMObjectiveFermGroups, free = NA, IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"),
                     p = ParmsFermGroups[,1], 
                     lower = ParmsFermGroups[, 2], 
                     upper = ParmsFermGroups[, 3], niter = 30000)
summary(nullGuess)
#Estimate
null_p <- abc_optim(fn = subMObjectiveFermGroups, free = NA, IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"),
                    par = ParmsFermGroups[,1], 
                    lb = ParmsFermGroups[,2], 
                    ub = ParmsFermGroups[,3],
                    maxCycle = 3000, FoodNumber = 50, criter = 100)
round(null_p$par, 6)
null_model <- subMFitFermGroups(null_p$par, free = NA, IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"),
                                 Soil = 'Certovo', Treatment = 'Aerobic')
null_model$errors
null_model$R2all

#subMObjectiveFermGroups(null_p$par, free = NA, IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"))

#Alternative models with microbial group-specific parameters
vpars0 <- as.matrix(expand.grid(c("Im", "yA", "k", "eta"), 
                                c("Im", "yA", "k", "eta"),
                                c("Im", "yA", "k", "eta"),
                                c("Im", "yA", "k", "eta")))
#==================Running in parallel
registerDoParallel(150)

vpars <- foreach(i=1:nrow(vpars0), .combine=rbind, .multicombine = TRUE) %dopar% {
  
  matrix(c(setdiff(c("Im", "yA", "k", "eta"), 
                   setdiff(c("Im", "yA", "k", "eta"), unique(vpars0[i, ]))), 
           rep(NA, times = (4-length(setdiff(c("Im", "yA", "k", "eta"), 
                                             setdiff(c("Im", "yA", "k", "eta"), unique(vpars0[i, ]))))))),
         nrow = 1, ncol = 4)
  
}

stopImplicitCluster()
#=================Running sequentially (takes ages)
#Keep the unique rows only
vpars <- matrix(nrow = 0, ncol = 4)
for(i in 1:nrow(vpars0)){
  vpars <- rbind(vpars, matrix(c(setdiff(c("Im", "yA", "k", "eta"), 
                                         setdiff(c("Im", "yA", "k", "eta"), unique(vpars0[i, ]))), 
                                 rep(NA, times = (4-length(setdiff(c("Im", "yA", "k", "eta"), 
                                                                   setdiff(c("Im", "yA", "k", "eta"), unique(vpars0[i, ]))))))),
                               nrow = 1, ncol = 4))
}
#====================================
vpars <- vpars[!duplicated(vpars), ]
rm(vpars0)
#====================================
source("Python/GroupsParallel.R")

registerDoParallel(nrow(vpars))

res <- foreach(i=1:nrow(vpars), .combine=rbind, .multicombine = TRUE,
               .packages=c("FME", "dplyr", "deSolve", "reticulate", "ABCoptim")) %dopar% {
                 
                 GroupsParallel(i)
                 
               }

stopImplicitCluster()

summary(res)

resD <- as.data.frame(rbind(c(rep("NULL", 4), null_model$errors, null_model$R2all), res))

for(i in 5:17){
  resD[, i] <- as.numeric(resD[, i])
}

resD$pr <- pf(q=(resD[1, "Fnorm"] - resD[, "Fnorm"])*(resD[, "n"] - resD[, "p"])/resD[, "Fnorm"]/(resD[, "p"] - resD[1, "p"]), 
             df1=(resD[, "p"] - resD[1, "p"]), 
             df2=(resD[, "n"] - resD[, "p"]), 
             lower.tail=F)
resD <- resD[order(resD$pr), ]
#==============================fitting the best model - Im, and k differs between bacteria and fungi
ParmsFermFinal = rbind(Im, Im, k, k, yA, eta, a)
#Best model - emax and Im differs between bacteria and fungi
##First guess by MCMC 
# CAFinalGuess <- modMCMC(subMObjectiveFermGroups, free = c("Im", "yA", "k", "eta"), 
#                         IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"),
#                         p = as.numeric(ParmsFermFinal[,1]), 
#                         lower = as.numeric(ParmsFermFinal[, 2]), 
#                         upper = as.numeric(ParmsFermFinal[, 3]), niter = 30000, updatecov = 100, burninlength = 500)
# summary(CAFinalGuess)
#Estimate
# CAFinal_p <- abc_optim(fn = subMObjectiveFermGroups, free = c("Im", "yA", "k", "eta"),
#                     par = as.numeric(ParmsFermFinal[,1]), #ParmsFermFinal[,1], 
#                     lb = as.numeric(ParmsFermFinal[,2]), #ParmsFermFinal[,2],
#                     ub = as.numeric(ParmsFermFinal[,3]), #ParmsFermFinal[,3], 
#                     IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"),
#                     maxCycle = 3000, FoodNumber = 50, criter = 100)
# names(CAFinal_p$par) <- c("ImF", "ImB", "yAF", "yAB", "kF", "kB", "etaF", "etaB", "a") #"nb", "ng", "a"
# round(CAFinal_p$par, 6)
# SimCTAFinal <- subMFitFermGroups(CAFinal_p$par, c("Im", "yA", "k", "eta"),
#                                  IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"),
#                                  Soil = "Certovo", Treatment = "Aerobic")
# SimCTAFinal$errors
# SimCTAFinal$R2all
#Making 100 parameter estimates to get the best parameters and calculate error of estimate
##Run on server
n_iter = 100
source("Python/finalParallel.R")
registerDoParallel(n_iter)

CTfinalPars <- foreach(i=1:n_iter, .combine=rbind, .multicombine = TRUE,
               .packages=c("FME", "dplyr", "deSolve", "reticulate", "ABCoptim")) %dopar% {
                 
                 finalParallel(i)
                 
               }

stopImplicitCluster()

CTfinalPars <- read.csv("CTfinalPars.csv", row.names = NULL)[, -1]
summary(CTfinalPars)
#adding goodness of fit
CTGF <- data.frame(R2 = numeric(length = 100), R2adj = numeric(length = 100), ll = numeric(length = 100),
                   AIC = numeric(length = 100), Form = numeric(length = 100), Gl = numeric(length = 100),
                   CO2 = numeric(length = 100), Fungi = numeric(length = 100), Bacteria = numeric(length = 100),
                   Cflush = numeric(length = 100), pH = numeric(length = 100)) 

for(i in 1:100){
  SimOut <- subMFitFermGroups(as.numeric(CTfinalPars[i, ]), c("Im", "k"),
                                   IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"),
                                   Soil = "Certovo", Treatment = "Aerobic")
  CTGF[i, ] <- c(SimOut$errors[1:5], SimOut$R2all)
}
CTbestP <- which(CTGF$ll == max(CTGF$ll))
CTfinalPars[CTbestP, ]
CTGF[CTbestP, ]

CTfinalBest <- subMFitFermGroups(as.numeric(CTfinalPars[CTbestP, ]), c("Im", "k"),
                                 IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"),
                                 Soil = "Certovo", Treatment = "Aerobic")
#Glucose concentration
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Certovo") %>% 
  summarise(y = mean(Gl, na.rm = T), ySD = sd(Gl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Glucose (", mu, "mol C g ", DW^{-1}, ")"))) +
  geom_line(data = CTfinalBest$Simulation, aes(x = Time, y = Gl))# +
  #geom_line(data = null_model$Simulation, aes(x = Time, y = Gl), color = "red")
#Cumulative respiration from glucose
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Certovo") %>% 
  summarise(y = mean(CumulativeRg, na.rm = T), ySD = sd(CumulativeRg, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Cumulative ", CO[2], " production (", mu, "mol C ", g~(DW)^{-1}, ")"))) +
  geom_line(data = CTfinalBest$Simulation, aes(x = Time, y = CO2))#+
  #geom_line(data = null_model$Simulation, aes(x = Time, y = CO2), color = "red")
#Chloroform flush
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Certovo") %>% 
  summarise(y = mean(CFlushgl, na.rm = T), ySD = sd(CFlushgl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in ",CHCl[3]~flush, " (", mu, "mol C ",g~(DW)^{-1}, ")" ))) +
  geom_line(data = CTfinalBest$Simulation, aes(x = Time, y = Cflush))# +
  #geom_line(data = null_model$Simulation, aes(x = Time, y = Cflush), color = "red")
#pH
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Certovo") %>% 
  summarise(y = mean(pH, na.rm = T), ySD = sd(pH, na.rm = T),
            ExtraH = mean(ExtraH, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab("pH") +
  geom_line(data = CTfinalBest$Simulation, aes(x = Time, y = pH))# +
  #geom_line(data = null_model$Simulation, aes(x = Time, y = pH), color = "red")+
  #geom_point(aes(x = Time, y = -log10(10^-4.3+ExtraH)), color = "green")


#PLFA
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Certovo") %>% 
  summarise(y = mean(Fungl, na.rm = T), ySD = sd(Fungl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab("Fungi") +
  geom_line(data = CTfinalBest$Simulation, aes(x = Time, y = Fungi))# +
  #geom_line(data = null_model$Simulation, aes(x = Time, y = Fungi), color = "red")

IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Certovo") %>% 
  summarise(y = mean(Bacgl, na.rm = T), ySD = sd(Bacgl, na.rm = T),
            y2 = mean(Bacgl, na.rm = T), ySD = sd(Bacgl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab("Bacteria") +
  geom_line(data = CTfinalBest$Simulation, aes(x = Time, y = Bacteria))# +
  #geom_line(data = null_model$Simulation, aes(x = Time, y = Bacteria), color = "red")

#Calculating residuals for all variables
CTmeans <- as.matrix(IE %>% filter(Treatment == "Aerobic" & Soil == "Certovo") %>% 
  group_by(Time) %>% summarize(Glucose = mean(Gl, na.rm = T), 
                               CO2 = mean(CumulativeRg, na.rm = T),
                               Fungi = mean(Fungl, na.rm = T),
                               Bacteria = mean(Bacgl, na.rm = T),
                               Cflush = mean(CFlushgl, na.rm = T),
                               pH = mean(pH, na.rm = T)))
CTresiduals <- ((CTfinalBest$Yhat - CTmeans[, -1])/CTfinalBest$W)^2      
CTresiduals <- as.data.frame(cbind(CTmeans[, 1], CTresiduals))
colnames(CTresiduals)[1] <- c("Time")  
write.csv(CTresiduals, "CTresiduals.csv", row.names = F)

#Comparing DNA a PLFA based fungi to bacteria ratio again applying estimated conversion factors
InitialMBCCT = mean(as.numeric(IE[(IE$Time == 0 & IE$Soil == "Certovo" & IE$Treatment == "Aerobic"), "Cflush"]), na.rm = T)/0.24
InitialFCT = InitialMBCCT*tail(as.numeric(CTfinalPars[43, ]), 1)
InitialBCT = InitialMBCCT*(1 - tail(as.numeric(CTfinalPars[43, ]), 1))
kfCT <- mean(as.numeric(IE[(IE$Time == 0 & IE$Soil == "Certovo" & IE$Treatment == "Aerobic"), "PLFAf"]), na.rm = T)/InitialFCT
kbCT <- mean(as.numeric(IE[(IE$Time == 0 & IE$Soil == "Certovo" & IE$Treatment == "Aerobic"), "PLFAb"]), na.rm = T)/InitialBCT

IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(PLFAf/kfCT/PLFAb*kbCT, na.rm = T), ySD = sd(PLFAf/kfCT/PLFAb*kbCT, na.rm = T),
            y2 = mean(DNAf/DNAb, na.rm = T), ySD2 = sd(PLFAf/PLFAb, na.rm = T)) %>% 
  ggplot(aes(y2, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  #geom_point(cex = 6, pch = 21, col = "red", aes(Time, y2)) + 
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + #geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  #geom_errorbar(aes(ymin = y2 - ySD2, ymax = y2 + ySD2), col = "red") + 
  xlab("Time (days)") + ylab(expression(paste("Fungi to bacteria ratio (unitless)"))) +
  stat_smooth(method = lm)

#==================================================For AI
#Matrix of parameters for 100 different simulations (Xp) - 100 rows (100 simulations) and 7 columns - 7 parameters to optimize
Xp <- matrix(c(sample(seq(from = Im[2], to = Im[3], length.out = 999), size = 100), #ImF
        sample(seq(from = Im[2], to = Im[3], length.out = 999), size = 100), #ImF
        sample(seq(from = k[2], to = k[3], length.out = 999), size = 100), #kF
        sample(seq(from = k[2], to = k[3], length.out = 999), size = 100), #kB
        sample(seq(from = yA[2], to = yA[3], length.out = 999), size = 100), #yA
        sample(seq(from = eta[2], to = eta[3], length.out = 999), size = 100), #eta
        sample(seq(from = a[2], to = a[3], length.out = 999), size = 100) #a
        ), ncol = 7, nrow = 100)
Xp[1, ] <- as.numeric(CTfinalPars[43, ])

#Matrix of simulations Xs - 100 rows (100 simulations) and 17*6 columns - 17 times and 6 measured variables
Xs <- matrix(nrow = 100, ncol = 17*6)
for(i in 1:100){
  Xs[i, ] <- as.numeric(subMFitFermGroups(as.numeric(Xp[i, ]), c("Im", "k"),
                                          IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"),
                                          Soil = "Certovo", Treatment = "Aerobic")$Yhat)
}

#Matrix of data (Xd) - 1 row and 17*6 columns - 17 times and 6 measured variables
Xd <- matrix(nrow = 1, ncol = 17*6)
Xd[1, ] <- as.numeric(as.matrix(IE %>% filter(Treatment == "Aerobic" & Soil == "Certovo") %>% 
                       group_by(Time) %>% summarize(Glucose = mean(Gl, na.rm = T), 
                                                    CO2 = mean(CumulativeRg, na.rm = T),
                                                    Fungi = mean(Fungl, na.rm = T),
                                                    Bacteria = mean(Bacgl, na.rm = T),
                                                    Cflush = mean(CFlushgl, na.rm = T),
                                                    pH = mean(pH, na.rm = T)))[, -1])
write.csv(Xp, "Xp.csv", row.names = FALSE)
write.csv(Xs, "Xs.csv", row.names = FALSE)
write.csv(Xd, "Xd.csv", row.names = FALSE)
#==================================================Plesne - aerobic conditions
#==================================================#
#~~~~~~~~~~~~~Parameters estimation~~~~~~~~~~~~~~~~#
#==================================================#
#=======================fermentation products and PLFA data
#===================Basic model and respective solver
source_python("Python/subMFermGroups.py")
#source_python("Python/subMSolverFermGroups.py")
#===================================================
#Objective function and initial parameter guess
source("Python/subMObjectiveFermGroups.R")
##Optimization
###Model parameters (Initial guess, lower and upper bound) 
Im = c(6, 0.1, 15)
Km = c(221, 1, 400)
Ki = c(100, 10, 1000)
yA = c(0.9, 0.24, 1)
Gm = c(3, 0.5, 10)
m = c(1e-3, 1e-8, 1)
edemand = c(1.6, 0.8, 5)
emax = c(0.19, 0.01, 1)
etaf = c(2.55, 0.1, 20)
etar = c(1.33, 0.1, 20)
pb = c(0.71, 0, 1)
pg = c(0.61, 0, 1)
ng = c(0.4, 1e-3, 1)
nb = c(0.2, 1e-3, 1)
k = c(1e-3, 1e-8, 1)
kpf = c(0.05, 0, 1)
kpb = kpf
kp = c(0.14, 1e-6, 1)
s = c(0.04, 0, 1)
phm = c(0.03, 0, 0.2)
eta = c(0.5, 0.1, 3)
a = c(0.5, 0.001, 0.999)

ParmsFermGroups = rbind(Im, yA, k, eta, a)
#=============================Testing physiological differences between bacteria and fungi
source_python("Python/subMFermGroups.py")
#source_python("Python/subMSolverFermGroups.py")
source("Python/subMObjectiveFermGroups.R")
source("Python/subMFitFermGroups.R")

#Null model - all parameters fixed for bacteria and fungi
# ##First guess by MCMC 
# nullGuess <- modMCMC(subMObjectiveFermGroups, free = NA, IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"),
#                      p = ParmsFermGroups[,1], 
#                      lower = ParmsFermGroups[, 2], 
#                      upper = ParmsFermGroups[, 3], niter = 30000)
# summary(nullGuess)
#Estimate
null_pPL <- abc_optim(fn = subMObjectiveFermGroups, free = NA, IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Plesne"),
                    par = ParmsFermGroups[,1], 
                    lb = ParmsFermGroups[,2], 
                    ub = ParmsFermGroups[,3],
                    maxCycle = 3000, FoodNumber = 50, criter = 100)
round(null_pPL$par, 6)
null_modelPL <- subMFitFermGroups(null_pPL$par, free = NA, IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Plesne"),
                                Soil = 'Plesne', Treatment = 'Aerobic')
null_modelPL$errors
null_modelPL$R2all

#================================Variable parameters
source("Python/GroupsParallel.R")

registerDoParallel(nrow(vpars))

resPL <- foreach(i=1:nrow(vpars), .combine=rbind, .multicombine = TRUE,
               .packages=c("FME", "dplyr", "deSolve", "reticulate", "ABCoptim")) %dopar% {
                 
                 GroupsParallel(i)
                 
               }

stopImplicitCluster()

summary(resPL)

resDPL <- as.data.frame(rbind(c(rep("NULL", 4), null_modelPL$errors, null_modelPL$R2all), resPL))

for(i in 5:17){
  resDPL[, i] <- as.numeric(resDPL[, i])
}

resDPL$pr <- pf(q=(resDPL[1, "Fnorm"] - resDPL[, "Fnorm"])*(resDPL[, "n"] - resDPL[, "p"])/resDPL[, "Fnorm"]/(resDPL[, "p"] - resDPL[1, "p"]), 
              df1=(resDPL[, "p"] - resDPL[1, "p"]), 
              df2=(resDPL[, "n"] - resDPL[, "p"]), 
              lower.tail=F)
resDPL <- resDPL[order(resDPL$pr), ]
#==============================fitting the best model - Im, and k differs between bacteria and fungi
ParmsFermFinal = rbind(Im, Im, k, k, yA, eta, a)
#Making 100 parameter estimates to get the best parameters and calculate error of estimate
##Run on server
n_iter = 100
source("Python/finalParallel.R")
registerDoParallel(n_iter)

PLfinalPars <- foreach(i=1:n_iter, .combine=rbind, .multicombine = TRUE,
                       .packages=c("FME", "dplyr", "deSolve", "reticulate", "ABCoptim")) %dopar% {
                         
                         finalParallel(i)
                         
                       }


PLfinalPars <- read.csv("PLfinalPars.csv", row.names = NULL)[, -1]
#adding goodness of fit
PLGF <- data.frame(R2 = numeric(length = 100), R2adj = numeric(length = 100), ll = numeric(length = 100),
                   AIC = numeric(length = 100), Form = numeric(length = 100), Gl = numeric(length = 100),
                   CO2 = numeric(length = 100), Fungi = numeric(length = 100), Bacteria = numeric(length = 100),
                   Cflush = numeric(length = 100), pH = numeric(length = 100)) 

for(i in 1:100){
  SimOut <- subMFitFermGroups(as.numeric(PLfinalPars[i, ]), c("Im", "k"),
                              IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Plesne"),
                              Soil = "Plese", Treatment = "Aerobic")
  PLGF[i, ] <- c(SimOut$errors[1:5], SimOut$R2all)
}
PLbestP <- which(PLGF$ll == max(PLGF$ll))
PLfinalPars[PLbestP, ]
PLGF[PLbestP, ]

PLfinalBest <- subMFitFermGroups(as.numeric(PLfinalPars[PLbestP, ]), c("Im", "k"),
                                 IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Plesne"),
                                 Soil = "Plesne", Treatment = "Aerobic")
#===========================================================================
#Glucose concentration
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Plesne") %>% 
  summarise(y = mean(Gl, na.rm = T), ySD = sd(Gl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Glucose (", mu, "mol C g ", DW^{-1}, ")"))) +
  geom_line(data = PLfinalBest$Simulation, aes(x = Time, y = Gl)) #+
  #geom_line(data = null_modelPL$Simulation, aes(x = Time, y = Gl), color = "red")
#Cumulative respiration from glucose
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Plesne") %>% 
  summarise(y = mean(CumulativeRg, na.rm = T), ySD = sd(CumulativeRg, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Cumulative ", CO[2], " production (", mu, "mol C ", g~(DW)^{-1}, ")"))) +
  geom_line(data = PLfinalBest$Simulation, aes(x = Time, y = CO2))#+
  #geom_line(data = null_modelPL$Simulation, aes(x = Time, y = CO2), color = "red")
#Chloroform flush
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Plesne" & CFlushgl<100) %>% 
  summarise(y = mean(CFlushgl, na.rm = T), ySD = sd(CFlushgl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in ",CHCl[3]~flush, " (", mu, "mol C ",g~(DW)^{-1}, ")" ))) +
  geom_line(data = PLfinalBest$Simulation, aes(x = Time, y = Cflush))# +
  #geom_line(data = null_modelPL$Simulation, aes(x = Time, y = Cflush), color = "red")
#pH
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Plesne") %>% 
  summarise(y = mean(pH, na.rm = T), ySD = sd(pH, na.rm = T),
            ExtraH = mean(ExtraH, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab("pH") +
  geom_line(data = PLfinalBest$Simulation, aes(x = Time, y = pH))# +
  #geom_line(data = null_modelPL$Simulation, aes(x = Time, y = pH), color = "red")

#PLFA
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Plesne") %>% 
  summarise(y = mean(Fungl, na.rm = T), ySD = sd(Fungl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab("Fungi") +
  geom_line(data = PLfinalBest$Simulation, aes(x = Time, y = Fungi))# +
  #geom_line(data = null_modelPL$Simulation, aes(x = Time, y = Fungi), color = "red")

IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Plesne") %>% 
  summarise(y = mean(Bacgl, na.rm = T), ySD = sd(Bacgl, na.rm = T),
            y2 = mean(Bacgl, na.rm = T), ySD = sd(Bacgl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab("Bacteria") +
  geom_line(data = PLfinalBest$Simulation, aes(x = Time, y = Bacteria))# +
  #geom_line(data = null_modelPL$Simulation, aes(x = Time, y = Bacteria), color = "red")

#Calculating residuals for all variables
PLmeans <- as.matrix(IE %>% filter(Treatment == "Aerobic" & Soil == "Plesne") %>% 
                       group_by(Time) %>% summarize(Glucose = mean(Gl, na.rm = T), 
                                                    CO2 = mean(CumulativeRg, na.rm = T),
                                                    Fungi = mean(Fungl, na.rm = T),
                                                    Bacteria = mean(Bacgl, na.rm = T),
                                                    Cflush = mean(CFlushgl, na.rm = T),
                                                    pH = mean(pH, na.rm = T)))
PLmeans[4, 6] <- NA
PLresiduals <- ((PLfinalBest$Yhat - PLmeans[, -1])/PLfinalBest$W)^2      
PLresiduals <- as.data.frame(cbind(PLmeans[, 1], PLresiduals))
colnames(PLresiduals)[1] <- c("Time")  
write.csv(PLresiduals, "PLresiduals.csv", row.names = F)

#Comparing DNA a PLFA based fungi to bacteria ratio again applying estimated conversion factors
InitialMBCPL = mean(as.numeric(IE[(IE$Time == 0 & IE$Soil == "Plesne" & IE$Treatment == "Aerobic"), "Cflush"]), na.rm = T)/0.24
InitialFPL = InitialMBCPL*tail(as.numeric(PLfinalPars[PLbestP, ]), 1)
InitialBPL = InitialMBCPL*(1 - tail(as.numeric(PLfinalPars[PLbestP, ]), 1))
kfPL <- mean(as.numeric(IE[(IE$Time == 0 & IE$Soil == "Plesne" & IE$Treatment == "Aerobic"), "PLFAf"]), na.rm = T)/InitialFPL
kbPL <- mean(as.numeric(IE[(IE$Time == 0 & IE$Soil == "Plesne" & IE$Treatment == "Aerobic"), "PLFAb"]), na.rm = T)/InitialBPL

PLresiduals$Soil <- c("Plesne")
CTresiduals$Soil <- c("Certovo")

residualsAll <- rbind(CTresiduals, PLresiduals)
write.csv(residualsAll, "residualsAll.csv", row.names = F)
