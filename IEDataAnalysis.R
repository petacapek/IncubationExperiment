#==============================Libraries
library(dplyr)
library(ggplot2)
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
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(500-Glraw, na.rm = T), ySD = sd(500-Glraw, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Glucose (", mu, "mol C g ", DW^{-1}, ")"))) +
  geom_hline(yintercept = 500)
#pH
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(pH, na.rm = T), ySD = sd(pH, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("pH")))

IE %>% group_by(Time, Treatment, Soil) %>% filter(Time == 0) %>% 
  summarise(y = mean(pH, na.rm = T), ySD = sd(pH, na.rm = T))
IE %>% group_by(Time, Treatment, Soil) %>% filter(Time == 3) %>% 
  summarise(y = mean(pH, na.rm = T))
IE %>% group_by(Time, Treatment, Soil) %>% filter(Time == 58) %>% 
  summarise(y = mean(pH, na.rm = T))
#Amount of acetic acid in soil to decrease pH in:
#Certovo - aerobic from 4.3 to 4.04 is
((10^(-4.04) - 10^(-4.3))^2/1.8e-5)*((40*(1-0.28)/(40*0.28))/1000)*1e6*2 #in umol C/g DW
#Plesne - aerobic from 4.12 to 4.03 is
((10^(-4.03) - 10^(-4.12))^2/1.8e-5)*((40*(1-0.28)/(40*0.28))/1000)*1e6*2 #in umol C/g DW
#Certovo - anaerobic from 4.3 to 3.83 is
((10^(-3.83) - 10^(-4.3))^2/1.8e-5)*((40*(1-0.28)/(40*0.28))/1000)*1e6*2 #in umol C/g DW
#Plesne - aerobic from 4.12 to 3.70 is
((10^(-3.70) - 10^(-4.12))^2/1.8e-5)*((40*(1-0.28)/(40*0.28))/1000)*1e6*2 #in umol C/g DW

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
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
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
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(DOPw, na.rm = T), ySD = sd(DOPw, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(DOP, " in water (", mu, "mol P g ", DW^{-1}, ")")))
#SRPw
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(PO4w, na.rm = T), ySD = sd(PO4w, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(SRP, " in water (", mu, "mol P g ", DW^{-1}, ")")))
#Al
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
  summarise(y = mean(Al, na.rm = T), ySD = sd(Al, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(Al, " in water (", mu, "mol Al g ", DW^{-1}, ")")))
#Fe
IE %>% group_by(Time, Treatment, Soil) %>% #filter(Treatment == "Aerobic") %>% 
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

ggplot(IE, aes(PO4w, VmaxP)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment)) + theme_min

#===========Inhibition of enzyme activity by changing product concentration due to microbial uptake===========#
#=====================Script that calculate product concentration is uploaded=============================#
source("/mnt/580CBE2464C5F83D/pracovni/data_statistika/Junior/Inkubace/aerobni/Enzymes_raw/ECalculatorJunior.R")
#==========================================Additional libraries===============================================#
library(openxlsx)
library(reshape2)
library(bbmle)
library(FME)
library(ABCoptim)
library(foreach)
library(doParallel)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#===========================Test for phosphatase activity and samples on M1 plate
#==============M1 - sample IDs 2, 4, 9, 10 , 11, 12, 25, 32, 33, 34, 35, 37, 38, 39, 40, 41
Samples=rbind(c(2,4,9,10), c(11,12,25,32), c(33, 34,35,37), c(38:41))
DW=cbind(c(0.28, 0.29,0.3,0.28), rep(0.28, 4), rep(0.28, 4), rep(0.28, 4))
rownames(DW)<-c("A", "B", "C", "D")
#=====Calculation
M1<-ECalculatorJunior(dataset = "/mnt/580CBE2464C5F83D/pracovni/data_statistika/Junior/Inkubace/aerobni/Enzymes_raw/MUF M1.xlsx", #zdrojovy soubor
                      Kalibracni = c(1,5,10), 
                      Substraty = c(224/20, 224/2, 224),
                      Nmeasure = 37, empty = 1,
                      Samples = Samples,
                      DW = DW)
#=====Data
M1d<-M1$data
#=====Quick vizualization
ggplot(subset(M1d, Enzyme=="Phosph"), aes(Time, Pcorr2)) + geom_point(cex=6, aes(color=ConcEnzyme2)) +
  facet_wrap(.~Sample, scales="free") + xlab("Time (h)")

#=====Specifying initial concentration of product (Pcorr2 at t = 0) and substrate (ConcEnz2 - Pcorr2 at t = 0)
M1d$P0<-NA
M1d$S0<-NA
for(i in unique(M1d$Sample)){
  for(n in unique(M1d$ConcEnzyme2)){
    M1d[(M1d$Sample == i & M1d$ConcEnzyme2 == n), "P0"] <- 
      M1d[(M1d$Sample == i & M1d$ConcEnzyme2 == n & M1d$Time == 0), "Pcorr2"]
    M1d[(M1d$Sample == i & M1d$ConcEnzyme2 == n), "S0"] <- 
      M1d[(M1d$Sample == i & M1d$ConcEnzyme2 == n & M1d$Time == 0), "ConcEnzyme2"] - M1d[(M1d$Sample == i & M1d$ConcEnzyme2 == n & M1d$Time == 0), "Pcorr2"]
  }
}

Phosph <- subset(M1d, Enzyme == "Phosph")
#Add measured concentration of SRP in water extract and test kinetic equations' fit again
names(Phosph)

Phosph$P0all <- numeric(length = nrow(Phosph))
for(i in unique(Phosph$Sample)){
  Phosph[Phosph$Sample == i, "P0all"] <- Phosph[Phosph$Sample == i, "P0"] + IE[(IE$Treatment == "Aerobic" & IE$Time != 0.17 & IE$Time != 0.19 & IE$Time != 4), "PO4w"][i]
}

ggplot(Phosph, aes(Sample, P0)) + geom_point(cex = 6) + geom_point(cex = 6, color = "red", aes(Sample, P0all))

WI2 <- nls(Time ~ -(Km*log((S0 - Pcorr2)/S0) - Pcorr2)/Vmax, data = Phosph, subset = Sample == 2, start = list(Vmax = 1, Km = 10))
summary(WI2)
CI2 <- nls(Time ~ -(Km*(S0/Ki + P0/Ki + 1)*))
#==================Biological parameters
#Cflush
IE$CflushOut <- c("False")
IE[(IE$Soil == "Plesne" & IE$Treatment == "Aerobic" & IE$Time == 2), "CflushOut"] <- c("True")

IE %>% group_by(Time, Treatment, Soil) %>% filter(CflushOut == "False") %>% 
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
#====================================Functional groups
#Gram+
PLFAIdent <- read.csv("MicrobialGroupsPLFAIdentities.csv")
PLFAIdent[PLFAIdent$Group == "Gram positive", "PLFA"]
GplusCols <- numeric()
for(i in 1:length(PLFAIdent[PLFAIdent$Group == "Gram positive", "PLFA"])){
  GplusCols <- append(GplusCols, which(names(IE) == PLFAIdent[PLFAIdent$Group == "Gram positive", "PLFA"][i]))
}

GplusCols <- GplusCols
IE$GPlus <- numeric(length = nrow(IE))

for(i in 1:nrow(IE)){
  IE$GPlus[i] <- sum(IE[i, GplusCols], na.rm = T)
}

#Gram-
GminusCols <- numeric()
for(i in 1:length(PLFAIdent[PLFAIdent$Group == "Gram negative", "PLFA"])){
  GminusCols <- append(GminusCols, which(names(IE) == PLFAIdent[PLFAIdent$Group == "Gram negative", "PLFA"][i]))
}

IE$GMinus <- numeric(length = nrow(IE))

for(i in 1:nrow(IE)){
  IE$GMinus[i] <- sum(IE[i, GminusCols], na.rm = T)
}

#Fungi
FunCols <- numeric()
for(i in 1:length(PLFAIdent[PLFAIdent$Group == "Fungi", "PLFA"])){
  FunCols <- append(FunCols, which(names(IE) == PLFAIdent[PLFAIdent$Group == "Fungi", "PLFA"][i]))
}

IE$Fun <- numeric(length = nrow(IE))

for(i in 1:nrow(IE)){
  IE$Fun[i] <- sum(IE[i, FunCols], na.rm = T)
}

#Protozoa
ProCols <- numeric()
for(i in 1:length(PLFAIdent[PLFAIdent$Group == "Protozoa", "PLFA"])){
  ProCols <- append(ProCols, which(names(IE) == PLFAIdent[PLFAIdent$Group == "Protozoa", "PLFA"][i]))
}

IE$Prot <- numeric(length = nrow(IE))

for(i in 1:nrow(IE)){
  IE$Prot[i] <- sum(IE[i, ProCols], na.rm = T)
}

#Actinobacteria
ActCols <- numeric()
for(i in 1:length(PLFAIdent[PLFAIdent$Group == "Gram positive or Actinobacteria", "PLFA"])){
  ActCols <- append(ActCols, which(names(IE) == PLFAIdent[PLFAIdent$Group == "Gram positive or Actinobacteria", "PLFA"][i]))
}

IE$Act <- numeric(length = nrow(IE))

for(i in 1:nrow(IE)){
  IE$Act[i] <- sum(IE[i, ActCols], na.rm = T)
}

#Non-Specific
NSCols <- numeric()
for(i in 1:length(PLFAIdent[PLFAIdent$Group == "Non-specific bacterial marker", "PLFA"])){
  NSCols <- append(NSCols, which(names(IE) == PLFAIdent[PLFAIdent$Group == "Non-specific bacterial marker", "PLFA"][i]))
}

IE$NS <- numeric(length = nrow(IE))

for(i in 1:nrow(IE)){
  IE$NS[i] <- sum(IE[i, NSCols], na.rm = T)
}

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

IE %>% group_by(Time, Treatment, Soil) %>% filter(CflushOut == "False") %>% 
  summarise(y = mean(CFlushgl, na.rm = T), ySD = sd(CFlushgl, na.rm = T),
            y2 = mean(Cflush, na.rm = T), y2SD = sd(Cflush, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  #geom_point(cex = 6, pch = 21, fill = "green", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "green") +
  #geom_point(cex = 6, pch = 21, fill = "blue", alpha = 0.5, aes(x = Time, y = y3)) + geom_errorbar(aes(ymin = y3 - y3SD, ymax = y3 + y3SD), color = "blue") +
  #geom_point(cex = 6, pch = 21, fill = "red", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "red") +
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in ",CHCl[3]~flush, " (", mu, "mol C ",g~(DW)^{-1}, ")" )))

#Comparing glucose in Cflush with residual glucose C (everything that was not respired)
IE %>% group_by(Time, Treatment, Soil) %>% filter(CflushOut == "False") %>% 
  summarise(y = mean(CFlushgl, na.rm = T), ySD = sd(CFlushgl, na.rm = T),
            y2 = mean((Gl0 - Gl - CumulativeRg), na.rm = T), y2SD = sd((Gl0 - Gl - CumulativeRg), na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  #geom_point(cex = 6, pch = 21, fill = "red", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "red") +
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in ",CHCl[3]~flush, " (", mu, "mol C ",g~(DW)^{-1}, ")" )))

#===========================kec calculation
IE$kec <- with(IE, CFlushgl/(Gl0 - Gl - CumulativeRg))
IE$kec <- ifelse(IE$kec <= 0, NA, IE$kec)
IE$kec <- ifelse(IE$kec >= 1, NA, IE$kec)

IE %>% group_by(Time, Treatment, Soil) %>% filter(CflushOut == "False" & Treatment == "Aerobic") %>% 
  summarise(y = mean(kec, na.rm = T), ySD = sd(kec, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab("kec factor")

#Glucose C in PLFAs in umol C/g DW
NC <- names(IE)[52:75]
for(i in 1:length(NC)){NC[i] <- paste0(NC[i], "gl")}
NCframe <- as.data.frame(matrix(nrow = nrow(IE), ncol = length(NC)))
colnames(NCframe) <- NC
IE <- cbind(IE, NCframe)

NC2 <- names(IE)[52:75]
NC3 <- names(IE)[77:100]

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
  IE$PLFAgl[i] <- sum(IE[i, c(112:135)], na.rm = T)
}

IE[(IE$PLFAgl == 0 & IE$Time > 0), "PLFAgl"] <- NA
#Gram+
PLFAIdent <- read.csv("MicrobialGroupsPLFAIdentities.csv")
PLFAIdent[PLFAIdent$Group == "Gram positive", "PLFA"]
GplusCols <- numeric()
for(i in 1:length(PLFAIdent[PLFAIdent$Group == "Gram positive", "PLFA"])){
  GplusCols <- append(GplusCols, which(names(IE[, 112:135]) == paste0(PLFAIdent[PLFAIdent$Group == "Gram positive", "PLFA"][i], "gl")))
}

GplusCols <- GplusCols+111
IE$GPlusgl <- numeric(length = nrow(IE))

for(i in 1:nrow(IE)){
  IE$GPlusgl[i] <- sum(IE[i, GplusCols], na.rm = T)
}

IE[(IE$GPlusgl == 0 & IE$Time >0), "GPlusgl"] <- NA

#Gram-
GminusCols <- numeric()
for(i in 1:length(PLFAIdent[PLFAIdent$Group == "Gram negative", "PLFA"])){
  GminusCols <- append(GminusCols, which(names(IE[, 112:135]) == paste0(PLFAIdent[PLFAIdent$Group == "Gram negative", "PLFA"][i], "gl")))
}

GminusCols <- GminusCols+111
IE$GMinusgl <- numeric(length = nrow(IE))

for(i in 1:nrow(IE)){
  IE$GMinusgl[i] <- sum(IE[i, GminusCols], na.rm = T)
}

IE[(IE$GMinusgl == 0 & IE$Time >0), "GMinusgl"] <- NA

#Fungi
FunCols <- numeric()
for(i in 1:length(PLFAIdent[PLFAIdent$Group == "Fungi", "PLFA"])){
  FunCols <- append(FunCols, which(names(IE[, 112:135]) == paste0(PLFAIdent[PLFAIdent$Group == "Fungi", "PLFA"][i], "gl")))
}

FunCols <- FunCols+111
IE$Fungl <- numeric(length = nrow(IE))

for(i in 1:nrow(IE)){
  IE$Fungl[i] <- sum(IE[i, FunCols], na.rm = T)
}

IE[(IE$Fungl == 0 & IE$Time >0), "Fungl"] <- NA

#Protozoa
ProCols <- numeric()
for(i in 1:length(PLFAIdent[PLFAIdent$Group == "Protozoa", "PLFA"])){
  ProCols <- append(ProCols, which(names(IE[, 112:135]) == paste0(PLFAIdent[PLFAIdent$Group == "Protozoa", "PLFA"][i], "gl")))
}

ProCols <- ProCols+111
IE$Protgl <- numeric(length = nrow(IE))

for(i in 1:nrow(IE)){
  IE$Protgl[i] <- sum(IE[i, ProCols], na.rm = T)
}

IE[(IE$Protgl == 0 & IE$Time >0), "Protgl"] <- NA

#Actinobacteria
ActCols <- numeric()
for(i in 1:length(PLFAIdent[PLFAIdent$Group == "Gram positive or Actinobacteria", "PLFA"])){
  ActCols <- append(ActCols, which(names(IE[, 112:135]) == paste0(PLFAIdent[PLFAIdent$Group == "Gram positive or Actinobacteria", "PLFA"][i], "gl")))
}

ActCols <- ActCols+111
IE$Actgl <- numeric(length = nrow(IE))

for(i in 1:nrow(IE)){
  IE$Actgl[i] <- sum(IE[i, ActCols], na.rm = T)
}

IE[(IE$Actgl == 0 & IE$Time >0), "Actgl"] <- NA

#Non-Specific
NSCols <- numeric()
for(i in 1:length(PLFAIdent[PLFAIdent$Group == "Non-specific bacterial marker", "PLFA"])){
  NSCols <- append(NSCols, which(names(IE[, 112:135]) == paste0(PLFAIdent[PLFAIdent$Group == "Non-specific bacterial marker", "PLFA"][i], "gl")))
}

NSCols <- NSCols+111
IE$NSgl <- numeric(length = nrow(IE))

for(i in 1:nrow(IE)){
  IE$NSgl[i] <- sum(IE[i, NSCols], na.rm = T)
}

IE[(IE$NSgl == 0 & IE$Time >0), "NSgl"] <- NA
#===================Total PLFA
IE %>% group_by(Time, Treatment, Soil) %>% filter(CflushOut == "False") %>% 
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

View(IE %>% group_by(Time, Treatment, Soil) %>% filter(Soil == "Plesne" & Treatment == "Anaerobic") %>% 
       summarise(y = mean(PLFAgl, na.rm = T), ySD = sd(PLFAgl, na.rm = T)))
print(IE[(IE$Soil == "Plesne" & IE$Treatment == "Anaerobic" & IE$Time == 3), "PLFAgl"])
which(IE$PLFAgl == (IE[(IE$Soil == "Plesne" & IE$Treatment == "Anaerobic" & IE$Time == 3), "PLFAgl"])[1])
IE$PLFAgl[169] <- NA
IE$NSgl[169] <- NA
print(IE[(IE$Soil == "Plesne" & IE$Treatment == "Anaerobic" & IE$Time == 35), "PLFAgl"])
which(IE$PLFAgl == (IE[(IE$Soil == "Plesne" & IE$Treatment == "Anaerobic" & IE$Time == 35), "PLFAgl"])[3])
IE$PLFAgl[243] <- NA
IE$Fungl[243] <- NA
#======================G+ bacteria
IE %>% group_by(Time, Treatment, Soil) %>% filter(CflushOut == "False") %>% 
  summarise(y = mean(GPlusgl, na.rm = T), ySD = sd(GPlusgl, na.rm = T),
            y2 = mean(PLFAgl, na.rm = T), y2SD = sd(PLFAgl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  #geom_point(cex = 6, pch = 21, fill = "green", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "green") +
  #geom_point(cex = 6, pch = 21, fill = "blue", alpha = 0.5, aes(x = Time, y = y3)) + geom_errorbar(aes(ymin = y3 - y3SD, ymax = y3 + y3SD), color = "blue") +
  geom_point(cex = 6, pch = 21, fill = "red", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "red") +
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in PLFA (", mu, "mol C ",g~(DW)^{-1}, ")" )))
#======================G- bacteria
IE %>% group_by(Time, Treatment, Soil) %>% filter(CflushOut == "False") %>% 
  summarise(y = mean(GMinusgl, na.rm = T), ySD = sd(GMinusgl, na.rm = T),
            y2 = mean(PLFAgl, na.rm = T), y2SD = sd(PLFAgl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  #geom_point(cex = 6, pch = 21, fill = "green", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "green") +
  #geom_point(cex = 6, pch = 21, fill = "blue", alpha = 0.5, aes(x = Time, y = y3)) + geom_errorbar(aes(ymin = y3 - y3SD, ymax = y3 + y3SD), color = "blue") +
  geom_point(cex = 6, pch = 21, fill = "red", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "red") +
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in PLFA (", mu, "mol C ",g~(DW)^{-1}, ")" )))
#======================Fungi
IE %>% group_by(Time, Treatment, Soil) %>% filter(CflushOut == "False") %>% 
  summarise(y = mean(Fungl, na.rm = T), ySD = sd(Fungl, na.rm = T),
            y2 = mean(PLFAgl, na.rm = T), y2SD = sd(PLFAgl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  #geom_point(cex = 6, pch = 21, fill = "green", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "green") +
  #geom_point(cex = 6, pch = 21, fill = "blue", alpha = 0.5, aes(x = Time, y = y3)) + geom_errorbar(aes(ymin = y3 - y3SD, ymax = y3 + y3SD), color = "blue") +
  geom_point(cex = 6, pch = 21, fill = "red", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "red") +
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in PLFA (", mu, "mol C ",g~(DW)^{-1}, ")" )))
#======================Actinobacteria
IE %>% group_by(Time, Treatment, Soil) %>% filter(CflushOut == "False") %>% 
  summarise(y = mean(Actgl, na.rm = T), ySD = sd(Actgl, na.rm = T),
            y2 = mean(PLFAgl, na.rm = T), y2SD = sd(PLFAgl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  #geom_point(cex = 6, pch = 21, fill = "green", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "green") +
  #geom_point(cex = 6, pch = 21, fill = "blue", alpha = 0.5, aes(x = Time, y = y3)) + geom_errorbar(aes(ymin = y3 - y3SD, ymax = y3 + y3SD), color = "blue") +
  geom_point(cex = 6, pch = 21, fill = "red", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "red") +
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in PLFA (", mu, "mol C ",g~(DW)^{-1}, ")" )))
#======================Non-specific
IE %>% group_by(Time, Treatment, Soil) %>% filter(CflushOut == "False") %>% 
  summarise(y = mean(NSgl, na.rm = T), ySD = sd(NSgl, na.rm = T),
            y2 = mean(PLFAgl, na.rm = T), y2SD = sd(PLFAgl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  #geom_point(cex = 6, pch = 21, fill = "green", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "green") +
  #geom_point(cex = 6, pch = 21, fill = "blue", alpha = 0.5, aes(x = Time, y = y3)) + geom_errorbar(aes(ymin = y3 - y3SD, ymax = y3 + y3SD), color = "blue") +
  geom_point(cex = 6, pch = 21, fill = "red", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "red") +
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in PLFA (", mu, "mol C ",g~(DW)^{-1}, ")" )))
#======================Protozoa
IE %>% group_by(Time, Treatment, Soil) %>% filter(CflushOut == "False") %>% 
  summarise(y = mean(Protgl, na.rm = T), ySD = sd(Protgl, na.rm = T),
            y2 = mean(PLFAgl, na.rm = T), y2SD = sd(PLFAgl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  #geom_point(cex = 6, pch = 21, fill = "green", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "green") +
  #geom_point(cex = 6, pch = 21, fill = "blue", alpha = 0.5, aes(x = Time, y = y3)) + geom_errorbar(aes(ymin = y3 - y3SD, ymax = y3 + y3SD), color = "blue") +
  geom_point(cex = 6, pch = 21, fill = "red", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "red") +
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in PLFA (", mu, "mol C ",g~(DW)^{-1}, ")" )))

#======================Bac vs Fungi
IE %>% group_by(Time, Treatment, Soil) %>% filter(CflushOut == "False") %>% 
  summarise(y = mean(Fungl/(GPlusgl + GMinusgl + NSgl + Actgl), na.rm = T), 
            ySD = sd(Fungl/(GPlusgl + GMinusgl + NSgl + Actgl), na.rm = T),
            y2 = mean(Fun/(GPlus + GMinus + NS + Act), na.rm = T), 
            y2SD = sd(Fun/(GPlus + GMinus + NS + Act), na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  #geom_point(cex = 6, pch = 21, fill = "green", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "green") +
  #geom_point(cex = 6, pch = 21, fill = "blue", alpha = 0.5, aes(x = Time, y = y3)) + geom_errorbar(aes(ymin = y3 - y3SD, ymax = y3 + y3SD), color = "blue") +
  geom_point(cex = 6, pch = 21, fill = "red", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "red") +
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in PLFA (", mu, "mol C ",g~(DW)^{-1}, ")" )))

#======================G+ vs G-
IE %>% group_by(Time, Treatment, Soil) %>% filter(CflushOut == "False") %>% 
  summarise(y = mean(GPlusgl/GMinusgl, na.rm = T), 
            ySD = sd(GPlusgl/GMinusgl, na.rm = T),
            y2 = mean(GPlus/GMinus, na.rm = T), 
            y2SD = sd(GPlus/GMinus, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  #geom_point(cex = 6, pch = 21, fill = "green", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "green") +
  #geom_point(cex = 6, pch = 21, fill = "blue", alpha = 0.5, aes(x = Time, y = y3)) + geom_errorbar(aes(ymin = y3 - y3SD, ymax = y3 + y3SD), color = "blue") +
  geom_point(cex = 6, pch = 21, fill = "red", alpha = 0.5, aes(x = Time, y = y2)) + geom_errorbar(aes(ymin = y2 - y2SD, ymax = y2 + y2SD), color = "red") +
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in PLFA (", mu, "mol C ",g~(DW)^{-1}, ")" )))

#===========================kPLFA calculation
IE$kPLFA <- with(IE, PLFAgl/(Gl0 - Gl - CumulativeRg))
IE$kPLFA <- ifelse(IE$kPLFA <= 0, NA, IE$kPLFA)
IE$kPLFA <- ifelse(IE$kPLFA >= 1, NA, IE$kPLFA)

IE %>% group_by(Time, Treatment, Soil) %>% filter(CflushOut == "False") %>% 
  summarise(y = mean(kPLFA, na.rm = T), ySD = sd(kPLFA, na.rm = T),
            y2 = mean(kec, na.rm = T), y2SD = sd(kec, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_wrap(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  #geom_point(aes(x = Time, y = y2), cex = 6, pch = 21, fill = "red", alpha = 0.5, show.legend = F) +
  xlab("Time (days)") + ylab("kplfa factor")

#==============================================SubMicrobial model==============================================
library(deSolve)
library(reticulate)
library(FME)
library(ABCoptim)
library(caRamel)
#============================Basic model without death and explicit PLFA======================================#
#===================Basic model and respective solver
source_python("Python/subM0.py")
source_python("Python/subMSolver0.py")
#===================================================
#Objective function and initial parameter guess
source("Python/subMObjective0.R")
p0 <- c(6.5, 221, 0.95, 3, 4e-4, 2, 0.41, 0.24)#, 0.05
pl <- c(0.01, 0.1, 0.1, 0.5, 1e-7, 1e-7, 0.3, 0.04)#, 0.003
pu <- c(20, 400, 1, 25, 1e-2, 5, 1, 0.6)#, 0.1
#==================================================Plesne - aerobic conditions
IEactive <- subset(IE, Treatment == "Aerobic" & Soil == "Plesne")
#==================================================#
#~~~~~~~~~~~~~Parameters estimation~~~~~~~~~~~~~~~~#
#==================================================#
##First guess by MCMC 
PAguess <- modMCMC(subMObjective0, p = p0, lower = pl, upper = pu, niter = 30000)
summary(PAguess)

#Estimate
PAOpt <- abc_optim(fn = subMObjective0,
                   par = as.numeric(apply(PAguess$pars, 2, quantile, 0.5)), 
                   lb = as.numeric(apply(PAguess$pars, 2, quantile, 0.05)), 
                   ub = as.numeric(apply(PAguess$pars, 2, quantile, 0.95)))
PAOpt[6]
round(PAOpt$par, 3)
#=============================================Model fit
source("Python/subMFit0.R")
PAout <- subMFit0(PAOpt$par, "Plesne", "Aerobic")
PAout$errors
PAout$R2all

#Glucose concentration
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Plesne") %>% 
  summarise(y = mean(Gl, na.rm = T), ySD = sd(Gl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil) + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Glucose (", mu, "mol C g ", DW^{-1}, ")"))) +
  geom_line(data = PAout$Simulation, aes(x = Time, y = Gl))
#Cumulative respiration from glucose
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Plesne") %>% 
  summarise(y = mean(CumulativeRg, na.rm = T), ySD = sd(CumulativeRg, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Cumulative ", CO[2]^{"Glucose"}, " production (", mu, "mol C ", g~(DW)^{-1}, ")"))) +
  geom_line(data = PAout$Simulation, aes(x = Time, y = CO2))
#kec
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Plesne") %>% 
  summarise(y = mean(kec, na.rm = T), ySD = sd(kec, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(italic(k[ec])))) +
  geom_line(data = PAout$Simulation, aes(x = Time, y = kec))
#Chloroform flush
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Plesne" & CFlushgl < 100) %>% 
  summarise(y = mean(CFlushgl, na.rm = T), ySD = sd(CFlushgl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in ",CHCl[3]~flush, " (", mu, "mol C ",g~(DW)^{-1}, ")" ))) +
  geom_line(data = PAout$Simulation, aes(x = Time, y = Cflush))
#kplfa
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Certovo") %>% 
  summarise(y = mean(kPLFA, na.rm = T), ySD = sd(kPLFA, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste(italic(k[PLFA])))) +
  geom_line(data = PAout$Simulation, aes(x = Time, y = kPLFA))
#PLFA
IE %>% group_by(Time, Treatment, Soil) %>% filter(Treatment == "Aerobic" & Soil == "Certovo") %>% 
  summarise(y = mean(PLFAgl, na.rm = T), ySD = sd(PLFAgl, na.rm = T)) %>% 
  ggplot(aes(Time, y)) + geom_point(cex = 6, pch = 21, aes(fill = Treatment), show.legend = F) +
  scale_fill_manual(values = c("white", "grey")) + theme_min + theme(legend.title = element_blank(),
                                                                     legend.position = c(0.85, 0.85)) +
  facet_grid(Treatment~Soil, scales = "free") + geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) + 
  xlab("Time (days)") + ylab(expression(paste("Glucose-C in PLFA (", mu, "mol C ",g~(DW)^{-1}, ")" ))) +
  geom_line(data = PAout$Simulation, aes(x = Time, y = PLFA))
#=========================Multi-objective optimization
# Number of objectives
nobj <- 6
# Number of variables
nvar <- length(p0)
# All the objectives are to be minimized
minmax <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
# Define the bound constraints
bounds <- matrix(data = 1, nrow = nvar, ncol = 2)
bounds[, 1] <- as.numeric(apply(CAguess$pars, 2, quantile, 0.05))
bounds[, 2] <- as.numeric(apply(CAguess$pars, 2, quantile, 0.95))
# Caramel optimization
CaOptMulti <-
  caRamel(nobj = nobj,
          nvar = nvar,
          minmax = minmax,
          bounds = bounds,
          func = subMObjective,
          popsize = 100,
          archsize = 300,
          maxrun = 1500,
          prec = matrix(1.e-3, nrow = 1, ncol = nobj),
          carallel = 0, graph = FALSE)
CaOptMulti$success
apply(CaOptMulti$objectives, 1, sum)

