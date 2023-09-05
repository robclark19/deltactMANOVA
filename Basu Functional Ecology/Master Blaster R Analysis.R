############# Packages ############
#install.packages("pacman")
#install.packages("lme4","car","multcomp", "mulcompView", "ggplot2","dplyr", "emmeans")
#install.packages("dplyr")
#install.packages("emmeans","car","ggplot2")
library(dplyr)
library(emmeans)
library(lme4)


#install.packages("xtable")
########### Data ###################
rm(list = ls())

setwd("C:/Users/Myrmica/Desktop/Master Blaster")

saumik1 <- read.csv("Master Blaster Data.csv", header= TRUE)
str(saumik1)

#drop non-experimental observations for 2x3 style design
#this should always drop A+,W-
saumik4 <- subset(saumik1, Drop.for.2x3 == "No")

#get summary statistics (mean) across all technical replicates
#pool by technical replicate
exp.dat2 <- saumik4 %>%
  group_by(Weevil.Timing,Aphid.Type,BioSample) %>%
  summarize_at(.vars = vars(Chitenase:lectin), # everything between Chitenase and lectin
               .funs = mean, na.rm = TRUE)

# 
Chitenase.mod <- lm(Chitenase ~ Weevil.Timing:Aphid.Type, data=exp.dat2)
summary(Chitenase.mod)

#fine estimated marginal means for Weevil*Aphid
chit.lsm <- emmeans(Chitenase.mod, ~Weevil.Timing:Aphid.Type, sort=FALSE, type="response")
chit.lsm

# ways to change confidence intervals if necessary
confint(chit.lsm, adjust = "none", level = 0.95)

#convert to compact letter display for group letters
#currently no tukey adjustment

chit.cld <- cld(chit.lsm, adjust="none", sort=FALSE, type="response")
chit.cld

#quickly plotting cld output with dispersion (conf interval)
plot(chit.cld, by = "disp")

#convert emmeans object to dataframe for transformation (delta delta ct and 2^-x)
chitenase.dat <- as.data.frame(chit.cld)
chitenase.dat

# subtract the control from all the values and transform parameter estimates to 2^-x
# also fine upper and lower bounds of standard error
#probably jsue SE works with transformation

chitenase.final <- chitenase.dat %>%
  mutate(emmean = emmean - chitenase.dat[5, "emmean"],
         emmean = 2^(emmean),
         SE = 2^(-SE)
  )

chitenase.final

#ggplot code for a nice figure

chit.plot <-  ggplot(chitenase.final, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")),y=emmean, fill=factor(Aphid.Type,levels=c("Sham","Infective")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), position=position_dodge(0.9), width=0.5) +
  #geom_linerange(aes(ymin = emmean - lower.CL, ymax = upper.CL), position=position_dodge(0.9))
  theme_bw(base_size = 16) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Relative Expression", x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
    #geom_text(aes(label=c("","","","","Control","")), position=position_dodge(width=0.9)) +
  scale_y_continuous(breaks=c(0,1, 10, 20))
chit.plot

############ Final Figures ###############
#import data from steps above (double checked with pilot chitenase figures)
saumik5 <- exp.dat2

#Clare suggested 4 types of MANOVA combinations
#JA, SA, Defensin-related, others

################## JA MANOVA #########################
ja.manova <- manova(cbind(OPR1,LOX.2) ~ Aphid.Type*Weevil.Timing, data=saumik5)

#make AOV table for univariate statistics (each response variable)
ja.summary <- summary.aov(ja.manova)
ja.summary

#make MANOVA table, default is Pillai test
#manova table for linear combination of responses
ja.pillai <- summary(ja.manova)
ja.pillai

#partition results from AOV.SUMMARY to each response
str(ja.summary)
opr1.table <- ja.summary$` Response OPR1`
lox2.table <- ja.summary$` Response LOX2`

#view to verify its a single table
str(opr1.table)
#export html table to working directory
print(xtable(opr1.table), type="html", file="opr1.html")
#view html in R
tableHTML(lox2.table)

str(lox2.table)
#export html table to working directory
print(xtable(lox2.table), type="html", file="lox2.html")
#html table in R
lox.html <- tableHTML(lox2.table)
write_tableHTML(lox.html, "out.html")
#another way to see it
lox.k <- kable(lox2.table)
write.html(lox.k, "out.html")

#estimates for figure construction
ja.lsm <- emmeans(ja.manova, ~ Aphid.Type*Weevil.Timing|Gene, mult.name = "Gene", adjust="none")
ja.cld <- cld(ja.lsm, sort = FALSE, adjust="none")
ja.cld$tukey <- as.numeric(ja.cld$.group)
ja.cld

tableHTML(ja.cld)

#ja figure
#how to change order of levels of a factor
#aes(x=factor(Weevil.Treatment,levels=c("None","Before Aphids","After Aphids")),

#how to get delta and beta
#hist(h, main = expression(paste("Sampled values, ", mu, "=5, ", sigma,"=1")))

#example figure
ja.fig <- ggplot(ja.cld, aes(x=Weevil.Timing, y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")), label=.group)) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=emmean-(SE), ymax=emmean+(SE)), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 16) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y=expression(paste(Delta,"CT  " , beta, "-tubulin")), x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  facet_wrap( ~ Gene)
ja.fig

#how to change facet wrap title
#labeller=label_bquote(c("My New Title", "Group 2"))

################## SA MANOVA #########################

sa.manova <- manova(cbind(ICS1,PR1,Chitenase) ~ Aphid.Type*Weevil.Timing, data=saumik5)
summary.aov(sa.manova)
#estimates
sa.lsm <- emmeans(sa.manova, ~ Aphid.Type*Weevil.Timing|Gene, mult.name = "Gene", adjust="none")
sa.cld <- cld(sa.lsm, sort = FALSE, adjust="none")
sa.cld$tukey <- as.numeric(sa.cld$.group)
sa.cld

tableHTML(sa.cld)

#sa figure

# USE THIS TO SUBSET FOR FIGURES!!!!!!!!!!!
sa.small <- subset (sa.cld, Gene == c("PR1","Chitenase"))
sa.small

#example figure
sa.fig <- ggplot(sa.small, aes(x=Weevil.Timing, y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")), label=.group)) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=emmean-(SE), ymax=emmean+(SE)), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y=expression(paste(Delta,"CT  " , beta, "-tubulin")), x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  facet_wrap( ~ Gene)
sa.fig

################## Defensin MANOVA #########################

defensin.manova <- manova(cbind(defensin,ACS2) ~ Aphid.Type*Weevil.Timing, data=saumik5)
summary.aov(defensin.manova)
#estimates
defensin.lsm <- emmeans(defensin.manova, ~ Aphid.Type*Weevil.Timing|Gene, mult.name = "Gene", adjust="none")
defensin.lsm

defensin.cld <- cld(defensin.lsm, sort = FALSE, adjust="none")
defensin.cld$tukey <- as.numeric(defensin.cld$.group)
defensin.cld

tableHTML(defensin.cld)

#defensin figure

defensin.fig <- ggplot(defensin.cld, aes(x=Weevil.Timing, y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")), label=.group)) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=emmean-(SE), ymax=emmean+(SE)), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 

  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y=expression(paste(Delta,"CT  " , beta, "-tubulin")), x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
                axis.line.y = element_line(color="black", size = 0.5)) +
  #geom_text(nudge_y=(1)) +
facet_wrap(~ Gene)
defensin.fig

###################### JACKPOT JUST USE ABOVE #################

facet_wrap( ~ c("defensin","ACS2"))

################## Other MANOVA #########################

others.manova <- manova(cbind(PsOX.11,GA,AOX.3,NLectin,ChlAbBp) ~ Aphid.Type*Weevil.Timing, data=saumik5)
summary.aov(others.manova)
#estimates
others.lsm <- emmeans(others.manova, ~ Aphid.Type*Weevil.Timing|Gene, mult.name = "Gene", adjust="none")
others.cld <- cld(others.lsm, sort = FALSE, adjust="none")
others.cld$tukey <- as.numeric(others.cld$.group)
others.cld

tableHTML(others.cld)

#others figure

others.fig <- ggplot(others.cld, aes(x=Weevil.Timing, y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")), label=.group)) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=emmean-(SE), ymax=emmean+(SE)), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y=expression(paste(Delta,"CT  " , beta, "-tubulin")), x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  facet_wrap( ~ Gene, nrow=2)
others.fig

########## JA/SA expression #####################




############ Amino Acid NMDS ####################
amino.dat <- read.csv("amino nmds.csv", header=TRUE)
response.dat <- select(amino.dat, -"Sample.name", - "Asp.Asn", - "Aphids", - "Weevils")
response.dat <- amino.dat[c(5:15)] 
predictor.dat <- amino.dat[c(1:2)]


vegdist(response.dat, method="raup")
amino.mds <- metaMDS(response.dat, k=2, plot = TRUE)
plot(amino.mds)
ordiplot(amino.mds)

ord.fit <- envfit(amino.mds ~ Weevils, data=predictor.dat, perm=999)
ord.fit

###############amino acid MANOVA

amino.manova <- manova(cbind(Ser.Gln,Glu,Gly,His,Hser,Arg,Thr,Ala,Pro,Tyr,Val,Lys,Ile,Leu,Phe) ~ Aphids*Weevils, data=amino.dat)
summary(amino.manova)
summary.aov(amino.manova)
#estimates
others.lsm <- emmeans(others.manova, ~ Aphid.Type*Weevil.Timing|Gene, mult.name = "Gene", adjust="none")
others.cld <- cld(others.lsm, sort = FALSE, adjust="none")
others.cld$tukey <- as.numeric(others.cld$.group)
others.cld

tableHTML(others.cld)

#others figure

others.fig <- ggplot(others.cld, aes(x=Weevil.Timing, y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")), label=.group)) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=emmean-(SE), ymax=emmean+(SE)), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y=expression(paste(Delta,"CT  " , beta, "-tubulin")), x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  facet_wrap( ~ Gene, nrow=2)
others.fig

