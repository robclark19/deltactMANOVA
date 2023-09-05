# Statistical analysis of real-time PCR data (DDCT) using multivariate and generalized linear models in R.

#packages
library("tidyverse")
library("car")
library("emmeans")
library("multcomp")
library("lme4")
library("ggplot2")
library("ggpubr")


# Part 0 ######
# Basic example for different types of data selection for GLM/MANOVA

culex.dat <- read.csv("goodman culex.csv")
culex.dat

# Have known difference in CT and Delta CT with a given mean and variance defined ahead of time
# We know true difference in viral titer, run simulation with multiple experiments with 10 replicates (pull out 10 samples out of 1000)
# does the 90%CI overlap with true difference

# do this for each of the 3 things below

#A#####
#CT values (not normalized with delta ct) goes into model

a.glm <- manova(cbind(Int.Ct,Ref.Ct) ~ Insulin:Virus, data=culex.dat)
summary(a.glm)
summary.aov(a.glm)


#B####
#Delta CT raw values (not transformed) goes into model

b.glm <- glm(Raw.Delta.Ct ~ Insulin:Virus, data=culex.dat)
summary(b.glm)
Anova(b.glm)

#C####
#Transformed Delta Delta CT values goes into model

c.glm <- glm(Transformed.Delta.Delta.CT ~ Insulin:Virus, data=culex.dat)
summary(c.glm)
Anova(c.glm)



# Part 1a  ####
# Data Clark et al 2019 Ecology
# verification of PCR presence absence
plate.dat <- read.csv("pemv titer data.csv")
head(plate.dat)
str(plate.dat)

# model evaluting whether gel data is predicted by rtpcr data
pcr.lm <- glm(PEMV ~ PEMV.ct, family="binomial", data=plate.dat)
summary(pcr.lm)

#reversed model to see if PEMV+ hit was associated with higher viral titer (actually what we did in the manuscript)
pcr.lm2 <- glm(PEMV.ct ~ PEMV, data=plate.dat)
summary(pcr.lm2)

# new method using delta delta ct calculations from raw data




# Part 1a: GLM for mRNA ####
# from proc b, may be removed due to issues with data being redundant with a much better example in part 2


#defensive genes model (has 4 treatment levels)
pr.genes <- read.csv("pr genes all.csv")

# Data from Chisholm et al. in review in Proceedings of the Royal Society B
# modified from this paper (I did a MANOVA for two genes, but this should be simplified to a glm so its an easy example)


str(pr.genes)

pr.genes.2 <- pr.genes %>%
  group_by(location,aphid.treatment,bio.rep) %>% #treatment combinations you want to save
  summarize_at(.vars = vars(pr.1:opr.1), # everything for both pr genes
               .funs = mean, na.rm = TRUE)
pr.genes.2


# something is odd here and the transformation worked better BEFORE doing the statistics
pr.genes.2$pr.1.trans <- 2^-(pr.genes.2$pr.1)
pr.genes.2$opr.1.trans <- 2^-(pr.genes.2$opr.1)

opr.mod <- manova(cbind(pr.1.trans,opr.1.trans) ~ location*aphid.treatment, data=pr.genes.2)
summary.aov(opr.mod)
# summary(opr.mod)
# hist(pr.genes.2$pr.1.trans)

#get parameter estimates at each treatment, then present as compact letter display for posthoc tests
opr.mod.lsm <- emmeans(opr.mod, ~ location|aphid.treatment|virus, mult.name = "virus", adjust="none", response=TRUE)
opr.mod.lsm

opr.mod.cld <- cld(opr.mod.lsm, sort = FALSE, adjust="none", response=TRUE, Letters=c("abcdef"))
opr.mod.cld$.group=gsub(" ", "", opr.mod.cld$.group)

opr.mod.cld <- as.data.frame(opr.mod.cld)
opr.mod.cld <- subset(opr.mod.cld, aphid.treatment != "Control")

labels <- c(pr.1.trans = "PR1", opr.1.trans = "OPR1")
#labels <- c("italic('PR1')", "italic('OPR1')")

# updated figure with "Weevils only" treatment now labeled as "Sham and Weevils"

opr.mod.fig.2 <- ggplot(opr.mod.cld, aes(x=factor(aphid.treatment), 
                                         y=emmean, fill=factor(location,levels=c("Top", "Middle","Bottom")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Relative expression", x="Treatment") + 
  labs(fill="Location on Plant") +
  scale_x_discrete(labels=c("Infective", "Infective and Weevils", "Sham and Weevils")) +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  #scale_y_log10()  +
  geom_text(aes(x = aphid.treatment, y = (emmean+SE+1), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ virus, nrow=2, labeller=labeller(virus = labels)) +
  theme(strip.text = element_text(face = "italic"))
opr.mod.fig.2

# ultimately this should just show PR1 or OPR1



# Part 2: MANOVA ###########
# code from "master blaster" paper: plant gene expression in response to weevils, aphids, and aphid-vectored viruses
# this code only looks at mRNA transcript accumulation
# cited as Basu et al 2019/2020 in prep


# Data organization
# Load raw data and check structure
saumik1 <- read.csv("Master Blaster Data.csv", header=TRUE)
str(saumik1)

#Modify data to only use data for 2x3 design
# Data for MANOVA ####
saumik4 <- subset(saumik1, Drop.for.2x3 == "No")

#get summary statistics (mean) across all technical replicates
#pool by technical replicate
exp.dat2 <- saumik4 %>%
  group_by(Weevil.Timing,Aphid.Type,BioSample) %>%
  summarize_at(.vars = vars(Chitenase:lectin), # everything between Chitenase and lectin
               .funs = mean, na.rm = TRUE)

str(exp.dat2)

# MANOVA Model for all genes #
# Model includes 13 plant genes known to respond to herbivory
# Endogenous control gene for this experiment was CT for beta-tubulin in the same plant tissue sample

genes.manova <- manova(cbind(Chitenase,OPR1,PR1,PsOX.11,GA,AOX.3,ACS2,LOX.2,ICS1,defensin,UBQ,ChlAbBp,lectin) ~ Aphid.Type*Weevil.Timing, data=exp.dat2)

#MANOVA TABLE (response of all genes with Pillai index)

summary(genes.manova)

#AOV TABLE (response of each gene)
manova.table <- summary.aov(genes.manova)
manova.table



# Estimates for DCT ########
# Generates parameter estimates for average, untransformed CT value of each gene in response to aphid and weevil treatments

genes.lsm <- emmeans(genes.manova, ~ Aphid.Type*Weevil.Timing|Gene, mult.name = "Gene", adjust="none", response=TRUE)

# Forma values as compact letter display and run posthoc tests
genes.cld <- cld(genes.lsm, sort = FALSE, adjust="none", response=TRUE, Letters=c("abcdef"))
genes.cld$.group=gsub(" ", "", genes.cld$.group)
genes.cld

#for very small values, Use / options(scipen = 999) / command in console to disable scientific notation

genes.table <- as.data.frame(genes.cld)
genes.table

# To write output to a csv for supplementary materials, use the following command
# write.csv(genes.table, file = "all genes emmeans.csv")

# Subtract each control (sham,none) within each gene from the values. control should = 1 in each gene.
expression.table <- genes.table

expression.table <- genes.table %>%
  group_by(Gene) %>% # each gene considered a batch
  mutate(control = expression.table[which(Aphid.Type == "Sham" & Weevil.Timing == "None"), "emmean"],
         emmean = 2^(-emmean), # will operate with new emmean column
         SE = 2^(-SE))

# Make new names and columes to get the control for comparison
new.data <- genes.table %>%
  filter(Aphid.Type == "Sham", Weevil.Timing == "None")

# Make new table for "control" treatment condition
colnames(new.data) = c("Aphid.Type", "Weevil.Timing", "Gene", "cemmean","cSE","cdf","clower.CL","cupper.CL","c.group")
new.data = new.data %>% dplyr::select(Gene, cemmean)

# Left-join table now adds control as a new column for caculating the second delta in delta delta ct
total.data <- left_join(genes.table, new.data, by=c("Gene"))
head(total.data)

#Calculate delta delta ct with left_joined table 

final.table <- total.data %>% 
  mutate(emmean = cemmean - emmean,
         upper.SE = 2^(-(emmean-SE)),
         lower.SE = 2^(-(emmean+SE)),
         emmean = 2^(-emmean), # will operate with new emmean column
         SE = 2^(-SE)
  )


######## Figure R1 ############
str(final.table)

# Subset figure into multiple groups (this does not change analysis, but makes visualization easier)
# group 1 includes only three genes with related functionality

R1.table <- subset(final.table, Gene == "Chitenase" | Gene == "PR1" | Gene == "defensin")

# This command makes a string of labels that will appear on the header of each facet (gene)
labels.R1 <- c(Chitenase = "Chitenase", PR1 = "PR1", defensin = "Defensin")

# Plotting uses ggplot2 package entirely
r1.fig <- ggplot(R1.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=emmean-(SE), ymax=emmean+(SE)), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Fold Change in Expression", x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  #scale_y_continuous(breaks=c(0,1, 10, 25)) + 
  #the above code allows manipulation of the y axis scale
  #theme(legend.position = c(0.2, 0.7), legend.text = element_text(size=8), legend.title = element_text(size=8)) + 
  #the above code will move legend to top right of figure
  theme(legend.position="none") + #this one is just easier, and is not inset
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+3), label = .group), position = position_dodge(width = 0.8)) +
  # the above code adds tukey letters to the figure. Edit the value after SE (e.g. +3) if letters are not positioned clearly
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.R1)) +
  theme(strip.text = element_text(face = "italic"))
r1.fig

#Figure R2 ####

R2.table <- subset(final.table, Gene == "ACS2" | Gene == "OPR1" | Gene == "PsOX.11")
labels.R2 <- c(ACS2 = "ACS2", OPR1 = "OPR1", PsOX.11 = "POX11")

r2.fig <- ggplot(R2.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=emmean-(SE), ymax=emmean+(SE)), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Fold Change in Expression", x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  theme(legend.position="none") + #this drops the legend entirely so its only shown in r1
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+3), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.R2)) +
  theme(strip.text = element_text(face = "italic"))
r2.fig

#Figure R3 ####


R3.table <- subset(final.table, Gene == "GA" | Gene == "lectin" | Gene == "LOX.2")
labels.R3 <- c(GA = "GA2ox", lectin = "Lectin", LOX.2 = "LOX2")

r3.fig <- ggplot(R3.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=emmean-(SE), ymax=emmean+(SE)), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Fold Change in Expression", x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(legend.position="none") + #this drops the legend entirely so its only shown in r1
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+0.75), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.R3)) +
  theme(strip.text = element_text(face = "italic"))
r3.fig


#Figure R4 ####


R4.table <- subset(final.table, Gene == "ICS1" | Gene == "ChlAbBp" | Gene == "AOX.3")
labels.R4 <- c(ICS1 = "ICS1", ChlAbBp = "Chla/bBP", AOX.3 = "AO3")

r4.fig <- ggplot(R4.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=emmean-(SE), ymax=emmean+(SE)), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Fold Change in Expression", x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(legend.position="bottom") + #this slaps it down at the bottom
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+0.4), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.R4)) +
  theme(strip.text = element_text(face = "italic"))
r4.fig

# Arrange all four figures into one large panel using ggpubr package (ggarrange function)

figure.r.all <- ggarrange(r1.fig, r2.fig, r3.fig, r4.fig, labels = c("", "", "", ""), ncol = 1, nrow = 4)
figure.r.all




####### Part 3: Viral Titer ############

# Infective aphids with PEMV were placed on pea plants. Tissue was collected at 4 days and 10 days.
# PEMV is an RNA virus, so ddct can be calculated for viral coat protein that should be replicating during infection
# cited as Basu et al 2019/2020 in prep, but a seperate paper from the MANOVA data 

# Data for PEMV transcript accumulation ####
titer.dat <- read.csv("pemv titer.csv", header=TRUE)


########## 10 DPI Model ##############
# PEMV titer 10dpi, relative amount of viral titer, this is ten days post infection


# exclude day four observations
titer.ten <- subset(titer.dat, DPI != "Four")
str(titer.ten)

# 10 day duration GLM
c.ten.glm <- glm(DeltaCQ ~ Treatment, data=titer.ten)
Anova(c.ten.glm)
dct <- as.data.frame(cld(emmeans(c.ten.glm, ~ Treatment), sort=FALSE))

# group together
dct$.group=gsub(" ", "", dct$.group)

#delta ct for 10.dpi #
dct$Delta.CT <- 2^(-dct$emmean)
dct$upper.SEM <- 2^(-(dct$emmean-pre.log$SE))
dct$lower.SEM <- 2^(-(dct$emmean+pre.log$SE))
dct$tukey <- as.factor(dct$.group)

# 10 DPI Fig ####
raw.dct <- ggplot(dct, aes(x=Treatment, y=Delta.CT)) +
  geom_bar(stat="identity", width=0.5, position="dodge") +
  geom_errorbar(aes(ymin=lower.SEM, ymax=upper.SEM), position=position_dodge(0.5), width=0.2) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Relative Viral Titer (10 DPI)") + 
  scale_x_discrete(labels=c("Autoclave", "Rhizobia", "Nitrogen","Control")) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.title.x=element_blank()) +
  geom_text(aes(x = Treatment, y = upper.SEM+10, label = .group), position = position_dodge(width = 0.8))
raw.dct




########## 4 DPI Model ############

titer.four <- subset(titer.dat, DPI == "Four")
c.four.glm <- glm(DeltaCQ ~ Treatment, data=titer.four)
Anova(c.four.glm)
dct.four <- as.data.frame(cld(emmeans(c.four.glm, ~ Treatment), sort=FALSE))
dct.four

#fix group labels so posthoc test looks nice
dct.four$.group=gsub(" ", "", dct.four$.group)
dct.four

# delta ct for 4dpi 
dct.four$Delta.CT <- 2^(-dct.four$emmean)
dct.four$upper.SEM <- 2^(-(dct.four$emmean-pre.log$SE))
dct.four$lower.SEM <- 2^(-(dct.four$emmean+pre.log$SE))
dct.four$tukey <- as.factor(dct.four$.group)

# 4 DPI model ###
raw.dct.four <- ggplot(dct.four, aes(x=Treatment, y=Delta.CT)) +
  geom_bar(stat="identity", width=0.5, position="dodge") +
  geom_errorbar(aes(ymin=lower.SEM, ymax=upper.SEM), position=position_dodge(0.5), width=0.2) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Relative Viral Titer (4 DPI)") + 
  scale_x_discrete(labels=c("Autoclave", "Rhizobia", "Nitrogen","Control")) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.title.x=element_blank()) +
  geom_text(aes(x = Treatment, y = upper.SEM+10, label = .group), position = position_dodge(width = 0.8))
#ggtitle("PEMV Titer at 4 days post-infection")
raw.dct.four

# 10 days DDct ####
# estimate for control delta cq for treatment is -7.4831
# this uses control.noa as the baseline so its relative to 1 (control levels)

#test is with 10 days post infection

titer.ten$dd.ct <- titer.ten$DeltaCQ - titer.ten$ControlDeltaCQ

ten.glm <- glm(dd.ct ~ Treatment, data=titer.ten)
summary(ten.glm)

pre.log <- as.data.frame(cld(emmeans(ten.glm, ~ Treatment), sort=FALSE))
pre.log
pre.log$Delta.Delta.CT <- 2^(-pre.log$emmean)
pre.log$upper.SEM <- 2^(-(pre.log$emmean-pre.log$SE))
pre.log$lower.SEM <- 2^(-(pre.log$emmean+pre.log$SE))

pre.log$tukey <- as.factor(pre.log$.group)

pre.log

ten.titer.fig <- ggplot(pre.log, aes(x=Treatment, y=Delta.Delta.CT)) +
  geom_bar(stat="identity", width=0.5, position="dodge") +
  geom_errorbar(aes(ymin=lower.SEM, ymax=upper.SEM), position=position_dodge(0.5), width=0.2) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Fold change in expression compared to control") + 
  #scale_x_discrete(labels=c("Rhizobia", "Rhizobia + Pea Aphids", "Rhizobia + PEMV")) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.title.x=element_blank()) +
  geom_text(aes(x = Treatment, y = upper.SEM+0.1, label = tukey), position = position_dodge(width = 0.8))
ten.titer.fig



