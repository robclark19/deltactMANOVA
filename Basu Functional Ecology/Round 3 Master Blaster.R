# This is updated Master Blaster as of 2 19 2021 #
# Round 1 and 2 Master Blaster are deprecated #

# packages ############
library(dplyr)
library(emmeans)
library(multcomp)
library(lme4)
library(vegan)
library(ggplot2)
library(ggpubr)
library(car)
library("xlsx")
library("svglite")

# Gene Data ####
# original data
# saumik1 <- read.csv("Master Blaster Data.csv", header= TRUE)
# str(saumik1)

# here is the new file I made, it is just raw delta ct from machine output
# to make it fit within the same data frame, i had to trip 2 replicates of lox, defensin and ics2
# otherwise the data would be very, very unbalanced (5 reps vs 3 reps for all other genes)

saumik1 <- read.xlsx("Saumiks Experiment 1 Data.xlsx", 1, header=TRUE)
str(saumik1)

#Modify data to only use data for 2x3 design

saumik4 <- subset(saumik1, Drop.for.2x3 == "No")

# get summary statistics (mean) across all technical replicates
# this means we are pooling by technical replicate, and not including technical rep as a random effect!
exp.dat2 <- saumik4 %>%
  group_by(Weevil.Timing,Aphid.Type,BioSample) %>%
  summarize_at(.vars = vars(Chitenase:lectin), # everything between Chitenase and lectin
               .funs = mean, na.rm = TRUE)


# Full MANOVA model with raw cycle threshold values of each gene predicted by treatment
# Gene-MANOVA ###########
# original pre feb 2021 manova
# UBQ removed because it has singularity warning

# genes.manova <- manova(cbind(Chitenase,OPR1,PR1,PsOX.11,GA,AOX.3,ACS2,LOX.2,ICS1,defensin,ChlAbBp,lectin) ~ Aphid.Type*Weevil.Timing, data=exp.dat2)

# new MANOVA with only genes of interest for feb 2021 revisions
# ICS2, AOX3 (aox3), lox2, ga2ox (ga), pr1, defensin (ddr230), lectin

genes.manova <- manova(cbind(PR1,GA,AOX.3, LOX.2,ICS1,defensin,lectin) ~ Aphid.Type*Weevil.Timing, data=exp.dat2)


# MANOVA table (response of all genes with pillai index)
summary(genes.manova)
#AOV TABLE (response of each gene with a separate f-statistic)
manova.table <- summary.aov(genes.manova)
manova.table

# get parameter estimates for all genes based on aphids and weevils
# then present as compact letter display for posthoc tests

genes.lsm <- emmeans(genes.manova, ~ Aphid.Type*Weevil.Timing|Gene, mult.name = "Gene", adjust="none", response=TRUE)
genes.lsm

genes.cld <- cld(genes.lsm, sort = FALSE, adjust="none", response=TRUE, Letters=c("abcdef"))
genes.cld$.group=gsub(" ", "", genes.cld$.group)
genes.cld

# make parameter estimates into a data.table
# for very small values, I used / options(scipen = 999)...
# this is a command in console to disable scientific notation from displaying

genes.table <- as.data.frame(genes.cld)
genes.table

# write output to a csv for supplementary materials
# write.csv(genes.table, file = "all genes emmeans.csv")

#subtract each control (sham,none) within each gene from the values. 
# control should = 1 in each gene for delta delta ct method
expression.table <- genes.table

# I've calculated this in two ways, but "final.table" is used in the analyses
expression.table <- genes.table %>%
  group_by(Gene) %>% # each gene considered a batch
  mutate(control = expression.table[which(Aphid.Type == "Sham" & Weevil.Timing == "None"), "emmean"],
         emmean = 2^(-emmean), # will operate with new emmean column
         SE = 2^(-SE),
         lower.CL = 2^(-lower.CL), # this conserves the 95% CI if you want to use them later
         upper.CL = 2^(-upper.CL)
         )

# make new names and columns to get the control for comparison
# This is effectively the "control" set of data for delta delta ct comparison
new.data <- genes.table %>%
  filter(Aphid.Type == "Sham", Weevil.Timing == "None")

# make a new frame with values plus the values of the controls
# for example, cemmean is the estimated marginal mean for sham aphids and no weevils
colnames(new.data) = c("Aphid.Type", "Weevil.Timing", "Gene", "cemmean","cSE","cdf","clower.CL","cupper.CL","c.group")
new.data = new.data %>% dplyr::select(Gene, cemmean)

#left join table now adds control as a new column for finally calculating the second delta in delta delta ct
total.data <- left_join(genes.table, new.data, by=c("Gene"))
head(total.data)

# do delta delta ct with left_joined table 

final.table <- total.data %>% 
  mutate(emmean = emmean - cemmean,
         upper.SE = 2^(-(emmean-SE)),
         lower.SE = 2^(-(emmean+SE)),
         emmean = 2^(-emmean),
         SE = 2^(-SE),
         lower.CL = 2^(-lower.CL), 
         upper.CL = 2^(-upper.CL))

# check to make sure everything is included. "Final table" is used for analyses from the full anova


# 2018 Plots ####
# plot figures by categories determined ahead of time 
# (sa related genes, ja related genes, defensin and others)

### Pr1,chit,ics1 ########

sa.table <- subset(final.table, Gene == "PR1" | Gene == "Chitenase" | Gene == "ICS1")

labels.sa <- c(Chitenase = "Chitenase", PR1 = "PR1", ICS1 = "ICS1")

sa.fig <- ggplot(sa.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=lower.SE, ymax=upper.SE), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Fold Change in Expression", x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  scale_y_continuous(breaks=c(0,1, 10, 25)) +
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+1), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.sa)) +
  theme(strip.text = element_text(face = "italic"))
sa.fig

### opr1,lox2 #####
ja.table <- subset(final.table, Gene == "OPR1" | Gene == "LOX.2")
labels.ja <- c(OPR1 = "OPR1", LOX.2 = "LOX 2")

ja.fig <- ggplot(ja.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=lower.SE, ymax=upper.SE), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Fold Change in Expression", x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  scale_y_continuous(breaks=c(0,1, 10, 25)) +
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+1), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.ja)) +
  theme(strip.text = element_text(face = "italic"))
ja.fig

### defensin, acs2 ########

def.table <- subset(final.table, Gene == "defensin" | Gene == "ACS2")
as.data.frame(def.table)

labels.def <- c(defensin = "Defensin", ACS2 = "ACS2")

def.fig <- ggplot(def.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=lower.SE, ymax=upper.SE), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Fold Change in Expression", x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  scale_y_continuous(breaks=c(0,1, 10, 25)) +
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+2), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.def)) +
  theme(strip.text = element_text(face = "italic"))
def.fig

### other genes ##########
# GA, AOX, UBQ, chlapb, lectin

others.table <- subset(final.table, Gene == "GA" | Gene == "AOX.3" | Gene == "UBQ" | Gene == "ChlAbBp" | Gene == "lectin")
as.data.frame(others.table)

labels.others <- c(GA = "GA", AOX.3 = "AOX 3", UBQ = "UBC4", ChlAbBp = "ChlAbBp", lectin = "Lectin")

others.fig <- ggplot(others.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=lower.SE, ymax=upper.SE), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Fold Change in Expression", x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  scale_y_continuous(breaks=c(0,1, 10, 25)) +
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+1), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.others)) +
  theme(strip.text = element_text(face = "italic"))
others.fig

## psOx.11 ##########
ps.table <- subset(final.table, Gene == "PsOX.11")
as.data.frame(ps.table)

labels.ps <- c(PsOX.11 = "PsOX 11")

ps.fig <- ggplot(ps.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=lower.SE, ymax=upper.SE), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Fold Change in Expression", x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  scale_y_continuous(breaks=c(0,1, 10, 25, 50)) +
  geom_text(aes(x = Weevil.Timing, y = (upper.SE+5), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.ps)) +
  theme(strip.text = element_text(face = "italic"))
ps.fig

# 2019 Plots ####
#this is the updated layout based on manuscript revisions with dave july 18 2019
#row 1 is chitenase, pr1, defensin
#row 2 is asc2, opr1, psox11
#row 3 is GA, lectin, lox2
#row 4 is ICS1, ChalBP, aox3

R1.table <- subset(final.table, Gene == "Chitenase" | Gene == "PR1" | Gene == "defensin")
as.data.frame(R1.table)

labels.R1 <- c(Chitenase = "Chitenase", PR1 = "PR1", defensin = "Defensin")

r1.fig <- ggplot(R1.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=lower.SE, ymax=upper.SE), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Fold Change in Expression", x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  #scale_y_continuous(breaks=c(0,1, 10, 25)) +
  #theme(legend.position = c(0.2, 0.7), legend.text = element_text(size=8), legend.title = element_text(size=8)) + 
  #the above code will move legend to top right of figure
  theme(legend.position="none") + #this one is just easier, and is not inset
    geom_text(aes(x = Weevil.Timing, y = (emmean+SE+3), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.R1)) +
  theme(strip.text = element_text(face = "italic"))
r1.fig

R2.table <- subset(final.table, Gene == "ACS2" | Gene == "OPR1" | Gene == "PsOX.11")
labels.R2 <- c(ACS2 = "ACS2", OPR1 = "OPR1", PsOX.11 = "POX11")

r2.fig <- ggplot(R2.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=lower.SE, ymax=upper.SE), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Fold Change in Expression", x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  theme(legend.position="none") + #this drops the legend entirely so its only shown in r1
  #scale_y_continuous(breaks=c(0,1, 10, 25)) +
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+3), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.R2)) +
  theme(strip.text = element_text(face = "italic"))
r2.fig



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

R4.table <- subset(final.table, Gene == "ICS1" | Gene == "ChlAbBp" | Gene == "AOX.3")
labels.R4 <- c(ICS1 = "ICS1", ChlAbBp = "Chla/bBP", AOX.3 = "AO3")

r4.fig <- ggplot(R4.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=lower.SE, ymax=upper.SE), position=position_dodge(0.8), width=0.5) +
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

#arrange into one large panel using ggpubr package
figure.r.all <- ggarrange(r1.fig, r2.fig, r3.fig, r4.fig, 
                          labels = c("", "", "", ""), ncol = 1, nrow = 4)
figure.r.all

# 2020 Plots ####
# ics1, lox2, ao3, ga2ox

## 2020 Fig 2 ####
f2.table <- subset(final.table, Gene == "ICS1" | Gene == "LOX.2" | Gene == "AOX.3" | Gene == "GA")
labels.f2 <- c(ICS1 = "ICS1", LOX.2 = "LOX2", AOX.3 = "AO3", GA = "GA2ox")

# reorder everything in the original data to get the right values

neworder <- c("ICS1","AOX.3","LOX.2","GA")
f2.table <- arrange(transform(f2.table,
                          Gene=factor(Gene,levels=neworder)),Gene)


f2a.fig<- ggplot(f2.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=lower.SE, ymax=upper.SE), position=position_dodge(0.8), width=0.5) +
  #make sure to use upper and lower se, not raw SE!
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  #this formatting allows me to insert greek symbols
  labs(y=expression("Fold Change in Expression ( "~Delta*Delta~"CT)"), x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(legend.position="bottom") + #this slaps it down at the bottom
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  # facetwrap does some silly things so you have to make 2 sets of tukey letters and make one invisible...
  # just... just don't ask.
  geom_text(aes(x = Weevil.Timing, y = upper.SE+(upper.SE*0.2), label = .group), alpha=0, position = position_dodge(width = 0.8)) +
  geom_text(aes(x = Weevil.Timing, y = upper.SE, label = .group), vjust=-0.5, position = position_dodge(width = 0.8)) +
  #places tukey letters 20% higher than the upper bounds of the error bars
  facet_wrap( ~ Gene, scales="free_y", labeller=labeller(Gene = labels.f2), nrow=2) + 
  #free_y lets you vary the y axis for each facet
  theme(strip.text = element_text(face = "italic"))
  # plotting y on a log scale in ggplot makes error bars less confusing, but it will alter the control to be the axis
  # scale_y_continuous(trans='log10')
f2a.fig


# write figure 4 to folder, use arguments to modify size or file type!
ggsave(filename = "Fig 2 Master Blaster.svg", plot = f2a.fig, device = "svg",
       width = 6, height = 5, units = "in")


## 2020 Fig 3 ####

f3.table <- subset(final.table, Gene == "PR1" | Gene == "defensin" | Gene == "lectin")
labels.f3 <- c(PR1 = "PR1", defensin = "DDR230", lectin = "Lectin")

# reorder everything in the original data to get the right values

neworder <- c("PR1","defensin","lectin")
f3.table <- arrange(transform(f3.table,
                              Gene=factor(Gene,levels=neworder)),Gene)


f3.fig <- ggplot(f3.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=lower.SE, ymax=upper.SE), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y=expression("Fold Change in Expression ( "~Delta*Delta~"CT)"), x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(legend.position="bottom") + #this slaps it down at the bottom
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  geom_text(aes(x = Weevil.Timing, y = upper.SE+(upper.SE*0.2), label = .group), alpha=0, position = position_dodge(width = 0.8)) +
  geom_text(aes(x = Weevil.Timing, y = upper.SE, label = .group), vjust=-0.5, position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, scales="free_y", labeller=labeller(Gene = labels.f3), nrow=1) +
  theme(strip.text = element_text(face = "italic"))
f3.fig


# write figure 4 to folder, use arguments to modify size or file type!
ggsave(filename = "Fig 3 Master Blaster.svg", plot = f3.fig, device = "svg",
       width = 6, height = 4, units = "in")



# Amino NMDS ####################
## Amino Data ####
amino.dat <- read.csv("plant aminos.csv", header=TRUE)
# response.dat <- amino.dat[c(6:21)] 
# predictor.dat <- amino.dat[c(2:3)]

## NMDS1 ########
#make new object so commands below dont edit primary data
allmatrix.dat <- amino.dat

#drop env columns to make a species matrix
allmatrix.dat$Aphids <- NULL
allmatrix.dat$Weevils <- NULL
colSums(allmatrix.dat)


amino.distance <- vegdist(allmatrix.dat, method="bray")
amino.distance 

amino.mds <- metaMDS(allmatrix.dat, k=2, plot=TRUE)
stressplot(amino.mds, amino.distance)
amino.mds

attach(amino.dat)
detach(amino.dat)

ordiplot(amino.mds, type = "none")



plot(amino.mds)
ordiplot(amino.mds)
orditorp(amino.mds)
ordiellipse(amino.mds,amino.dat$Weevils, display="species")

ord.fit <- envfit(amino.mds ~ Weevils, data=predictor.dat, perm=999)
ord.fit

plot(ord.fit)

#ordination plot (starts with blank ordination that layers can be added on)
ordiplot(amino.mds, type="none")
#species plot (add insect species first)
orditorp(amino.mds,display="species",col="black",air=0.2,cex=1)
# draws a shape around it based on the environmental variable of interest
ordihull(amino.mds, groups=Weevils, draw="polygon",col="grey90",label=TRUE)

# ellipse might be better
ordiellipse(amino.mds,amino.dat$Weevils, display="species", label=TRUE)

## NMDS2 ########
#make new object so commands below dont edit primary data
allmatrix.dat <- all.dat

#drop env columns to make a species matrix
allmatrix.dat$Aphids <- NULL
allmatrix.dat$Weevils<- NULL

colSums(allmatrix.dat)


#check
colSums(allmatrix.dat)

#run nmds with k=2 dimensions
all.mds <- metaMDS(allmatrix.dat, k=2, plot=TRUE, na.rm=TRUE)
all.mds
# sites pooling removes some of the nasty warnings. the rest should be ok now.

attach(all.dat)
#just plot points (sites and species)
plot(all.mds)

#ordination plot (starts with blank ordination that layers can be added on)
ordiplot(all.mds, type="none")
#species plot (add insect species first)
orditorp(all.mds,display="species",col="black",air=0.2,cex=1)
# draws a shape around it based on the environmental variable of interest
ordihull(all.mds, groups=Weevils, draw="polygon", col="grey90",labels=F)

ordiarrows(ord.fit)


# ellipse might be better but this borks it
ordiellipse(amino.mds,amino.dat$Weevils, display="species")


ord.fit <- envfit(all.mds ~ Weevils, data=all.dat, perm=999)
ord.fit

#air determines how much space is between, <1 means it will let things overlap. tradeoff between neatness and accuracy.
#dropped overlapping text shown as a circle

#adding species codes (PE2,EE1, etc.) should help with overlapping points on ordination


#try to add arrows that indicate impact of a species on ordination
ordiarrows(all.mds,label=TRUE)

#site plot (no species data shown)
ordiplot(all.mds, type="n")
orditorp(all.mds,display="sites",col="black",air=1,cex=1)

ordiplot(all.mds, type="n")
ordiellipse(all.mds, site, label=T,air=0.01,cex=0.5)

# actual statistical test to for urban vs. rural has different community structure
ord.fit.all <- envfit(all.mds ~ Gradient, data=all.dat, perm=999)
ord.fit.all


## Amino GLM ####

#total amino acid content appears to be increased when weevils come late (second)
amino.dat.2 <- amino.dat
amino.dat.2$Aminosum <- amino.dat$Asp +
  amino.dat$Ser.Gln +
  amino.dat$Glu +
  amino.dat$Gly +
  amino.dat$His +
  amino.dat$Hser +
  amino.dat$Arg +
  amino.dat$Thr +
  amino.dat$Ala +
  amino.dat$Pro +
  amino.dat$Cys +
  amino.dat$Val +
  amino.dat$Lys +
  amino.dat$Ile +
  amino.dat$Leu +
  amino.dat$Phe

total.amino.glm <- glm(Aminosum ~ Aphids*Weevils, data=amino.dat.2)
summary(total.amino.glm)
Anova(total.amino.glm)

total.cld <- cld(emmeans(total.amino.glm, ~ Weevils|Aphids), 
                 adjust="tukey", type="response", Letters=c("baa"), sort=FALSE)
total.cld$.group=gsub(" ", "", total.cld$.group)

#raw pairs

total.cld$.group <- c("b","ab","a","b","b","b")

# 2020 Fig 5 #######
totals.fig.nap <- ggplot(total.cld, aes(x=factor(Weevils,levels=c("None","First","Second")), 
                                        y=emmean)) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=emmean-(SE), ymax=emmean+(SE)), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Amino acid concentration (nmol/mg DW)", x="Weevil Treatment") + 
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  geom_text(aes(x = Weevils, y = (emmean+SE+20), label = .group), 
            position = position_dodge(width = 0.8)) +
          facet_wrap( ~ Aphids, nrow=1)
totals.fig.nap


# No aphid facet_wrap
totals.fig <- ggplot(weevil.amino.cld, aes(x=factor(Weevils,levels=c("None","First","Second")), y=emmean, fill=factor(Aphids,levels=c("Sham", "Infective")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=emmean-(SE), ymax=emmean+(SE)), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Amino acid concentration (nmol/mg DW)", x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))
totals.fig

# Old weevil total amino acid anlayses

# the "totals" figure above did not partition the posthoc tests by aphid treatment, so the tukey correction (or lackthereof) is sketchy
# the figure below partitions into both aphid treatments, then shows some nice tukey letters

weevil.amino.cld <- cld(emmeans(total.amino.glm, ~ Weevils|Aphids), adjust="tukey", type="response",sort=FALSE)
weevil.amino.cld

##fix group labels so posthoc test looks nice
weevil.amino.cld$.group=gsub(" ", "", weevil.amino.cld$.group)

weevil.amino.fig <- ggplot(weevil.amino.cld, aes(x=Weevils, y=emmean, label=.group)) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=emmean-(SE), ymax=emmean+(SE)), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Total Amino Acids", x="Weevil Treatment") + 
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  geom_text(aes(x = Weevils, y = (emmean+SE+30), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Aphids)

weevil.amino.fig




# 2020 Hormone ################ 
## JA/SA Data ####
hormone.raw <- read.csv("master blaster hormones.csv", header=TRUE)

# Modify data to only use data for 2x3 design
hormone.dat <- subset(hormone.raw, Drop.for.2x3 == "no")
str(hormone.dat)

# this is not a delta-ct method, so no calcs required, just log-transformed hormone conc
hormone.mod <- manova(cbind(Log.SA,Log.JA,Log.ABA) ~ Aphid.Treatment*Weevil.Treatment, data=hormone.dat)
summary(hormone.mod)
summary.aov(hormone.mod)

# estimated marginal means for Weevil*Aphid
hormone.lsm <- cld(emmeans(hormone.mod, ~Weevil.Treatment*Aphid.Treatment|Hormone, mult.name = "Hormone"), Letters=c("abcd"), sort=FALSE)
hormone.lsm$.group=gsub(" ", "", hormone.lsm$.group)


#quickly plotting cld output with dispersion (conf interval)
plot(hormone.lsm, by = "disp")

#convert emmeans object to dataframe for transformation (delta delta ct and 2^-x)
hormone.lsm.dat <- as.data.frame(hormone.lsm)
hormone.lsm.dat

#change names of facets

facet_names <- c(
  `Log.SA` = "Salicylic Acid",
  `Log.JA` = "Jasmonic Acid",
  `Log.ABA` = "Abscisic acid"
)

#make a manual set tukey letters
# sketchy and will break if analysis changes

# 2020 Fig 4 #####

# I had to set these manually because ggplot is being impossible, be sure to update if analyses changes
hormone.lsm$poop <- c("a","a","b","c","c","c",     "ab", "ab","ab","b", "a", "ab",        "ab","ab","b","a","ab","a")
hormone.lsm

# Hormone figure from MANOVA
hormones.fig <- ggplot(hormone.lsm, aes(x=factor(Weevil.Treatment,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Treatment,levels=c("Sham", "Infective")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=emmean-(SE), ymax=emmean+(SE)), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Log of Hormone Content", x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  geom_text(aes(x = Weevil.Treatment, y = (emmean+SE+0.15), label = poop), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Hormone, labeller=as_labeller(facet_names))
hormones.fig

