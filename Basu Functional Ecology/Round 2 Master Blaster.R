########### packages ############
library(dplyr)
library(emmeans)
library(multcomp)
library(lme4)
library(vegan)
library(ggplot2)
library(ggpubr)
library(car)

# Set Working Directory
# setwd("C:/Users/Camponotus/Desktop/Master Blaster")

#Load raw data and check structure
saumik1 <- read.csv("Master Blaster Data.csv", header= TRUE)
str(saumik1)

#Modify data to only use data for 2x3 design
############## data for 2x3 ###############
saumik4 <- subset(saumik1, Drop.for.2x3 == "No")

#get summary statistics (mean) across all technical replicates
#pool by technical replicate
exp.dat2 <- saumik4 %>%
  group_by(Weevil.Timing,Aphid.Type,BioSample) %>%
  summarize_at(.vars = vars(Chitenase:lectin), # everything between Chitenase and lectin
               .funs = mean, na.rm = TRUE)



str(exp.dat2)
#full manova model with raw delta ct of each gene predicted by treatment (or cycles)
############ all gene MANOVA ###########
genes.manova <- manova(cbind(Chitenase,OPR1,PR1,PsOX.11,GA,AOX.3,ACS2,LOX.2,ICS1,defensin,UBQ,ChlAbBp,lectin) ~ Aphid.Type*Weevil.Timing, data=exp.dat2)

#MANOVA TABLE (response of all genes)
summary(genes.manova)
#AOV TABLE (response of each gene)
manova.table <- summary.aov(genes.manova)
manova.table

#get parameter estimates for all genes based on aphids and weevils, then present as compact letter display for posthoc tests

genes.lsm <- emmeans(genes.manova, ~ Aphid.Type*Weevil.Timing|Gene, mult.name = "Gene", adjust="none", response=TRUE)
genes.lsm

genes.cld <- cld(genes.lsm, sort = FALSE, adjust="none", response=TRUE, Letters=c("abcdef"))
genes.cld$.group=gsub(" ", "", genes.cld$.group)
genes.cld

#make parameter estimates into a data.table
#for very small values, i used / options(scipen = 999) / command in console to disable scientific notation

genes.table <- as.data.frame(genes.cld)
genes.table

# write output to a csv for supplementary materials
# write.csv(genes.table, file = "all genes emmeans.csv")

#subtract each control (sham,none) within each gene from the values. control should = 1 in each gene.
expression.table <- genes.table

expression.table <- genes.table %>%
  group_by(Gene) %>% # each gene considered a batch
  mutate(control = expression.table[which(Aphid.Type == "Sham" & Weevil.Timing == "None"), "emmean"],
         emmean = 2^(-emmean), # will operate with new emmean column
         SE = 2^(-SE))

#make new names and columes to get the control for comparison
new.data <- genes.table %>%
  filter(Aphid.Type == "Sham", Weevil.Timing == "None")

colnames(new.data) = c("Aphid.Type", "Weevil.Timing", "Gene", "cemmean","cSE","cdf","clower.CL","cupper.CL","c.group")
new.data = new.data %>% dplyr::select(Gene, cemmean)

#left join table now adds control as a new column for caculating the second delta in delta delta ct
total.data <- left_join(genes.table, new.data, by=c("Gene"))
head(total.data)

#do delta delta ct with left_joined table 

final.table <- total.data %>% 
  mutate(emmean = cemmean - emmean,
         upper.SE = 2^(-(emmean-SE)),
         lower.SE = 2^(-(emmean+SE)),
         emmean = 2^(-emmean), # will operate with new emmean column
         SE = 2^(-SE)
  )
head(final.table)

#write.csv(final.table, file = "all genes emmeans.csv")

#transform mean and SE into 2^-(delta delta ct) for all values in the emmean output

#plot figures by categories determined ahead of time (sa related genes, ja related genes, defensin and others)

#### SA group plot ########

sa.table <- subset(final.table, Gene == "PR1" | Gene == "Chitenase" | Gene == "ICS1")
# sa.table$LSD <- as.numeric(sa.table$.group)
as.data.frame(sa.table)

labels.sa <- c(Chitenase = "Chitenase", PR1 = "PR1", ICS1 = "ICS1")

sa.fig <- ggplot(sa.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
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
  scale_y_continuous(breaks=c(0,1, 10, 25)) +
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+1), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.sa)) +
  theme(strip.text = element_text(face = "italic"))
sa.fig

### ja group plot #####

ja.table <- subset(final.table, Gene == "OPR1" | Gene == "LOX.2")
as.data.frame(ja.table)

labels.ja <- c(OPR1 = "OPR1", LOX.2 = "LOX 2")

ja.fig <- ggplot(ja.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
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
  scale_y_continuous(breaks=c(0,1, 10, 25)) +
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+1), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.ja)) +
  theme(strip.text = element_text(face = "italic"))
ja.fig

######### defensin group ########

def.table <- subset(final.table, Gene == "defensin" | Gene == "ACS2")
as.data.frame(def.table)

labels.def <- c(defensin = "Defensin", ACS2 = "ACS2")

def.fig <- ggplot(def.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
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
  scale_y_continuous(breaks=c(0,1, 10, 25)) +
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+2), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.def)) +
  theme(strip.text = element_text(face = "italic"))
def.fig

# "others" plot ##########
#Chitenase,OPR1,PR1,PsOX.11,GA,AOX.3,ACS2,LOX.2,ICS1,defensin,UBQ,ChlAbBp,lectin

others.table <- subset(final.table, Gene == "GA" | Gene == "AOX.3" | Gene == "UBQ" | Gene == "ChlAbBp" | Gene == "lectin")
as.data.frame(others.table)

labels.others <- c(GA = "GA", AOX.3 = "AOX 3", UBQ = "UBC4", ChlAbBp = "ChlAbBp", lectin = "Lectin")

others.fig <- ggplot(others.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
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
  scale_y_continuous(breaks=c(0,1, 10, 25)) +
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+1), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.others)) +
  theme(strip.text = element_text(face = "italic"))
others.fig

#psOx.11 plot ##########
ps.table <- subset(final.table, Gene == "PsOX.11")
as.data.frame(ps.table)

labels.ps <- c(PsOX.11 = "PsOX 11")

ps.fig <- ggplot(ps.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
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
  scale_y_continuous(breaks=c(0,1, 10, 25, 50)) +
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+5), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.ps)) +
  theme(strip.text = element_text(face = "italic"))
ps.fig

############ july 2019 12 panel gene plot ########
#this is the updated layout based on manuscript revisions with dave july 18 2019
#row 1 is chitenase, pr1, defensin
#row 2 is asc2, opr1, psox11
#row 3 is GA, lectin, lox2
#row 4 is ICS1, ChalBP, aox3

######## R1: Chit,PR1,def ############
str(final.table)


R1.table <- subset(final.table, Gene == "Chitenase" | Gene == "PR1" | Gene == "defensin")
as.data.frame(R1.table)

labels.R1 <- c(Chitenase = "Chitenase", PR1 = "PR1", defensin = "Defensin")

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
  #theme(legend.position = c(0.2, 0.7), legend.text = element_text(size=8), legend.title = element_text(size=8)) + 
  #the above code will move legend to top right of figure
  theme(legend.position="none") + #this one is just easier, and is not inset
    geom_text(aes(x = Weevil.Timing, y = (emmean+SE+3), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.R1)) +
  theme(strip.text = element_text(face = "italic"))
r1.fig

#### R2 #########

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
  #scale_y_continuous(breaks=c(0,1, 10, 25)) +
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+3), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.R2)) +
  theme(strip.text = element_text(face = "italic"))
r2.fig

# R3 ####################


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


########### R4 #############


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

####### R-all together ##########
r1.fig
r2.fig
r3.fig
r4.fig

#arrange into one large panel using ggpubr package
figure.r.all <- ggarrange(r1.fig, r2.fig, r3.fig, r4.fig, labels = c("", "", "", ""), ncol = 1, nrow = 4)
figure.r.all

# Fig 2 final ####
# ics1, lox2, ao3, ga2ox


# fix y scale so they look better #####
 
#add a line to the ggplot of --- facet_wrap( ~ Gene, scales="free_y")



f2.table <- subset(final.table, Gene == "ICS1" | Gene == "LOX.2" | Gene == "AOX.3" | Gene == "GA")
labels.f2 <- c(ICS1 = "ICS1", LOX.2 = "LOX2", AOX.3 = "AO3", GA = "GA2ox")

# reorder everything in the original data to get the right values

neworder <- c("ICS1","AOX.3","LOX.2","GA")
f2.table <- arrange(transform(f2.table,
                          Gene=factor(Gene,levels=neworder)),Gene)


f2.fig <- ggplot(f2.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
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
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+0.8), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.f2), nrow=2) +
  theme(strip.text = element_text(face = "italic"))
f2.fig

# Fig 3 final ######

f3.table <- subset(final.table, Gene == "PR1" | Gene == "defensin" | Gene == "lectin")
labels.f3 <- c(PR1 = "PR1", defensin = "DDR230", lectin = "Lectin")

# reorder everything in the original data to get the right values

neworder <- c("PR1","defensin","lectin")
f3.table <- arrange(transform(f3.table,
                              Gene=factor(Gene,levels=neworder)),Gene)


f3.fig <- ggplot(f3.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
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
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+1.5), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.f3), nrow=1) +
  theme(strip.text = element_text(face = "italic"))
f3.fig


##### Redo all ###########

# D1 #####
D1.table <- subset(final.table, Gene == "Chitenase" | Gene == "PR1" | Gene == "defensin")


labels.D1 <- c(Chitenase = "Chitenase", PR1 = "PR1", defensin = "Defensin")


d1.fig <- ggplot(D1.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
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
  #theme(legend.position = c(0.2, 0.7), legend.text = element_text(size=8), legend.title = element_text(size=8)) + 
  #the above code will move legend to top right of figure
  theme(legend.position="none") + #this one is just easier, and is not inset
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+3), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.D1)) +
  theme(strip.text = element_text(face = "italic"))
d1.fig


# D2 #####
D2.table <- subset(final.table, Gene == "PsOX.11" | Gene == "GA" | Gene == "LOX.2")
D2.table

labels.D2 <- c(PsOX.11 = "POX11", GA = "GA2ox", LOX.2 = "LOX2")
labels.D2

d2.fig <- ggplot(D2.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
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
  #theme(legend.position = c(0.2, 0.7), legend.text = element_text(size=8), legend.title = element_text(size=8)) + 
  #the above code will move legend to top right of figure
  theme(legend.position="none") + #this one is just easier, and is not inset
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+3), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.D2)) +
  theme(strip.text = element_text(face = "italic"))
d2.fig


# D3 #####

D3.table <- subset(final.table, Gene == "lectin" | Gene == "AOX.3" | Gene == "ICS1")
D3.table

#the only way to fix this is to find a way to manually re-arrange the d3 table (they currently keep going in alphabetical order)

labels.D3 <- c(lectin = "Lectin", AOX.3 = "AO3", ICS1 = "ICS1")
labels.D3

d3.fig <- ggplot(D3.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
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
  #theme(legend.position = c(0.2, 0.7), legend.text = element_text(size=8), legend.title = element_text(size=8)) + 
  #the above code will move legend to top right of figure
  theme(legend.position="none") + #this one is just easier, and is not inset
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+0.75), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.D3)) +
  theme(strip.text = element_text(face = "italic"))
d3.fig



# i have no idea why its putting them in alphabetical order

# D4 #####

D4.table <- subset(final.table, Gene == "OPR1" | Gene == "ACS2" | Gene == "ChlAbBp")
D4.table

labels.D4 <- c(OPR1 = "OPR1",ACS2 = "ACS2",ChlAbBp = "Chla/bBP")
labels.D4

d4.fig <- ggplot(D4.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
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
  #theme(legend.position = c(0.2, 0.7), legend.text = element_text(size=8), legend.title = element_text(size=8)) + 
  #the above code will move legend to top right of figure
  theme(legend.position="none") + #this one is just easier, and is not inset
  geom_text(aes(x = Weevil.Timing, y = (emmean+SE+1), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ Gene, labeller=labeller(Gene = labels.D4)) +
  theme(strip.text = element_text(face = "italic"))
d4.fig






####### D-all together ##########
d1.fig
d2.fig
d3.fig
d4.fig

#arrange into one large panel using ggpubr package
figure.d.all <- ggarrange(d1.fig, d2.fig, d3.fig, d4.fig, labels = c("", "", "", ""), ncol = 1, nrow = 4)
figure.d.all










###################### junk MANOVA code below #######################

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

# tableHTML(defensin.cld)

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


############ Amino Acid NMDS ####################
# Fig 5b ####
amino.dat <- read.csv("plant aminos.csv", header=TRUE)
# response.dat <- amino.dat[c(6:21)] 
# predictor.dat <- amino.dat[c(2:3)]

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
ordihull(amino.mds, groups=Weevils, draw="polygon",col="grey90",label=T)

# ellipse might be better
ordiellipse(amino.mds,amino.dat$Weevils, display="species")

# new amino acid ########


all.dat <- read.csv("plant aminos.csv")
str(all.dat)

#make new object so commands below dont edit primary data
allmatrix.dat <- all.dat

#drop env columns to make a species matrix
allmatrix.dat$Aphids <- NULL
allmatrix.dat$Weevils<- NULL

colSums(allmatrix.dat)


#nmds full
#drop mostly zeroes (species with <5 observations for now)
# matrix.dat$ERA <- NULL
# matrix.dat$CHR <- NULL
# matrix.dat$HEL <- NULL
# matrix.dat$SP1 <- NULL
# $MEE <- NULL

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
# ordiellipse(amino.mds,amino.dat$Weevils, display="species")


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



###############amino acid MANOVA
# probably not enough replication for such a large array of amino acids (or these responses are simply not impacted seperately)

amino.manova <- manova(cbind(Asp.Asn,	Ser.Gln,	Glu,	Gly,	His,	Hser,	Arg,	Thr,	Ala,	Pro,	Cys,	Val,	Lys,	Ile,	Leu,	Phe) ~ Aphids*Weevils, data=amino.dat)
summary(amino.manova)
summary.aov(amino.manova)
#estimates
amino.lsm <- emmeans(amino.manova, ~ Aphids*Weevils|Amino.Acid, mult.name = "Amino.Acid", adjust="none")
amino.cld <- cld(amino.lsm, sort = FALSE, adjust="none")
amino.cld
as.data.frame(amino.cld)

######## total amino acids #######
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


#total amino acids figure
# fig 5a #######
totals.fig.nap <- ggplot(total.cld, aes(x=factor(Weevils,levels=c("None","First","Second")), 
                                        y=emmean)) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=emmean-(SE), ymax=emmean+(SE)), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Total Amino Acids", x="Weevil Treatment") + 
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  geom_text(aes(x = Weevils, y = (emmean+SE+20), label = .group), 
            position = position_dodge(width = 0.8)) +
          facet_wrap( ~ Aphids, nrow=1)
totals.fig.nap


# original figure with aphids included
totals.fig <- ggplot(total.cld, aes(x=factor(Weevils,levels=c("None","First","Second")), y=emmean, fill=factor(Aphids,levels=c("Sham", "Infective")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=emmean-(SE), ymax=emmean+(SE)), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Total Amino Acids", x="Weevil Treatment") + 
  labs(fill="Aphid Type") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))
totals.fig

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




############# hormone data ################ 

hormone.raw <- read.csv("master blaster hormones.csv", header=TRUE)

#Modify data to only use data for 2x3 design
hormone.dat <- subset(hormone.raw, Drop.for.2x3 == "no")
str(hormone.dat)

hormone.mod <- manova(cbind(Log.SA,Log.JA,Log.ABA) ~ Aphid.Treatment*Weevil.Treatment, data=hormone.dat)
summary(hormone.mod)
summary.aov(hormone.mod)

#fine estimated marginal means for Weevil*Aphid
hormone.lsm <- cld(emmeans(hormone.mod, ~Weevil.Treatment*Aphid.Treatment|Hormone, mult.name = "Hormone"), Letters=c("abcd"), sort=FALSE)
hormone.lsm

hormone.lsm$.group=gsub(" ", "", hormone.lsm$.group)
hormone.lsm

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

#CRITICAL point if hormone data CHANGE ####

# fig 4 #####
hormone.lsm$poop <- c("a","a","b","c","c","c",     "ab", "ab","ab","b", "a", "ab",        "ab","ab","b","a","ab","a")
hormone.lsm

#example figure
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



#facet_grid(hospital ~ ., labeller = as_labeller(hospital_names))