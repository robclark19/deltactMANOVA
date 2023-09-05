# packages ############
library(dplyr)
library(emmeans)
library(multcomp)
library(lme4)
library(vegan)
library(ggplot2)
library(ggpubr)
library(car)

# Gene Data ####
saumik1 <- read.csv("Master Blaster Data.csv", header= TRUE)
str(saumik1)

#Modify data to only use data for 2x3 design

saumik4 <- subset(saumik1, Drop.for.2x3 == "No")

# Crit Notes ####
# two important notes: 

# 1 we are pooling by technical replicate, and not including technical rep as a random effect
# 2 we are using "raw" delta ct (gene of interest and beta-tubulin) values NOT raw ct values
# 1 and 2 are valid options given the data is there

# get summary statistics (mean) across all technical replicates
exp.dat2 <- saumik4 %>%
  group_by(Weevil.Timing,Aphid.Type,BioSample) %>%
  summarize_at(.vars = vars(Chitenase:lectin), # everything between Chitenase and lectin
               .funs = mean, na.rm = TRUE)


# Full MANOVA model with raw delta ct values of each gene predicted by treatment
# Gene-MANOVA ###########
genes.manova <- manova(cbind(Chitenase,OPR1,PR1,PsOX.11,GA,AOX.3,ACS2,LOX.2,ICS1,defensin,UBQ,ChlAbBp,lectin) ~ Aphid.Type*Weevil.Timing, data=exp.dat2)

# a roughly normal delta ct value would be nice i.e.
hist(exp.dat2$UBQ)
hist(exp.dat2$OPR1)

# MANOVA table (response of all genes with pillai index)
summary(genes.manova)
#AOV TABLE (response of each gene with a separate f-statistic)
manova.table <- summary.aov(genes.manova)
manova.table

##  Model estimates ####
# all genes based on aphids and weevils
# then present as compact letter display for posthoc tests

genes.lsm <- emmeans(genes.manova, ~ Aphid.Type*Weevil.Timing|Gene, mult.name = "Gene", adjust="none", response=TRUE)
genes.lsm

genes.cld <- cld(genes.lsm, sort = FALSE, adjust="none", response=TRUE, Letters=c("abcdef"))
genes.cld$.group=gsub(" ", "", genes.cld$.group)
genes.cld

# make cld () parameter estimates into a data.table so dplyer can manipulate
# for very small values, I used / options(scipen = 999)...
# this is a command in console to disable scientific notation from displaying
# you might have to do this with bonkers low gene expression levels (10^-8 or something silly)

genes.table <- as.data.frame(genes.cld)


## Set Control Treatment #### 
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
  mutate(emmean = cemmean - emmean,
         upper.SE = 2^(-(emmean-SE)),
         lower.SE = 2^(-(emmean+SE)), #upper and lower bound of model parameter (SE of the regression coefficient)
         emmean = 2^(-emmean),
         SE = 2^(-SE), # don't use this for plotting unless you want a migraine
         lower.CL = 2^(-lower.CL), 
         upper.CL = 2^(-upper.CL)) #95% confidence limits


## 2020 Fig 2 ####

# grab just the genes of interest from the larger MANOVA

f2.table <- subset(final.table, Gene == "ICS1" | Gene == "LOX.2" | Gene == "AOX.3" | Gene == "GA")

# make a labeling object to use later in ggplot
labels.f2 <- c(ICS1 = "ICS1", LOX.2 = "LOX2", AOX.3 = "AO3", GA = "GA2ox")

# use dplyer to reorder so the facets group correctly
neworder <- c("ICS1","AOX.3","LOX.2","GA")
f2.table <- arrange(transform(f2.table,
                              Gene=factor(Gene,levels=neworder)),Gene)

# make your plot!
f2a.fig<- ggplot(f2.table, aes(x=factor(Weevil.Timing,levels=c("None","First","Second")), y=emmean, fill=factor(Aphid.Type,levels=c("Sham", "Infective")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=lower.SE, ymax=upper.SE), position=position_dodge(0.8), width=0.5) +
  #make sure to use upper and lower se, not raw SE!
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  #this formatting allows me to insert greek symbols
  labs(y=expression("Fold Change in Expression ( "~ Delta*Delta~"CT)"), x="Weevil Treatment") + 
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
