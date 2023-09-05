library("lme4")
library("car")
library("multcomp")
library("ggplot2")
library("emmeans")
library("piecewiseSEM")
library("MASS")
library("dplyr")
library("effects")

################### Data and resulting path model #############################

#setwd("C:/Users/Camponotus/Desktop/pauls proc b paper")

#pooled data for 9 plants surrounding source plant
paul_data <- read.csv("Peaches2.csv", header = TRUE)
str(paul_data)

#run primary GLMs

######## Establishment #######
#factors that affect the number of adults that established

# This first link was never significant so it is not included in further models
link1 <- lm(Adult.established ~ Herbivory, data=paul_data)
summary(link1)

# what is distribtion of adult.established?
hist(paul_data$Adult.established)
# its a wierd quasipoisson, i'm not sure we should rely on this regression above

# attempting with mixed model
link1.mix <- lmer(Adult.established ~ Herbivory + (1|Block), data=paul_data)
summary(link1.mix)
#the variance component for block is zero (not surprising) so we can drop it from analyses
anova(link1.mix,link1)
#delta aic is so low the p-value for model comparison is 1; drop block


# attempting with negbin model
link1.nb <- glm.nb(Adult.established ~ Herbivory, data=paul_data)
summary(link1.nb)
#model does not converge AND isn't even remotely signfiicant if you force it

######### Nymph/Adult ############
# modeling nymph abundance with adults as offset
link.2.offset <- glm(Tot.nyms ~ Herbivory, offset=Adult.off, family=poisson, data=paul_data)
summary(link.2.offset)


#link 2_2 removes adult_established
link2_2 <- lm(Nym.adult ~ Herbivory, data=paul_data)
summary(link2_2)

# Modeling total nymphs as a function of herbivory and adult established...link 2_4 removes adult_established
link2_3 <- lm(Tot.nyms ~ Herbivory + Adult.established,data=paul_data)
summary(link2_3)

link2_4 <- lm(Tot.nyms ~ Herbivory,data=paul_data)
summary(link2_4)

# Creating links for movement variables, as a function of nymphs per adult and herbivory
# ratio is adult.off / adult.center
# ratio2 is adult.off / adult.lived

#these calculations work for getting very similar path models, but the offset/weighted binomial is superior
link3 <- lm(Ratio2 ~ Nym.adult+Herbivory,data=paul_data)
link3_2 <- lm(Ratio ~ Nym.adult+Herbivory,data=paul_data)
link3_3 <- lm(Adult.ratio ~ Nym.adult+Herbivory,data=paul_data)
link3_4 <- lm(Adult.off ~ Nym.adult+Herbivory,data=paul_data)
link3_5 <- lm(Infest ~ Nym.adult+Herbivory,data=paul_data)

#make a few vectors for frequency calculations for weighted binomial
paul_data$adult.totes <- paul_data$Adult.off + paul_data$Adult.lived
paul_data$ratio3 <- paul_data$Adult.off / paul_data$adult.totes

#run model that is a frequency of adults off weighted by the number of adults observed both on and off
#is now a binomial model rather than the normal model
link3.wt <- glm(ratio3 ~ Tot.nyms, data=paul_data, family=binomial, weights=adult.totes)
summary(link3.wt)

# Same models as above with herbivory effect removed (not significant in final path model)
link32 <- lm(Ratio2 ~ Nym.adult,data=paul_data)
link3_22 <- lm(Ratio ~ Nym.adult,data=paul_data)
link3_32 <- lm(Adult.ratio ~ Nym.adult,data=paul_data)
link3_42 <- lm(Adult.off ~ Nym.adult,data=paul_data)
link3_52 <- lm(Infest ~ Nym.adult,data=paul_data)

# Summary of these models to check output (not used in final path analysis)
summary(link3)
summary(link3_2)
summary(link3_3)
summary(link3_4)
summary(link3_5)

summary(link32)
summary(link3_22)
summary(link3_32)
summary(link3_42)
summary(link3_52)

# These models are the same as links 3 to 3_5, except nymphs per adult replaced by total nymphs
# not used in final path analysis

link3_6 <- lm(Ratio2 ~ Tot.nyms+Herbivory,data=paul_data)
link3_7 <- lm(Ratio ~ Tot.nyms+Herbivory,data=paul_data)
link3_8 <- lm(Adult.ratio ~ Tot.nyms+Herbivory,data=paul_data)
link3_9 <- lm(Adult.off ~ Tot.nyms+Herbivory,data=paul_data)
link3_10 <- lm(Infest ~ Tot.nyms+Herbivory,data=paul_data)

link3_62 <- lm(Ratio2 ~ Tot.nyms,data=paul_data)
link3_72 <- lm(Ratio ~ Tot.nyms,data=paul_data)
link3_82 <- lm(Adult.ratio ~ Tot.nyms,data=paul_data)
link3_92 <- lm(Adult.off ~ Tot.nyms,data=paul_data)
link3_102 <- lm(Infest ~ Tot.nyms,data=paul_data)

# Summary of these models
# not used in final path analysis
summary(link3_6)
summary(link3_7)
summary(link3_8)
summary(link3_9)
summary(link3_10)

summary(link3_62)
summary(link3_72)
summary(link3_82)
summary(link3_92)
summary(link3_102)

#  models where total plants infested is a function of adult movement and herbivory
# Change Plants.inf2 to Plants.inf if you want to model only the ring (not all 9 plants)

link4 <- glm(cbind(Plants.inf2,Plants.healthy)~Ratio2+Herbivory,family="binomial",data=paul_data)
link4_2 <- glm(cbind(Plants.inf2,Plants.healthy)~Ratio+Herbivory,family="binomial",data=paul_data)
link4_3 <- glm(cbind(Plants.inf2,Plants.healthy)~Adult.ratio+Herbivory,family="binomial",data=paul_data)
link4_4 <- glm(cbind(Plants.inf2,Plants.healthy)~Adult.off+Herbivory,family="binomial",data=paul_data)
link4_5 <- glm(cbind(Plants.inf2,Plants.healthy)~Infest+Herbivory,family="binomial",data=paul_data)

# from interwebs: Note that when using the frequency form of a binomial glm, 
# you should supply the number of observations per trial in the weights argument.
# It would look like: glm(events/n ~ x, data=*, weights=n, ...)

#calculate ratio of infected based on number of "trials"
paul_data$inf.ratio <- paul_data$Plants.inf2 / 9
paul_data$inf.ratio

#make a vector that is the number of plants for each trial (in case this may have varied due to plant death)
paul_data$w.nine <- paul_data$Plants.inf2 + paul_data$Plants.healthy

#run new model without cbind statement (cbind confuses SEM and parameter estimation)
link4.w <- glm(inf.ratio ~ ratio3 + Herbivory, family="binomial", data=paul_data, weights=w.nine)
summary(link4.w)
#model output should be similar, if not exactly equivalent, to the cbind ratio model
summary(link4)

# Alternative models where total plants infected only affected by nymphs (these are a poor fit and dropped from analyses)
link4_6 <- glm(cbind(Plants.inf2,Plants.healthy)~Nym.adult,family="binomial",data=paul_data)
link4_7 <- glm(cbind(Plants.inf2,Plants.healthy)~Tot.nyms,family="binomial",data=paul_data)

# Summary of models (not used in path model)
summary(link4)
summary(link4_2)
summary(link4_3)
summary(link4_4)
summary(link4_5)
summary(link4_6)
summary(link4_7)

# Organize linear models into lists for series of equations using piecewiseSEM

# The five lists here involve different variables for adult movement. Right now these equations use
# link 2_2 and corresponding link 3(2) models, meaning that the link between adult establishment and
# nymphs, and the effects of herbivory on movement are not included

# This proved to be a better fit than models where nymphs were affected by adult establishment
# and where moved adults was a direct function of herbivory

paul.list <- list(link2_2,link32,link4)
paul.list2 <- list(link2_2,link3_22,link4_2)
paul.list3 <- list(link2_2,link3_32,link4_3)
paul.list4 <- list(link2_2,link3_42,link4_4)
paul.list5 <- list(link2_2, link3_52, link4_5)

#new list with updated binomial data for third model
list1.psem <- psem(link2_2,link32,link4.w, data=paul_data)
summary(list1.psem)
coefs(list1.psem)
dSep(list1.psem)

#new list with all updated models with weights and offsets using poisson glms
rob.psem <- psem(link.2.offset,link4.w,link3.wt, data=paul_data)
summary(rob.psem)

# Model where just nymphs affect infection (we assume this is not correct because nymphs not important)
# Nymphs are affected by herbivory

paul.list6 <- list(link2_2, link4_6)

# Same as list1 to 6, except total nymphs is the variable for fitness, not nymphs per adult

paul.list7 <- list(link2_4,link3_62,link4)
paul.list8 <- list(link2_4,link3_72,link4_2)
paul.list9 <- list(link2_4,link3_82,link4_3)
paul.list10 <- list(link2_4,link3_92,link4_4)
paul.list11 <- list(link2_4, link3_102, link4_5)

paul.list12 <- list(link2_4, link4_7)

#confirmatory path analysis pt1: are our path model valid?
sem.fit(paul.list, data=paul_data) # YES
sem.fit(paul.list2, data=paul_data) # YES
sem.fit(paul.list3, data=paul_data) # YES
sem.fit(paul.list4, data=paul_data) # YES
sem.fit(paul.list5, data=paul_data) # YES
sem.fit(paul.list6, data=paul_data) # NO

sem.fit(paul.list7, data=paul_data) # YES
sem.fit(paul.list8, data=paul_data) # YES
sem.fit(paul.list9, data=paul_data) # YES
sem.fit(paul.list10, data=paul_data) # YES
sem.fit(paul.list11, data=paul_data) # YES
sem.fit(paul.list12, data=paul_data) # NO

#confirmatory path analysis pt2: what is the contribution of each predictor?
#old code. sem.coefs is no longer updated in piecewiseSEM package
sem.coefs(paul.list, data=paul_data)
sem.coefs(paul.list2, data=paul_data)
sem.coefs(paul.list3, data=paul_data)
sem.coefs(paul.list4, data=paul_data)
sem.coefs(paul.list5, data=paul_data)
sem.coefs(paul.list6, data=paul_data)

sem.coefs(paul.list7, data=paul_data)
sem.coefs(paul.list8, data=paul_data)
sem.coefs(paul.list9, data=paul_data)
sem.coefs(paul.list10, data=paul_data)
sem.coefs(paul.list11, data=paul_data)
sem.coefs(paul.list12, data=paul_data)

################ Expression data ###########
#import positional data with delta ct raw values
paul.pemv <- read.csv("paul proc b.csv")
str(paul.pemv)

#pool by technical replicate with mean value
paul.pemv.2 <- paul.pemv %>%
  group_by(location,aphid.treatment,bio.rep) %>% #treatment combinations you want to save
  summarize_at(.vars = vars(pemv.1:pemv.2), # everything for pemv 1-2
               .funs = mean, na.rm = TRUE)

paul.pemv.2

#run univariate model with position * aphid treatment interaction to make sure everything chekcs out
ap.mod <- glm(pemv.1 ~ aphid.treatment*location, data=paul.pemv.2)
summary(ap.mod)
Anova(ap.mod)


#run manova model with pemv1 and pemv2
pemv.manova <- manova(cbind(pemv.1,pemv.2) ~ location*aphid.treatment, data=paul.pemv.2)
summary.aov(pemv.manova)


#get parameter estimates at each treatment, then present as compact letter display for posthoc tests
pemv.manova.lsm <- emmeans(pemv.manova, ~ location*aphid.treatment|virus, mult.name = "virus", adjust="none", response=TRUE)
pemv.manova.lsm

pemv.manova.cld <- cld(pemv.manova.lsm, sort = FALSE, adjust="none", response=TRUE, Letters=c("abcdef"))
pemv.manova.cld$.group=gsub(" ", "", pemv.manova.cld$.group)
pemv.manova.cld

# make a table where transformations are done on parameter estimates for upper and lower bounds of SE and mean
pemv.manova.table <- as.data.frame(pemv.manova.cld) %>% 
  mutate(upper.SE = 2^(-(emmean-SE)),
         lower.SE = 2^(-(emmean+SE)),
         emmean = 2^(-emmean)) # will operate with new emmean column
pemv.manova.table


#make preliminary figure in ggplot2

pemv.manova.fig <- ggplot(pemv.manova.table, aes(x=factor(aphid.treatment,levels=c("Control","Infective","Weevil and Infective")), 
                                                 y=emmean, fill=factor(location,levels=c("Top", "Middle","Bottom")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=emmean-(SE), ymax=emmean+(SE)), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Fold Change in PEMV expression", x="Aphid Treatment") + 
  labs(fill="Location on Plant") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  #scale_y_log10()  +
  #geom_text(aes(x = aphid.treatment, y = (emmean+SE+10), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ virus)
pemv.manova.fig






#defensive genes model (has 4 treatment levels)
pr.genes <- read.csv("pr genes all.csv")
str(pr.genes)

pr.genes.2 <- pr.genes %>%
  group_by(location,aphid.treatment,bio.rep) %>% #treatment combinations you want to save
  summarize_at(.vars = vars(pr.1:opr.1), # everything for both pr genes
               .funs = mean, na.rm = TRUE)
pr.genes.2

pr.genes.2$pr.1.trans <- 2^-(pr.genes.2$pr.1)
pr.genes.2$opr.1.trans <- 2^-(pr.genes.2$opr.1)


opr.mod <- manova(cbind(pr.1.trans,opr.1.trans) ~ location*aphid.treatment, data=pr.genes.2)
summary.aov(opr.mod)
summary(opr.mod)

hist(pr.genes.2$pr.1.trans)




#get parameter estimates at each treatment, then present as compact letter display for posthoc tests
opr.mod.lsm <- emmeans(opr.mod, ~ location|aphid.treatment|virus, mult.name = "virus", adjust="none", response=TRUE)
opr.mod.lsm

opr.mod.cld <- cld(opr.mod.lsm, sort = FALSE, adjust="none", response=TRUE, Letters=c("abcdef"))
opr.mod.cld$.group=gsub(" ", "", opr.mod.cld$.group)
opr.mod.cld

# make csv for messing with parameter estimates
#write.csv(as.data.frame(opr.mod.cld), file = "paul pemv para.csv")

# make a table where transformations are done on parameter estimates for upper and lower bounds of SE and mean
#opr.mod.table <- as.data.frame(opr.mod.cld) %>% 
 # mutate(upper.SE = 2^(-(emmean-SE)),
      #   lower.SE = 2^(-(emmean+SE)),
       #  emmean = 2^(-emmean)) # will operate with new emmean column
#opr.mod.table

#to make a control for reference, I used the mean expression level in the control treatment for all locations
#this should give the baseline expression level as an effective comparison

# make preliminary figure in ggplot2
# i dont think this will work without a designated 1=control setting. 
# 
opr.mod.cld <- as.data.frame(opr.mod.cld)
opr.mod.cld <- subset(opr.mod.cld, aphid.treatment != "Control")

labels <- c(pr.1.trans = "PR1", opr.1.trans = "OPR1")
#labels <- c("italic('PR1')", "italic('OPR1')")

# original figure
opr.mod.fig <- ggplot(opr.mod.cld, aes(x=factor(aphid.treatment), 
                                                 y=emmean, fill=factor(location,levels=c("Top", "Middle","Bottom")))) +
  geom_bar(stat="identity", width=0.8, position="dodge") +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), position=position_dodge(0.8), width=0.5) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Relative expression", x="Treatment") + 
  labs(fill="Location on Plant") +
  scale_fill_grey() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  #scale_y_log10()  +
  geom_text(aes(x = aphid.treatment, y = (emmean+SE+1), label = .group), position = position_dodge(width = 0.8)) +
  facet_wrap( ~ virus, labeller=labeller(virus = labels)) +
  theme(strip.text = element_text(face = "italic"))
opr.mod.fig

# updated figure July 29 with "Weevils only" now "Sham and Weevils"

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




###### new methods for expression data MANOVA?#######
# will do later if figures above didn't come out well


############niche part #############

#pooled data for 9 plants surrounding source plant
niche.dat <- read.csv("Weevilaphidnichepartitioning.csv", header = TRUE)
niche.dat

niche.dat$prop.notch <- niche.dat$Notches_at_node/niche.dat$Total_notches_on_plant

glm_1 <- glm(prop.notch ~ Node, weights=Total_notches_on_plant, family=binomial(link = probit), data=niche.dat)
Anova(glm_1)
summary(glm_1)

niche.dat$aphid.prop <- niche.dat$Total_aphids_at_node/niche.dat$Total_aphids_on_plant

glm_2 <- glm(aphid.prop ~ Notches_at_node, family=binomial(link=probit), weights = Total_aphids_on_plant, data=niche.dat)
summary(glm_2)
Anova(glm_2)


glm_3 <- glm(aphid.prop ~ Node*Weevil, weights=Total_aphids_on_plant, data=niche.dat, family=binomial(link=probit))
summary(glm_3)

#figure for glm1
# plot using effects package (nicer)
interaction.raw.1 <- effect("Node", glm_1, se=TRUE, xlevels=20)
# Data Frame (for ggplot2)
interaction.dat.1 <-as.data.frame(interaction.raw.1)
#Create plot
glm_1_plot <- ggplot(data=interaction.dat.1, aes(x=Node, y=fit)) +
  theme_bw() +
  geom_line(size=2) +
  #scale_color_manual(values=c("Red", "Blue"),labels=c("Present","Absent")) +
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se),alpha=.2) +
  #scale_fill_manual(values=c("Red", "Blue"),labels=c("Present","Absent")) +
  theme(text = element_text(size=16),
        legend.text = element_text(size=16),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="top") +
  #geom_point(data=pooled.dat, aes(x=(log.cats), y=post_treat, colour=ant_treatment), position=position_jitter()) +
  #facet_wrap( ~ sap_presence, nrow=2, labeller = as_labeller(sap.names)) +
  labs(x="Height (Node Number)", y="Proportion weevil notches on plant")
glm_1_plot 


#plot for glm2
# plot using effects package (nicer)
interaction.raw.2 <- effect("Notches_at_node", glm_2, se=TRUE, xlevels=20)
# Data Frame (for ggplot2)
interaction.dat.2 <-as.data.frame(interaction.raw.2)
#Create plot
glm_2_plot <- ggplot(data=interaction.dat.2, aes(x=Notches_at_node, y=fit)) +
  theme_bw() +
  geom_line(size=2) +
  #scale_color_manual(values=c("Red", "Blue"),labels=c("Present","Absent")) +
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se),alpha=.2) +
  #scale_fill_manual(values=c("Red", "Blue"),labels=c("Present","Absent")) +
  theme(text = element_text(size=16),
        legend.text = element_text(size=16),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="top") +
  #geom_point(data=pooled.dat, aes(x=(log.cats), y=post_treat, colour=ant_treatment), position=position_jitter()) +
  #facet_wrap( ~ sap_presence, nrow=2, labeller = as_labeller(sap.names)) +
  labs(x="Notches at node", y="Proportion of aphids")
glm_2_plot 






#figure for glm3
# plot using effects package (nicer)
interaction.raw.4 <- effect("Node*Weevil", glm_3, se=TRUE, xlevels=20)
interaction.raw.4 

# Data Frame
interaction.dat.4 <-as.data.frame(interaction.raw.4)


#Create plot
money.plot <- ggplot(data=interaction.dat.4, aes(x=Node, y=fit, group=Weevil)) +
  theme_bw() +
  geom_line(size=2, aes(color=Weevil)) +
  scale_color_manual(values=c("Red", "Blue"),labels=c("Present","Absent")) +
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se,fill=Weevil),alpha=.2) +
  scale_fill_manual(values=c("Red", "Blue"),labels=c("Present","Absent")) +
  theme(text = element_text(size=16),
        legend.text = element_text(size=16),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="top") +
  #geom_point(data=pooled.dat, aes(x=(log.cats), y=post_treat, colour=ant_treatment), position=position_jitter()) +
  #facet_wrap( ~ sap_presence, nrow=2, labeller = as_labeller(sap.names)) +
  labs(x="Node Number", y="Proportion of aphids at node", fill="Weevils", color="Weevils")

money.plot


