######################################################################
#
#  R code for Quantitative community ecology short course
#  ESA 2024, Long Beach, California
#  Jesse Miller, kawriver@gmail.com, August 2024
#  Partly based on materials created by Rob Smith for the 2023 workshop
#
##      GNU General Public License, Version 3.0    ###################


####    INSTALL DEPENDENCIES (uncomment these if running first time)   ####
# install.packages('vegan')     # community ecology
# install.packages('labdsv')    # community ecology
# install.packages('ade4')      # community ecology
# install.packages('ecodist')   # community ecology
# install.packages('fso')       # ordination
# install.packages('vegclust')  # clustering
# install.packages('ape')       # phylogenetics
# install.packages('picante')   # phylogenetics
# install.packages('mgcv')      # nonlinear regression
# install.packages('tidyr')		# data management
# install.packages('dplyr')		# data management
# install.packages('lme4')		# lienar mixed models
# install.packages('lmerTest')  # linear mixed models
# install.packages('effects')	# model effects


####    load packages     ####
require('vegan')
require('labdsv')
require('ade4')
require('ecodist')
require('fso')
require('vegclust')
require('ape')
require('picante')
require('mgcv')
require('tidyr')
require('dplyr')
require('lme4')
require('lmerTest')
require('effects')

### R

### file I/O

### load from RDA file (can contain multiple objects)
 load('./data/veg.rda')  # assumed to be in your working dir
xy  <- veg$xy    # spatial
spe <- veg$spe   # species
env <- veg$env   # environment
tra <- veg$tra   # traits
phy <- veg$phy   # phylogeny

###### data preparation ######
spe <- data.frame(log10(spe + 1))                           # transformation
env <- data.frame(decostand(scale(env, center=F), 'range')) # transformation
tra <- data.frame(decostand(tra, 'range'))                  # transformation
d   <- vegdist(spe, method='bray', binary=T)
D   <- stepacross(d, 'shortest', toolong=1, trace=F)
E   <- dist(xy)                         # spatial distance matrix
pc  <- pcnm(E, w=rowSums(spe)/sum(spe)) # PCNMs *weighted*  by abundances
cl  <- vegclustdist(D, mobileMemb=7, method='FCM', m=1.2, nstart=5)
grp <- defuzzify(cl, 'max')[[2]]
m1  <- metaMDS(D, k=2, try=200, trymax=500, trace=0)  # NMS
m2  <- cmdscale(D, k=2, add=T)                        # PCoA
m3  <- prcomp(spe)                                    # PCA
grp_k <- cut(env$k, quantile(env$k, probs=seq(0,1,by=0.25)), include.lowest=T,
             labels=c('lo','med','hi','veryhi'))
`get_palette` <- function() {
  pal <- c('#414487E6','#404688E6','#3F4889E6','#3E4989E6','#3E4C8AE6',
           '#3D4E8AE6','#3C508BE6','#3B528BE6','#3A548CE6','#39558CE6',
           '#38588CE6','#375A8CE6','#365C8DE6','#355E8DE6','#35608DE6',
           '#34618DE6','#33638DE6','#32658EE6','#31678EE6','#30698EE6',
           '#306A8EE6','#2F6C8EE6','#2E6E8EE6','#2D708EE6','#2C718EE6',
           '#2C738EE6','#2B748EE6','#2A768EE6','#2A788EE6','#297A8EE6',
           '#287C8EE6','#287D8EE6','#277F8EE6','#26818EE6','#26828EE6',
           '#25848EE6','#24868EE6','#24878EE6','#23898EE6','#228B8DE6',
           '#228D8DE6','#218F8DE6','#21908CE6','#20928CE6','#20938CE6',
           '#1F958BE6','#1F978BE6','#1F998AE6','#1F9A8AE6','#1E9C89E6',
           '#1F9E89E6','#1FA088E6','#1FA187E6','#20A386E6','#20A486E6',
           '#21A685E6','#22A884E6','#24AA83E6','#25AC82E6','#26AD81E6',
           '#28AE80E6','#2AB07FE6','#2DB27DE6','#2FB47CE6','#32B67AE6',
           '#34B679E6','#37B878E6','#3ABA76E6','#3DBC74E6','#40BD72E6',
           '#43BF71E6','#47C06FE6','#4AC16DE6','#4EC36BE6','#52C569E6',
           '#55C668E6','#59C864E6','#5DC863E6','#60CA60E6','#65CB5EE6',
           '#68CD5BE6','#6DCD59E6','#71CF57E6','#75D054E6','#7AD151E6',
           '#7FD34EE6','#83D44CE6','#87D549E6','#8CD646E6','#90D743E6',
           '#95D840E6','#9AD93CE6','#9FDA3AE6','#A3DA37E6','#A8DB34E6',
           '#ADDC30E6','#B2DD2DE6','#B7DE2AE6','#BBDF27E6')
  return(pal)
}
`colvec` <- function(x) {
  pal <- get_palette()
  return(pal[cut(as.numeric(x), breaks=length(pal), include.lowest=TRUE)])
}
`makecwm` <- function (spe, tra) {
  spe <- as.matrix(spe)
  tra <- as.matrix(tra)
  `standardize` <- function(x) {
    (x - min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE))
  }
  tra <- apply(tra, MARGIN = 2, FUN = standardize)
  awt <- spe %*% tra
  awt / rowSums(spe, na.rm=TRUE)
}
cwm <- data.frame(makecwm(spe, tra)) # make the CWM traits matrix

####    Data     ####

### Load data
# load('./veg.rda')  # assumed to be in your working dir
xy  <- veg$xy    # spatial
spe <- veg$spe   # species
env <- veg$env   # environment
tra <- veg$tra   # traits
phy <- veg$phy   # phylogeny
#rm(veg)          # cleanup
ls()             # objects now in this local environment

## Pre-analysis
### Check data structure
str(spe)     # community (species) abundance matrix
str(env)     # environmental matrix
str(tra)     # traits matrix
str(phy, 1)  # phylogeny
str(xy)      # spatial coordinates

######    Data transformations    ######

### load data
spe <- veg$spe   # species
env <- veg$env   # environment
tra <- veg$tra   # traits
### basic transformations
spe <- data.frame(log10(spe + 1))
env <- data.frame(vegan::decostand(scale(env, center=F), 'range'))
tra <- data.frame(vegan::decostand(tra, 'range'))

###### Visual data exploration  #########
# spatial
plot(xy, pch=19, col='grey', xlab='Eastness', ylab='Northness')
# species
vegan::tabasco(spe, col=get_palette())
# environment
vegan::tabasco(env, col=get_palette())
# traits
vegan::tabasco(tra, col=get_palette())
# phylogeny
plot(phy, cex=0.6, no.margin=TRUE)


####    Diversity measures    ####

### Gamma (regional) diversity
gamma <- sum(colSums(spe) > 0)
gamma
### Alpha (per-site) diversity
alpha <- rowSums(spe > 0)
alpha    # within-site
avgalpha <- mean(rowSums(spe > 0))
avgalpha # average within-site
### Beta (among-site) diversity: Whittaker's
gamma    <- sum(colSums(spe) > 0)
avgalpha <- mean(rowSums(spe > 0))
beta     <- gamma / avgalpha - 1
beta

### breakout group exercise: compare species richness among groups
par(mfrow=c(1,1))
dev.off()
grp <- cut(env$k, quantile(env$k, probs=seq(0,1,by=0.25)), include.lowest=T,
		   labels=c('lo','med','hi','veryhi'))   # group memberships
plot(grp, alpha, xlab="Potassium level", ylab="Plant species richness")
alpha.model <- lm(alpha ~ grp)
summary(alpha.model)



####    Dissimilarities    ####

### Dissimilarity matrix for species
### Using Bray-Curtis distance (abundance-based)
d <- vegdist(spe, method='bray', binary=T)
tabasco(as.matrix(d), col=get_palette())

### Using jaccard distance (presence-absence)
d2 <- vegdist(spe, method='jaccard', binary=T)
tabasco(as.matrix(d2), col=get_palette())

####    Ordination (unconstrained)    ####

# three kinds of unconstrained ordination
m1 <- metaMDS(D, k=2, try=50, trymax=51, trace=0)   # NMS
m2 <- cmdscale(D, k=2, add=T)                       # PCoA
m3 <- prcomp(spe)                                   # PCA
### Visualizing ordinations
# color vector for plotting
u <- get_palette()
u <- u[1:nrow(spe)]
# compare three kinds of ordination
par(mfrow=c(1,3), bty='l', las=1)
plot(m1$points, type='n', xlab='NMDS1', ylab='NMDS2')
text(m1$points, rownames(m1$points),cex=.8, col=u)
plot(m2$points, type='n', xlab='PCoA1', ylab='PCoA2')
text(m2$points, rownames(m2$points),cex=.8, col=u)
plot(m3$x, type='n', xlab='PCA1', ylab='PCA2')
text(m3$x, rownames(m3$x),cex=.8, col=u)
par(mfrow=c(1,1))

### overlay gradients on the NMS using point colors
par(mfrow=c(1,3), bty='l', las=1)
plot(vegan::scores(m1), pch=16, col=colvec(env$k), cex=2)
plot(vegan::scores(m2), pch=16, col=colvec(env$k), cex=2)
plot(vegan::scores(m3), pch=16, col=colvec(env$k), cex=2)
par(mfrow=c(1,1))
### overlay gradients on the NMS by fitting a GAM surface
dev.off()
f1 <- ordisurf(m1 ~ env$k, plot = FALSE)
plot(vegan::scores(m1), pch = 16, col = colvec(env$k))
plot(f1, add = TRUE, col = 1, lwd = 2)

### Dimensionality selection in NMS
# define screeplot function, running NMS for varying 'k' dimensions
`scree_nms` <- function(D, k=5, ...) {
  stress <- rep(NA, k)
  for (i in 1:k) {
    cat('calculating', i, 'of', k, 'dimensions...\n')
    stress[i] <- metaMDS(D, k=i, trace=0, ...)$stress
  }
  plot(1:k, stress, main='', xlab='Dimension', ylab='Stress',
       ylim=c(0, max(stress)*1.05), pch=16, las=1, bty='l')
  lines(1:k, stress)
  abline(0.20, 0, col='red', lty = 2)
  data.matrix(stress)
}
scree <- scree_nms(D, k=5, trymax=10)
scree



#####    Ordination (constrained)    ######

### RDA: redundancy analysis (assumes linear species responses)
m4 <- rda(spe ~ k + sand + conductivity, data=env)
### CCA: constrained correspondence analysis (assumes unimodal species responses)
m5 <- cca(spe ~ k + sand + conductivity, data=env)
### dbRDA: distance-based RDA
m6 <- dbrda(D ~ k + sand + conductivity,
            data = env,
            comm = spe,
            add  = 'lingoes')

### constrained ordination summaries (explore on your own) ####
print(m4)
print(m5)
print(m6)

### compare three kinds of constrained ordination ####
u <- colvec(env$k)         # color by soil potassium
par(mfrow=c(1,3), bty='l', las=1)
plot(m4, disp='wa')
points(m4, pch=16, col=u)
plot(m5, disp='wa')
points(m5, pch=16, col=u)
plot(m6, disp='wa')
points(m6, pch=16, col=u)


####    Group differences  (PERMANOVA, PERMDISP)  ####

### define and visualize four groups
par(mfrow=c(1,1))
brk     <- quantile(env$k, probs=seq(0,1,by=0.25))   # define breaks
grp_k <- cut(env$k, brk, include.lowest=T,
             labels=c('lo','med','hi','veryhi'))   # group memberships
table(grp_k, useNA='always')                       # group tally
tapply(env$k, INDEX = grp_k, FUN = mean)           # group means
plot(m1, cex=0.01)                              # visualize on the NMS
text(m1, labels=grp_k, col=as.numeric(grp_k)) # group memberships
ordispider(m1, groups=grp_k, col=1:4)              # group centroids

### Test for difference in community compositions
### permanova: test for differences in multivariate *centroid*
a1 <- adonis(D ~ grp_k, permu=999)
a1$aov.tab

### Test for homogeneity of community compositions
### permdisp: test for differences in multivariate *dispersion*
b1 <- betadisper(D, grp_k)
b1
permutest(b1, pairwise=TRUE, permu=999)
boxplot(b1)

##### PERMANOVA breakout activity


####    Indicator species analysis   ####

### group plots based on quartiles of potassium levels
grp <- cut(env$k, quantile(env$k, probs=seq(0,1,by=0.25)), include.lowest=T,
           labels=c('lo','med','hi','veryhi'))   # group memberships

iv  <- labdsv::indval(spe, grp) # indicator species analysis for *real* groups
summary(iv)                     # IndVal observed

### random groups
rnd <- sample(grp, length(grp), replace=T) # define random groups by bootstrapping
ivr <- labdsv::indval(spe, rnd) # indicator species analysis for *random* groups
summary(ivr)                    # IndVal expected at random
### null expectation setting alpha = 5%
ceiling( ncol(spe) * 0.05 )


####    Community traits    ####

names(tra)
head(tra)
hist(tra$lfp, breaks=11)

#### Create new data frame for environmental variables plus two trait CWMs
mast <- cbind(veg$env, cwm$spikiness, cwm$lfp)
mast$spikiness <- mast$`cwm$spikiness` #rename column
mast$lfp <- mast$`cwm$lfp` #rename column

#### simple regression with CWM (this is an example of what not to do) ######

#### question 1: Relationship between lfp and k
plot(mast$k, mast$lfp)
c1 <- lm(lfp ~ k, data= mast)
summary(c1)
plot(allEffects(c1))


#### question 2: Relationship between spikiness and k
plot(mast$k, mast$spikiness)
c2 <- lm(spikiness ~ k, data= mast)
summary(c2)
plot(allEffects(c2))

####### alternate approach using linear mixed models (this is a better way to do it) #####

### set up data in long form
spe <- veg$spe
env$site <- rownames(env)
spe$site <- rownames(spe)
tra$species <- rownames(tra)
long_spe <- spe %>% pivot_longer(cols=arisarum_vulgare:hedysarum_coronarium, names_to = "species", values_to = "pres")
long_spe$k <- env$k [match(long_spe$site, env$site)]  ## add in data for k trait
long_spe$spikiness <- tra$spikiness[match(long_spe$species, tra$species)]
long_spe$lfp <- tra$lfp[match(long_spe$species, tra$species)]

l1 <- lmer (pres ~ k #+ I(k^2)
			* lfp + (k|species) #+ (lfp|site)
			, data = long_spe)
summary(l1)
plot(allEffects(l1))


l2 <- lmer (pres ~ k #+ I(k^2)
			* spikiness + (k|species) #+ (lfp|site)
			, data = long_spe)
summary(l2)
plot(allEffects(l2))


####### Trait-environment breakout group ######



####    Community phylogenetics    ####

### Plot phylogenies
### basic plotting
dev.off()
par(mfrow=c(1,3), mar=c(0,0,0,0), oma=c(0,0,0,0)) # plotting parameters
phy                                               # basic structure
plot(phy, cex=0.75, no.margin=T)                  # basic plotting
### color tip labels by trait values
plot(phy, cex=0.75, tip.color=colvec(tra$lfp), no.margin=T)
# axisPhylo(cex=0.6)
### symbolize trait values at tips
plot(phy, cex=0.75, label.offset = 48, no.margin=T)
tiplabels(pch = 21, bg = c(tra$annual), adj = 0)
tiplabels(pch = 21, bg = c(tra$biennial), adj = 15)
tiplabels(pch = 21, bg = c(tra$perennial), adj = 30)

### Phylogenetic diversity - alpha ####
### phylogenetic diversity metrics
spe <- veg$spe
Dp  <- cophenetic(phy) # phylogenetic distances
### Faith's PD phylogenetic diversity = total branch length connecting all species in a site
fpd <- pd(spe, phy)    # Faith's PD
### mean pairwise distance
mpd <- ses.mpd(spe, Dp, null.model='independentswap')
### mean nearest taxon distance
mnd <- ses.mntd(spe, Dp, null.model='independentswap')
### bring all together
phy_div <- cbind(fpd, mpd = mpd$mpd.obs.z, mntd = mnd$mntd.obs.z)
head(phy_div)
pairs(phy_div)


####    Community spatial analysis    ####

### Mantel test
E <- dist(xy)                            # spatial distance matrix
vegan::mantel(D, E, method='spearman')   # spearman *rank* correlation
# in `ecodist` package, we also get useful bootstrap CIs:
ecodist::mantel(D ~ E, mrank=T)          # spearman *rank* correlation
brk <- seq(0, round(max(E),-1), by=10)
plot(ecodist::mgram(D, E, breaks=brk, nperm=99, mrank=T, nboot=500))
abline(h=0, col='red')
plot(vegan::mantel.correlog(D, E, break.pts=brk, cutoff=F, r.type='spearman',
                            nperm=99, mult='holm', progressive=T))

####    END    ####

