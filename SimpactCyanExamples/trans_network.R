# Before running the code below, make sure you have followed the instructions for installing all the required software,
# as explained in the README file.
setwd("/home/acer/SimpactCyanExamples")

# Make sure compiled tools (Seq-Gen and FastTree) are in same working directory
Sys.setenv(PATH=paste("/home/acer/simpact-cyan-0.21.0/build",Sys.getenv("PATH"),sep=":"))
Sys.setenv(PYTHONPATH="/home/acer/simpact-cyan-0.21.0/python")
Sys.setenv(SIMPACT_DATA_DIR="/home/acer/simpact-cyan-0.21.0/data/")

# Enable installing from bioconductor
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("ggtree")
## Load required packages
library(RSimpactCyan)
library(devtools)
#install_github("wdelva/RSimpactHelp")
library(RSimpactHelper)
library(tidyverse)
library(ggplot2)
library(geomnet)


# Sourcing the IDs.Seq.Random.skew function and renaming IDS in transmission table 
# to match tips names in the tree
source("util.seq.cov.weight.R")

#######################
# Step 1: Run Simpact # 
#######################

EAAA.SciRep.revision.Seq.loaded.object <- load(file = paste0(getwd(), "/EAAA.SciRep.revision2.Seq.RData"))
# load(file = "/Users/delvaw/Google Drive/SimpactCyanPaper/Scientific Reports/Revision/EAAA.SciRep.revision2.Seq.RData")

Seq.object <- EAAA.SciRep.revision2.Seq 
# The vector of target features was:
features.pop.growth <- exp(0.015)
features.hiv.prev <- c(0.143, 0.008, 0.315, 0.066, 0.467, 0.213, 0.538, 0.366, 0.491, 0.47, 0.397, 0.455, 0.316, 0.425)
features.hiv.inc <- exp(c(0.038, 0.008, 0.043, 0.016, 0.02, 0.026, 0.027, 0.031, 0.04, 0.004, 0.021, 0.012, 0.012, 0))
features.art.cov <- c(0.37, 0.40, 0.44, 0.49, 0.58, 0.67, 0.76, 0.85) # c(0.33, 0.38, 0.45, 0.51, 0.61, 0.7, 0.8)
features.vl.suppr <- 0.74  # 0.68
features.unaids.prev <- c(0.017, 0.033, 0.057, 0.088, 0.124, 0.161, 0.194, 0.220, 0.239, 0.251, 0.258, 0.261, 0.261, 0.259, 0.257, 0.255, 0.256, 0.259, 0.263, 0.268, 0.274, 0.278, 0.282, 0.284, 0.284, 0.283, 0.279, 0.274)
target.features.EAAA <- c(features.pop.growth, features.hiv.prev, features.hiv.inc, features.art.cov, features.vl.suppr, features.unaids.prev)

# To find the best fitting model, we compute the RMSE for the 250 parameter combinations of the posterior.
posterior.size <- nrow(Seq.object$stats)
RMSE.vect <- rep(NA, posterior.size)


for (i.posterior in 1:posterior.size){
  fitting.features.EAAA <- Seq.object$stats[i.posterior, ]
  # The RMSE compared to the target statistics was:
  RMSE.vect[i.posterior] <- sqrt(sum(((fitting.features.EAAA - target.features.EAAA)/target.features.EAAA)^2) / length(target.features.EAAA))
}
index.bestfit <- which(RMSE.vect == min(RMSE.vect))
bestfitting.features.EAAA <- Seq.object$stats[index.bestfit, ]


# Using the best-fitting model from the calibration
inputvector.EAAA.example <- Seq.object$param[index.bestfit, ]


inputvector <- c(0, inputvector.EAAA.example)

age.distr <- agedistr.creator(shape = 5, scale = 65)

cfg.list <- input.params.creator(population.eyecap.fraction = 1,
                                 population.simtime = 52,
                                 population.nummen = 500,
                                 population.numwomen = 500,
                                 population.msm = "no",
                                 hivseed.time = 8.9 ,
                                 hivseed.type = 'amount',
                                 hivseed.amount = 20, #30,
                                 hivseed.age.min = 20,
                                 hivseed.age.max = 50,
                                 hivtransmission.param.a = -1,
                                 hivtransmission.param.b = -90,
                                 hivtransmission.param.c = 0.5,
                                 hivtransmission.param.f1 = log(2),
                                 hivtransmission.param.f2 = log(log(sqrt(2)) / log(2)) / 5,
                                 formation.hazard.agegapry.gap_factor_man_age = -0.01,
                                 formation.hazard.agegapry.gap_factor_woman_age = -0.01,
                                 formation.hazard.agegapry.meanage = -0.025,
                                 formation.hazard.agegapry.gap_factor_man_const = 0,
                                 formation.hazard.agegapry.gap_factor_woman_const = 0,
                                 formation.hazard.agegapry.gap_factor_man_exp = -1,
                                 formation.hazard.agegapry.gap_factor_woman_exp = -1,
                                 formation.hazard.agegapry.gap_agescale_man = 0.25,
                                 formation.hazard.agegapry.gap_agescale_woman = 0.25,
                                 dissolution.alpha_4 = -0.05,
                                 debut.debutage = 15,
                                 conception.alpha_base = -2.7,
                                 dropout.interval.dist.type = "uniform")

#standard deviation of 200 CD4 cells
#mu = ln(mean / sqrt(1 + variance/mean^2))
#sigma^2 = ln(1 + variance/mean^2)
#Here, we say mean = 825 and variance = 200^2
mu.cd4 <- 800
var.cd4 <- 200^2
mu.cd4.end <- 20
var.cd4.end <- 5
cfg.list["person.cd4.start.dist.type"] <- "lognormal"
cfg.list["person.cd4.start.dist.lognormal.zeta"] <- log(mu.cd4/sqrt(1+var.cd4/mu.cd4^2))
cfg.list["person.cd4.start.dist.lognormal.sigma"] <- sqrt(log(1+var.cd4/mu.cd4^2))
cfg.list["person.cd4.end.dist.type"] <- "lognormal"
cfg.list["person.cd4.end.dist.lognormal.zeta"] <- log(mu.cd4.end/sqrt(1+var.cd4.end/mu.cd4.end^2))
cfg.list["person.cd4.end.dist.lognormal.sigma"] <- sqrt(log(1+var.cd4.end/mu.cd4.end^2))

cfg.list["formation.hazard.agegapry.baseline"] <- 2
cfg.list["mortality.aids.survtime.C"] <- 65
cfg.list["mortality.aids.survtime.k"] <- -0.2
cfg.list["monitoring.fraction.log_viralload"] <- 0 #0.3

cfg.list["dropout.interval.dist.uniform.min"] <- 1000
cfg.list["dropout.interval.dist.uniform.max"] <- 2000

cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1

cfg.list["person.agegap.man.dist.type"] <- "normal"
cfg.list["person.agegap.woman.dist.type"] <- "normal"

cfg.list["monitoring.cd4.threshold"] <- 1 # 0
cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.5 # 1 # 0.9 # 0.75 # 0.5
cfg.list["diagnosis.baseline"] <- -99999 # -2
cfg.list["periodiclogging.interval"] <- 0.25
# cfg.list["dropout.interval.dist.exponential.lambda"] <- 0.1


cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 6

cfg.list["person.vsp.toacute.x"] <- 5 # See Bellan PLoS Medicine

seedid <- inputvector[1]
cfg.list["hivtransmission.param.f1"] = log(inputvector[2])
cfg.list["hivtransmission.param.f2"] = log(log(sqrt(inputvector[2])) / log(inputvector[2])) / 5
cfg.list["formation.hazard.agegapry.gap_agescale_man"] = inputvector[3]
cfg.list["formation.hazard.agegapry.gap_agescale_woman"] = inputvector[3]
cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[4]
cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[4]
cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[5]
cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[5]
cfg.list["person.eagerness.man.dist.gamma.a"] <- inputvector[6]
cfg.list["person.eagerness.woman.dist.gamma.a"] <- inputvector[7]
cfg.list["person.eagerness.man.dist.gamma.b"] <- inputvector[8]
cfg.list["person.eagerness.woman.dist.gamma.b"] <- inputvector[9]
cfg.list["formation.hazard.agegapry.gap_factor_man_exp"] <- inputvector[10]
cfg.list["formation.hazard.agegapry.gap_factor_woman_exp"] <- inputvector[10]
cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[11]
cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[12]
cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[13]
cfg.list["conception.alpha_base"] <- inputvector[14]
cfg.list["dissolution.alpha_0"] <- inputvector[15]

# Introducing ART
art.intro <- list()
art.intro["time"] <- 20
art.intro["diagnosis.baseline"] <- inputvector[16] # prior [-4 , 0] # -2
art.intro["monitoring.cd4.threshold"] <- 100

art.intro1 <- list()
art.intro1["time"] <- 22
art.intro1["diagnosis.baseline"] <- inputvector[16]  #+ inputvector[17] # prior [0, 2] # -1.8
art.intro1["monitoring.cd4.threshold"] <- 500

art.intro2 <- list()
art.intro2["time"] <- 23
art.intro2["diagnosis.baseline"] <- inputvector[16] + inputvector[17] + inputvector[18] # prior [0, 2] # -1.5
art.intro2["monitoring.cd4.threshold"] <- 200
# 
art.intro3 <- list()
art.intro3["time"] <- 30
art.intro3["diagnosis.baseline"] <- inputvector[16] + inputvector[17] + inputvector[18] + inputvector[19] # prior [0, 2] #-1
art.intro3["monitoring.cd4.threshold"] <- 350
  
art.intro4 <- list()
art.intro4["time"] <- 33.5
art.intro4["diagnosis.baseline"] <- inputvector[16] #+ inputvector[17] + inputvector[18] + inputvector[19] + inputvector[20] # prior [0, 2]
art.intro4["monitoring.cd4.threshold"] <- 500

art.intro5 <- list()
art.intro5["time"] <- 36.75
art.intro5["monitoring.cd4.threshold"] <- 6000


ART.factual <- list(art.intro, art.intro1, art.intro2, art.intro3, art.intro4, art.intro5)
#art.intro,art.intro1, art.intro2, art.intro3, art.intro4, art.intro5

identifier <- paste0(seedid)
destDir <- paste0(getwd(), "/temp")

results <- tryCatch(simpact.run(configParams = cfg.list,
                                destDir = destDir,
                                agedist = age.distr,
                                intervention = list(art.intro1),
                                seed = seedid),
                    error = simpact.errFunction)
datalist.phylo <- readthedata(results)

###########################################
# Step 2: Construct transmission networks #
###########################################
# Produce a list of transmission networks in epi object format
simpact.trans.net <- transmission.network.builder(datalist = datalist.phylo, endpoint = 40)

net.sizes <- purrr::map(simpact.trans.net, 1) %>%
  lapply(., length) %>%
  unlist()
max.net.size <- max(net.sizes)
max.net.size.index <- which(net.sizes %in% max.net.size) # So we know which network and tree to visualise

revised.network.df <- data.frame(from_id = factor(simpact.trans.net[[max.net.size.index]]$parent),
                                 to_id = factor(simpact.trans.net[[max.net.size.index]]$id), stringsAsFactors = FALSE)
revised.network.vertices.df <- data.frame(vertices = as.character(unlist(revised.network.df)), stringsAsFactors = FALSE)
revised.network.fortified <- fortify(as.edgedf(revised.network.df), revised.network.vertices.df, stringsAsFactors = FALSE)

revised.network.fortified_na = na.omit(revised.network.fortified)
library(igraph)
library(data.tree)
edges1 = matrix(c(revised.network.fortified$from_id, revised.network.fortified$to_id),byrow = TRUE, ncol = 2)

#G = graph_from_data_frame(revised.network.fortified_na, directed = TRUE, vertices = NULL)


# A. Transmission network
# The transmission event from the "ghost infector with ID "-1" to the seed individual with ID "858" is excluded from the transmission network

transmissionnetwork.plot <- ggplot(data = revised.network.fortified[2:nrow(revised.network.fortified), ]) +
  geomnet::geom_net(data = revised.network.fortified[2:nrow(revised.network.fortified), ], na.rm = TRUE, 
                    aes(from_id = from_id,
                        to_id = to_id),
                    directed = TRUE,
                    size = 1, # 2.5
                    layout.alg = "kamadakawai",
                    alpha = 1,
                    arrowsize = 0.25, # 0.5
                    arrowgap = 0.004,
                    ecolour = "darkgrey",
                    colour = "black") +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        panel.background = element_rect(fill = "grey97")) +
  ylab("")
print(transmissionnetwork.plot)

  ggsave(filename = "network_revised.pdf",
       plot = transmissionnetwork.plot,
       path = paste0(getwd(), "/plots"),
       width = 10, height = 10, units = "cm")

plot_degree_distribution <-function(graph,a) {
  d = igraph::degree(graph, mode = a)
  dd = igraph::degree.distribution(graph, mode = a, cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  
  TYPE<-paste0("Degree_",a)
  TITLErev<-paste0("Degree Distribution_",a)
  plotdata<-cbind(probability,degree)
  plotdata<-as.data.frame(plotdata)
  colnames(plotdata)<-c("probability","degree")
  #graphics::plot(probability ~ degree, log = "xy", xlab = TYPE, ylab = "Probability (log)",
  #     col = 1, main = TITLErev)
  ggplot2::ggplot(plotdata, ggplot2::aes(x=degree, y=probability)) +
    ggplot2::geom_point(color="darkred")+
    ggplot2::geom_jitter(position = "jitter",color="darkred")+
    ggplot2::labs(title=TITLErev,
                  x=TYPE, y = "Probability")+
    ggplot2::theme_gray()
}

#plot_degree_distribution(G, 'out')

#similarity_metric <- similarity(graph = G, vids = V(G), method = "invlogweighted")




edges1 <- na.omit(edges1)
inds <- which(duplicated(edges1) == TRUE)
edges_1 <- edges1[-c(inds),]
Edges <- data.frame(Parent=c(edges_1[,1]), Child=c(edges_1[,2]))

# length(edges_1)
G <- graph(edges=t(Edges), directed = TRUE)
write_graph(graph = G, file = '/home/acer/Transmission_Networks/mixed//transmission_network_15_500.csv', format = 'edgelist')

  # distance matrix 
dist_trans_network = dist(Edges)

# hierarchical clustering analysis
clus_trans_network = hclust(dist_trans_network, method = "average")

# plot dendrogram
plot(clus_trans_network, hang = 0.1)

library(ape)

# convert 'hclust' to 'phylo' object

phylo_tree = as.phylo(clus_trans_network)

plot(phylo_tree)

library(phyloTop)
library(ggtree)

#ggtree(phylo_tree)

#topo_sort(G, mode='out')

