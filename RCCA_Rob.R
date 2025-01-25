###LOAD LIBRARIES####

pacman::p_load(igraph, mixOmics, tidyverse, data.table) 
select<-dplyr::select

###IMPORT AND EDIT INPUT TABLES####

#READ NORMALISED GENE COUNTS AND ANNOTATED PRE-TREATED METABOLOMICS MATRIC, BOTH CONTAINING ONLY DEGS/DAMS

met0 <- read.csv("rob_DAM.csv", header = TRUE)
gen0 <- read.csv("rob_DEG.csv", header = TRUE)


###EDIT METABOLOMICS (Y) MATRIX####
#SUBSET METALOMICS MATRIX TO SAME TIMEPOINTS OF RNA-SEQ (0, 48, 96 HPI)

met <- met0 %>% select(grep("X", names(met0)), grep("T48", names(met0)), grep("T96", names(met0)))

#CALCULATE THE MEANS OF 4 REPLICATES, SO FROM 12 REPLICATES, WE'LL HAVE 3 REPS IN EACH CONDITION, TO MATCH RNA-SEQ MATRIX                   
#USE COLUMN NAMES THAT MATCH THE NAMES IN THE RNA-SEQ DATA

met$R29T48R1 <- rowMeans(subset(met, select = c(Rob_I_T48, Rob_I_T48.1, Rob_I_T48.2, Rob_I_T48.3)), na.rm=T)
met$R29T48R2 <- rowMeans(subset(met, select = c(Rob_I_T48.4, Rob_I_T48.5, Rob_I_T48.6, Rob_I_T48.7)), na.rm=T)
met$R29T48R3 <- rowMeans(subset(met, select = c(Rob_I_T48.8, Rob_I_T48.9, Rob_I_T48.10, Rob_I_T48.11)), na.rm=T)

met$RMT48R1 <- rowMeans(subset(met, select = c(Rob_M_T48, Rob_M_T48.1, Rob_M_T48.2, Rob_M_T48.3)), na.rm=T)
met$RMT48R2 <- rowMeans(subset(met, select = c(Rob_M_T48.4, Rob_M_T48.5, Rob_M_T48.6, Rob_M_T48.7)), na.rm=T)
met$RMT48R3 <- rowMeans(subset(met, select = c(Rob_M_T48.8, Rob_M_T48.9, Rob_M_T48.10, Rob_M_T48.11)), na.rm=T)

met$R29T96R1 <- rowMeans(subset(met, select = c(Rob_I_T96, Rob_I_T96.1, Rob_I_T96.2, Rob_I_T96.3)), na.rm=T)
met$R29T96R2 <- rowMeans(subset(met, select = c(Rob_I_T96.4, Rob_I_T96.5, Rob_I_T96.6, Rob_I_T96.7)), na.rm=T)
met$R29T96R3 <- rowMeans(subset(met, select = c(Rob_I_T96.8, Rob_I_T96.9, Rob_I_T96.10, Rob_I_T96.11)), na.rm=T)

met$RMT96R1 <- rowMeans(subset(met, select = c(Rob_M_T96, Rob_M_T96.1, Rob_M_T96.2, Rob_M_T96.3)), na.rm=T)
met$RMT96R2 <- rowMeans(subset(met, select = c(Rob_M_T96.4, Rob_M_T96.5, Rob_M_T96.6, Rob_M_T96.7)), na.rm=T)
met$RMT96R3 <- rowMeans(subset(met, select = c(Rob_M_T96.8, Rob_M_T96.9, Rob_M_T96.10, Rob_M_T96.11)), na.rm=T)


#SUBSET THE COLUMNS WE JUST CREATED

Y0 <- met %>% select(grep("X", names(met)), grep("RM", names(met)), grep("R29", names(met)))

#COLUMN "X" CONTAINS THE M/Z ANNOTATIONS

Y0<-Y0[!duplicated(Y0$X), ]

rownames(Y0) <- Y0$X

Y0 <- Y0[,-1]

write.csv(Y0, "rob_metabolites_avgintensity.csv", row.names = T)

#METABOLITES SHOULD BE IN THE COLUMNS, AND SAMPLES SHOULD BE IN THE ROWS. BOTH MATRICES SHOULD HAVE THE SAME AMOUNT OF ROWS, AND SAME ORDER

Yt <- as.data.frame(t(Y0))

#ROWS SHOULD BE IN THE SAME ORDER AS RNA-SEQ DATA

Yt$names <- paste0(rownames(Yt))
Yt <- arrange(Yt, Yt$names)
Y <- Yt

#MATRIX SHOULD BE NUMERIC. REMOVE ROWNAMES, AND COLUMN WE USED TO SORT ROWS

rownames(Y) <- c()
Y <- Y[, -which(names(Y) %in% "names")]


###EDIT RNA-SEQ COUNTS MATRIX (X)####

rownames(gen0) <- gen0$X

#ARRANGE IN THE SAME MANNER AS METABOLOMICS MATRIX. HERE I'M USING ALPHABETIC ORDER
gen0 <- arrange(gen0, gen0$X)
gen <- gen0[,-1]
X <- gen
rownames(X) <- c()

###CREATE METADATA####

#PASTE METABOLOMICS NAMES

metadata <- as.data.frame(paste0(rownames(Yt)))

#PASTE RNA-SEQ NAMES
metadata$RNAnames <- paste0(rownames(gen))

#CHECK IF THEY MATCH. THEIR POSITION/NAMES SHOULD BE THE SAME

#REMOVE REPLICATE INFO TO CREATE CLASS
metadata$Group <- as.factor(paste0(substr(metadata$RNAnames, 1, nchar(metadata$RNAnames)-2)))

#EXTRACT TREATMENT, GENOTYPE AND TIMEPOINT INFO. THESE SHOULD BE FACTORS SO YOU CAN USE TO PLOT LATTER
metadata$Treatment <- as.factor(ifelse(metadata$Group %like% "29", "I", "M"))
metadata$Genotype <- as.factor(ifelse(metadata$Group %like% "H", "Her", "Rob"))
metadata$Background <- as.factor(paste0(metadata$Genotype, "_", metadata$Treatment))
metadata$hpi <- as.factor(paste0(sub(".*T", "", metadata$Group)))


###STANDARDISE MATRICES: ZERO MEAN/UNIT VARIENCE####

#CHECK IF MATRIX IS STANDARDISED. SD SHOULD BE 1

sapply(Y, sd, na.rm=T)
sapply(X, sd, na.rm=T)

#IF NOT, USE SCALE FUNCTION TO STANDARDISE DATA

Ysc <- as.data.frame(scale(Y))
sapply(Ysc, sd, na.rm=T)
Xsc <- as.data.frame(scale(X))
sapply(Xsc, sd, na.rm=T)

#X AND Y MATRICES ARE FINALLY READY TO BE INTEGRATED

Y <- Xsc
X <- Ysc

dim(X); dim(Y)


###RCCA####
shrink <- rcc(X, Y, ncomp = 3, method = 'shrinkage')

plot(shrink, type = "barplot")

corr <- shrink$cor

#Calculate correlation matrix

bisect = shrink$variates$X[, 1:3] + shrink$variates$Y[, 1:3]
cord.X = cor(shrink$X, bisect, use = "pairwise")
cord.Y = cor(shrink$Y, bisect, use = "pairwise")
simMat = as.matrix(cord.X %*% t(cord.Y))

write.csv(simMat, "rcc_matrix.csv")

#EXTRACT THE CORRELATION MATRIX BY USING SUMMARY FUNCTION

more <- summary(shrink, cutoff = 0.6)

Cm.X <- more$Cm.X
Cm.Y <- more$Cm.Y

write.csv(Cm.X, "communalities_DEGs.csv")
write.csv(Cm.Y, "communalities_DAMs.csv")

###EXPORT NETWORK IN GML

myNetwork.rcc <- network(shrink, comp = 1:3, cutoff = 0.6)

write_graph(myNetwork.rcc$gR, file = "robrccnetwork06.gml", format = "gml")

