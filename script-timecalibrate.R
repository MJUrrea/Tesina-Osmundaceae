#Load packages
library(phytools)
library(ape)
library(treeio)
library(TreeTools)
library(Claddis)
library(ggtree)
library(tidytree)
library(tidyverse)
library(stringr)
library(ggplot2)
library(deeptime)
library(dplyr)
library(castor)
library(paleotree)
library(svglite)
library(phangorn)
library(tibble)
library(ggplotify)
library(tidyr)
library(plyr)
library(scales)
library(patchwork)
library(zoo)



#00. FUNCTIONS ####
  vectoreame <- function(dated_multiphylo_list, dated_multiphylo_tibble_list, reference_phylo, reference_tibble, A){
  
  if(!inherits(dated_multiphylo_list, "multiPhylo")){ stop("Not multiPhylo list object given")}
  if(!is.list(dated_multiphylo_tibble_list)){ stop("Not dated trees tibble list object given")}
  
  if(!inherits(reference_phylo, "phylo")){ stop("Not reference Phylo object given")}
  if(!inherits(reference_tibble, "tbl_tree")){ stop("Not reference tibble tree object given")}
  
  #Actual code starts here
  eqterminals <- vector(mode = "logical", length = length(reference_phylo$tip.label))
  eqterminals[] <- NA
  tr2labels <- matchLabels(reference_phylo, dated_multiphylo_list[[A]])
  tr2labels <- as.vector(tr2labels[,2])
  eqterminals[1:length(tr2labels)] <- tr2labels
  
  eqnodes <- vector(mode = "logical", length = (nrow(reference_tibble)-length(reference_phylo$tip.label)) )
  eqnodes[] <- NA
  tr2nodes <- matchNodes(reference_phylo, dated_multiphylo_list[[A]], method = "descendants")
  tr2nodes <- as.vector(tr2nodes[,2])
  eqnodes[1:length(tr2nodes)] <- tr2nodes
  
  eqgroups <- c(eqterminals,eqnodes)
  eqages <- vector(mode = "integer", length = length(eqgroups))
  eqages[] <- NA
  
  eqbranchlen <- vector(mode = "integer", length = length(eqgroups))
  eqbranchlen[] <- NA
  
    
  for (zz in 1:length(eqgroups)) {
    if(!is.na(eqgroups[zz])){
      eqages[zz] <- dated_multiphylo_tibble_list[[A]]$nodeages[eqgroups[zz]]
      eqbranchlen[zz] <- dated_multiphylo_tibble_list[[A]]$branch.length[eqgroups[zz]]
    }
  }
  #end of actual code
  
  stringA <- paste("node", A, sep = ".")
  stringB <- paste("ages", A, sep = ".")
  stringC <- paste("branch.length", A, sep = ".")
  
  output_matrix <- matrix(data = c(eqgroups, eqages, eqbranchlen), ncol = 3)
  colnames(output_matrix) <- c(stringA, stringB, stringC)
  
  return(output_matrix)
  }
  
  #To process LTT and paleoclimate data. It uses Judd et al climate data
  #As compared to traditional LTT plots, this function "smooths" the LTT profile
  ltt_and_clime <- function(mytrees, clime="PhanDA_GMSTandCO2_percentiles.csv"){
    #Output: 
     #A list with:
     #'lttclime' tibble that stores median ltt and median temperature, and their respective time points.
     #'ltt_ribbon_data' dataframe with 95% limits and time points for creating a ribbon around median ltt
     #'clim_ribbon_data' dataframe with 95% limits and time points for creating a ribbon around median temperature   
    
    #LTT from mytrees
    if(inherits(mytrees, "multiPhylo")){
      cat("Handling more than a single input tree - they are assumed to be the same tree but with different inferred ages\n")
      flag_multiphylo = 1
      mytree_ltt <- ltt.plot.coords(mytrees[[1]])
      mytree_ltt[,1] <- mytree_ltt[,1] * (-1)
      colnames(mytree_ltt)[1] <- paste("_","1",".", colnames(mytree_ltt)[1], sep = "")
      colnames(mytree_ltt)[2] <- paste("_","1",".", colnames(mytree_ltt)[2], sep = "")
      for (i in 2:length(mytrees)) {
        toadd <- ltt.plot.coords(mytrees[[i]])
        toadd[,1] <- toadd[,1] * (-1)
        colnames(toadd)[1] <- paste("_",i,".", colnames(toadd)[1], sep = "")
        colnames(toadd)[2] <- paste("_",i,".", colnames(toadd)[2], sep = "")
        mytree_ltt <- cbind(mytree_ltt, toadd)
      }
      
      mytree_ltt <- as.data.frame(mytree_ltt)
      
      #Select LTT columns
      median_ltt <- mytree_ltt %>% select(contains(".N"))
      #Calculate median + 95% CI (2.5% and 97.5% quantiles) for each row
      ltt.row_stats <- t(apply(median_ltt, 1, function(x) {
        quantile(x, probs = c(0.5, 0.025, 0.975), na.rm = TRUE)
      }))
      #Rename columns for clarity
      colnames(ltt.row_stats) <- c("median.ltt", "ltt.lower_95", "ltt.upper_95")
      
      #Select time columns
      median_time <- mytree_ltt %>% select(contains(".time"))
      #Calculate median + 95% CI (2.5% and 97.5% quantiles) for each row
      time.row_stats <- t(apply(median_time, 1, function(x) {
        quantile(x, probs = c(0.5, 0.025, 0.975), na.rm = TRUE)
      }))
      #Rename columns for clarity
      colnames(time.row_stats) <- c("median.time", "time.lower_95", "time.upper_95")
      
      summary_ltt <- cbind(ltt.row_stats, time.row_stats)
      summary_ltt <- as.data.frame(summary_ltt)
      
      # Create a separate dataset for geom_ribbon
      ltt_ribbon_data <- summary_ltt %>%
        select(median.time, ltt.lower_95, median.ltt, ltt.upper_95) %>%
        drop_na()
      # Where the lower limits have no difference with median, I create an artificiial lower limit of 1
      difflow = ltt_ribbon_data[[3]] - ltt_ribbon_data[[2]]
      difflow_index = which(diff == 0)
      ltt_ribbon_data[[2]][difflow_index] <- ltt_ribbon_data[[3]][difflow_index] - 1
      # Where the upper limits have no difference with median, I create an artificiial lower limit of 1
      diffup = ltt_ribbon_data[[4]] - ltt_ribbon_data[[3]]
      diffup_index = which(diff == 0)
      ltt_ribbon_data[[4]][difflow_index] <- ltt_ribbon_data[[3]][difflow_index] + 1
      
      
    } #end of multiPhylo
    else{
      cat("There is a single input tree\n")
      flag_phylo = 1
      mytree_ltt <- ltt.plot.coords(mytrees)
      mytree_ltt[,1] <- mytree_ltt[,1] * (-1)
      colnames(mytree_ltt)[1] <- paste("median",".", colnames(mytree_ltt)[1], sep = "")
      colnames(mytree_ltt)[2] <- paste("median.", colnames(mytree_ltt)[2], sep = "")
      
      summary_ltt <- mytree_ltt
    }
    
    #Paleoclimate dataset
    
    # Read paleoclimate data from Judd et al. (2024): 10.1126/science.adk3705
    paleodata <- read.csv(clime, header = TRUE, sep = ",")
    paleodata <- paleodata[, (names(paleodata) %in% c("AverageAge", "GMST_05", "GMST_50", "GMST_95"))]
    colnames(paleodata)[1] <- "clim_time"
    
    # Create a separate dataset for geom_ribbon
    paleo_ribbon_data <- combined_paleodata_ltt %>%
      select(clim_time, GMST_05, GMST_50, GMST_95) %>%
      drop_na()
    
    # Combining with summary_ltt
    combined_paleodata_ltt <- merge(summary_ltt, paleodata, all = TRUE) %>% 
      drop_na()
    
    # Reshape data for being plotted
    paleoltt_long_data <- combined_paleodata_ltt %>%
      pivot_longer(
        cols = c(median.ltt, GMST_50),
        names_to = "variable",
        values_to = "value"
      ) %>%
      mutate(time = case_when(
        variable == "median.ltt" ~ median.time,
        variable == "GMST_50" ~ clim_time
      ))
    
    
    
    return(list(lttclime = paleoltt_long_data, clim_ribbon_data = paleo_ribbon_data, ltt_ribbon_data = ltt_ribbon_data))
  }
  
  #Function to rescale GMST values relative to LTT values
  gmst_rescaler <- function(ltt_data, gmst_data) {
    function(x) {
      (x - min(ltt_data$value, na.rm = TRUE)) /
        (max(ltt_data$value, na.rm = TRUE) - min(ltt_data$value, na.rm = TRUE)) *
        (max(gmst_data$GMST_50, na.rm = TRUE) - min(gmst_data$GMST_50, na.rm = TRUE)) +
        min(gmst_data$GMST_50, na.rm = TRUE)
    }
  }
  
  

#0.  STEP - INITIALIZING SOME VARIABLES FOR TIME SCALING ####
  vartimer = 30
  sampletime = 1000
  
  inferior_limit = (-325)
  rgbfill <- rgb(red = 185, green = 239, blue = 249, maxColorValue = 255)
  breaks_vec <- c(inferior_limit,
                  -193, -175, #Pleinsbachian-Toarcian, Toarcian-Aalenian (PTo-E and T-OAE events)
                  -49, -16    #Opening of the Drake Passage and thermal isolation of Antarctica
                              #Even though its opening it is assumed to have concluded by 16-17 Ma, its beginning is blurry but estimated
                              #to be around the Eocene-Oligocene transition.
                              #See:
                              #Scher & Martin (2006): 10.1126/science.1120044
                              #Vincze et al. (2021): https://doi.org/10.1038/s41598-021-99123-0
                  )
  
#1.  STEP - PREPARING DATA ####
  #Open calibrated non-parametric bootstrap replicates from Sterli et al.'s TNT script
  #branch lengths here actually are discarded
  optimal_tree <- ReadTntTree("arbol1_pi5.tre")
  suboptimal_trees <- ReadTntTree("arbol1_pi5_subopt.tre")
  suboptimal_trees <- sample(suboptimal_trees, 10, replace = F)

  # checks whether I have one or multiple trees in "optimal_tree"
  if (inherits(optimal_tree, "multiPhylo")) {
    cat("The object is of class 'multiPhylo' -- a strict consensus is computed \n")
    
    #Establish a consensus tree, to function as reference
    reference_tree <- consensus(optimal_tree, rooted = TRUE)
    
    } else if (inherits(optimal_tree, "phylo"))
      { cat("The object is of class 'phylo' -- no strict consensus computed \n")
        reference_tree <- optimal_tree} else { stop("The object is neither 'phylo' nor 'multiPhylo'.\n") }
  
  # all individual trees
    allreplicates <- c(optimal_tree, suboptimal_trees)
  
  #The max-minFAD (or FAD-LAD) matrix
  #check that tiplabels and names in FAD match
    FAD <- read.csv("FAD.csv", header=T, row.names=1)
    macheado <- match(reference_tree$tip.label, rownames(FAD)) 
    excludeme <- NULL
    for(x in 1:length(macheado)) {if(is.na(macheado[x])) excludeme <- append(excludeme, x, after = length(excludeme))}
  
    if(!is.null(excludeme)){macheado <- macheado[-c(excludeme)]}
    
    FAD <- FAD[macheado,]
  
    for(x in 1:nrow(FAD)){
      if(FAD$FAD[x] < 6) FAD$FAD[x] <- 0
      if(FAD$LAD[x] < 6) FAD$LAD[x] <- 0
      }
    rangeCont <- as.matrix(FAD)
  
  #Plot reference
  reference_tree_tibble <- as_tibble(reference_tree)
  reference_tree_tibble$median_age <- vector(mode = "logical", length = nrow(reference_tree_tibble))
  reference_tree_tibble$min_age <- vector(mode = "logical", length = nrow(reference_tree_tibble))
  reference_tree_tibble$max_age <- vector(mode = "logical", length = nrow(reference_tree_tibble))
  reference_tree_tibble$branch.length <- vector(mode = "logical", length = nrow(reference_tree_tibble))
  ggtree(as.phylo(reference_tree_tibble), branch.length = "none") + geom_tiplab(size=2)


#2.  STEP - TIME-SCALING INDIVIDUAL TREES ####
  #Time-scale with timePaleoPhy each tree found in "optimal_tree" and "suboptimal_trees"
  #This will return optimal and suboptimal trees (replicates) with root time and branch lengths.
  #Such replicates are then used to estimate node ages in the reference tree
  #using "dateNodes" and the respective root time.
  dated_trees <- NULL
  for (i in 1:length(allreplicates)) {
      tree_i <- allreplicates[[i]]
      tree_i_dated <- timePaleoPhy(tree = tree_i,
                                   ntrees = sampletime, 
                                   randres = F, 
                                   vartime = vartimer,
                                   dateTreatment = "minMax",
                                   timeData = rangeCont, 
                                   type = "equal", 
                                   plot = F)
  
      dated_trees <- append(dated_trees, tree_i_dated, after = length(dated_trees))
  }
  #Create a list of tibbles, where each tibble is a dated tree
  dated_trees_tibbled <- vector(mode = "list", length = length(dated_trees))
  for (z in 1:length(dated_trees_tibbled)) {
    dated_trees_tibbled[[z]] <- as_tibble(dated_trees[[z]])
    dated_trees_tibbled[[z]][["nodeages"]] <- dateNodes(tree = dated_trees[[z]])
    dated_trees_tibbled[[z]][["nodeages"]] <- round(dated_trees_tibbled[[z]][["nodeages"]], digits = 3)
    dated_trees_tibbled[[z]][["nodeages"]] <- as.vector(dated_trees_tibbled[[z]][["nodeages"]])
  }
  
  #Now time-scale only the optimal trees to compute consensus edge.lengths in the reference tree
  optimal_dated <- NULL
  if(inherits(optimal_tree, "phylo")){
    
    tree_i_optimal <- optimal_tree
    tree_i_dated_optimal <- timePaleoPhy(tree = tree_i_optimal,
                                         ntrees = sampletime, 
                                         randres = F, 
                                         vartime = vartimer,
                                         dateTreatment = "minMax",
                                         timeData = rangeCont, 
                                         type = "equal", 
                                         plot = F)
    optimal_dated <- append(optimal_dated, tree_i_dated_optimal, after = length(optimal_dated))
    
  } else if(inherits(optimal_tree, "multiPhylo")) {
    
    for (i in 1:length(optimal_tree)) {
      tree_i_optimal <- optimal_tree[[i]]
      tree_i_dated_optimal <- timePaleoPhy(tree = tree_i_optimal,
                                           ntrees = sampletime, 
                                           randres = F, 
                                           vartime = vartimer,
                                           dateTreatment = "minMax",
                                           timeData = rangeCont, 
                                           type = "equal", 
                                           plot = F)
      optimal_dated <- append(optimal_dated, tree_i_dated_optimal, after = length(optimal_dated))}
    
  }
  consedges <- consensus.edges(optimal_dated)
  
  
  
  
  final_mat <- NULL
  for (j in 1:length(dated_trees)) {
    temporal_mat <- vectoreame(dated_multiphylo_list = dated_trees, dated_multiphylo_tibble = dated_trees_tibbled,
                               reference_phylo = reference_tree, reference_tibble = reference_tree_tibble, A = j)
    final_mat <- cbind(final_mat, temporal_mat)
  }
  
  # Remove columns with names matching "node.X" and "ages.X"
  final_branch.lengths <- final_mat[, !grepl("^node\\..*$", colnames(final_mat))]
  final_branch.lengths <- final_branch.lengths[, !grepl("^ages\\..*$", colnames(final_branch.lengths))]
  # Estimate median of branch lengths as final columns
  final_branch.lengths <- cbind(final_branch.lengths, NA)
  colnames(final_branch.lengths)[ncol(final_branch.lengths)] <- "branch.length"
  for (yy in 1:nrow(final_branch.lengths)) {
    final_branch.lengths[yy,ncol(final_branch.lengths)] <- median(unname(unlist(final_branch.lengths[yy,1:ncol(final_branch.lengths)])),na.rm = T) 
  }
  
  reference_tree_tibble <- cbind(reference_tree_tibble, final_mat)
  
  # Remove columns with names matching "node.X" and "branch.length.X"
  reference_tree_tibble <- reference_tree_tibble[, !grepl("^node\\..*$", colnames(reference_tree_tibble))]
  reference_tree_tibble <- reference_tree_tibble[, !grepl("^branch\\.length\\..*$", colnames(reference_tree_tibble))]
  
  # Estimate median values from ages, minimum and maximum ages and median branch lengths
  for (ii in 1:nrow(reference_tree_tibble)) {
    
    #Ages:
    reference_tree_tibble$median_age[ii] <- median(unname(unlist(reference_tree_tibble[ii,8:ncol(reference_tree_tibble)])),
                                                   na.rm = T)
    reference_tree_tibble$min_age[ii] <- min(unname(unlist(reference_tree_tibble[ii,8:ncol(reference_tree_tibble)])),
                                             na.rm = T)
    reference_tree_tibble$max_age[ii] <- max(unname(unlist(reference_tree_tibble[ii,8:ncol(reference_tree_tibble)])),
                                             na.rm = T)
  
  }
  
  
  
  #Tibble to plot later on
  main_tibble <- reference_tree_tibble[,1:6]
  main_tibble <- cbind(main_tibble,final_branch.lengths[,ncol(final_branch.lengths)])
  colnames(main_tibble)[ncol(main_tibble)] <- "branch.length"
  #Create ranges from min and max ages
  range_string <-  list()
  for (x in 1:nrow(main_tibble)) {
    
      range_string[[x]] <- c(main_tibble$min_age[[x]], main_tibble$max_age[[x]])
  
  }
  cols_to_delete <- c("min_age", "max_age")
  main_tibble <- main_tibble[,!colnames(main_tibble) %in% cols_to_delete]
  main_tibble$range <- range_string
  
  main_tibble <- as_tibble(main_tibble)
  write.nexus(consedges, file = "temporal.temp")
  main_tree <- as_tibble(TNTOrder(read.nexus("temporal.temp")))
  main_tree$median_ages <- main_tibble$median_age
  main_tree$range <- main_tibble$range
  for (i in 1:length(main_tree$label)) { 
    if(main_tree$label[i] == "1") {
      main_tree$label[i] <- paste("HTU_", i, sep = "")
    }
  }
  
  main_tree <- as.treedata(main_tree)

#2B. STEP - PLOTTING TREES ####
  # Finding ancestors to highlight clade
  mrca_I <- getMRCA(main_tree@phylo, tip = c("Osmundopsis_sturii","Todea_papuana"))
  mrca_lep <- getMRCA(main_tree@phylo, tip = c("Todites_cacereii","Todea_papuana"))
  # Replace "_" by blank
  main_tree@phylo$tip.label <- gsub("_", " ", main_tree@phylo$tip.label)
  main_tree@phylo$edge.length[26] <- 65.36209
  #Start plotting
  f3 <- ggtree(main_tree, layout = "rectangular", 
               ladderize=TRUE,
               right=TRUE,
               size=1.75,
               color="darkgrey")+
               #aes(color=body_mass))+
    geom_rect(aes(xmin = breaks_vec[3], xmax = breaks_vec[2], ymin = -2, ymax = Ntip(main_tree)+2.5), 
              fill = rgbfill, alpha = 0.20, inherit.aes = FALSE, linewidth = 0) +
    geom_rect(aes(xmin = breaks_vec[5], xmax = breaks_vec[4], ymin = -2, ymax = Ntip(main_tree)+2.5), 
              fill = rgbfill, alpha = 0.20, inherit.aes = FALSE, linewidth = 0) +
    theme_tree2() +
    geom_tree(layout = "rectangular", 
              #ladderize=TRUE,
              right=TRUE,
              size=1.75,
              color="darkgrey")+
    geom_range(range = 'range', colour='red', size=1.75, alpha=0.22,
               position = position_nudge(x = -.5))+
    #geom_text(aes(x=branch, label=range), size = 2)+ #check range limits if needed
    #geom_text(aes(x=branch, label=node), size = 2)+ #check node numbers if needed
    #scale_color_continuous(low='blue', high='red')+ #Body mass
    geom_tiplab(size=3.0, offset = 0.5, fontface = "italic")+
    coord_geo(xlim=c(inferior_limit,Ntip(main_tree)+1), ylim=c(.5,Ntip(main_tree)), expand=T,
              dat = list("eras", "periods", "epochs"), abbrv = list(F, F, T), 
              skip = c("Pliocene", "Pleistocene", "Holocene","Quaternary"),
              pos = list("top", "top", "top"),
              alpha = 1, height = unit(0.75, "line"),
              size = "auto", fittext_args = list(size=22),
              rot = 0, neg = TRUE) +
    scale_x_continuous(breaks = breaks_vec,
                       labels = abs(breaks_vec),
                       expand = expansion(mult = c(0.009, 0.14)))
  # Reverse axis
  f3_rev <- revts(f3)
  # Extract data from ggtree reversed object
  tree_data <- f3_rev$data  
  # Filter non-tip nodes to plot median_ages
  non_tip_labels <- tree_data[!tree_data$isTip & tree_data$median_ages >= 6, ]
  # Plot
  f3_final <- f3_rev + 
    geom_text(data = non_tip_labels, 
              aes(x = branch, 
                  label = round(median_ages, 1)), 
              size = 2.75, 
              vjust = 1.5, 
              hjust = .5) +
    geom_cladelab(node = mrca_I,
                  label = "clade I",
                  align = T,
                  offset = 70,
                  offset.text = 2,
                  angle = -90,
                  fontsize = 3) +
    geom_cladelab(node = mrca_lep,
                  label = "clade leptopteroid",
                  align = F,
                  offset = 63,
                  offset.text = 2,
                  angle = -90,
                  fontsize = 3);f3_final
    
  # Save as SVG
  svglite("MAIN_DATEDTREE.svg", width = 15, height = 10)  # width/height in inches
  grid::grid.draw(ggplotGrob( f3_final ))
  dev.off()
  
  
  
    
  #my nice tree
  trial <- ggplot(main_tree, ladderize=F) +
    # Rectangles (unchanged)
    geom_rect(aes(xmin = breaks_vec[3], xmax = breaks_vec[2], ymin = -2, ymax = Ntip(main_tree)+2.5), 
              fill = rgbfill, alpha = 0.10, inherit.aes = FALSE, linewidth = 0) +
    geom_rect(aes(xmin = breaks_vec[5], xmax = breaks_vec[4], ymin = -2, ymax = Ntip(main_tree)+2.5), 
              fill = rgbfill, alpha = 0.10, inherit.aes = FALSE, linewidth = 0) +
    theme_tree2() +
    coord_geo(xlim=c(inferior_limit, Ntip(main_tree)), 
              ylim=c(.15, Ntip(main_tree)), 
              expand = T,  # Changed to FALSE to prevent automatic expansion
              dat = list("eras", "periods", "epochs"), 
              abbrv = list(F, F, T),
              skip = c("Pliocene", "Pleistocene", "Holocene","Quaternary"),
              pos = list("top", "top", "top"),
              alpha = 1, height = unit(0.75, "line"),
              size = "auto", fittext_args = list(size=22),
              rot = 0, neg = TRUE) +
    geom_tree(layout = "rectangular", 
              #ladderize=TRUE,
              right=TRUE,
              size=1.15,
              color="darkgrey")+
    geom_tiplab(size=2.8, offset = 0.5, fontface = "italic") +  # Reduced offset from 3 to 0.5
    geom_range(range = 'range', colour='red', size=1, alpha=0.12) +
    scale_x_continuous(breaks = breaks_vec,
                       labels = abs(breaks_vec),
                       expand = expansion(mult = c(0.02, 0.17))); revts(trial)

  aposttree <- main_tree
 main_tree@phylo$node.label <- NULL
 write.tree(main_tree@phylo, file = "exportado.nex") #This is to be used in PhyGeo

#2C. STEP - SENSITIVITY PLOT ####
 
 # Read comparisons
 comptax_k5  <- read.csv("~/gTNT/mydata/majo/res_jun/comptax_k5.csv",  comment.char = "#")
 comptax_k10 <- read.csv("~/gTNT/mydata/majo/res_jun/comptax_k10.csv", comment.char = "#")
 comptax_k15 <- read.csv("~/gTNT/mydata/majo/res_jun/comptax_k15.csv", comment.char = "#")
 comptax_ew  <- read.csv("~/gTNT/mydata/majo/res_jun/comptax_ew.csv",  comment.char = "#")
 
 # Rename first column to RowID
 colnames(comptax_k5)[1] <- "RowID"
 colnames(comptax_k10)[1] <- "RowID"
 colnames(comptax_k15)[1] <- "RowID"
 colnames(comptax_ew)[1] <- "RowID"
 
 # Combine all by RowID
 allcomptax <- comptax_k5 %>%
   left_join(comptax_k10, by = "RowID") %>%
   left_join(comptax_k15, by = "RowID") %>%
   left_join(comptax_ew,  by = "RowID")
 
 # Identify the columns that contain "A2"
 a2_cols <- grep("A2", colnames(allcomptax), value = TRUE)
 
 # Replace A2 values for "clade_I"
 allcomptax[allcomptax$RowID == "clade_I", a2_cols] <-
   allcomptax[allcomptax$RowID == "clade_I_No_Ossturii", a2_cols]
 
 # Replace A2 values for "Osmundaceae"
 allcomptax[allcomptax$RowID == "Osmundaceae", a2_cols] <-
   allcomptax[allcomptax$RowID == "Osmundaceae_No_Ossturii", a2_cols]
 
 # Replace A2 values for "Osmundopsis"
 allcomptax[allcomptax$RowID == "Osmundopsis", a2_cols] <-
   allcomptax[allcomptax$RowID == "Os_rafaelii_Os_zunigai", a2_cols]
 
 # Discard useless rows 
 leftme <- c("Os_rafaelii_Os_zunigai", "Osmundaceae_No_Ossturii", "clade_I_No_Ossturii")
 allcomptax <- allcomptax %>%
   filter(!RowID %in% leftme)
 
 #Rename
 allcomptax$RowID[1] <- "clade I"
 allcomptax$RowID[2] <- "clade leptopteroid"
 
 # Reshape to long format
 df_long_allcomptax <- allcomptax %>%
   pivot_longer(
     cols = -RowID,
     names_to = "Variable",
     values_to = "Value"
   ) %>%
   mutate(
     Category = str_extract(Variable, "k5|k10|k15|ew"),
     Condition = str_extract(Variable, "A\\d+")
   )
 
 # Plot heatmap
 sensitivy_plot <- ggplot(df_long_allcomptax, aes(x = Condition, y = fct_rev(RowID), fill = factor(Value))) +
   geom_tile(color = "white") +
   facet_wrap(~Category, ncol = 2) +
   scale_fill_manual(values = c("0" = "white", "1" = "steelblue")) +
   guides(fill = "none") +
   labs(fill = "Value", x = "", y = "", title = "Sensitivity plot") +
   theme(
     plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
     axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
     axis.text.y = element_text(size = 9, face = "italic"),
     strip.text = element_text(size = 8, face = "bold")
   )
 
#3.  STEP (OLD RUDIMENTARY VERSION IT WORKS BUT IT IS ODD) - NOW LTT - ####
  
  # Now extracting the number of lineages through time from the dated phylogenies #
  lttmcmctree <- ltt.plot.coords(aposttree@phylo)
  lttapost <- ltt.plot.coords(aposttree@phylo)
  
  # Set the time to a positive scale #
  lttmcmctree[,1] <- lttmcmctree[,1] * (-1)  
  lttapost[,1] <- lttapost[,1] * (-1)
  
  # Remove rows where the first column (time) is > 390 or < 12 Ma
  lttmcmctree <- lttmcmctree[lttmcmctree[,1] <= 390 & lttmcmctree[,1] > 0, ]
  lttapost <- lttapost[lttapost[,1] <= 390 & lttapost[,1] > 0, ]
  
  # Find the maximum number of rows
  max_rows <- max(nrow(lttmcmctree), nrow(lttapost))
  
  # Extend both matrices with NA to match row lengths
  lttmcmctree_extended <- rbind(lttmcmctree, matrix(NA, nrow = max_rows - nrow(lttmcmctree), ncol = ncol(lttmcmctree)))
  lttapost_extended <- rbind(lttapost, matrix(NA, nrow = max_rows - nrow(lttapost), ncol = ncol(lttapost)))
  
  # Rename columns of each matrix to avoid conflicts
  colnames(lttmcmctree_extended)[1] <- "mcmctree_time"
  colnames(lttmcmctree_extended)[2] <- "mcmctree_ltt"
  colnames(lttapost_extended)[1] <- "apost_time"
  colnames(lttapost_extended)[2] <- "apost_ltt"
  
  # Merge the matrices column-wise
  merged_ltt <- cbind(lttmcmctree_extended, lttapost_extended)
  
  # Read paleoclimate data from Judd et al. (2024): 10.1126/science.adk3705
  paleodata <- read.csv("PhanDA_GMSTandCO2_percentiles.csv", header = TRUE, sep = ",")
  paleodata <- paleodata[, (names(paleodata) %in% c("AverageAge", "GMST_05", "GMST_50", "GMST_95"))]
  colnames(paleodata)[1] <- "clim_time"
  
  # Combining with merged_ltt
  #This works: combined_paleodata_ltt <- NULL
  #This works: combined_paleodata_ltt <- merge(merged_ltt, paleodata, all = TRUE)
  combined_paleodata_ltt <- merge(merged_ltt, paleodata, all = TRUE) %>% 
    drop_na()
  
  # Reshape data for being plotted
  paleoltt_long_data <- combined_paleodata_ltt %>%
    pivot_longer(
      cols = c(mcmctree_ltt, apost_ltt, GMST_50),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(time = case_when(
      variable == "mcmctree_ltt" ~ mcmctree_time,
      variable == "apost_ltt" ~ apost_time,
      variable == "GMST_50" ~ clim_time
    ))
  
  # Create a separate dataset for geom_ribbon
  paleo_ribbon_data <- combined_paleodata_ltt %>%
    select(clim_time, GMST_05, GMST_50, GMST_95) %>%
    drop_na()
  
  # Define a function to rescale GMST values relative to LTT values
  gmst_rescale <- function(x) {(x - min(paleoltt_long_data$value, na.rm = TRUE)) / 
      (max(paleoltt_long_data$value, na.rm = TRUE) - min(paleoltt_long_data$value, na.rm = TRUE)) *
      (max(paleodata$GMST_50, na.rm = TRUE) - min(paleodata$GMST_50, na.rm = TRUE)) + 
      min(paleodata$GMST_50, na.rm = TRUE)}

#3B. STEP (OLD RUDIMENTARY VERSION IT WORKS BUT IT IS ODD) - LTT PLOTS - ####
  # Plot the data
  LTT <- ggplot(paleoltt_long_data, aes(x = time, y = value, color = variable)) +
    geom_line(linewidth = 1) + 
    geom_ribbon(data = paleo_ribbon_data, aes(x = clim_time, ymin = GMST_05, ymax = GMST_95), 
                fill = "#ffb09c", alpha = 0.45, inherit.aes = FALSE) +  # Adjust ribbon color
    scale_color_manual(values = c("GMST_50" = "#ee2400",  # Red for GMST_50
                                  "GMST_05" = "#ffb09c",  # Light red for GMST_05
                                  "GMST_95" = "#ffb09c",  # Light red for GMST_95
                                  "mcmctree_ltt" = "darkgreen",  # Green for MCMC LTT
                                  "apost_ltt" = "#00008B")) +  # Blue for Apost LTT
    labs(x = "Geological Time (Ma)",
         y = "Lineages Through Time") +
    scale_y_continuous(sec.axis = sec_axis(~ gmst_rescale(.), name = "Global Mean Surface Temperature (°C)")) +
    theme(legend.position = "none")+
    xlim(5, 350) +  
    coord_geo(xlim = c(350, 5), 
              size = "auto", fittext_args = list(size=15),
              dat = list("eras", "periods", "epoch"),
              pos = list("b", "b", "b"), abbrv = list(F, F, F),
              height = unit(0.8, "line")) + 
    scale_x_continuous(breaks = seq(350, 0, by = -25))+
    scale_x_reverse(); LTT
  
  
  # Alternative plot
  Alt_LTT <- ggplot(paleoltt_long_data, aes(x = time, y = value, color = variable)) +
    geom_rect(xmin = 25, xmax = 65, ymin = -Inf, ymax = Inf, fill="lightgrey")+  
    geom_line(linewidth = 1) + 
    geom_ribbon(data = paleo_ribbon_data, aes(x = clim_time, ymin = GMST_05, ymax = GMST_95), 
                fill = "#ffb09c", alpha = 0.45, inherit.aes = FALSE) +  # Adjust ribbon color
    scale_color_manual(values = c("GMST_50" = "#ee2400",  # Red for GMST_50
                                  "GMST_05" = "#ffb09c",  # Light red for GMST_05
                                  "GMST_95" = "#ffb09c",  # Light red for GMST_95
                                  "mcmctree_ltt" = "darkgreen",  # Green for MCMC LTT
                                  "apost_ltt" = "#00008B")) +  # Blue for Apost LTT
    labs(x = "Geological Time (Ma)",
         y = "Lineages Through Time") +
    scale_y_continuous(
      sec.axis = sec_axis(~ gmst_rescale(.), name = "Global Mean Surface Temperature (°C)")
    ) + 
    scale_x_reverse(limits = c(300, 5), breaks = seq(300, 5, by = -25)) + # Set limits and breaks
    theme_minimal() +
    #coord_cartesian(xlim = c(125, 5)) +
    #scale_x_continuous(breaks = seq(125, 0, by = -25)) +
    theme(legend.position = "none");Alt_LTT
  
  
  # Alternative plot
  Alt_LTT <- ggplot(paleoltt_long_data, aes(x = time, y = value, color = variable)) +
    geom_rect(aes(xmin = (breaks_vec[3]*(-1)), xmax = (breaks_vec[2]*(-1)), ymin = -Inf, ymax = Inf), 
              fill = rgbfill, alpha = 0.10, inherit.aes = FALSE, linewidth = 0) +
    geom_rect(aes(xmin = (breaks_vec[5]*(-1)), xmax = (breaks_vec[4]*(-1)), ymin = -Inf, ymax = Inf), 
              fill = rgbfill, alpha = 0.10, inherit.aes = FALSE, linewidth = 0) +  
    geom_line(linewidth = 1) + 
    geom_ribbon(data = paleo_ribbon_data, aes(x = clim_time, ymin = GMST_05, ymax = GMST_95), 
                fill = "#ffb09c", alpha = 0.45, inherit.aes = FALSE) +
    scale_color_manual(values = c("GMST_50" = "#ee2400",
                                  "GMST_05" = "#ffb09c",
                                  "GMST_95" = "#ffb09c",
                                  "mcmctree_ltt" = "darkgreen",
                                  "apost_ltt" = "#00008B")) +
    labs(x = "Geological Time (Ma)",
         y = "Lineages Through Time") +
    scale_y_continuous(
      sec.axis = sec_axis(~ gmst_rescale(.), name = "GMST (°C)")
    ) + 
    scale_x_reverse(limits = c(325, 2),
                    breaks = c(325, 193, 175, 34, 30, 2)) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 8)  #  adjust this value
    );Alt_LTT
  
  

#4.  STEP NEW LTT SECTION (see function above)####
  
 my_ltt <- ltt_and_clime(dated_trees)
 
 #Data is in my_ltt$lttclime and my_ltt$clim_ribbon_data
 rescale_fun <- gmst_rescaler(my_ltt$lttclime, my_ltt$clim_ribbon_data)
 
 #Plot LTT and paleotemperature
 ltt_plot <- ggplot(my_ltt$lttclime, aes(x = time, y = value, color = variable)) +
   # coord_geo(size = "auto", fittext_args = list(size=15),
   #           dat = list("eras", "periods", "epoch"),
   #           pos = list("t", "t", "t"),
   #           abbrv = list(F, F, T),
   #           height = unit(0.8, "line")
   # ) +
   geom_rect(aes(xmin = (breaks_vec[3]*(-1)), xmax = (breaks_vec[2]*(-1)), ymin = -Inf, ymax = Inf), alpha=0.25, 
             fill = rgbfill, alpha = 0.5, inherit.aes = FALSE, linewidth = 0) +
   geom_rect(aes(xmin = (breaks_vec[5]*(-1)), xmax = (breaks_vec[4]*(-1)), ymin = -Inf, ymax = Inf), alpha=0.25, 
             fill = rgbfill, alpha = 0.5, inherit.aes = FALSE, linewidth = 0) +
   geom_line(linewidth = 1) + 
   geom_ribbon(data = my_ltt$clim_ribbon_data, aes(x = clim_time, ymin = GMST_05, ymax = GMST_95), 
               fill = "#ffb09c", alpha = 0.45, inherit.aes = FALSE) +  # Adjust ribbon color
   geom_ribbon(data = my_ltt$ltt_ribbon_data, aes(x = median.time, ymin = ltt.lower_95, ymax = ltt.upper_95), 
               fill = "#6f7ab7", alpha = 0.45, inherit.aes = FALSE) +  # Adjust ribbon color
   scale_color_manual(values = c("GMST_50" = "#ee2400",  # Red for GMST_50
                                 "GMST_05" = "#ffb09c",  # Light red for GMST_05
                                 "GMST_95" = "#ffb09c",  # Light red for GMST_95
                                 "median.ltt" = "blue"  # Green for MCMC LTT
                                 )) + 
   labs(x = "Geological Time (Ma)",
        y = "LTT",
        title = "Lineages Through Time") +
   scale_y_continuous( sec.axis = sec_axis(~ rescale_fun(.), name = "GMST (°C)")) +
   theme_minimal()+
   theme(legend.position = "none",
         plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
         axis.title.x = element_text(size = 9),
         axis.title.y = element_text(size = 9),
         axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
         axis.text.y = element_text(size = 9))+
   scale_x_reverse(limits = c(350, 0),
                   breaks = c(350, 193, 175, 49, 16, 0)
                   )
 # Save as SVG
 svglite("ltt_clime.svg", width = 10, height = 6)  # width/height in inches
 grid::grid.draw(ggplotGrob(ltt_plot))
 dev.off()

#5.  STEP (OPTIONAL) Alternative with ltt phytools ####
    secltt <- ltt(main_tree@phylo)
    plot(secltt,show.tree=TRUE,lwd=2, 
         log.lineages=F,log="y",bty="n",cex.lab=0.9, 
         transparency=0.1,axes=FALSE, 
         xlab="millions of year bp")
    h<-max(nodeHeights(main_tree@phylo)) 
    axis(1,at=h-seq(0,350,by=5),labels=seq(0,350,by=5),las=1, 
         cex.axis=0.8) 
    axis(2,las=1,cex.axis=0.8)

#6.  STEP (OPTIONAL) PLOTTING DATED TREE AND SENS_PLOT AND LTT PLOT TOGETHER ####
    sensitivity_plot <- sensitivity_plot +
      theme(
        plot.title = element_text(margin = margin(b = 5)),  # reduce bottom margin
        axis.title.x = element_text(margin = margin(t = 5)),  # reduce top margin
        axis.title.y = element_text(margin = margin(r = 5))
      )
    
    ltt_plot <- ltt_plot +
      theme(
        plot.title = element_text(margin = margin(b = 5)),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.title.y = element_text(margin = margin(r = 5))
      )
    
    left_column <- sensitivity_plot / ltt_plot + plot_layout(heights = c(.6,.4))
    final_plot <- (left_column | f3_final) + 
      plot_layout(widths = c(0.4, 1.25)); final_plot
    
    ggsave("Figure_2.svg", plot = final_plot, width = 16, height = 8, device = svglite)
    
#7.  STEP (OPTIONAL) PLOTTING TREES FOR SUPPLEMENTARY MATERIAL ####
    #fsize
    fsize <- 3.15
    #Read all trees
    pi5 <- c("pi5_Arbroles1.tre", "pi5_Arbroles2.tre", "pi5_Arbroles3.tre",
             "pi5_Arbroles4.tre", "pi5_Arbroles5.tre", "pi5_Arbroles6.tre")
    pi10 <- c("pi10_Arbroles1.tre", "pi10_Arbroles2.tre", "pi10_Arbroles3.tre",
             "pi10_Arbroles4.tre", "pi10_Arbroles5.tre", "pi10_Arbroles6.tre")
    pi15 <- c("pi15_Arbroles1.tre", "pi15_Arbroles2.tre", "pi15_Arbroles3.tre",
             "pi15_Arbroles4.tre", "pi15_Arbroles5.tre", "pi15_Arbroles6.tre")
    ew <- c("Arbroles1.tre", "Arbroles2.tre", "Arbroles3.tre",
             "Arbroles4.tre", "Arbroles5.tre", "Arbroles6.tre")
    
    pi5_trees  <- vector("list", 6) ; pi5_ggtrees  <- vector("list", 6)
      for (i in 1:6) {
      titleme <- paste("A",i,sep = "")
      pi5_trees[[i]] <- ReadTntTree(pi5[[i]])
      if(inherits(pi5_trees[[i]], "multiPhylo")) { pi5_trees[[i]] <- consensus(pi5_trees[[i]]) } 
      pi5_trees[[i]] <- root.phylo(pi5_trees[[i]], outgroup = "Dipteris_conjugata")
      
      # Abbreviate tip labels
      pi5_trees[[i]]$tip.label <- gsub("^([A-Za-z])[a-z]+_", "\\1_", pi5_trees[[i]]$tip.label)
      
      # Now plot
      pi5_ggtrees[[i]] <- ggtree(pi5_trees[[i]], ladderize = F, size = 1.4, color="darkgrey") +
        geom_tiplab(size = fsize, offset = 0.0, fontface = "italic") +
        xlim(0, 20) +
        labs(title = titleme) +
        theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5))
      }
    
    upperrow <- pi5_ggtrees[[1]] | pi5_ggtrees[[2]] | pi5_ggtrees[[3]]
    lowerrow <- pi5_ggtrees[[4]] | pi5_ggtrees[[5]] | pi5_ggtrees[[6]]
    patchtitle <- "k5"
    FigS1_plot <- (upperrow / lowerrow)  + plot_annotation(title = patchtitle) & theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    FigS1_plot
    ggsave("FigS1.svg", plot = FigS1_plot, width = 10, height = 15, device = svglite)
    
    pi10_trees <- vector("list", 6); pi10_ggtrees <- vector("list", 6)
      for (i in 1:6) {
        titleme <- paste("A",i,sep = "")
        pi10_trees[[i]] <- ReadTntTree(pi10[[i]])
        if(inherits(pi10_trees[[i]], "multiPhylo")) { pi10_trees[[i]] <- consensus(pi10_trees[[i]]) } 
        pi10_trees[[i]] <- root.phylo(pi10_trees[[i]], outgroup = "Dipteris_conjugata")

        # Abbreviate tip labels
        pi10_trees[[i]]$tip.label <- gsub("^([A-Za-z])[a-z]+_", "\\1_", pi10_trees[[i]]$tip.label)
        
        # Now plot
        pi10_ggtrees[[i]] <- ggtree(pi10_trees[[i]], ladderize = F, size = 1.4, color="darkgrey") +
          geom_tiplab(size = fsize, offset = 0.0, fontface = "italic") +
          xlim(0, 20) +
          labs(title = titleme) +
          theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5))
        
      }
    upperrow <- pi10_ggtrees[[1]] | pi10_ggtrees[[2]] | pi10_ggtrees[[3]]
    lowerrow <- pi10_ggtrees[[4]] | pi10_ggtrees[[5]] | pi10_ggtrees[[6]]
    patchtitle <- "k10"
    FigS2_plot <- (upperrow / lowerrow) + plot_annotation(title = patchtitle) & theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    FigS2_plot
    ggsave("FigS2.svg", plot = FigS2_plot, width = 10, height = 15, device = svglite)
    
    
    pi15_trees <- vector("list", 6); pi15_ggtrees <- vector("list", 6)
      for (i in 1:6) {
        titleme <- paste("A",i,sep = "")
        pi15_trees[[i]] <- ReadTntTree(pi15[[i]])
        if(inherits(pi15_trees[[i]], "multiPhylo")) { pi15_trees[[i]] <- consensus(pi15_trees[[i]]) } 
        pi15_trees[[i]] <- root.phylo(pi15_trees[[i]], outgroup = "Dipteris_conjugata")
        
        # Abbreviate tip labels
        pi15_trees[[i]]$tip.label <- gsub("^([A-Za-z])[a-z]+_", "\\1_", pi15_trees[[i]]$tip.label)
        
        # Now plot
        pi15_ggtrees[[i]] <- ggtree(pi15_trees[[i]], ladderize = F, size = 1.4, color="darkgrey") +
          geom_tiplab(size = fsize, offset = 0.0, fontface = "italic") +
          xlim(0, 20) +
          labs(title = titleme) +
          theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5))
      }
    upperrow <- pi15_ggtrees[[1]] | pi15_ggtrees[[2]] | pi15_ggtrees[[3]]
    lowerrow <- pi15_ggtrees[[4]] | pi15_ggtrees[[5]] | pi15_ggtrees[[6]]
    patchtitle <- "k15"
    FigS3_plot <- (upperrow / lowerrow) + plot_annotation(title = patchtitle) & theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    FigS3_plot
    ggsave("FigS3.svg", plot = FigS3_plot, width = 10, height = 15, device = svglite)
    
        
    ew_trees   <- vector("list", 6); ew_ggtrees <- vector("list", 6)
      for (i in 1:6) {
        titleme <- paste("A",i,sep = "")
        ew_trees[[i]] <- ReadTntTree(ew[[i]])
        if(inherits(ew_trees[[i]], "multiPhylo")) { ew_trees[[i]] <- consensus(ew_trees[[i]]) } 
        ew_trees[[i]] <- root.phylo(ew_trees[[i]], outgroup = "Dipteris_conjugata")
        
        # Abbreviate tip labels
        ew_trees[[i]]$tip.label <- gsub("^([A-Za-z])[a-z]+_", "\\1_", ew_trees[[i]]$tip.label)
        
        # plot
        ew_ggtrees[[i]] <- ggtree(ew_trees[[i]], ladderize = F, size = 1.4, color="darkgrey") +
          geom_tiplab(size = fsize, offset = 0.0, fontface = "italic") +
          xlim(0, 20) +
          labs(title = titleme) +
          theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5))
        
      }
    upperrow <- ew_ggtrees[[1]] | ew_ggtrees[[2]] | ew_ggtrees[[3]]
    lowerrow <- ew_ggtrees[[4]] | ew_ggtrees[[5]] | ew_ggtrees[[6]]
    patchtitle <- "ew"
    FigS4_plot <- (upperrow / lowerrow) + plot_annotation(title = patchtitle) & theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    FigS4_plot
    ggsave("FigS4.svg", plot = FigS4_plot, width = 10, height = 15, device = svglite)
