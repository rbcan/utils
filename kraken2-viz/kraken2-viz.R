#' ---
#' title: "kraken2 Report"
#' author: "Andrea Telatin & Rebecca Ansorge"
#' date: "15/01/2021"
#' params:
#'     inbracken: ""
#'     inkraken: ""
#'     topspecies: ""
#'     outdir: ""
#' output: 
#'     html_document:
#'         fig_width: 10
#'         fig_height: 5
#'         toc: true
#'         number_sections: false
#'         toc_float: true
#'         toc_depth: 3
#' editor_options: 
#'     chunk_output_type: console
#' ---

#+ include=FALSE
####### FIRST THINGS FIRST #######
library("tidyverse")
library("plotly") 
library("RColorBrewer")
library("vegan") #
library("ggdendro") # 
library("htmltools")
library("rmarkdown")

#if (pandoc_available())
#  cat("pandoc", as.character(pandoc_version()), "is available!\n")
#Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")

theme_set(theme_bw())

# input arguments
k2In<-params$inkraken
bIn<-params$inbracken
topS<-params$topspecies
outDir<-params$outdir

#k2In<-"kraken2-merged.report"
#bIn<-"bracken-merged.report"
#topS<-10
#outDir<-"./"

##################################



### FILE IMPORT AND FORMATTING ###

# import files
b_data <- read_delim(bIn, delim="\t", comment= "#") %>%
  mutate(name = gsub("^\\s+","",name)) %>%
  filter(rowSums(.[,1:(ncol(.)-3)])>0) # filtering out rows where counts are 0 for all samples

k_data <- read_delim(k2In, delim="\t", comment= "#") %>%
  mutate(name = gsub("^\\s+","",name)) %>%
  filter(lvl_type == "R" | lvl_type == "U") %>%
  gather(key="sample",value="reads", -c("lvl_type","taxid","name"))


# taxa level categories
tax_lvl <- c("D","P","C","O","F","G","S")
tax_lvl2 <- c("Domain","Phylum","Class","Order","Family","Genus","Species")

# split bracken table into one dataframe per taxonomic level
taxlist<- lapply(tax_lvl, function(i) {
  tmp <- b_data %>%
    subset(lvl_type==i)
  write_delim(tmp,paste(outDir,"/",i,"_abundances.tsv",sep=""),delim="\t")
  tmp<- tmp %>%
    gather(key="sample",value="reads", -c("lvl_type","taxid","name"))
  tmp 
})

##################################

##### TOP TAXA DETERMINATION #####
# pre-formatting
l <- length(names(b_data))
samples <- names(b_data)[1:(l-3)]
header <- names(b_data)[(l-2):l]

# create table for top taxa
# df<-b_data[header] %>%
#   subset(lvl_type=="S")
# tmp<-data.frame("S","0000","a_other")
# names(tmp)<-header

# make top 10 taxa dataframe and sum all 'others'
topS_taxids<-data.frame()
for (i in samples) {
x<-b_data %>%
  subset(lvl_type=="S") %>%
  select(c(i,"taxid")) %>%
  arrange(desc(.[,1])) %>%
  slice(1:topS) %>%
  na_if(0) %>%
  select("taxid")
topS_taxids<-rbind(topS_taxids,x)
}
topS_taxids<-unique(topS_taxids)

toptaxa<-b_data %>% 
  subset(lvl_type=="S") %>%
  right_join(topS_taxids, by="taxid")

other<-b_data %>% 
  subset(lvl_type=="S") %>%
  anti_join(topS_taxids, by="taxid") 
otherS<-colSums(other[,1:length(samples)]) %>%
  as.data.frame() %>%
  t() %>%
  as_tibble() %>%
  mutate(lvl_type="S") %>%
  mutate(taxid=0) %>%
  mutate(name="a_other")

df<-rbind(toptaxa,otherS)

# df<-rbind(df,tmp)
# for (i in samples) {
#   x <- b_data %>%
#     subset(lvl_type=="S") %>%
#     select(c(i,"taxid")) %>%
#     arrange(desc(.[,1])) %>%
#     slice(1:topS) %>%
#     na_if(0)
#   y <- b_data %>%
#     subset(lvl_type=="S") %>%
#     select(c(i,"taxid")) %>%
#     arrange(desc(.[,1])) %>%
#     slice((topS+1):n()) %>%
#     summarize(sum=sum(.[,1])) %>%
#     rename_with(.fn = ~paste0(i)) %>%
#     mutate(taxid="0000") 
#   print(head(x))
#   print(head(y))
#   z<-rbind(x,y)
#   print(tail(z))
#   df<-full_join(df,z,by="taxid")
# }

# remove species that were not in top taxa in any sample - entire row has NA
topspec <- df %>% filter_at(vars(all_of(samples)),any_vars(!is.na(.)))

# make top species relative 
#topspec2<-topspec[,(4:length(topspec))]
topspec2<-topspec[,1:length(samples)]
topspec2[is.na(topspec2)] <- 0
topspec2 <- apply(topspec2,2,function(x){x/sum(x)})
topspec3 <- cbind((topspec[(length(samples)+1):ncol(topspec)]),topspec2)
tops2 <-topspec3 %>% 
  gather(key="sample",value="reads", -c("lvl_type","taxid","name"))
##################################


###### PREPARE FOR PLOTTING ######

# set dimentsions of plots
panel<-((nrow(topspec)*25))*1
if (panel<350){panel<-350}
##################################


######## PLOT UNCLASSIFIED #######

#' ##  1. Unclassified reads {.tabset}

#' Count of classified and unclassified reads for raw **kraken2** output

#+ echo=FALSE
uncl <- ggplot(k_data) +
  geom_bar(stat="identity",aes(x=sample,y=reads, fill=name)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  scale_fill_brewer(palette="Set2") +
  ggtitle("Number of unclassified reads (raw counts)") 

uncl2 <- ggplot(k_data) +
  geom_bar(stat="identity",aes(x=sample,y=reads, fill=name), position="fill") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  scale_fill_brewer(palette="Set2") +
  ggtitle("Number of unclassified reads (relative coutns)")

#' ### absolute counts
#+ echo=FALSE
p_uncl<-ggplotly(uncl)
p_uncl

#' ### relative
#+ echo=FALSE
p_uncl2<-ggplotly(uncl2)
p_uncl2
##################################


######## PLOT TOP SPECIES ########

#' ## 2. Top species

#' Retrieved from **bracken** re-estimated read counts 

#' The top 10 species per sample are shown, the remaining species in a sample were summed into the "a_other" category
#'

#' ### Bar plot {.tabset}

#' #### absolute counts
#+ echo=FALSE
tops <-topspec %>% 
  gather(key="sample",value="reads", -c("lvl_type","taxid","name"))

ptop <- ggplot(tops) +
  geom_bar(stat="identity",aes(x=sample,y=reads, fill=name)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  ggtitle(paste("Top", topS,"species per sample"))
p_ptop<-ggplotly(ptop)
p_ptop

#' #### relative
#+ echo=FALSE
ptop2 <- ggplot(tops) +
  geom_bar(stat="identity",aes(x=sample,y=reads, fill=name), position="fill") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  ggtitle(paste("Top", topS,"species per sample"))
p_ptop2<-ggplotly(ptop2)
p_ptop2

#' ### Bubble plot {.tabset}

#' #### absolute counts
#+ echo=FALSE
ptop3 <- ggplot(tops) +
  geom_point(aes(x=sample,y=name, size=reads, color=reads), alpha = 0.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_size_continuous(limits = c(0.00001,max(tops$reads)))+
  ggtitle(paste("Top", topS,"species per sample"))
p_ptop3<-ggplotly(ptop3,height = panel)
p_ptop3

#' #### relative
#+ echo=FALSE
ptop4 <- ggplot(tops2) +
  geom_point(aes(x=sample,y=name, size=reads, color=reads), alpha = 0.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_size_continuous(limits = c(0.00001,max(tops2$reads)))+
  ggtitle(paste("Top", topS,"species per sample"))
p_ptop4<-ggplotly(ptop4,height = panel)
p_ptop4
##################################


######## PLOT NMDS #######

#' ## 3. Beta-diversity

#' NMDS on bray-curtis dissimilarites of relative reads counts. 
#' 
#' **NOTE: The stress value is displayed below - when stress is > 0.2, the clustering might be unreliable**

#+ echo=FALSE, results='hide'
set.seed(200)

speclist <- b_data %>%
  subset(lvl_type=="S") %>%
  select(-c("lvl_type","taxid","name"))
speclist <- apply(speclist,2,function(x){x/sum(x)})

o_nmds <- NULL
try(invisible(o_nmds <-metaMDS(speclist, distance="bray", k=3)), silent=TRUE)

#+ echo=FALSE
if (is(o_nmds,"metaMDS")) {
  o_nmds_data<-as.data.frame(o_nmds$species) %>%
    rownames_to_column("sample")
  
  # #+ echo=FALSE
  paste("NMDS stress:",round(o_nmds$stress,4))
  
  p_nmds <- plot_ly(o_nmds_data, x = ~MDS1, y = ~MDS2, z = ~MDS3,
                    type="scatter3d", mode="text", 
                    text = ~sample)
  p_nmds <- p_nmds %>% add_markers()
  p_nmds <- p_nmds %>% layout(scene = list(xaxis = list(title = 'MDS1'),
                                           yaxis = list(title = 'MDS2'),
                                           zaxis = list(title = 'MDS3')))
  p_nmds
} else {
  print("Caught an error during NMDS generation - skipping ordination")
}
##################################


######## PLOT DENDROGRAM #########

#' ## 4. Dendrogram

#' Relative species abundances used to perform hierarchical clustering of samples using 'euclidean' distances

#+ echo=FALSE
dhc <- hclust(dist(t(speclist)))

# Rectangular lines
ddata <- dendro_data(dhc, type = "rectangle")
p_dendro <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = ddata$labels, aes(x, y, label =label), hjust=0,size = 3)+
  ylim(max(ddata$segments$y),-((max(ddata$segments$y)/5))) +
  coord_flip() +
  theme(panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
p_dendro
##################################


### PLOT TAX LEVEK ABUNDANCES ####

#' ## 5. Taxa abundances
#' 
#' Read count are shown for each taxonomic level. Kraken2 classification was re-estimated using **bracken**.

#+ echo=FALSE
plotlist<- htmltools::tagList(lapply(seq_along(taxlist), function(ii) {
  tmp2 <- ggplot(taxlist[[ii]]) +
    geom_bar(stat="identity",aes(x=sample,y=reads, fill=name)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank()) +
    ggtitle(tax_lvl2[ii])
  as_widget(ggplotly(tmp2))
}))

plotlist2<- htmltools::tagList(lapply(seq_along(taxlist), function(ii) {
  tmp3 <- ggplot(taxlist[[ii]]) +
    geom_bar(stat="identity",aes(x=sample,y=reads, fill=name), position="fill") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank()) +
    ggtitle(tax_lvl2[ii])
  as_widget(ggplotly(tmp3))
}))

#' ### Domain {.tabset}
#' #### absolute counts
#+ echo=FALSE
plotlist[[1]][[1]]
#' #### relative
#+ echo=FALSE
plotlist2[[1]][[1]]

#' ### Phylum {.tabset}
#' #### absolute counts
#+ echo=FALSE
plotlist[[1]][[2]]
#' #### relative
#+ echo=FALSE
plotlist2[[1]][[2]]

#' ### Class {.tabset}
#' #### absolute counts
#+ echo=FALSE
plotlist[[1]][[3]]
#' #### relative
#+ echo=FALSE
plotlist2[[1]][[3]]

#' ### Order {.tabset}
#' #### absolute counts
#+ echo=FALSE
plotlist[[1]][[4]]
#' #### relative
#+ echo=FALSE
plotlist2[[1]][[4]]

#' ### Family {.tabset}
#' #### absolute counts
#+ echo=FALSE
plotlist[[1]][[5]]
#' #### relative
#+ echo=FALSE
plotlist2[[1]][[5]]

#' ### Genus {.tabset}
#' #### absolute counts
#+ echo=FALSE
plotlist[[1]][[6]]
#' #### relative
#+ echo=FALSE
plotlist2[[1]][[6]]

#' ### Species {.tabset}
#' #### absolute counts
#+ echo=FALSE
plotlist[[1]][[7]]
#' #### relative
#+ echo=FALSE
plotlist2[[1]][[7]]
##################################


##### EXPORT PLOT HTML FILES #####

#+ echo=FALSE

htmlwidgets::saveWidget(as_widget(p_uncl),paste(outDir,"/p_uncl.html",sep=""))
htmlwidgets::saveWidget(as_widget(p_uncl2),paste(outDir,"/p_uncl2.html",sep=""))
htmlwidgets::saveWidget(as_widget(p_ptop),paste(outDir,"/p_ptop.html",sep=""))
htmlwidgets::saveWidget(as_widget(p_ptop2),paste(outDir,"/p_ptop2.html",sep=""))
htmlwidgets::saveWidget(as_widget(p_ptop3),paste(outDir,"/p_ptop3.html",sep=""))
htmlwidgets::saveWidget(as_widget(p_ptop4),paste(outDir,"/p_ptop4.html",sep=""))
try(htmlwidgets::saveWidget(as_widget(p_nmds),paste(outDir,"/p_nmds.html",sep="")),silent=TRUE)
ggsave(paste(outDir,"/p_dendro.png",sep=""),p_dendro)
##################################


#### EXPORT PLOT PDF FILES ####

pdf(paste(outDir,"/high_res_figs.pdf",sep=""))
uncl
uncl2
ptop
ptop2
ptop3
ptop4
p_dendro
dev.off()
##################################

####### EXPORT DATA TABLES #######
# top 10 species per sample (counts)
write_delim(topspec,paste(outDir,"/top_species.tsv",sep=""),delim="\t")
# top 10 species per sample (relative)
write_delim(as.data.frame(topspec2),paste(outDir,"/top_species_relative.tsv",sep=""),delim="\t")
# relative abundances of species matrix (basis for NMDS and dendrogram)
write_delim(as.data.frame(speclist),paste(outDir,"/species_relabundance_matrix.tsv",sep=""),delim="\t")
# NMDS coordinates
try(write_delim(o_nmds_data,paste(outDir,"/NMDS_data.tsv",sep=""),delim="\t"),silent=TRUE)

##################################


