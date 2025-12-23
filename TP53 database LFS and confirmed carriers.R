#download CSV of LFS patients and mutations from TP53 database to computer = The TP53 Database (R20, July 2019): https://tp53.isb-cgc.org/view_data?bq_view_name=GermlineDownload with filters as below
#how to get to download page --> https://tp53.isb-cgc.org/ --> Germline Variants --> Data Downloads Germline Variants --> TP53 germline variants and family history (data file) --> click Preview icon to get to filterable data table
#type in filters into columns before download - LFS class, confirmed Germline_carrier; 654 entries left after filtering & were downloaded

#set seed for reproducibility
set.seed(1212)

#load (and if necessary install) important packages that will be used henceforth
pacman::p_load_gh("jokergoo/ComplexHeatmap")
pacman::p_load(patchwork, dplyr, data.table, scales, grid, gridExtra, cowplot, ggplot2, janitor, pivottabler, tibble, plotly, RColorBrewer, ComplexHeatmap, ggalluvial)

#read in downloaded CSV file into R as dataframe
germlinedataorig <- read.csv("C:/Users/nwali/Downloads/GermlineDownload_r20.csv", sep = ",", header = TRUE) %>% remove_rownames()

#separate data into LFS patients without or with cancers at follow-up based on unaffected column --> false means patient has cancer and true means patient does not have cancer
germlinenocancers <- germlinedataorig[germlinedataorig$Unaffected == "true",] %>% remove_rownames()
germlinejustcancers <- germlinedataorig[germlinedataorig$Unaffected == "false",] %>% remove_rownames()

#remove data with unknown site or other sites short topology for cancer
germlinejustcancers <- germlinejustcancers %>% filter(!Short_topo %in% c("OTHER SITES", "UNKNOWN SITE")) %>% remove_rownames()
germlinenocancers <- germlinenocancers %>% filter(!Short_topo %in% c("OTHER SITES", "UNKNOWN SITE")) %>% remove_rownames()

#concatenate codon number and mut type to p.? mutations which have unspecified residues
germlinejustcancers <- germlinejustcancers %>% mutate(ProtDescription = ifelse(ProtDescription %in% c("p.?"), paste(ProtDescription, Codon_number, Effect, sep = "-"), no = ProtDescription))
germlinenocancers <- germlinenocancers %>% mutate(ProtDescription = ifelse(ProtDescription %in% c("p.?"), paste(ProtDescription, Codon_number, Effect, sep = "-"), no = ProtDescription))

#remove p. prefix from all mutations to clean up visually
germlinejustcancers$ProtDescription <- substring(germlinejustcancers$ProtDescription, 3)
germlinenocancers$ProtDescription <- substring(germlinenocancers$ProtDescription, 3)

#bear in mind, everyone in this dataset is already p53-mutated! no p53-unaltered patients or samples here, so all calculations are WITHIN an already mutated population! cannot do prevalence across altered and unaltered cancers then, at best can do mutation type proportion within mutated samples; so cannot do exact comparisons with cbioportal calculations

#create dataframes of unique (MUTATED) patient IDs of LFS patients with and without cancer and combine into 1
#germlinejustcancersuniquepatientIDs <- distinct(germlinejustcancers, Individual_code) 
#germlinenocancersuniquepatientIDs <- distinct(germlinenocancers, Individual_code)
#germlinealluniquepatientIDs <- rbind(germlinejustcancersuniquepatientIDs, germlinenocancersuniquepatientIDs)

#sort short topology by unique tumor ID and descending counts (calculating number of MUTATED samples)
germlinecancersunique <- germlinejustcancers %>%
  group_by(Short_topo) %>%
  summarise(count = n_distinct(Tumor_ID))
germlinecancersunique <- germlinecancersunique[order(-germlinecancersunique$count),] %>% remove_rownames()

#keep only short topologies with at least 10 unique (MUTATED) tumor IDs 
germlinecancersover10 <- germlinecancersunique %>% 
  subset(count >= 10) 
germlinecancersover10 <- germlinecancersover10[order(-germlinecancersover10$count),] %>% remove_rownames()

#add total row at bottom of short topologies with at least 10 unique (MUTATED) tumor IDs
germlinecancersover10withtotal <- germlinecancersover10 %>% adorn_totals() %>% remove_rownames()

#graph 10+ (MUTATED) IDs topologies as plot and save
germlinecancersover10grid <- tableGrob(germlinecancersover10withtotal, rows = NULL, theme = ttheme_minimal(base_size = 10, core = list(padding=unit(c(2, 2), "mm"))))
germlinecancersgrid <- grid.arrange(germlinecancersover10grid)
#save_plot("FILEPATH.svg", germlinecancersgrid)
write_csv(germlinecancersover10withtotal, file = "C:/Users/nwali/Downloads/germlinecancersover10.csv")

#make subset of germline data which only have topologies with at least 10 tumor IDs
germlinejustcancersover10 <- germlinejustcancers %>% 
  subset(germlinejustcancers$Short_topo %in% germlinecancersover10$Short_topo) %>% 
  remove_rownames()
#check that cancers in new subset of germline data are only those with at least 10 tumor IDs
all(germlinejustcancersover10$Short_topo %in% germlinecancersover10$Short_topo) #should be TRUE
#identify which of 10+ samples cancers don't have p53 muts and will not be included in pivot tables (should be 0 because by definition everything assessed is already mutated!)
germlinecancersover10$Short_topo[!germlinecancersover10$Short_topo %in% germlinejustcancersover10$Short_topo]

#make pivot tables based on mutation type, sorted by descending counts for rows and columns; need 2 pivot tables because one is for patients with cancers and one is for non-cancer patients because they still have mutant p53
muttypepivotgermcancer <- PivotTable$new()
muttypepivotgermcancer$addData(germlinejustcancersover10)
muttypepivotgermcancer$addRowDataGroups("Short_topo")
muttypepivotgermcancer$addColumnDataGroups("Effect")
muttypepivotgermcancer$defineCalculation(calculationName = "Count of mutations", summariseExpression = "n()")
muttypepivotgermcancer$sortColumnDataGroups(levelNumber = 1, orderBy = "calculation", sortOrder = "desc")
muttypepivotgermcancer$sortRowDataGroups(levelNumber = 1, orderBy = "calculation", sortOrder = "desc")
muttypepivotgermcancer$renderPivot()

muttypepivotgermnocancer <- PivotTable$new()
muttypepivotgermnocancer$addData(germlinenocancers)
muttypepivotgermnocancer$addRowDataGroups("Short_topo")
muttypepivotgermnocancer$addColumnDataGroups("Effect")
muttypepivotgermnocancer$defineCalculation(calculationName = "Count of mutations", summariseExpression = "n()")
muttypepivotgermnocancer$sortColumnDataGroups(levelNumber = 1, orderBy = "calculation", sortOrder = "desc")
muttypepivotgermnocancer$sortRowDataGroups(levelNumber = 1, orderBy = "calculation", sortOrder = "desc")
muttypepivotgermnocancer$renderPivot()

#make pivot tables as dataframes for further calculations
pivotdfgermcancer <- muttypepivotgermcancer$asDataFrame(rowGroupsAsColumns = TRUE) %>% remove_rownames()
pivotdfgermnocancer <- muttypepivotgermnocancer$asDataFrame(rowGroupsAsColumns = TRUE) %>% remove_rownames()

#make 100% stacked barchart based on topologies with at least 10 muts
#remove bottom total row again
pivotdfgermcancernobottomtotal10muts <- pivotdfgermcancer %>% subset(Short_topo != "Total") %>% remove_rownames()
pivotdfgermnocancernobottomtotal10muts <- pivotdfgermnocancer %>% subset(Short_topo != "Total") %>% remove_rownames()

#subset dataframes to cancers with at least 10 muts (should not cause any changes as all are over 10 muts)
pivotdfgermcancernobottomtotal10muts <- pivotdfgermcancernobottomtotal10muts %>% 
  subset(Total >= 10) %>% remove_rownames()
pivotdfgermnocancernobottomtotal10muts <- pivotdfgermnocancernobottomtotal10muts %>% 
  subset(Total >= 10) %>% remove_rownames()

#divide each column by total # muts
pivotdfgermcancernobottomtotal10muts <- pivotdfgermcancernobottomtotal10muts %>%
  mutate_at(vars(missense:Total), .funs = ~./Total) %>% remove_rownames()
pivotdfgermnocancernobottomtotal10muts <- pivotdfgermnocancernobottomtotal10muts %>%
  mutate_at(vars(missense:Total), .funs = ~./Total) %>% remove_rownames()

#remove total column so it's not included in barchart
pivotdfgermcancernobottomtotal10muts <- pivotdfgermcancernobottomtotal10muts %>% 
  subset(select = -c(Total)) %>% remove_rownames()
pivotdfgermnocancernobottomtotal10muts <- pivotdfgermnocancernobottomtotal10muts %>% 
  subset(select = -c(Total)) %>% remove_rownames()

#melt dataframes and rename columns to use for ggplot2 100% stacked barcharts
germcancer_data_barchart <- reshape2::melt(pivotdfgermcancernobottomtotal10muts) %>% dplyr::rename(Mutation = variable, Cancer = Short_topo, Proportion = value)
head(germcancer_data_barchart)

germnocancer_data_barchart <- reshape2::melt(pivotdfgermnocancernobottomtotal10muts) %>% dplyr::rename(Mutation = variable, Cancer = Short_topo, Proportion = value)
head(germnocancer_data_barchart)

#say "No cancer at follow-up" for all no cancer rows which are otherwise empty
germnocancer_data_barchart$Cancer[germnocancer_data_barchart$Cancer == ''] <- 'No cancer at follow-up'
head(germnocancer_data_barchart)

#combine cancer and no cancer data into 1 for stacked barchart, can rbind because same columns
combinedgerm_data_barchart <- rbind(germcancer_data_barchart, germnocancer_data_barchart) %>% remove_rownames()

#sort cancers by lowest to highest missense mutation proportion for stacked barcharts
germlinebarchartorder = combinedgerm_data_barchart[combinedgerm_data_barchart$Mutation == 'missense',]
combinedgerm_data_barchart$Cancer = factor(combinedgerm_data_barchart$Cancer, levels = germlinebarchartorder$Cancer[order(germlinebarchartorder$Proportion)])

#missensemutordergermcancerbarchart = germcancer_data_barchart[germcancer_data_barchart$Mutation == 'missense',]
#germcancer_data_barchart$Cancer = factor(germcancer_data_barchart$Cancer, levels = missensemutordergermcancerbarchart$Cancer[order(missensemutordergermcancerbarchart$Proportion)])

#missensemutordergermnocancerbarchart = germnocancer_data_barchart[germnocancer_data_barchart$Mutation == 'missense',]
#germnocancer_data_barchart$Cancer = factor(germnocancer_data_barchart$Cancer, levels = missensemutordergermnocancerbarchart$Cancer[order(missensemutordergermnocancerbarchart$Proportion)])

#make 100% stacked barcharts with ggplot2 with reversed missense order i.e. highest to lowest
ggplotstackedgermlinebarchart <- ggplot(combinedgerm_data_barchart, aes(x = Cancer, y = Proportion, fill = Mutation)) + geom_bar(position = "fill", stat = "identity") + scale_y_continuous(labels = scales::percent_format(), expand = c(0,0), breaks = scales::pretty_breaks(n = 6)) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + scale_x_discrete(limits = rev(levels(combinedgerm_data_barchart$Cancer)), expand = c(0,0)) + ggtitle("# muts/# total muts") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text = element_text(color = "black"), axis.title = element_text(color = "black"))
ggplotstackedgermlinebarchart
save_plot(file = "C:/Users/nwali/Downloads/germlinebarchart.svg", ggplotstackedgermlinebarchart, base_width = 6, base_height = 4.5)

# if needed, can also make stacked barchart horizontal
origdb_ggplotstackedbarchart_horiz <- ggplot(origdbpivottype_melted, 
                                             aes(x = Cancer, 
                                                 y = Proportion, 
                                                 fill = Mutation)) + 
  geom_bar(position = position_fill(reverse = TRUE), 
           stat = "identity") + 
  ggtitle("# muts/# total muts") + 
  scale_y_continuous(labels = scales::percent_format(), 
                     expand = c(0,0), 
                     breaks = scales::pretty_breaks(n = 6), 
                     position = "right") +
  scale_x_discrete(expand = c(0,0)) + 
  theme_classic() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5), 
        axis.text = element_text(color = "black",
                                 size = 12), 
        axis.title = element_text(color = "black",
                                  face = "bold",
                                  size = 14),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  color = "black",
                                  size = 16),
        legend.text = element_text(color = "black",
                                   size = 12),
        legend.title = element_text(color = "black",
                                    face = "bold",
                                    size = 14),
        axis.text.y = element_text(vjust = 0.5))
origdb_ggplotstackedbarchart_horiz
# save_plot(file = "C:/Users/nwali/Downloads/origdbbarchart_horiz_germline.svg", 
#           origdb_ggplotstackedbarchart_horiz, 
#           base_width = 5,
#           base_height = 5)

#ggplotstackedgermnocancerbarchart <- ggplot(germnocancer_data_barchart, aes(x = Cancer, y = Proportion, fill = Mutation)) + geom_bar(position = "fill", stat = "identity") + scale_y_continuous(labels = scales::percent_format(), expand = c(0,0), breaks = scales::pretty_breaks(n = 6)) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + scale_x_discrete(limits = rev(levels(germnocancer_data_barchart$Cancer)), expand = c(0,0)) + ggtitle("# muts/# total muts") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text = element_text(color = "black"), axis.title = element_text(color = "black"), axis.title.x = element_blank())
#ggplotstackedgermnocancerbarchart
#save_plot(file = "C:/Users/nwali/Downloads/germnocancerbarchart.svg", ggplotstackedgermnocancerbarchart, base_width = 3, base_height = 4.5)

#convert proportions to percentages rounded to 1 decimal point in new column for no cancer data so that can show percentages on upcoming pie chart
#germnocancer_data_barchart <- germnocancer_data_barchart %>% mutate(Labels = scales::percent(Proportion, accuracy = 0.1))

#can make pie chart based on no cancer data since only 1 column in stacked barchart
#ggplotpiechartgermnocancer <- ggplot(germnocancer_data_barchart, aes(x = Cancer, y = Proportion, fill = Mutation)) + geom_bar(position = "fill", stat = "identity") + ggtitle("# muts/# total muts, No cancer at follow up") + geom_label(aes(label = Labels, x = 1.32), position = position_stack(vjust = 0.5), show.legend = FALSE) + coord_polar(theta = "y", start = 0, direction = -1) + theme_void() + theme(plot.title = element_text(hjust = 0.5))
#ggplotpiechartgermnocancer
#save_plot(file = "C:/Users/nwali/Downloads/germnocancerpiechart.svg", ggplotpiechartgermnocancer, base_width = 12, base_height = 8)

#make pivot tables based on mutation position, sorted by descending counts for rows and columns; need 2 pivot tables because one is for patients with cancers and one is for non-cancer patients because they still have mutant p53
mutpospivotgermcancer <- PivotTable$new()
mutpospivotgermcancer$addData(germlinejustcancersover10)
mutpospivotgermcancer$addRowDataGroups("Short_topo")
#mutpospivotgermcancer$addColumnDataGroups("Codon_number")
mutpospivotgermcancer$addColumnDataGroups("ProtDescription") #, addTotal = FALSE)
mutpospivotgermcancer$defineCalculation(calculationName = "Count of mutations", summariseExpression = "n()")
#mutpospivotgermcancer$sortColumnDataGroups(levelNumber = 1, orderBy = "value", sortOrder = "asc")
mutpospivotgermcancer$sortColumnDataGroups(levelNumber = 1, orderBy = "calculation", sortOrder = "desc")
mutpospivotgermcancer$sortRowDataGroups(levelNumber = 1, orderBy = "calculation", sortOrder = "desc")
mutpospivotgermcancer$renderPivot()

mutpospivotgermnocancer <- PivotTable$new()
mutpospivotgermnocancer$addData(germlinenocancers)
mutpospivotgermnocancer$addRowDataGroups("Short_topo")
#mutpospivotgermnocancer$addColumnDataGroups("Codon_number")
mutpospivotgermnocancer$addColumnDataGroups("ProtDescription") #, addTotal = FALSE)
mutpospivotgermnocancer$defineCalculation(calculationName = "Count of mutations", summariseExpression = "n()")
#mutpospivotgermnocancer$sortColumnDataGroups(levelNumber = 1, orderBy = "value", sortOrder = "asc")
mutpospivotgermnocancer$sortColumnDataGroups(levelNumber = 1, orderBy = "calculation", sortOrder = "desc")
mutpospivotgermnocancer$sortRowDataGroups(levelNumber = 1, orderBy = "calculation", sortOrder = "desc")
mutpospivotgermnocancer$renderPivot()

#make pivot tables as dataframes for further calculations
pivotdfposgermcancer <- mutpospivotgermcancer$asDataFrame(rowGroupsAsColumns = TRUE) %>% remove_rownames()
pivotdfposgermnocancer <- mutpospivotgermnocancer$asDataFrame(rowGroupsAsColumns = TRUE) %>% remove_rownames()

#remove initial codon numbers from dataframe columns to clean up visually
#names(pivotdfposgermcancer) <- trimws(colnames(pivotdfposgermcancer), whitespace = "\\d+\\s") %>% str_trim("right") %>% remove_rownames()
#names(pivotdfposgermnocancer) <- trimws(colnames(pivotdfposgermnocancer), whitespace = "\\d+\\s") %>% str_trim("right") %>% remove_rownames()

#make joint heatmap based on mutation positions with at least 10 positions
#remove bottom total row again
pivotdfposgermcancernobottomtotal <- pivotdfposgermcancer %>% subset(Short_topo != "Total") %>% remove_rownames()
pivotdfposgermnocancernobottomtotal <- pivotdfposgermnocancer %>% subset(Short_topo != "Total") %>% remove_rownames()

#subset dataframes to cancers with at least 10 muts (should not cause any changes as all are over 10 muts)
pivotdfposgermcancernobottomtotal <- pivotdfposgermcancernobottomtotal %>% 
  subset(Total >= 10) %>% remove_rownames()
pivotdfposgermnocancernobottomtotal <- pivotdfposgermnocancernobottomtotal %>% 
  subset(Total >= 10) %>% remove_rownames()

#divide each column by total # muts
pivotdfposgermcancernobottomtotal <- pivotdfposgermcancernobottomtotal %>%
  mutate_at(vars(2:Total), .funs = ~./Total) %>% remove_rownames()
pivotdfposgermnocancernobottomtotal <- pivotdfposgermnocancernobottomtotal %>%
  mutate_at(vars(2:Total), .funs = ~./Total) %>% remove_rownames()

#remove total column so it's not included in heatmap
pivotdfposgermcancernobottomtotal <- pivotdfposgermcancernobottomtotal %>% 
  subset(select = -c(Total)) %>% remove_rownames()
pivotdfposgermnocancernobottomtotal <- pivotdfposgermnocancernobottomtotal %>% 
  subset(select = -c(Total)) %>% remove_rownames()

#melt dataframes and rename columns to use for ggplot2 heatmap
germcancer_data_heatmap <- reshape2::melt(pivotdfposgermcancernobottomtotal) %>% dplyr::rename(Mutation = variable, Cancer = Short_topo, Proportion = value)
head(germcancer_data_heatmap)

germnocancer_data_heatmap <- reshape2::melt(pivotdfposgermnocancernobottomtotal) %>% dplyr::rename(Mutation = variable, Cancer = Short_topo, Proportion = value)
head(germnocancer_data_heatmap)

#say "No cancer at follow-up" for all no cancer rows which are otherwise empty
germnocancer_data_heatmap$Cancer[germnocancer_data_heatmap$Cancer == ''] <- 'No cancer at follow-up'
head(germnocancer_data_heatmap)

#combine melted data from cancer and no cancer patients to make 1 heatmap, can rbind because same column names
combinedgerm_data_heatmap <- rbind(germcancer_data_heatmap, germnocancer_data_heatmap) %>% remove_rownames()

#remove mutations with less than 1% proportion in total muts
combinedgerm_data_heatmap_over0.01 <- combinedgerm_data_heatmap %>% subset(Proportion >= 0.01)

#sort data by lowest to highest 0-splice mutation proportion for heatmap (all cancers must have this mutation at 1% prevalence otherwise it will make cancers NA)
# germlineheatmaporder = combinedgerm_data_heatmap_over0.01[combinedgerm_data_heatmap_over0.01$Mutation == '?-0-splice',]
# combinedgerm_data_heatmap_over0.01$Cancer = factor(combinedgerm_data_heatmap_over0.01$Cancer, levels = germlineheatmaporder$Cancer[order(germlineheatmaporder$Proportion, na.last = FALSE)])

#make ggplot heatmap
ggplotheatmapgermlinepos <- ggplot(combinedgerm_data_heatmap_over0.01, aes(x = Mutation, y = Cancer)) + geom_tile(aes(fill = Proportion)) + scale_fill_gradient(low = "white", high = "red", na.value = "white") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 0), axis.text = element_text(color = "black"), axis.title = element_text(color = "black")) + scale_x_discrete(position = "top", expand = c(0,0)) + scale_fill_continuous(breaks = seq(0, 0.3, by = 0.05), low = "white", high = "red", na.value = "white") + ggtitle("# muts/# total muts, ≥ 1% proportion") + theme(plot.title = element_text(hjust = 0.5)) + scale_y_discrete(expand = c(0,0))
ggplotheatmapgermlinepos
save_plot(file = "C:/Users/nwali/Downloads/germlineposhm.svg", ggplotheatmapgermlinepos, base_width = 20, base_height = 3.5)

#try making better formatted heatmap with complexheatmap
#create new dataframe of positions with 1% proportion and keep distinct positions
posgermline_totalmuts <- data.frame(combinedgerm_data_heatmap_over0.01$Mutation, stringsAsFactors = FALSE) %>% dplyr::rename(Mutation = combinedgerm_data_heatmap_over0.01.Mutation)
head(posgermline_totalmuts)
posgermline_totalmuts <- distinct(posgermline_totalmuts, Mutation) %>% remove_rownames()
head(posgermline_totalmuts)

#extract columns of position, codon, domain, and buried or not in protein from original germline data and only keep completely distinct rows of information
germlinemutinfo <- germlinedataorig %>% subset(select = c("Codon_number", "ProtDescription", "Domain_function", "Residue_function", "Effect"))
head(germlinemutinfo)
germlinemutinfo <- germlinemutinfo %>% distinct()
head(germlinemutinfo)

#concatenate codon number and mut type to p.? mutations which have unspecified residues
germlinemutinfo <- germlinemutinfo %>% mutate(ProtDescription = ifelse(ProtDescription %in% c("p.?"), paste(ProtDescription, Codon_number, Effect, sep = "-"), no = ProtDescription))
head(germlinemutinfo)

#remove p. prefix from all mutations to clean up visually
germlinemutinfo$ProtDescription <- substring(germlinemutinfo$ProtDescription, 3)
head(germlinemutinfo)

#only keep mutations that are in the 1% proportion list
germlinemutinfo <- germlinemutinfo %>% 
  subset(germlinemutinfo$ProtDescription %in% posgermline_totalmuts$Mutation)
#check that all mutations in new subset of germline data are only those with at least 1% proportion
all(germlinemutinfo$ProtDescription %in% posgermline_totalmuts$Mutation) #should be TRUE

#domains and codons already assigned for the most part, so just sort position dataframe by ascending codon
germlinemutinfo <- germlinemutinfo[order(germlinemutinfo$Codon_number),] %>% remove_rownames()
head(germlinemutinfo)

#need rownames to be things to match for annotation in complexheatmap, so add positions to rownames
rownames(germlinemutinfo) <- germlinemutinfo$ProtDescription
head(germlinemutinfo)

#as some annotations are NA and will not show up in heatmap annotations, need to actually write cells as NA text
germlinemutinfo <- germlinemutinfo %>% replace_na(list(Domain_function = 'NA', Residue_function = 'NA'))

#don't need to make annotations dataframe for cancers because already by tissue essentially, so will just annotate positions by p53 domain and if buried/exposed, etc.

#convert melted dataframe to matrix for complexheatmap, filling anything below 0.01 as 0 although not 0 in reality, just placeholder to not have errors in making heatmap later on
posgermline_matrix_totalmuts <- combinedgerm_data_heatmap_over0.01 %>% reshape2::acast(Cancer ~ Mutation, value.var = 'Proportion', fill = '0')
class(posgermline_matrix_totalmuts) <- "numeric"

#set up annotation dataframes to just have columns used for annotation, otherwise all columns used for annotation in heatmap
germlinemutinfo_annotate <- germlinemutinfo %>% select(-c(Codon_number, ProtDescription, Effect))

#order domain by p53 domain as a factor and buried/exposed as factor so that they are maintained in legend too
unique(germlinemutinfo_annotate$Domain_function) #to find list of domains in order to use as factors
germlinemutinfo_annotate$Domain_function <- factor(germlinemutinfo_annotate$Domain_function, levels = c("Transactivation TAD1", "Transactivation TAD2", "SH3-like/Pro-rich", "DNA binding", "Tetramerisation", "Tetramerisation/NES", "NA"))
sort(as.vector(unique(germlinemutinfo_annotate$Residue_function))) #to find list of domains in order to use as factors
germlinemutinfo_annotate$Residue_function <- factor(germlinemutinfo_annotate$Residue_function, levels = c("ADP-ribosylation site", "Buried", "DNA binding", "Exposed", "Partially exposed", "Phosphorylation site", "S-glutathionylation site", "Tetramerisation", "Tetramerisation/Methylation site", "Zn binding", "NA"))

#show all color-blind friendly palettes in RColorBrewer
display.brewer.all(colorblindFriendly = TRUE)

#create color range for proportions on complexheatmap
color_germline <- colorRampPalette((c("white", "red")))(50)

#assign color hex codes to annotations
annotationcolorgermline_totalmuts <- list(Domain_function = RColorBrewer::brewer.pal(n = 7, name = "Set2"), Residue_function = RColorBrewer::brewer.pal(n = 11, name = "Paired"))
annotationcolornamedgermline_totalmuts <- list(
  Domain_function = c(`Transactivation TAD1` = "#E78AC3", 
                      `Transactivation TAD2` = "#FFD92F",
                      `SH3-like/Pro-rich` = "#FC8D62", 
                      `DNA binding` = "#8DA0CB", 
                      `Tetramerisation` = "brown",
                      `Tetramerisation/NES` = "#A6D854", 
                      `NA` = "black"), 
  Residue_function = c(`ADP-ribosylation site` = "#CAB2D6",
                        Buried = "#FF10F0", 
                       `DNA binding` = "#6A3D9A", 
                       Exposed = "#1F78B4", 
                       `Partially exposed` = "#A6CEE3", 
                       `Phosphorylation site` = "#33A02C", 
                       `S-glutathionylation site` = "#B2DF8A", 
                       Tetramerisation = "#8B8000", 
                       `Tetramerisation/Methylation site` = "#FF7F00", 
                       `Zn binding` = "#E30B5C", 
                        `NA` = "#5A5A5A"))

#order matrix based on codon order and alphabetical cancer so that heatmap organized by ascending codon and alphabetical cancer
posgermline_matrix_ordered_totalmuts <- posgermline_matrix_totalmuts[order(row.names(posgermline_matrix_totalmuts)),rownames(germlinemutinfo_annotate)]

#need to recreate annotations as annotation class for complexheatmap
column_ha_germline <- ComplexHeatmap::columnAnnotation(df = germlinemutinfo_annotate, col = annotationcolornamedgermline_totalmuts, annotation_label = c("Domain", "Feature"), show_annotation_name = FALSE)

#draw complexheatmaps, a pair without and with cancer clustering for total muts
svg(file = "C:/Users/nwali/Downloads/complexheatmap_germline.svg", width = 18, height = 4.5)
complexhm_germline <- ComplexHeatmap::Heatmap(posgermline_matrix_ordered_totalmuts,
                                     col = color_germline,
                                     name = "Proportion",
                                     na_col = "white",
                                     row_names_side = "left",
                                     column_names_side = "top",
                                     cluster_columns = FALSE,
                                     row_dend_side = "left",
                                     column_dend_side = "top",
                                     column_title = "# muts/# total muts, ≥ 1% proportion",
                                     column_title_side = "top",
                                     column_title_gp = gpar(fontface = "bold"),
                                     row_names_gp = gpar(fontsize = 10),
                                     row_names_max_width = max_text_width(
                                       rownames(posgermline_matrix_ordered_totalmuts), 
                                       gp = gpar(fontsize = 10)),
                                     column_names_gp = gpar(fontsize = 10),
                                     top_annotation = column_ha_germline,
                                     heatmap_legend_param = list(at = seq(from = 0, to = 0.3, by = 0.05))
)
draw(complexhm_germline)
dev.off() #run multiple times until following error shows up: Error in dev.off() : cannot shut down device 1 (the null device)

svg(file = "C:/Users/nwali/Downloads/complexheatmap_germline_ordered.svg", width = 18, height = 4.5)
complexhm_ordered_germline <- ComplexHeatmap::Heatmap(posgermline_matrix_ordered_totalmuts,
                                             col = color_germline,
                                             name = "Proportion",
                                             na_col = "white",
                                             row_names_side = "left",
                                             column_names_side = "top",
                                             cluster_rows = FALSE,
                                             cluster_columns = FALSE,
                                             row_dend_side = "left",
                                             column_dend_side = "top",
                                             column_title = "# muts/# total muts, ≥ 1% proportion",
                                             column_title_side = "top",
                                             column_title_gp = gpar(fontface = "bold"),
                                             row_names_gp = gpar(fontsize = 10),
                                             row_names_max_width = max_text_width(
                                               rownames(posgermline_matrix_ordered_totalmuts), 
                                               gp = gpar(fontsize = 10)),
                                             column_names_gp = gpar(fontsize = 10),
                                             top_annotation = column_ha_germline,
                                             heatmap_legend_param = list(at = seq(from = 0, to = 0.3, by = 0.05))
)
draw(complexhm_ordered_germline)
dev.off() #run multiple times until following error shows up: Error in dev.off() : cannot shut down device 1 (the null device)

# also can make alluvial plots to show flow of proportions from tissues
# need to combine all strata labels of alluvial plot (i.e. heatmap annotation labels) into melted dfs so that everything assigned to a line, convert things you want ordered into factors so that they have defined order in alluvial plot, and rename columns to better names so that if in legend for alluvial plot automatically correct name
alluvial_type <- merge(origdbpivottype_melted,
                       origdbcancersinfo,
                       by = "Cancer")
alluvial_type$Germ_Layer <- factor(alluvial_type$Germ_Layer, 
                                   levels = unique(origdbcancersinfo$Germ_Layer))
alluvial_type$Cancer <- factor(alluvial_type$Cancer, 
                               levels = unique(origdbcancersinfo$Cancer))
alluvial_type <- alluvial_type %>% 
  rename(Tissue = Cancer,
         `Germ Layer` = Germ_Layer,
         `Mutation Type` = Mutation)

alluvial_pos <- merge(origdbpivotpos_melted_over0.01,
                      origdbcancersinfo,
                      by = "Cancer")
alluvial_pos <- merge(alluvial_pos,
                      origdbmutinfo,
                      by.x = "Mutation",
                      by.y = "ProtDescription")
alluvial_pos <- alluvial_pos[order(alluvial_pos$Codon_number),]
alluvial_pos$Germ_Layer <- factor(alluvial_pos$Germ_Layer, 
                                  levels = unique(origdbcancersinfo$Germ_Layer))
alluvial_pos$Cancer <- factor(alluvial_pos$Cancer, 
                              levels = unique(origdbcancersinfo$Cancer))
alluvial_pos$Domain <- factor(alluvial_pos$Domain, 
                              levels = levels(origdbmutinfo_annotate$Domain))
alluvial_pos$Mutation <- factor(alluvial_pos$Mutation, 
                                levels = unique(alluvial_pos$Mutation))
alluvial_pos <- alluvial_pos %>% 
  rename(Tissue = Cancer,
         `Germ Layer` = Germ_Layer,
         `p53 Domain` = Domain,
         Function = Residue_function)

# draw alluvial plots with tissue as coloring
alluvialplot_type <- ggplot(data = alluvial_type,
                            aes(axis1 = `Germ Layer`,   # First variable on the X-axis
                                axis2 = Tissue, # Second variable on the X-axis
                                axis3 = `Mutation Type`,   # Third variable on the X-axis
                                y = Proportion)) +
  geom_alluvium(aes(fill = Tissue),
                aes.bind = "alluvia") +
  geom_stratum(alpha = 0.2) +     # makes strata somewhat transparent so can see colors feeding into them
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            min.y = 0.05) + # so that very small categories are not labeled which are getting too tight to see
  ggtitle("# muts/# total muts") +
  scale_x_continuous(breaks = 1:3, 
                     labels = c("Germ Layer", "Tissue", "Mutation Type"),
                     expand = c(0,0)) +
  theme_classic() + 
  scale_y_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        axis.text = element_text(color = "black",
                                 size = 12,
                                 face = "bold"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(face = "bold",
                                    size = 14,
                                    color = "black"),
        legend.text = element_text(size = 12,
                                   color = "black"))
alluvialplot_type

alluvialplot_pos <- ggplot(data = alluvial_pos,
                           aes(axis1 = `Germ Layer`,   # First variable on the X-axis
                               axis2 = Tissue, # Second variable on the X-axis
                               axis3 = `p53 Domain`, # Third variable on the X-axis
                               axis4 = Mutation,   # Fourth variable on the X-axis
                               y = Proportion)) +
  geom_alluvium(aes(fill = Tissue),
                aes.bind = "alluvia") +
  geom_stratum(alpha = 0.2) +     # makes strata somewhat transparent so can see colors feeding into them
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            min.y = 0.05) + # so that very small categories are not labeled which are getting too tight to see
  ggtitle("# muts/# total muts, ≥ 1% proportion") +
  scale_x_continuous(breaks = 1:4, 
                     labels = c("Germ Layer", "Tissue", "p53 Domain", "Mutation"),
                     expand = c(0,0)) +
  theme_classic() + 
  scale_y_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        axis.text = element_text(color = "black",
                                 size = 12,
                                 face = "bold"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(face = "bold",
                                    size = 14,
                                    color = "black"),
        legend.text = element_text(size = 12,
                                   color = "black"))
alluvialplot_pos

# draw alluvial plots with germ layer as coloring
alluvialplot_type_side <- ggplot(data = alluvial_type,
                                 aes(axis1 = `Germ Layer`,   # First variable on the X-axis
                                     axis2 = Tissue, # Second variable on the X-axis
                                     axis3 = `Mutation Type`,   # Third variable on the X-axis
                                     y = Proportion)) +
  geom_alluvium(aes(fill = `Germ Layer`),
                aes.bind = "alluvia") +
  geom_stratum(alpha = 0.2) +     # makes strata somewhat transparent so can see colors feeding into them
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            min.y = 0.05) + # so that very small categories are not labeled which are getting too tight to see
  ggtitle("# muts/# total muts") +
  scale_x_continuous(breaks = 1:3, 
                     labels = c("Germ Layer", "Tissue", "Mutation Type"),
                     expand = c(0,0)) +
  scale_fill_manual(values = annotationcolornamedorigdb[["Germ_Layer"]]) +
  theme_classic() + 
  scale_y_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        axis.text = element_text(color = "black",
                                 size = 12,
                                 face = "bold"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(face = "bold",
                                    size = 14,
                                    color = "black"),
        legend.text = element_text(size = 12,
                                   color = "black"))
alluvialplot_type_side

alluvialplot_pos_side <- ggplot(data = alluvial_pos,
                                aes(axis1 = `Germ Layer`,   # First variable on the X-axis
                                    axis2 = Tissue, # Second variable on the X-axis
                                    axis3 = `p53 Domain`, # Third variable on the X-axis
                                    axis4 = Mutation,   # Fourth variable on the X-axis
                                    y = Proportion)) +
  geom_alluvium(aes(fill = `Germ Layer`),
                aes.bind = "alluvia") +
  geom_stratum(alpha = 0.2) +     # makes strata somewhat transparent so can see colors feeding into them
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            min.y = 0.05) + # so that very small categories are not labeled which are getting too tight to see
  ggtitle("# muts/# total muts, ≥ 1% proportion") +
  scale_x_continuous(breaks = 1:4, 
                     labels = c("Germ Layer", "Tissue", "p53 Domain", "Mutation"),
                     expand = c(0,0)) +
  scale_fill_manual(values = annotationcolornamedorigdb[["Germ_Layer"]]) +
  theme_classic() + 
  scale_y_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        axis.text = element_text(color = "black",
                                 size = 12,
                                 face = "bold"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(face = "bold",
                                    size = 14,
                                    color = "black"),
        legend.text = element_text(size = 12,
                                   color = "black"))
alluvialplot_pos_side

# draw flipped alluvial plots
alluvialplot_type_horiz <- alluvialplot_type +
  coord_flip() + 
  scale_y_reverse() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
alluvialplot_type_horiz

alluvialplot_pos_horiz <- alluvialplot_pos +
  coord_flip() + 
  scale_y_reverse() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
alluvialplot_pos_horiz

alluvialplot_type_horiz_side <- alluvialplot_type_side +
  coord_flip() + 
  scale_y_reverse() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
alluvialplot_type_horiz_side

alluvialplot_pos_horiz_side <- alluvialplot_pos_side +
  coord_flip() + 
  scale_y_reverse() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
alluvialplot_pos_horiz_side









