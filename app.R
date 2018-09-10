##all
library(shiny)
library(kBET)
library(ggplot2)
library(reshape2)
library(matrixStats)
library(markdown)
#source('BatchMixingEntropy.R')
library(plotly)
library(kBET)
library(akima)
library(reshape2)
library(crosstalk)
library(tidyr)
library(dplyr)
library (randomForest)
library(biomaRt)
library(org.Mm.eg.db)
library(topGO)
library (tree)
library(gridExtra)
library(htmlwidgets)
library(htmltools)
library(jsonlite)
#library(kableExtra)
library(grDevices)
library(wgeesel)
library(MASS)
library(MuMIn) 
library(dplyr)
library(rmarkdown)
library(data.table)
library(ggplot2)
library(grid)
library(gee)
library(gridExtra)
library(lme4)
library(lmtest)
library(DT)
library(Rtsne)
library(plotly)
library(htmlwidgets)
require(TSCAN)
library(shinycssloaders)

load("rnaseq.rda")
load("raw.rda")
load("info.rda")
load("mnn.rda")
load("comb.rda")
load("iter.rda")
#load('ave_dist.rda')

load('true_promoter.rda')
load('pred_promoter.rda')

load('true_exon.rda')
load('pred_exon.rda')

load('true_intron.rda')
load('pred_intron.rda')

load('true_enhancer.rda')
load('pred_enhancer.rda')

load('true_tf.rda')
load('pred_tf.rda')



rna <- readRDS("RNA_data_norm_hg19_all_name.rds")
DHScluster <- as.numeric(readLines("DH_cluster_2000_hg19_all.txt"))

ui <-  shinyUI(navbarPage("APP",
                          tabPanel("APP", fluidPage(
                            #headerPanel("GrowthAnalyst"),
                            sidebarLayout(
                              sidebarPanel(
                                h3('Data processing'),
                                # checkboxInput('upload', "Upload dataset"),
                                # conditionalPanel(
                                #   condition = "input.cate== true",
                                #   fileInput("file1", "Upload DNase CSV File")
                                # ),
                                # checkboxInput("new", label = 'Show new data'),
                                # #selectInput('search','Search',multiple = T),
                                uiOutput('search'),
                                h3('Plotting option'),
                                checkboxInput("cate", label = 'Plot by category'),
                                conditionalPanel(
                                  condition = "input.cate== true",
                                  selectInput("cate_color", "Select category (color)", 
                                              list('Assay' = 'assay',"Biosample type"="bio", "Cell type"="cell","Germ layer"="layer",'Cancer' = 'cancer')
                                  ),
                                  uiOutput('cate_shape')
                                ),
                                h3('Method selection'),
                                checkboxInput("method", 'Correction methods'),
                                conditionalPanel(
                                  condition = "input.method== true",
                                  selectInput("corr", "Batch correction method", 
                                              list("uncorrected"="uncorrected", "combat"="comb","mnn"="mnn","iterated mnn"="iter")
                                  )),
                                h3("Filtering"),
                                checkboxInput("filter", 'Filter data'),
                                conditionalPanel(
                                  condition = "input.filter== true",
                                  checkboxInput('cancer','Cancer'),
                                  conditionalPanel(
                                    condition = "input.cancer == true",
                                    uiOutput('sel_cancer')
                                  ),
                                  checkboxInput('assay','Assay'),
                                  conditionalPanel(
                                    condition = "input.assay == true",
                                    uiOutput('sel_assay')
                                  ),
                                  checkboxInput('bio','Biosample type '),
                                  conditionalPanel(
                                    condition = "input.bio == true",
                                    uiOutput('sel_bio')
                                  ),
                                  checkboxInput('cell','Cell type'),
                                  conditionalPanel(
                                    condition = "input.cell == true",
                                    uiOutput('sel_cell')
                                  ),
                                  checkboxInput('layer','Germ layer'),
                                  conditionalPanel(
                                    condition = "input.layer == true",
                                    uiOutput('sel_layer')
                                  ),
                                  checkboxInput('var','View low variability cell type'),
                                  conditionalPanel(
                                    condition = "input.var == true",
                                    sliderInput("quantile",'View cell type with distance variance below quantile:' ,
                                                min = 0.1, max = 1, value = 0.5,step =0.05)
                                  )
                                ),
                                
                                
                                h3("Heatmap"),
                                checkboxInput("heat", 'Heatmap'),
                                conditionalPanel(
                                  condition = "input.heat == true",
                                  checkboxInput("geneexpre", 'Gene expression heatmap'),
                                    uiOutput('gene'),
                                    uiOutput('gene_pro'),
                                    uiOutput('gene_exon'),
                                    uiOutput('gene_intron'),
                                    uiOutput('gene_enhancer'),
                                    uiOutput('gene_tf')
                                ), 
                                
                                h3("Curves"),
                                checkboxInput("curve", 'Curves'),
                                conditionalPanel(
                                  condition = "input.curve == true",
                                  checkboxInput("genecurve", 'Gene expression curve'),
                                   conditionalPanel(
                                    condition = "input.genecurve == true",
                                    radioButtons(inputId="data_type_gene", "Data type:",
                                               choices = list("RNA-seq" = "gene_rna"),inline=TRUE),
                                  
                                  uiOutput('gene2'),
                                  uiOutput('celltype')
                                  ),
                                checkboxInput("procurve", 'Promoter curve'),
                                conditionalPanel(
                                  condition = "input.procurve == true",
                                  radioButtons(inputId="data_type_pro", "Data type:",
                                               choices = list("Predicted DNase" = "pro_pred","True DNase" = "pro_true",
                                                              "Predicted DNase + True DNase-seq" = "pro_both"),inline=TRUE),
                                  uiOutput('gene_pro2'),
                                  uiOutput('celltype_pro')
                                  ),
                                checkboxInput("exoncurve", 'Exon curve'),
                                conditionalPanel(
                                  condition = "input.exoncurve == true",
                                  radioButtons(inputId="data_type_exon", "Data type:",
                                               choices = list("Predicted DNase" = "exon_pred","True DNase" = "exon_true",
                                                              "Predicted DNase + True DNase-seq" = "exon_both"),inline=TRUE),
                                  uiOutput('gene_exon2'),
                                  uiOutput('celltype_exon')
                                  ),
                                  checkboxInput("introncurve", 'Intron curve'),
                                conditionalPanel(
                                  condition = "input.introncurve == true",
                                  radioButtons(inputId="data_type_intron", "Data type:",
                                               choices = list("Predicted DNase" = "intron_pred","True DNase" = "intron_true",
                                                              "Predicted DNase + True DNase-seq" = "intron_both"),inline=TRUE),
                                  uiOutput('gene_intron2'),
                                  uiOutput('celltype_intron')
                                  ),
                                checkboxInput("enhancercurve", 'Enhancer curve'),
                                conditionalPanel(
                                  condition = "input.enhancercurve == true",
                                  radioButtons(inputId="data_type_enhancer", "Data type:",
                                               choices = list("Predicted DNase" = "enhancer_pred","True DNase" = "enhancer_true",
                                                              "Predicted DNase + True DNase-seq" = "enhancer_both"),inline=TRUE),
                                  uiOutput('gene_enhancer2'),
                                  uiOutput('celltype_enhancer')
                                  ),
                                checkboxInput("tfcurve", 'Transcription factor curve'),
                                conditionalPanel(
                                  condition = "input.tfcurve == true",
                                  radioButtons(inputId="data_type_tf", "Data type:",
                                               choices = list("Predicted DNase" = "tf_pred", "True DNase" = "tf_true",
                                                              "Predicted DNase + True DNase-seq" = "tf_both"),inline=TRUE),
                                  uiOutput('gene_tf2'),
                                  uiOutput('celltype_tf')
                                  )
                                ),
                                h3("Enrichment"),
                                checkboxInput("enri", 'Enrichment'),
                                conditionalPanel(
                                  condition = "input.enri == true",
                                    uiOutput('celltype_enri'),
                                    selectizeInput("clus",
                                                 label = "Cluster of interest",
                                                 choices = c(1:20),
                                                 multiple = F),
                                  uiOutput('gene_clu')
                                )
                              ),
                              # Show a plot of the generated distribution
                              mainPanel(
                                tags$style(type="text/css",
                                           ".shiny-output-error { visibility: hidden; }",
                                           ".shiny-output-error:before { visibility: hidden; }"),
                                
                                tabsetPanel(
                                  tabPanel("PCA", 
                                           h4('PCA plot'),
                                           plotlyOutput('pca1')%>% withSpinner()
                                  ), 
                                  tabPanel("Tsne", 
                                           h4('Tsne plot'),
                                           plotlyOutput('plot1')%>% withSpinner()
                                  ), 
                                  
                                  tabPanel("Heatmaps", 
                                           h4('Gene expression heatmap'),
                                           plotlyOutput('heat1')%>% withSpinner(),
                                           plotlyOutput('heat_pca1')%>% withSpinner(),
                                           

                                           h4('Promoter heatmap'),
                                           plotlyOutput('heat_promoter')%>% withSpinner(),
                                           plotlyOutput('heat_promoter_pca')%>% withSpinner(),
                                           h4('Exon heatmap'),
                                           plotlyOutput('heat_exon')%>% withSpinner(),
                                           plotlyOutput('heat_exon_pca')%>% withSpinner(),
                                           h4('Intron heatmap'),
                                           plotlyOutput('heat_intron')%>% withSpinner(),
                                           plotlyOutput('heat_intron_pca')%>% withSpinner(),
                                           h4('Enhancer heatmap'),
                                           plotlyOutput('heat_enhancer')%>% withSpinner(),
                                           plotlyOutput('heat_enhancer_pca')%>% withSpinner(),
                                           h4('Transcription factor heatmap'),
                                           plotlyOutput('heat_tf')%>% withSpinner(),
                                           plotlyOutput('heat_tf_pca')%>% withSpinner()
                                  ),
                                  tabPanel("Gene expression curve", 
                                           h4('Gene expression curve'),
                                         
                                           plotlyOutput('curve1')%>% withSpinner(),
                                           
                                           plotOutput('pca_curve')%>% withSpinner(),
                                           verbatimTextOutput("impo"),
                                           verbatimTextOutput('curve_text'),
                                           verbatimTextOutput('exon_curve_text'),
                                           verbatimTextOutput('intron_curve_text'),
                                           h4('Promoter curve'),
                                           plotlyOutput('promoter_curve')%>% withSpinner(),
                                           plotOutput('promoter_pca_curve')%>% withSpinner(),

                                           h4('Exon curve'),
                                           plotlyOutput('exon_curve')%>% withSpinner(),
                                           plotOutput('exon_pca_curve')%>% withSpinner(),

                                           h4('Intron curve'),
                                           plotlyOutput('intron_curve')%>% withSpinner(),
                                           plotOutput('intron_pca_curve')%>% withSpinner(),

                                           h4('Enhancer curve'),
                                           plotlyOutput('enhancer_curve')%>% withSpinner(),
                                           plotOutput('enhancer_pca_curve')%>% withSpinner(),

                                           h4('Transcription factor curve'),
                                           plotlyOutput('tf_curve')%>% withSpinner(),
                                           plotOutput('tf_pca_curve')%>% withSpinner()
                                           

                                  ),
                                  tabPanel("Enrichment", 
                                           h4('Clustered RNA gene expression curve'),
                                           
                                           plotlyOutput('rna_curve')%>% withSpinner(),
                                           verbatimTextOutput('gene_cluster'),
                                           
                                           h4('GO analysis of chosen cluster'),
                                           
                                           verbatimTextOutput("go"),
                                           h4('Clustered promoter expression curve'),
                                           plotlyOutput('promoter_enri')%>% withSpinner(),
                                           h4('Clustered exon expression curve'),
                                           plotlyOutput('exon_enri')%>% withSpinner(),
                                           h4('Clustered intron expression curve'),
                                           plotlyOutput('intron_enri')%>% withSpinner()
                                  )
                                           
                                  # tabPanel("Correction methods evaluation",
                                  #          # 
                                  # 
                                  #         # plotOutput('dist'),
                                  #          plotOutput('silhouette'),
                                  #         # plotOutput('entropies'),
                                  #          plotOutput('k_unc'),
                                  #          plotOutput('k_comb'),
                                  #          plotOutput('k_mnn'),
                                  #          plotOutput('k_iter')
                                  )
                                  
                                )# end of tabset panel
                              )# end of main panel
                              
                            )
                          )
                          ))





server <- function(input, output) {
    cell=info[,1]
    cancer = info[,2]
    batch=info[,3]
    batch[info[,3]=='match']='Predicted DNase'
    batch[info[,3]=='predicted']='Predicted DNase'
    batch[info[,3]=='true']='True DNase' 
    bat = ifelse(batch=='True DNase',2,1)
    link= gsub("@.*","",info[,4])
    layer=info[,5]
    bio = info[,6]
    
    
    ###input 
    output$search <- renderUI({
      selectizeInput("search",
                     label = "Sample of Interest",
                     choices = colnames(raw),
                     multiple = T,
                     options = list(maxItems = nrow(raw), placeholder = 'Select a sample')
      )
    })
    
    
    ########### gene expression 
    output$gene <- renderUI({
      selectizeInput("gene",
                     label = "Gene of Interest",
                     choices = rownames(rna),
                     multiple = F,
                     options = list(placeholder = 'Select genes')
      )
    })
    output$gene2 <- renderUI({
      selectizeInput("gene2",
                     label = "Gene of Interest",
                     choices = rownames(rna),
                     multiple = T,
                     options = list(maxItems = nrow(rna),placeholder = 'Select genes')
      )
    })

    output$celltype <- renderUI({
      selectizeInput("celltype",
                     label = "Cell types of Interest",
                     choices = cell_predicted(),
                     multiple = T,
                     options = list(placeholder = 'Select cell types')
      )
    })
    
    output$celltype_enri <- renderUI({
      selectizeInput("celltype_enri",
                     label = "Cell types of Interest",
                     choices = cell_predicted(),
                     multiple = T,
                     options = list(placeholder = 'Select cell types')
      )
    })
    
    ######### promoter
      output$gene_pro <- renderUI({
        selectizeInput("gene_pro",
                       label = "Promoter of Interest",
                       choices = rownames(true_promoter),
                       multiple = F,
                       #selected = input$gene,
                       options = list(placeholder = 'Select promoters')
             
        )
      })
      output$gene_pro2 <- renderUI({
        selectizeInput("gene_pro2",
                       label = "Promoter of Interest",
                       choices = rownames(true_promoter),
                       multiple = T,
                       #selected = input$gene2,
                       options = list(maxItems = nrow(true_promoter),placeholder = 'Select promoters')
        )
      })

      output$celltype_pro <- renderUI({
        selectizeInput("celltype_pro",
                       label = "Cell types of Interest",
                       choices = cell_pro(),
                       multiple = T,
                    #selected = input$celltype,
                       options = list(placeholder = 'Select cell types')
        )
      })
      output$gene_clu <- renderUI({
        selectizeInput("gene_clu",
                       label = "Gene of Interest",
                       choices = rownames(rna),
                       multiple = T,
                       options = list(placeholder = 'Select genes')

        )
      })

    
      #######exon
      output$gene_exon <- renderUI({
        selectizeInput("gene_exon",
                       label = "Exon of Interest",
                       choices = rownames(true_exon),
                       multiple = F,
                    #selected = input$gene,
                       options = list(placeholder = 'Select exons')
        )
      })
      output$gene_exon2 <- renderUI({
        selectizeInput("gene_exon2",
                       label = "Exon of Interest",
                       choices = rownames(true_exon),
                       multiple = T,
                      #selected = input$gene2,
                       
                       options = list(maxItems = nrow(true_exon),placeholder = 'Select exons')
        )
      })

      output$celltype_exon <- renderUI({
        selectizeInput("celltype_exon",
                       label = "Cell types of Interest",
                       choices = cell_exon(),#change
                       multiple = T,
                       #selected = input$celltype,
                       options = list(placeholder = 'Select cell types')

        )
      })

      #######Intron
      output$gene_intron <- renderUI({
        selectizeInput("gene_intron",
                       label = "Intron of Interest",
                       choices = rownames(true_intron),
                       multiple = F,
                       #selected = input$gene,
                       options = list(placeholder = 'Select introns')
        )
      })
      output$gene_intron2 <- renderUI({
        selectizeInput("gene_intron2",
                       label = "Intron of Interest",
                       choices = rownames(true_intron),
                       multiple = T,
                       #selected = input$gene2,
                       
                       options = list(maxItems = nrow(true_intron),placeholder = 'Select introns')
        )
      })

      output$celltype_intron <- renderUI({
        selectizeInput("celltype_intron",
                       label = "Cell types of Interest",
                       choices = cell_intron(),
                       multiple = T,
                       #selected = input$celltype,
                       options = list(placeholder = 'Select cell types')
        )
      })

      #######Transcription factor
      output$gene_tf <- renderUI({
        selectizeInput("gene_tf",
                       label = "Transcription factor of Interest",
                       choices = rownames(true_tf),
                       multiple = F,
                       #selected = input$gene,
                       options = list(placeholder = 'Select transcription factors')
        )
      })

      output$gene_tf2 <- renderUI({
        selectizeInput("gene_tf2",
                       label = "Transcription factor of Interest",
                       choices = rownames(true_tf),
                       multiple = T,
                       #selected = input$gene2,
                       options = list(maxItems = nrow(true_tf),placeholder = 'Select transcription factors')
        )
      })
      output$celltype_tf <- renderUI({
        selectizeInput("celltype_tf",
                       label = "Cell types of Interest",
                       choices = cell_tf(),
                       multiple = T,
                       #selected = input$celltype,
                       options = list(placeholder = 'Select cell types'))
      })

      ##enhancer
      output$gene_enhancer <- renderUI({
        selectizeInput("gene_enhancer",
                       label = "Enhancer of Interest",
                       choices = rownames(true_enhancer),
                       multiple = F,
                       #selected = input$gene,
                       
                       options = list(placeholder = 'Select Enhancers')
        )
      })
      
      output$gene_enhancer2 <- renderUI({
        selectizeInput("gene_enhancer2",
                       label = "Enhancer of Interest",
                       choices = rownames(true_enhancer),
                       multiple = T,
                       #selected = input$gene2,
                       options = list(maxItems = nrow(true_enhancer),placeholder = 'Select Enhancers')
        )
      })
      
      output$celltype_enhancer <- renderUI({
        selectizeInput("celltype_enhancer",
                       label = "Cell types of Interest",
                       choices = cell_enhancer(),
                       multiple = T,
                       #selected = input$celltype,
                       options = list(placeholder = 'Select cell types')
        )
      })
      
    
    ## category 
    output$sel_assay <- renderUI({
      selectizeInput("sel_assay",
                     label = "Assay of Interest",
                     choices = batch,
                     multiple = T,
                     options = list(maxItems = length(unique(batch)))
      )
    })
    output$sel_bio <- renderUI({
      selectizeInput("sel_bio",
                     label = "Biosample type of Interest",
                     choices = bio,
                     multiple = T,
                     options = list(maxItems = length(unique(bio)))
      )
    })
    output$sel_cell <- renderUI({
      selectizeInput("sel_cell",
                     label = "Cell type of Interest",
                     choices = cell,
                     multiple = T,
                     options = list(maxItems = length(unique(cell)))
      )
    })
    
    output$sel_layer <- renderUI({
      selectizeInput("sel_layer",
                     label = "Germ Layer of Interest",
                     choices = layer,
                     multiple = T,
                     options = list(maxItems = length(unique(layer)))
      )
    })  
    
    output$sel_cancer <- renderUI({
      selectizeInput("sel_cancer",
                     label = "Tumor & normal cell",
                     choices = cancer,
                     multiple = T,
                     options = list(maxItems = length(unique(cancer)))
      )
    })
    
    output$cate_shape  <- renderUI({
      selectInput("cate_shape", "Select category (shape)",
                  list('Assay' = 'assay',"Biosample type"="bio", "Cell type"="cell","Germ layer"="layer",'cancer' = 'Cancer'),
                  selected = input$cate_color)
    })
    
    pos_sel <- reactive ({
      if (length(input$sel_cell)==0){
        cell_pos = c(1:length(cell))
      }else{
        cell_pos = which(cell%in%input$sel_cell)
      }
      
      if (length(input$sel_assay)==0){
        assay_pos = c(1:length(cell))
      }else{
        assay_pos = which(batch%in%input$sel_assay)
      }
      
      if (length(input$sel_layer)==0){
        layer_pos = c(1:length(cell))
      }else{
        layer_pos = which(layer%in%input$sel_layer)
      }
      if (length(input$sel_cancer)==0){
        cancer_pos = c(1:length(cell))
      }else{
        cancer_pos = which(cancer%in%input$sel_cancer)
      }
      if (length(input$sel_bio)==0){
        bio_pos = c(1:length(cell))
      }else{
        bio_pos = which(bio%in%input$sel_bio)
      }
      Reduce(intersect, list(cell_pos,assay_pos,layer_pos,cancer_pos,bio_pos))
    })
    
    
    ### selected whole data (correction)
    dat_sel = reactive({
      if (input$corr=='uncorrected'){
        d = raw[,pos_sel()]
        colnames(d) = colnames(raw)[pos_sel()]
        d
      }else if (input$corr=='comb'){
        comb[,pos_sel()]
      }else if (input$corr=='mnn'){
        mnn[,pos_sel()]
      }else{
        iter[,pos_sel()]
      }
    })
    
    cell_sel = reactive({
      cell[pos_sel()]
    })
    
    ######### cv cutoff
    keep = reactive({
      if (input$var==FALSE){
        cutoff=1
      }else{
        cutoff = input$quantile
      }
      set.seed(10)
      rtsne=Rtsne(prcomp(t((dat_sel())),scale=T)$x[,1:50])$Y  
      rownames(rtsne)=colnames(dat_sel())
      tsne_raw=as.data.frame(rtsne)
      type = unique(cell_sel()) 
      type=type[type!='stem-cell']
      distance= lapply(c(1:length(type)), function(i) mean(as.matrix(dist(tsne_raw[which(cell_sel()==type[i]),]))) )
      distance=unlist(distance)
      dd=distance[distance!=0]
      cut = quantile(dd,cutoff)
      high_var=which(distance>cut)
      type_del=type[high_var]
      keep_type = which(is.na(match(cell_sel(),type_del)))
      keep_type
      
    })
    
    pos=reactive({
      if (input$var==F){
        pos_sel()  
      }else{
        pos_sel()[keep()]
      }
    })
    
    dat = reactive({
      dat_sel()[,keep()]
    })
    #rna 
    rna_dat = reactive({
      rnaseq[,na.omit(match(pos(),which(batch=='Predicted DNase')))] 
    })
    
    
    ####promoter

      pro_pred_dat = reactive({
          d =  pred_promoter[,na.omit(match(pos(),which(batch=='Predicted DNase')))]
          colnames(d) = colnames(pred_promoter)[na.omit(match(pos(),which(batch=='Predicted DNase')))]
          d
        })

      pro_true_dat = reactive({
        d =true_promoter[,na.omit(match(pos(),which(batch=='True DNase')))]
        colnames(d) = colnames(true_promoter)[na.omit(match(pos(),which(batch=='True DNase')))]
        d
      })

      ##########exon
      exon_pred_dat = reactive({
        d =  pred_exon[,na.omit(match(pos(),which(batch=='Predicted DNase')))]
        colnames(d) = colnames(pred_exon)[na.omit(match(pos(),which(batch=='Predicted DNase')))]
        d
      })

      exon_true_dat = reactive({
        d =true_exon[,na.omit(match(pos(),which(batch=='True DNase')))]
        colnames(d) = colnames(true_exon)[na.omit(match(pos(),which(batch=='True DNase')))]
        d
      })

      ##########intron
      intron_pred_dat = reactive({
        d =  pred_intron[,na.omit(match(pos(),which(batch=='Predicted DNase')))]
        colnames(d) = colnames(pred_intron)[na.omit(match(pos(),which(batch=='Predicted DNase')))]
        d
      })

      intron_true_dat = reactive({
        d =true_intron[,na.omit(match(pos(),which(batch=='True DNase')))]
        colnames(d) = colnames(true_intron)[na.omit(match(pos(),which(batch=='True DNase')))]
        d
      })

      #####################
      ##########enhancer
      enhancer_pred_dat = reactive({
        d =  pred_enhancer[,na.omit(match(pos(),which(batch=='Predicted DNase')))]
        colnames(d) = colnames(pred_enhancer)[na.omit(match(pos(),which(batch=='Predicted DNase')))]
        d
      })

      enhancer_true_dat = reactive({
        d =true_enhancer[,na.omit(match(pos(),which(batch=='True DNase')))]
        colnames(d) = colnames(true_enhancer)[na.omit(match(pos(),which(batch=='True DNase')))]
        d
      })

    ####################
    ##########tf
    tf_pred_dat = reactive({
      d =  pred_tf[,na.omit(match(pos(),which(batch=='Predicted DNase')))]
      colnames(d) = colnames(pred_tf)[na.omit(match(pos(),which(batch=='Predicted DNase')))]
      d
    })

    tf_true_dat = reactive({
      d =true_tf[,na.omit(match(pos(),which(batch=='True DNase')))]
      colnames(d) = colnames(true_tf)[na.omit(match(pos(),which(batch=='True DNase')))]
      d
    })
    # 
    #####################
    #####################
    
    info_var= reactive({
      info[pos(),]
    })
    
    cell_var = reactive({
      info_var()[,1]
    })
    
    layer_var = reactive({
      info_var()[,5]
    }) 
    bio_var = reactive({
      info_var()[,6]
    }) 
    cancer_var = reactive({
      info_var()[,2]
    })
    batch_var = reactive({
      batch[pos()]
    })
    
    bat_var = reactive({
      bat[pos()]
    })
    
    cate_color <- reactive ({
      switch(input$cate_color,
             "assay" = batch_var(),
             "bio" = bio_var(),
             'cell'=cell_var(),
             'layer' = layer_var(), 
             'cancer' = cancer_var())
    })
    
    cate_shape <- reactive ({
      switch(input$cate_shape,
             "assay" = batch_var(),
             "bio" = bio_var(),
             'cell'=cell_var(),
             'layer' = layer_var(), 
             'cancer' = cancer_var())
    })
    
    shape_var=reactive({
      allshape=NULL
      for (j in 1: length(unique(cate_shape()))){
        for (i in 1:length(cate_shape())){
          if (cate_shape()[i] == sort((unique(cate_shape())))[j]){
            allshape[i] = j
          }
        }
      }
      allshape=unique(allshape)
      shape=ifelse(allshape<=18,allshape,allshape-18)
      shape
    })
    

    
    del_cv =  reactive({
      dd = dat()
      cv = apply(dd,1,sd)/rowMeans(dd)
      which(cv>0)
    })
    
    dat_cv =  reactive({
      dd = dat()
      dd[del_cv(),]
    })

    
    pca = reactive({
      pca_unc=prcomp(t((dat_cv())),scale=T)$x
      pca_unc = as.data.frame(pca_unc)
      pca_unc
    })
    
    tsne = reactive({
      set.seed(10)
      rtsne=Rtsne(prcomp(t((dat_cv())),scale=T)$x[,1:50])$Y
      rownames(rtsne)=colnames(dat_cv())
      tsne=as.data.frame(rtsne)
      colnames(tsne)=c('Tsne1','Tsne2')
      tsne
    })
    ## promoter data  
    
    rep_pro = reactive({
      na.omit(match(colnames(pro_true_dat()),colnames(pro_pred_dat())))
    })
    
    pro_dat = reactive({
      if (input$data_type_pro == 'pro_true'){
        pro_true_dat()
      }else if (input$data_type_pro == 'pro_pred'){
        pro_pred_dat()
      }else if (input$data_type_pro == 'pro_both'){
        ### replicate sample, use true data 
        cbind(pro_pred_dat()[,-rep_pro()],pro_true_dat())
      }
    })
    

    ########exon
    rep_exon = reactive({
      na.omit(match(colnames(exon_true_dat()),colnames(exon_pred_dat())))
    })

    exon_dat = reactive({
      if (input$data_type_exon == 'exon_true'){
        exon_true_dat()
      }else if (input$data_type_exon == 'exon_pred'){
        exon_pred_dat()
      }else if (input$data_type_exon == 'exon_both'){
        ### replicate sample, use true data
        cbind(exon_pred_dat()[,-rep_exon()],exon_true_dat())
      }
    })

    ########tf
    rep_tf = reactive({
      na.omit(match(colnames(tf_true_dat()),colnames(tf_pred_dat())))
    })

    tf_dat = reactive({
      if (input$data_type_tf == 'tf_true'){
        tf_true_dat()
      }else if (input$data_type_tf == 'tf_pred'){
        tf_pred_dat()
      }else if (input$data_type_tf == 'tf_both'){
        ### replicate sample, use true data
        cbind(tf_pred_dat()[,-rep_tf()],tf_true_dat())
      }
    })

    ########intron
    rep_intron = reactive({
      na.omit(match(colnames(intron_true_dat()),colnames(intron_pred_dat())))
    })

    intron_dat = reactive({
      if (input$data_type_intron == 'intron_true'){
        intron_true_dat()
      }else if (input$data_type_intron == 'intron_pred'){
        intron_pred_dat()
      }else if (input$data_type_intron == 'intron_both'){
        ### replicate sample, use true data
        cbind(intron_pred_dat()[,-rep_intron()],intron_true_dat())
      }
    })

    ########enhancer
    rep_enhancer = reactive({
      na.omit(match(colnames(enhancer_true_dat()),colnames(enhancer_pred_dat())))
    })

    enhancer_dat = reactive({
      if (input$data_type_enhancer == 'enhancer_true'){
        enhancer_true_dat()
      }else if (input$data_type_enhancer == 'enhancer_pred'){
        enhancer_pred_dat()
      }else if (input$data_type_enhancer == 'enhancer_both'){
        ### replicate sample, use true data
        cbind(enhancer_pred_dat()[,-rep_enhancer()],enhancer_true_dat())
      }
    })
 
    
    pos_sample <- reactive ({
      match(input$search,colnames(dat()))
    })
    
    
    pos_gene <- reactive ({
      match(input$gene,rownames(rna))
    })
    
    gene_curve = reactive({
      match(input$gene2,rownames(rna))
    })
    
    dat_predicted = reactive({
      dat_cv()[,which(batch_var()=='Predicted DNase')]
    })
    
    dat_true = reactive({
      dat_cv()[,which(batch_var()=='True DNase')]
    })
    
    cell_predicted = reactive({
      cell_var()[which(batch_var()=='Predicted DNase')]
    })
    
    cell_true = reactive({
      cell_var()[which(batch_var()=='True DNase')]
    })
    
    
    ####### true dnase data   
      pos_new = reactive({
      c(1:ncol(DNase_new()))
    })

    output$plot1 <- renderPlotly({
      if (!is.null(tsne())){
          if (length(input$search) == 0) {
          plot_ly(tsne(), x = ~Tsne1, y = ~Tsne2,
                  color = cate_color(),colors='Set1',
                  symbol=~cate_shape(),symbols =shape_var(),  marker=list(size=10, opacity=0.7)) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                        '</br> Cell: ', cell_var(),
                                                        '</br> Batch: ', batch_var()))
          }else{
          group = rep('unselected',ncol(dat()))
          group[pos_sample()] = cate_color()[pos_sample()]
          plot_ly(tsne(), x = ~Tsne1, y = ~Tsne2,
                  color = ~group,colors='Set1',
                  symbol= ~cate_shape() ,symbols =shape_var(),  marker=list(size=10, opacity=0.7)) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                        '</br> Cell: ', cell_var(),
                                                        '</br> Batch: ', batch_var()))
        }
      }
    })
    
    
    output$pca1 <- renderPlotly({
      if (!is.null(pca())){
        if (length(input$search) == 0) {
          plot_ly(pca(), x = ~PC1, y = ~PC2,
                  color = cate_color(),colors='Set1',
                  symbol= ~cate_shape() ,symbols =shape_var(),  marker=list(size=10, opacity=0.7)) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                        '</br> Cell: ', cell_var(),
                                                        '</br> Batch: ', batch_var()))
        }else{
          group = rep('unselected',ncol(dat()))
          group[pos_sample()] = cate_color()[pos_sample()]
          plot_ly(pca(), x = ~PC1, y = ~PC2,
                  color = ~group,colors='Set1',
                  symbol= ~cate_shape() ,symbols =shape_var(),  marker=list(size=10, opacity=0.7)) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                        '</br> Cell: ', cell_var(),
                                                        '</br> Batch: ', batch_var()))
        }
      }
    })
    
    mc = reactive({
      cell_mc = cell_predicted() 
      data_mc = dat_predicted()
      dd = data_mc[,which(cell_mc%in%input$celltype)]
      colnames(dd) = colnames(data_mc)[which(cell_mc%in%input$celltype)]
      set.seed(10)
      exprmclust(dd)
    })
    
    order <- reactive({
      TSCANorder(mc(),listbranch = F)
    })
    
    time = reactive({
      data_time = rna_dat()
      if (length(input$gene2) == 1){
        type = rna[gene_curve(),which(colnames(data_time)%in%order())]
        order = na.omit(match(order(),names(type)))
        type[order]
      }else{
        type = rna[gene_curve(),which(colnames(data_time)%in%order())]
        order = na.omit(match(order(),colnames(type)))
        type[,order]
      }
    })
    
    
    long = reactive({
      long = melt(time(), value.name = "gene")
      long
    })
    
    line = reactive({
      line = t(sapply(1:length(input$gene2),function(i) fitted(loess(time()[i,] ~ c(1:ncol(time()))))))
      melt(line)
    })
    # 


    output$curve1 =renderPlotly({
      if (length(input$gene2) == 1){
        line1 = fitted(loess(time() ~ c(1:length(time()))))
        line = data.frame(cbind(Var1 = 1,Var2 = c(1:length(line1)) , value = line1))
        long= data.frame(cbind(Var1 = input$gene2,Var2 = names(time()) , gene = time()))
        rownames(long)=NULL
        long[,2]= line[,2]
        plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='', name= input$celltype,
                marker = list(opacity=0.5,width = 2)) %>% 
          add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers'
          ) %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_predicted()[which(cell_predicted()%in%input$celltype)],
                                                      '</br> Gene: ', long$Var1,
                                                      '</br> Sample: ',long$Var2))%>%
          layout(title = paste0( 'Promoter expression curve of ',input$celltype_pro),
                 xaxis = list(
                   title = "Pseudotime",
                   showticklabels = FALSE,
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')),
                 yaxis = list(
                   title = 'Gene expression',
                   titlefont = list(
                     size = 16,
                     color = 'rgb(107, 107, 107)'),
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')))
        
      }else{
        long=long()
        line=line()
        plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                color = ~Var1, 
                marker = list(opacity=0.5,width = 2)) %>%
          add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',showlegend = F,
                    color = ~long$Var1)  %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_predicted()[which(cell_predicted()%in%input$celltype)],
                                                      '</br> Gene: ', long$Var1,
                                                      '</br> Sample: ',long$Var2))%>%
          layout(title = paste0( 'Gene expression curve of ',input$celltype),
                 xaxis = list(
                   title = "Pseudotime",
                   showticklabels = FALSE,
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')),
                 yaxis = list(
                   title = 'Gene expression',
                   titlefont = list(
                     size = 16,
                     color = 'rgb(107, 107, 107)'),
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')))
      }
    })

    
    output$pca_gene <- renderPlot({
      if (!is.null(input$celltype)){
      dd = dat_predicted()[,which(cell_predicted()%in%input$celltype)]
      colnames(dd) = colnames(data_mc)[which(cell_mc%in%input$celltype)]
      cell = cell_predicted()[which(cell_predicted()%in%input$celltype)]
      pca = as.data.frame(prcomp(t(dd),scale=T)$x)
      batch = (batch_var()[[which(batch_var()=='Predicted DNase')]])[which(cell_predicted()%in%input$celltype)]
      plot_ly(pca, x = ~PC1, y = ~PC2,
              color = cell,colors='Set1',
              symbol= ~cell ,symbols = cell,  marker=list(size=10, opacity=0.7)) %>%
        add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dd),
                                                    '</br> Cell: ', cell,
                                                    '</br> Batch: ', batch))
      }
    })
 
    output$pca_curve <- renderPlot({
      if (!is.null(input$celltype)){
        cell_time = cell_predicted()[which(cell_predicted()%in%input$celltype)]
        mc=mc()
        cell_line=unique(cell_time)
        n.clust = c(1:max(unique(mc[[3]])))
        clu=mc[[3]]
        
        type_clu=lapply(n.clust, function(i) cell_time[which(clu==i)])  # find cell type in each cluster
        n = lapply(n.clust,function(i) length(type_clu[[i]]))
        x = lapply(n.clust,function(i) 
          unlist(lapply(cell_line,function(j) sum(type_clu[[i]]==j)))
        )
        pvalue = lapply(n.clust,function(i) 
          unlist(lapply(c(1:length(cell_line)),function(j) poisson.test(unlist(x[[i]][j]),unlist(n[i]))$p.value))
        )
        pvalue_k=lapply(n.clust, function(i) pvalue[[i]][which(x[[i]]!=0)])
        per = lapply(n.clust,function(i) 
          unlist(lapply(c(1:length(cell_line)),function(j) paste(cell_line[j],':',round(100*unlist(x[[i]][[j]])/unlist(n[i]),2), "%")))
        )
        per_k = lapply(n.clust, function(i) per[[i]][which(x[[i]]!=0)])
        enrich_gene_per = lapply(n.clust,function(i) 
          toString(per_k[[i]][sort(unlist(pvalue_k[[i]]), decreasing = T,index=TRUE)$ix[1:min(3,length(per_k[[i]]))]]))
        #label 
        
        enrich=as.factor(unlist(lapply(c(1:length(cell_time)), function(i)  
          enrich_gene_per[[clu[i]]])))
        
        names(enrich)=names(mc$clusterid)
        mc$clusterid=enrich
        plotmclust(mc,show_cell_names = F)+
          theme(legend.position = "right", legend.key.size = unit(0.3, "in"),
                legend.text = element_text(size = 10),legend.title=element_text(size = 10)) + 
          theme(legend.key = element_blank())+
          labs(title = "PCA pseudo time")
      }
      
    })
    output$heat1 =renderPlotly({
      if (!is.null(tsne())){
        tsne=tsne()[which(batch_var()=='Predicted DNase'),]
        pos = pos()[pos() <= ncol(rna)]
        RNAseq = rna[pos_gene(),pos]
        plot_ly(tsne, x = ~tsne$Tsne1, y = ~tsne$Tsne2, 
                color = ~RNAseq,size=~RNAseq) %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()[which(batch_var()=='Predicted DNase')]),
                                                      '</br> Cell: ', cell_predicted(),
                                                      '</br> Batch: ', batch_var()[which(batch_var()=='Predicted DNase')])) %>%
          layout(title = paste0('Gene expression Tsne heatmap of ',input$gene))
      }
    })

    output$heat_pca1 =renderPlotly({
      if (!is.null(pca())){
        pca = pca()[which(batch_var()=='Predicted DNase'),]
        pos = pos()[pos() <= ncol(rna)]
        RNAseq = rna[pos_gene(),pos]
        plot_ly(pca, x = ~pca$PC1, y = ~pca$PC2, 
                color = ~RNAseq,size=~RNAseq) %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()[which(batch_var()=='Predicted DNase')]),
                                                      '</br> Cell: ', cell_predicted(),
                                                      '</br> Batch: ', batch_var()[which(batch_var()=='Predicted DNase')])) %>%
          layout(title = paste0('Gene expression PCA heatmap of ',input$gene))
      }
    })
    
    ############promoter
      pos_gene_pro <- reactive ({
        match(input$gene_pro,rownames(true_promoter))
      })

      gene_curve_pro = reactive({
        match(input$gene_pro2,rownames(true_promoter))
      })

    cell_pro = reactive({
      if (input$data_type_pro == 'pro_pred'){
        cell_mc = cell_predicted()
      }else if(input$data_type_pro == 'pro_true'){
        cell_mc = cell_true()
      }else{
        cell_mc = cell_var()[-rep_pro()]
      }
      cell_mc
    })

      mc_pro = reactive({
       if (input$data_type_pro == 'pro_true'){
         data_mc = dat_true()
        }else if (input$data_type_pro == 'pro_pred'){
         data_mc = dat_predicted()
        }else if (input$data_type_pro == 'pro_both'){
          ### replicate sample, use true data 
          data_mc = cbind(dat_predicted()[,-rep_pro()],dat_true())
        }
        cell_mc = cell_pro()
        dd = data_mc[,which(cell_mc%in%input$celltype_pro)]
        colnames(dd) = colnames(data_mc)[which(cell_mc%in%input$celltype_pro)]
        set.seed(10)
        exprmclust(dd)
      })

      order_pro <- reactive({
        TSCANorder(mc_pro(),listbranch = F)
      })

      time_pro = reactive({
          data_time = pro_dat()
       if (length(input$gene_pro2) == 1){
          type = pro_dat()[gene_curve_pro(),which(colnames(data_time)%in%order_pro())]
          order = na.omit(match(order_pro(),names(type)))
          type[order]
        }else{
          type = pro_dat()[gene_curve_pro(),which(colnames(data_time)%in%order_pro())]
          order = na.omit(match(order_pro(),colnames(type)))
          type[,order]
        }
      })

      long_pro = reactive({
          long = melt(time_pro(), value.name = "gene")
          long
        })

      line_pro = reactive({
        line = t(sapply(1:length(input$gene_pro2),function(i) fitted(loess(time_pro()[i,] ~ c(1:ncol(time_pro()))))))
        melt(line)
      })
# output$impo = renderPrint({
#   head(long_pro())
# })
# output$curve_text = renderPrint({
#   head(line_pro())
# })
# output$exon_curve_text = renderPrint({
#   dim(time_pro())
# })
# output$intron_curve_text = renderPrint({
#   order_pro()
# })
      output$promoter_curve =renderPlotly({
        if (length(input$gene_pro2) == 1){
          line1 = fitted(loess(time_pro() ~ c(1:length(time_pro()))))
          line = data.frame(cbind(Var1 = 1,Var2 = c(1:length(line1)) , value = line1))
          long= data.frame(cbind(Var1 = input$gene_pro2,Var2 = names(time_pro()) , gene = time_pro()))
          rownames(long)=NULL
          long[,2]= line[,2]
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='', 
                  marker = list(opacity=0.5,width = 2)) %>% 
            add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers',showlegend = F) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_pro()[which(cell_pro()%in%input$celltype_pro)],
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = paste0( 'Promoter expression curve of ',input$celltype_pro),
                   xaxis = list(
                     title = "Pseudotime",
                     showticklabels = FALSE,
                     tickfont = list(
                       size = 14,
                       color = 'rgb(107, 107, 107)')),
                   yaxis = list(
                     title = 'Promoter expression',
                     titlefont = list(
                       size = 16,
                       color = 'rgb(107, 107, 107)'),
                     tickfont = list(
                       size = 14,
                       color = 'rgb(107, 107, 107)')))

        }else{
          long=long_pro()
          line=line_pro()
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                  color = ~Var1,
                  marker = list(opacity=0.5,width = 2)) %>%
            add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',showlegend = F,
                      color = ~long$Var1)  %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_pro()[which(cell_pro()%in%input$celltype_pro)],
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = paste0( 'Promoter expression curve of ',input$celltype_pro),
                   xaxis = list(
                     title = "Pseudotime",
                     showticklabels = FALSE,
                     tickfont = list(
                       size = 14,
                       color = 'rgb(107, 107, 107)')),
                   yaxis = list(
                     title = 'Promoter expression',
                     titlefont = list(
                       size = 16,
                       color = 'rgb(107, 107, 107)'),
                     tickfont = list(
                       size = 14,
                       color = 'rgb(107, 107, 107)')))
        }
      })

      output$promoter_pca_curve <- renderPlot({
        if (!is.null(input$celltype_pro)){
          cell_time = cell_pro()[which(cell_pro()%in%input$celltype_pro)]
          mc=mc_pro()
          cell_line=unique(cell_time)
          n.clust = c(1:max(unique(mc[[3]])))
          clu=mc[[3]]

          type_clu=lapply(n.clust, function(i) cell_time[which(clu==i)])  # find cell type in each cluster
          n = lapply(n.clust,function(i) length(type_clu[[i]]))
          x = lapply(n.clust,function(i)
            unlist(lapply(cell_line,function(j) sum(type_clu[[i]]==j)))
          )
          pvalue = lapply(n.clust,function(i)
            unlist(lapply(c(1:length(cell_line)),function(j) poisson.test(unlist(x[[i]][j]),unlist(n[i]))$p.value))
          )
          pvalue_k=lapply(n.clust, function(i) pvalue[[i]][which(x[[i]]!=0)])
          per = lapply(n.clust,function(i)
            unlist(lapply(c(1:length(cell_line)),function(j) paste(cell_line[j],':',round(100*unlist(x[[i]][[j]])/unlist(n[i]),2), "%")))
          )
          per_k = lapply(n.clust, function(i) per[[i]][which(x[[i]]!=0)])
          enrich_gene_per = lapply(n.clust,function(i)
            toString(per_k[[i]][sort(unlist(pvalue_k[[i]]), decreasing = T,index=TRUE)$ix[1:min(3,length(per_k[[i]]))]]))
          #label

          enrich=as.factor(unlist(lapply(c(1:length(cell_time)), function(i)
            enrich_gene_per[[clu[i]]])))

          names(enrich)=names(mc$clusterid)
          mc$clusterid=enrich
          plotmclust(mc,show_cell_names = F)+
            theme(legend.position = "right", legend.key.size = unit(0.3, "in"),
                  legend.text = element_text(size = 10),legend.title=element_text(size = 10)) +
            theme(legend.key = element_blank())+
            labs(title = "PCA pseudo time")
        }

      })

      
 
      

      output$heat_promoter =renderPlotly({
        if (!is.null(tsne())){
          promoter = cbind(pro_pred_dat(),pro_true_dat())[pos_gene_pro(),]
          plot_ly(tsne(), x = ~Tsne1, y = ~Tsne2,
                  color = ~promoter,size=~promoter) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                        '</br> Cell: ', cell_var(),
                                                        '</br> Batch: ', batch_var())) %>%
          layout(title = paste0('Promoter TSNE heatmap of ',input$gene_pro))
        }
      })

      output$heat_promoter_pca =renderPlotly({
        if (!is.null(pca())){
          promoter = cbind(pro_pred_dat(),pro_true_dat())[pos_gene_pro(),]
          plot_ly(pca(), x = ~PC1, y = ~PC2,
                  color = ~promoter,size=~promoter) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                        '</br> Cell: ', cell_var(),
                                                        '</br> Batch: ', batch_var())) %>%
            layout(title = paste0('Promoter PCA heatmap of ',input$gene_pro))
        }
      })

      ###########exon
      gene_exon2 = reactive({
        input$gene_exon2
      })
      celltype_exon = reactive({
        input$celltype_exon
      })

      pos_gene_exon <- reactive ({
        match(input$gene_exon,rownames(true_exon))
      })

      gene_curve_exon = reactive({
        match(input$gene_exon2,rownames(true_exon))
      })

      cell_exon = reactive({
        if (input$data_type_exon == 'exon_pred'){
          cell_mc = cell_predicted()
        }else if(input$data_type_exon == 'exon_true'){
          cell_mc = cell_true()
        }else{
          cell_mc = cell_var()[-rep_exon()]
        }
        cell_mc
      })

      mc_exon = reactive({
        if (input$data_type_exon == 'exon_true'){
          data_mc = dat_true()
        }else if (input$data_type_exon == 'exon_pred'){
          data_mc = dat_predicted()
        }else if (input$data_type_exon == 'exon_both'){
          data_mc = cbind(dat_predicted()[,-rep_exon()],dat_true())
        }
        cell_mc = cell_exon()
        dd = data_mc[,which(cell_mc%in%input$celltype_exon)]
        colnames(dd) = colnames(data_mc)[which(cell_mc%in%input$celltype_exon)]
        set.seed(10)
        exprmclust(dd)
      })

      order_exon <- reactive({
        TSCANorder(mc_exon(),listbranch = F)
      })

      time_exon = reactive({
        data_time = exon_dat()
        if (length(input$gene_exon2) == 1){
          type = exon_dat()[gene_curve_exon(),which(colnames(data_time)%in%order_exon())]
          order = na.omit(match(order_exon(),names(type)))
          type[order]
        }else{
          type = exon_dat()[gene_curve_exon(),which(colnames(data_time)%in%order_exon())]
          order = na.omit(match(order_exon(),colnames(type)))
          type[,order]
        }
      })

      long_exon = reactive({
        long = melt(time_exon(), value.name = "gene")
        long
      })


      line_exon = reactive({
        line = t(sapply(1:length(input$gene_exon2),function(i) fitted(loess(time_exon()[i,] ~ c(1:ncol(time_exon()))))))
        #  rownames(line) = colnames(time())
        melt(line)
      })


      output$exon_curve =renderPlotly({
        if (length(input$gene_exon2) == 1){
          line1 = fitted(loess(time_exon() ~ c(1:length(time_exon()))))
          line = data.frame(cbind(Var1 = 1,Var2 = c(1:length(line1)) , value = line1))
          long= data.frame(cbind(Var1 = input$gene_exon2,Var2 = names(time_exon()) , gene = time_exon()))
          rownames(long)=NULL
          long[,2]= line[,2]
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='', 
                  marker = list(opacity=0.5,width = 2)) %>% 
            add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers'  , showlegend = F 
            ) %>%
            
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_exon()[which(cell_exon()%in%input$celltype_exon)],
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = paste0( 'Exon expression curve of ',input$celltype_exon),
                   xaxis = list(
                     title = "Pseudotime",
                     showticklabels = FALSE,
                     tickfont = list(
                       size = 14,
                       color = 'rgb(107, 107, 107)')),
                   yaxis = list(
                     title = 'Exon expression',
                     titlefont = list(
                       size = 16,
                       color = 'rgb(107, 107, 107)'),
                     tickfont = list(
                       size = 14,
                       color = 'rgb(107, 107, 107)')))

        }else{
          long=long_exon()
          line=line_exon()
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='', 
                  color = ~Var1,
                  marker = list(opacity=0.5,width = 2)) %>%
            add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',   showlegend = F,
                      color = ~long$Var1 )  %>%         
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_exon()[which(cell_exon()%in%input$celltype_exon)],
                                                                                                  '</br> Gene: ', long$Var1,
                                                                                                  '</br> Sample: ',long$Var2))%>%
            layout(title = paste0( 'Exon expression curve of ',input$celltype_exon),
                   xaxis = list(
                     title = "Pseudotime",
                     showticklabels = FALSE,
                     tickfont = list(
                       size = 14,
                       color = 'rgb(107, 107, 107)')),
                   yaxis = list(
                     title = 'Exon expression',
                     titlefont = list(
                       size = 16,
                       color = 'rgb(107, 107, 107)'),
                     tickfont = list(
                       size = 14,
                       color = 'rgb(107, 107, 107)')))
        }
      })

      output$exon_pca_curve <- renderPlot({
        if (!is.null(input$celltype_exon)){
          cell_time = cell_exon()[which(cell_exon()%in%input$celltype_exon)]
          mc=mc_exon()
          cell_line=unique(cell_time)
          n.clust = c(1:max(unique(mc[[3]])))
          clu=mc[[3]]

          type_clu=lapply(n.clust, function(i) cell_time[which(clu==i)])  # find cell type in each cluster
          n = lapply(n.clust,function(i) length(type_clu[[i]]))
          x = lapply(n.clust,function(i)
            unlist(lapply(cell_line,function(j) sum(type_clu[[i]]==j)))
          )
          pvalue = lapply(n.clust,function(i)
            unlist(lapply(c(1:length(cell_line)),function(j) poisson.test(unlist(x[[i]][j]),unlist(n[i]))$p.value))
          )
          pvalue_k=lapply(n.clust, function(i) pvalue[[i]][which(x[[i]]!=0)])
          per = lapply(n.clust,function(i)
            unlist(lapply(c(1:length(cell_line)),function(j) paste(cell_line[j],':',round(100*unlist(x[[i]][[j]])/unlist(n[i]),2), "%")))
          )
          per_k = lapply(n.clust, function(i) per[[i]][which(x[[i]]!=0)])
          enrich_gene_per = lapply(n.clust,function(i)
            toString(per_k[[i]][sort(unlist(pvalue_k[[i]]), decreasing = T,index=TRUE)$ix[1:min(3,length(per_k[[i]]))]]))
          #label

          enrich=as.factor(unlist(lapply(c(1:length(cell_time)), function(i)
            enrich_gene_per[[clu[i]]])))

          names(enrich)=names(mc$clusterid)
          mc$clusterid=enrich
          plotmclust(mc,show_cell_names = F)+
            theme(legend.position = "right", legend.key.size = unit(0.3, "in"),
                  legend.text = element_text(size = 10),legend.title=element_text(size = 10)) +
            theme(legend.key = element_blank())+
            labs(title = "PCA pseudo time")
        }

      })

      tsne_exon = reactive({
        if (input$data_type_exon == 'exon_true'){
          tsne()[which(batch_var()=='True DNase'),]
        }else if (input$data_type_exon == 'exon_pred'){
          tsne()[which(batch_var()=='Predicted DNase'),]
        }else if (input$data_type_exon == 'exon_both'){
          ### replicate sample, use true data
          tsne()[-rep_exon(),]
        }
      })

      output$heat_exon =renderPlotly({
        if (!is.null(tsne())){
          Exon = cbind(exon_pred_dat(),exon_true_dat())[pos_gene_exon(),]
          tsne=tsne()
          plot_ly(tsne, x = ~Tsne1, y = ~Tsne2,
                  color = ~Exon,size=~Exon) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                        '</br> Cell: ', cell_var(),
                                                        '</br> Batch: ', batch_var())) %>%
            layout(title = paste0('Exon Tsne heatmap of ',input$gene_exon))
        }
      })

      output$heat_exon_pca =renderPlotly({
        if (!is.null(pca())){
          Exon = cbind(exon_pred_dat(),exon_true_dat())[pos_gene_exon(),]
          pca=pca()
          plot_ly(pca, x = ~PC1, y = ~PC2,
                  color = ~Exon,size=~Exon) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                        '</br> Cell: ', cell_var(),
                                                        '</br> Batch: ', batch_var())) %>%
            layout(title = paste0('Exon PCA heatmap of ',input$gene_exon))
        }
      })


      #################intron
      ###########intron

      pos_gene_intron <- reactive ({
        match(input$gene_intron,rownames(true_intron))
      })

      gene_curve_intron = reactive({
        match(input$gene_intron2,rownames(true_intron))
      })

      cell_intron = reactive({
        if (input$data_type_intron == 'intron_pred'){
          cell_mc = cell_predicted()
        }else if(input$data_type_intron == 'intron_true'){
          cell_mc = cell_true()
        }else{
          cell_mc = cell_var()[-rep_intron()]
        }
        cell_mc
      })

      mc_intron = reactive({
        if (input$data_type_intron == 'intron_true'){
          data_mc = dat_true()
        }else if (input$data_type_intron == 'intron_pred'){
          data_mc = dat_predicted()
        }else if (input$data_type_intron == 'intron_both'){
          ### replicate sample, use true data 
          data_mc = cbind(dat_predicted()[,-rep_intron()],dat_true())
        }
        cell_mc = cell_intron()
        dd = data_mc[,which(cell_mc%in%input$celltype_intron)]
        colnames(dd) = colnames(data_mc)[which(cell_mc%in%input$celltype_intron)]
        set.seed(10)
        exprmclust(dd)
      })

      order_intron <- reactive({
        TSCANorder(mc_intron(),listbranch = F)
      })

      time_intron = reactive({
        data_time = intron_dat()
        if (length(input$gene_intron2) == 1){
          type = intron_dat()[gene_curve_intron(),which(colnames(data_time)%in%order_intron())]
          order = na.omit(match(order_intron(),names(type)))
          type[order]
        }else{
          type = intron_dat()[gene_curve_intron(),which(colnames(data_time)%in%order_intron())]
          order = na.omit(match(order_intron(),colnames(type)))
          type[,order]
        }
      })

      long_intron = reactive({
        long = melt(time_intron(), value.name = "gene")
        long
      })


      line_intron = reactive({
        line = t(sapply(1:length(input$gene_intron2),function(i) fitted(loess(time_intron()[i,] ~ c(1:ncol(time_intron()))))))
        #  rownames(line) = colnames(time())
        melt(line)
      })


      output$intron_curve =renderPlotly({
        if (length(input$gene_intron2) == 1){
          line1 = fitted(loess(time_intron() ~ c(1:length(time_intron()))))
          line = data.frame(cbind(Var1 = 1,Var2 = c(1:length(line1)) , value = line1))
          long= data.frame(cbind(Var1 = input$gene_intron2,Var2 = names(time_intron()) , gene = time_intron()))
          rownames(long)=NULL
          long[,2]= line[,2]
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',  
                  marker = list(opacity=0.5,width = 2)) %>% 
            add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers' , showlegend = F) %>%
            
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_intron()[which(cell_intron()%in%input$celltype_intron)],
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = paste0( 'Intron curve of ',input$celltype_intron),
                   xaxis = list(
                     title = "Pseudotime",
                     showticklabels = FALSE,
                     tickfont = list(
                       size = 14,
                       color = 'rgb(107, 107, 107)')),
                   yaxis = list(
                     title = 'Intron expression',
                     titlefont = list(
                       size = 16,
                       color = 'rgb(107, 107, 107)'),
                     tickfont = list(
                       size = 14,
                       color = 'rgb(107, 107, 107)')))

        }else{
          long=long_intron()
          line=line_intron()
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='', 
                  color = ~Var1,
                  marker = list(opacity=0.5,width = 2)) %>%
            add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',   showlegend = F,
                      color = ~long$Var1 )  %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_intron()[which(cell_intron()%in%input$celltype_intron)],
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
              layout(title = paste0( 'Intron expression curve of ',input$celltype_intron),
                   xaxis = list(
                     title = "Pseudotime",
                     showticklabels = FALSE,
                     tickfont = list(
                       size = 14,
                       color = 'rgb(107, 107, 107)')),
                   yaxis = list(
                     title = 'Gene expression',
                     titlefont = list(
                       size = 16,
                       color = 'rgb(107, 107, 107)'),
                     tickfont = list(
                       size = 14,
                       color = 'rgb(107, 107, 107)')))
        }
      })

      output$intron_pca_curve <- renderPlot({
        if (!is.null(input$celltype_intron)){
          cell_time = cell_intron()[which(cell_intron()%in%input$celltype_intron)]
          mc=mc_intron()
          cell_line=unique(cell_time)
          n.clust = c(1:max(unique(mc[[3]])))
          clu=mc[[3]]

          type_clu=lapply(n.clust, function(i) cell_time[which(clu==i)])  # find cell type in each cluster
          n = lapply(n.clust,function(i) length(type_clu[[i]]))
          x = lapply(n.clust,function(i)
            unlist(lapply(cell_line,function(j) sum(type_clu[[i]]==j)))
          )
          pvalue = lapply(n.clust,function(i)
            unlist(lapply(c(1:length(cell_line)),function(j) poisson.test(unlist(x[[i]][j]),unlist(n[i]))$p.value))
          )
          pvalue_k=lapply(n.clust, function(i) pvalue[[i]][which(x[[i]]!=0)])
          per = lapply(n.clust,function(i)
            unlist(lapply(c(1:length(cell_line)),function(j) paste(cell_line[j],':',round(100*unlist(x[[i]][[j]])/unlist(n[i]),2), "%")))
          )
          per_k = lapply(n.clust, function(i) per[[i]][which(x[[i]]!=0)])
          enrich_gene_per = lapply(n.clust,function(i)
            toString(per_k[[i]][sort(unlist(pvalue_k[[i]]), decreasing = T,index=TRUE)$ix[1:min(3,length(per_k[[i]]))]]))
          #label

          enrich=as.factor(unlist(lapply(c(1:length(cell_time)), function(i)
            enrich_gene_per[[clu[i]]])))

          names(enrich)=names(mc$clusterid)
          mc$clusterid=enrich
          plotmclust(mc,show_cell_names = F)+
            theme(legend.position = "right", legend.key.size = unit(0.3, "in"),
                  legend.text = element_text(size = 10),legend.title=element_text(size = 10)) +
            theme(legend.key = element_blank())+
            labs(title = "PCA pseudo time")
        }

      })

      tsne_intron = reactive({
        if (input$data_type_intron == 'intron_true'){
          tsne()[which(batch_var()=='True DNase'),]
        }else if (input$data_type_intron == 'intron_pred'){
          tsne()[which(batch_var()=='Predicted DNase'),]
        }else if (input$data_type_intron == 'intron_both'){
          ### replicate sample, use true data
          tsne()[-rep_intron(),]
        }
      })

      output$heat_intron =renderPlotly({
        if (!is.null(tsne())){
          tsne=tsne()
          Intron = cbind(intron_pred_dat(),intron_true_dat())[pos_gene_intron(),]
          plot_ly(tsne, x = ~Tsne1, y = ~Tsne2,
                  color = ~Intron,size=~Intron) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                        '</br> Cell: ', cell_var(),
                                                        '</br> Batch: ', batch_var())) %>%
            layout(title = paste0('Intron TSNE heatmap of ',input$gene_intron))
        }
      })

      output$heat_intron_pca =renderPlotly({
        if (!is.null(pca())){
          pca=pca()
          Intron = cbind(intron_pred_dat(),intron_true_dat())[pos_gene_intron(),]
          plot_ly(pca, x = ~PC1, y = ~PC2,
                  color = ~Intron,size=~Intron) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                        '</br> Cell: ', cell_var(),
                                                        '</br> Batch: ', batch_var())) %>%
            layout(title = paste0('Intron PCA heatmap of ',input$gene_intron))
        }
      })

      ###########enhancer

      pos_gene_enhancer <- reactive ({
        match(input$gene_enhancer,rownames(true_enhancer))
      })

      gene_curve_enhancer = reactive({
        match(input$gene_enhancer2,rownames(true_enhancer))
      })

      cell_enhancer = reactive({
        if (input$data_type_enhancer == 'enhancer_pred'){
          cell_mc = cell_predicted()
        }else if(input$data_type_enhancer == 'enhancer_true'){
          cell_mc = cell_true()
        }else{
          cell_mc = cell_var()[-rep_enhancer()]
        }
        cell_mc
      })

      mc_enhancer = reactive({
        if (input$data_type_enhancer == 'enhancer_true'){
          data_mc = dat_true()
        }else if (input$data_type_enhancer == 'enhancer_pred'){
          data_mc = dat_predicted()
        }else if (input$data_type_enhancer == 'enhancer_both'){
          ### replicate sample, use true data 
          data_mc = cbind(dat_predicted()[,-rep_enhancer()],dat_true())
        }
        cell_mc = cell_enhancer()
        dd = data_mc[,which(cell_mc%in%input$celltype_enhancer)]
        colnames(dd) = colnames(data_mc)[which(cell_mc%in%input$celltype_enhancer)]
        set.seed(10)
        exprmclust(dd)
      })

      order_enhancer <- reactive({
        TSCANorder(mc_enhancer(),listbranch = F)
      })

      time_enhancer = reactive({
        data_time = enhancer_dat()
        if (length(input$gene_enhancer2) == 1){
          type = enhancer_dat()[gene_curve_enhancer(),which(colnames(data_time)%in%order_enhancer())]
          order = na.omit(match(order_enhancer(),names(type)))
          type[order]
        }else{
          type = enhancer_dat()[gene_curve_enhancer(),which(colnames(data_time)%in%order_enhancer())]
          order = na.omit(match(order_enhancer(),colnames(type)))
          type[,order]
        }
      })

      long_enhancer = reactive({
        long = melt(time_enhancer(), value.name = "gene")
        long
      })


      line_enhancer = reactive({
        line = t(sapply(1:length(input$gene_enhancer2),function(i) fitted(loess(time_enhancer()[i,] ~ c(1:ncol(time_enhancer()))))))
        #  rownames(line) = colnames(time())
        melt(line)
      })


      output$enhancer_curve =renderPlotly({
        if (length(input$gene_enhancer2) == 1){
          line1 = fitted(loess(time_enhancer() ~ c(1:length(time_enhancer()))))
          line = data.frame(cbind(Var1 = 1,Var2 = c(1:length(line1)) , value = line1))
          long= data.frame(cbind(Var1 = input$gene_enhancer2,Var2 = names(time_enhancer()) , gene = time_enhancer()))
          rownames(long)=NULL
          long[,2]= line[,2]
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',showlegend = F,
                  marker = list(opacity=0.5,width = 2)) %>% 
            add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers'
            )  %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_enhancer()[which(cell_enhancer()%in%input$celltype_enhancer)],
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = paste0( 'Enhancer curve of ',input$celltype_enhancer),
                   xaxis = list(
                     title = "Pseudotime",
                     showticklabels = FALSE,
                     tickfont = list(
                       size = 14,
                       color = 'rgb(107, 107, 107)')),
                   yaxis = list(
                     title = 'Enhancer expression',
                     titlefont = list(
                       size = 16,
                       color = 'rgb(107, 107, 107)'),
                     tickfont = list(
                       size = 14,
                       color = 'rgb(107, 107, 107)')))

        }else{
          long=long_enhancer()
          line=line_enhancer()
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                  color = ~Var1,
               
                  marker = list(opacity=0.5,width = 2)) %>%
            add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',   showlegend = F,
                      color = ~long$Var1)  %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_enhancer()[which(cell_enhancer()%in%input$celltype_enhancer)],
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = paste0( 'Enhancer expression curve of ',input$celltype_enhancer),
                   xaxis = list(
                     title = "Pseudotime",
                     #showticklabels = FALSE,
                     tickfont = list(
                       size = 14,
                       color = 'rgb(107, 107, 107)')),
                   yaxis = list(
                     title = 'Enhancer expression',
                     titlefont = list(
                       size = 16,
                       color = 'rgb(107, 107, 107)'),
                     tickfont = list(
                       size = 14,
                       color = 'rgb(107, 107, 107)')))
        }
      })

      output$enhancer_pca_curve <- renderPlot({
        if (!is.null(input$celltype_enhancer)){
          cell_time = cell_enhancer()[which(cell_enhancer()%in%input$celltype_enhancer)]
          mc=mc_enhancer()
          cell_line=unique(cell_time)
          n.clust = c(1:max(unique(mc[[3]])))
          clu=mc[[3]]

          type_clu=lapply(n.clust, function(i) cell_time[which(clu==i)])  # find cell type in each cluster
          n = lapply(n.clust,function(i) length(type_clu[[i]]))
          x = lapply(n.clust,function(i)
            unlist(lapply(cell_line,function(j) sum(type_clu[[i]]==j)))
          )
          pvalue = lapply(n.clust,function(i)
            unlist(lapply(c(1:length(cell_line)),function(j) poisson.test(unlist(x[[i]][j]),unlist(n[i]))$p.value))
          )
          pvalue_k=lapply(n.clust, function(i) pvalue[[i]][which(x[[i]]!=0)])
          per = lapply(n.clust,function(i)
            unlist(lapply(c(1:length(cell_line)),function(j) paste(cell_line[j],':',round(100*unlist(x[[i]][[j]])/unlist(n[i]),2), "%")))
          )
          per_k = lapply(n.clust, function(i) per[[i]][which(x[[i]]!=0)])
          enrich_gene_per = lapply(n.clust,function(i)
            toString(per_k[[i]][sort(unlist(pvalue_k[[i]]), decreasing = T,index=TRUE)$ix[1:min(3,length(per_k[[i]]))]]))
          #label

          enrich=as.factor(unlist(lapply(c(1:length(cell_time)), function(i)
            enrich_gene_per[[clu[i]]])))

          names(enrich)=names(mc$clusterid)
          mc$clusterid=enrich
          plotmclust(mc,show_cell_names = F)+
            theme(legend.position = "right", legend.key.size = unit(0.3, "in"),
                  legend.text = element_text(size = 10),legend.title=element_text(size = 10)) +
            theme(legend.key = element_blank())+
            labs(title = "PCA pseudo time")
        }

      })

      tsne_enhancer = reactive({
        if (input$data_type_enhancer == 'enhancer_true'){
          tsne()[which(batch_var()=='True DNase'),]
        }else if (input$data_type_enhancer == 'enhancer_pred'){
          tsne()[which(batch_var()=='Predicted DNase'),]
        }else if (input$data_type_enhancer == 'enhancer_both'){
          ### replicate sample, use true data
          tsne()[-rep_enhancer(),]
        }
      })

      output$heat_enhancer =renderPlotly({
        if (!is.null(tsne())){
          tsne=tsne()
          Enhancer = cbind(enhancer_pred_dat(),enhancer_true_dat())[pos_gene_enhancer(),]
          plot_ly(tsne, x = ~Tsne1, y = ~Tsne2,
                  color = ~Enhancer,size=~Enhancer) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                        '</br> Cell: ', cell_var(),
                                                        '</br> Batch: ', batch_var())) %>%
            layout(title = paste0('Enhancer TSNE heatmap of ',input$gene_enhancer))
        }
      })

      output$heat_enhancer_pca =renderPlotly({
        if (!is.null(pca())){
          pca=pca()
          Enhancer = cbind(enhancer_pred_dat(),enhancer_true_dat())[pos_gene_enhancer(),]
          plot_ly(pca, x = ~PC1, y = ~PC2,
                  color = ~Enhancer,size=~Enhancer) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                        '</br> Cell: ', cell_var(),
                                                        '</br> Batch: ', batch_var())) %>%
            layout(title = paste0('Enhancer PCA heatmap of ',input$gene_enhancer))
        }
      })

    ###########Transcription factor
    pos_gene_tf <- reactive ({
      match(input$gene_tf,rownames(true_tf))
    })

    gene_curve_tf = reactive({
      match(input$gene_tf2,rownames(true_tf))
    })

    cell_tf = reactive({
      if (input$data_type_tf == 'tf_pred'){
        cell_mc = cell_predicted()
      }else if(input$data_type_tf == 'tf_true'){
        cell_mc = cell_true()
      }else{
        cell_mc = cell_var()[-rep_tf()]
      }
      cell_mc
    })

    mc_tf = reactive({
      if (input$data_type_tf == 'tf_true'){
        data_mc = dat_true()
      }else if (input$data_type_tf == 'tf_pred'){
        data_mc = dat_predicted()
      }else if (input$data_type_tf == 'tf_both'){
        ### replicate sample, use true data 
        data_mc = cbind(dat_predicted()[,-rep_tf()],dat_true())
      }
      cell_mc = cell_tf()
      dd = data_mc[,which(cell_mc%in%input$celltype_tf)]
      colnames(dd) = colnames(data_mc)[which(cell_mc%in%input$celltype_tf)]
      set.seed(10)
      exprmclust(dd)
    })

    order_tf <- reactive({
      TSCANorder(mc_tf(),listbranch = F)
    })

    time_tf = reactive({
      data_time = tf_dat()
      if (length(input$gene_tf2) == 1){
        type = tf_dat()[gene_curve_tf(),which(colnames(data_time)%in%order_tf())]
        order = na.omit(match(order_tf(),names(type)))
        type[order]
      }else{
        type = tf_dat()[gene_curve_tf(),which(colnames(data_time)%in%order_tf())]
        order = na.omit(match(order_tf(),colnames(type)))
        type[,order]
      }
    })

    long_tf = reactive({
      long = melt(time_tf(), value.name = "gene")
      long
    })


    line_tf = reactive({
      line = t(sapply(1:length(input$gene_tf2),function(i) fitted(loess(time_tf()[i,] ~ c(1:ncol(time_tf()))))))
      #  rownames(line) = colnames(time())
      melt(line)
    })


    output$tf_curve =renderPlotly({
      if (length(input$gene_tf2) == 1){
        line1 = fitted(loess(time_tf() ~ c(1:length(time_tf()))))
        line = data.frame(cbind(Var1 = 1,Var2 = c(1:length(line1)) , value = line1))
        long= data.frame(cbind(Var1 = input$gene_tf2,Var2 = names(time_tf()) , gene = time_tf()))
        rownames(long)=NULL
        long[,2]= line[,2]
        plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',showlegend = F,
                marker = list(opacity=0.5,width = 2)) %>% 
          add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers'
          ) %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_tf()[which(cell_tf()%in%input$celltype_tf)],
                                                      '</br> Gene: ', long$Var1,
                                                      '</br> Sample: ',long$Var2))%>%
          layout(title = paste0( 'Transcription factor curve of ',input$celltype_tf),
                 xaxis = list(
                   title = "Pseudotime",
                   showticklabels = FALSE,
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')),
                 yaxis = list(
                   title = 'Transcription factor expression',
                   titlefont = list(
                     size = 16,
                     color = 'rgb(107, 107, 107)'),
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')))

      }else{
        long=long_tf()
        line=line_tf()
        plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                color = ~Var1,
                marker = list(opacity=0.5,width = 2)) %>%
          add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',
                    color = ~long$Var1,
                    showlegend = F)  %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_tf()[which(cell_tf()%in%input$celltype_tf)],
                                                      '</br> Gene: ', long$Var1,
                                                      '</br> Sample: ',long$Var2))%>%
          layout(title = paste0( 'Transcription factor expression curve of ',input$celltype_tf),
                 xaxis = list(
                   title = "Pseudotime",
                   showticklabels = FALSE,
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')),
                 yaxis = list(
                   title = 'Transcription factor expression',
                   titlefont = list(
                     size = 16,
                     color = 'rgb(107, 107, 107)'),
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')))
      }
    })

    output$tf_pca_curve <- renderPlot({
      if (!is.null(input$celltype_tf)){
        cell_time = cell_tf()[which(cell_tf()%in%input$celltype_tf)]
        mc=mc_tf()
        cell_line=unique(cell_time)
        n.clust = c(1:max(unique(mc[[3]])))
        clu=mc[[3]]

        type_clu=lapply(n.clust, function(i) cell_time[which(clu==i)])  # find cell type in each cluster
        n = lapply(n.clust,function(i) length(type_clu[[i]]))
        x = lapply(n.clust,function(i)
          unlist(lapply(cell_line,function(j) sum(type_clu[[i]]==j)))
        )
        pvalue = lapply(n.clust,function(i)
          unlist(lapply(c(1:length(cell_line)),function(j) poisson.test(unlist(x[[i]][j]),unlist(n[i]))$p.value))
        )
        pvalue_k=lapply(n.clust, function(i) pvalue[[i]][which(x[[i]]!=0)])
        per = lapply(n.clust,function(i)
          unlist(lapply(c(1:length(cell_line)),function(j) paste(cell_line[j],':',round(100*unlist(x[[i]][[j]])/unlist(n[i]),2), "%")))
        )
        per_k = lapply(n.clust, function(i) per[[i]][which(x[[i]]!=0)])
        enrich_gene_per = lapply(n.clust,function(i)
          toString(per_k[[i]][sort(unlist(pvalue_k[[i]]), decreasing = T,index=TRUE)$ix[1:min(3,length(per_k[[i]]))]]))
        #label

        enrich=as.factor(unlist(lapply(c(1:length(cell_time)), function(i)
          enrich_gene_per[[clu[i]]])))

        names(enrich)=names(mc$clusterid)
        mc$clusterid=enrich
        plotmclust(mc,show_cell_names = F)+
          theme(legend.position = "right", legend.key.size = unit(0.3, "in"),
                legend.text = element_text(size = 10),legend.title=element_text(size = 10)) +
          theme(legend.key = element_blank())+
          labs(title = "PCA pseudo time")
      }

    })

    tsne_tf = reactive({
      if (input$data_type_tf == 'tf_true'){
        tsne()[which(batch_var()=='True DNase'),]
      }else if (input$data_type_tf == 'tf_pred'){
        tsne()[which(batch_var()=='Predicted DNase'),]
      }else if (input$data_type_tf == 'tf_both'){
        ### replicate sample, use true data
        tsne()[-rep_tf(),]
      }
    })

    output$heat_tf =renderPlotly({
      if (!is.null(tsne())){
        tsne=tsne()
        TF = cbind(tf_pred_dat(),tf_true_dat())[pos_gene_tf(),]
        plot_ly(tsne, x = ~Tsne1, y = ~Tsne2,
                color = ~TF,size=~TF) %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                      '</br> Cell: ', cell_var(),
                                                      '</br> Batch: ', batch_var())) %>%
          layout(title = paste0('Transcription factor  TSNE heatmap of ',input$gene_tf))

      }
    })

    output$heat_tf_pca =renderPlotly({
      if (!is.null(pca())){
        pca=pca()
        TF = cbind(tf_pred_dat(),tf_true_dat())[pos_gene_tf(),]
        plot_ly(pca, x = ~PC1, y = ~PC2,
                color = ~TF,size=~TF) %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                      '</br> Cell: ', cell_var(),
                                                      '</br> Batch: ', batch_var())) %>%
          layout(title = paste0('Transcription factor  PCA heatmap of ',input$gene_tf))
      }
    })
    

####################### enrichment 
    pos_rna = reactive({
      which(cell_predicted()%in%(input$celltype_enri))
    })
    cell_predicted_chose =  reactive({
      cell_predicted()[pos_rna()]
    })
    rna_chose =  reactive({
      rna[,pos_rna()]
    })
    
    raw_rna_chose = reactive({
      raw[,pos_rna()]
    })
    
    ####################### find order (use predicted dnase to do peudo time )
    rna_order = reactive({
      cell_mc = cell_predicted_chose()
      set.seed(10)
      mc = exprmclust(raw_rna_chose())
      order <- TSCANorder(mc,listbranch = F)
      pos_order = match(order,colnames(rna_chose()))
      rna_chose()[,pos_order]
    })
    
    fit = reactive({
      set.seed(1)
      kmeans(rna_order(), 20)     
    })
    
    clu = reactive({
      fit()$cluster
    })
    
    clu_mean =  reactive({
      t(sapply(c(1:max(clu())), function(i) colMeans(rna_order()[which(clu()==i),])))
    })
    

    line_enri = reactive({
      line= t(sapply(c(1:nrow(clu_mean())),function(i) 
                 fitted(loess(clu_mean()[i,] ~ c(1:ncol(clu_mean()))))))
      melt(line)
    })
    
    
    long_enri = reactive({
      melt(clu_mean(), value.name = "gene")
    })
    

    output$rna_curve =renderPlotly({
      if (length(input$celltype_enri) != 0){
        long=long_enri()
        line=line_enri()
        long[,2] = line[,2]
        plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                color = ~as.factor(Var1),showlegend = FALSE,
                marker = list(opacity=0.5,width = 2))  %>%
          add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',
                    color = ~as.factor(long$Var1))  %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_predicted()[which(cell_predicted()%in%input$celltype_enri)],
                                                      '</br> Cluster: ', long$Var1,
                                                      '</br> Sample: ',long_enri()$Var2))%>%
          layout(title = paste0( 'Clustered gene expression curve of ',input$celltype_enri),
                 xaxis = list(
                   title = "Pseudotime",
                   showticklabels = FALSE,
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')),
                 yaxis = list(
                   title = 'Clustered gene expression',
                   titlefont = list(
                     size = 16,
                     color = 'rgb(107, 107, 107)'),
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')))
      }
    })
    
    output$gene_cluster <- renderPrint({
      paste('Gene' , input$gene_clu ,'is in cluster ', clu()[which(names(clu())==input$gene_clu)])
    })
    
######## find corresponding cluster 
    ########## keep only same genes 
    gene_pos = reactive({
      na.omit(match(rownames(rna),rownames(true_promoter)))
    })
      
    del_enri = reactive({
      which(is.na(match(rownames(rna),rownames(true_promoter))))
    })
    
    rna_data_gene = reactive({
      rna[-del_enri(),]
    })
    
    clu_gene = reactive({
      clu()[-del_enri()]
    })
    
    pos_dnase = reactive({
      which(cell_enri()%in%input$celltype_enri)
    })
    

     
     enri_promoter = reactive({
         if (input$enri_type == 'true'){
           pro_true_dat()
         }else if (input$enri_type == 'pred'){
           pro_pred_dat()
         }else if (input$enri_type == 'both'){
           cbind(pro_pred_dat()[,-rep_pro()],pro_true_dat())
         }
       })
     
     dat_enri = reactive({
       if (input$enri_type == 'true'){
         dat_true()
       }else if (input$enri_type == 'pred'){
         dat_predicted()
       }else if (input$enri_type == 'both'){
         cbind(dat_predicted()[,-rep_tf()],dat_true())
       }
     })
     
    promoter_chose_gene = reactive({
      enri_promoter()[gene_pos(),pos_dnase()]
    })
    

    
    promoter_chose_raw = reactive({
      dat_enri()[,pos_dnase()]
    })
    
    ## construct peusdo time  
    promoter_order = reactive({
      set.seed(10)
      promoter_mc = exprmclust(promoter_chose_raw())
      order_d <- TSCANorder(promoter_mc,listbranch = F)
      pos_order_d = match(order_d,colnames(promoter_chose_gene()))
      promoter_chose_gene()[,pos_order_d]
     })

    clu_mean_promoter =  reactive({
      t(sapply(c(1:max(clu_gene())), function(i) colMeans(promoter_order()[which(clu_gene()==i),])))
    })
    
    
    line_d1 = reactive({
      line = t(sapply(c(1:20),function(i) 
          fitted(loess(clu_mean_promoter()[i,] ~ c(1:ncol(clu_mean_promoter()))))))
      melt(line)
    })
    
    long_d1 = reactive({
      melt(clu_mean_promoter(), value.name = "gene")
    })
    

    output$promoter_enri =renderPlotly({
      if (length(input$celltype_enri) != 0){
        long=long_d1()
        line=line_d1()
        long[,2]= line[,2]
        plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                color = ~as.factor(Var1),
                marker = list(opacity=0.5,width = 2))  %>%
          add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',showlegend = FALSE,
                    color = ~as.factor(long$Var1))  %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_predicted()[which(cell_predicted()%in%input$celltype_enri)],
                                                      '</br> Cluster: ', long$Var1,
                                                      '</br> Sample: ',long_d1()$Var2))%>%
          layout(title = paste0( 'Clustered promoter expression curve of ',input$celltype_enri),
                 xaxis = list(
                   title = "Pseudotime",
                   showticklabels = FALSE,
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')),
                 yaxis = list(
                   title = 'Clustered promoter expression',
                   titlefont = list(
                     size = 16,
                     color = 'rgb(107, 107, 107)'),
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')))
      }
    })
 
    ########## choose cluster 
    ########### do GO 
    output$go <- renderPrint({
      genelist = clu()
      GOdata <- new("topGOdata", ontology = "BP", allGenes = genelist, geneSel = function(p) p == input$clus, description = "Test", annot = annFUN.org, mapping = "org.Mm.eg.db", 
            ID = "Symbol")
      resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
      GenTable(GOdata, classicFisher = resultFisher)
    })
    
    ######## find corresponding cluster exon
    ########## keep only same genes 
    gene_pos_exon = reactive({
      na.omit(match(rownames(rna),rownames(true_exon)))
    })
    
    del_enri_exon = reactive({
      which(is.na(match(rownames(rna),rownames(true_exon))))
    })
    
    rna_exon_gene = reactive({
      rna[-del_enri_exon(),]
    })
    
    clu_gene_exon = reactive({
      clu()[-del_enri_exon()]
    })
    
    enri_exon = reactive({
      if (input$enri_type == 'true'){
        exon_true_dat()
      }else if (input$enri_type == 'pred'){
        exon_pred_dat()
      }else if (input$enri_type == 'both'){
        cbind(exon_pred_dat()[,-rep_exon()],exon_true_dat())
      }
    })
    
    exon_chose_gene= reactive({
      enri_exon()[gene_pos_exon(),pos_dnase()]
    })
    
    exon_chose_raw = reactive({
      dat_enri()[,pos_dnase()]
    })
    
    ## construct peusdo time  
    exon_order = reactive({
      set.seed(10)
      exon_mc = exprmclust(exon_chose_raw())
      order_d <- TSCANorder(exon_mc,listbranch = F)
      pos_order_d = match(order_d,colnames(exon_chose_gene()))
      exon_chose_gene()[,pos_order_d]
    })
    
    clu_mean_exon =  reactive({
      t(sapply(c(1:max(clu_gene_exon())), function(i) colMeans(exon_order()[which(clu_gene_exon()==i),])))
    })
    
    
    line_exon_clu = reactive({
      line = t(sapply(c(1:20),function(i) 
        fitted(loess(clu_mean_exon()[i,] ~ c(1:ncol(clu_mean_exon()))))))
      melt(line)
    })
    
    long_exon_clu = reactive({
      melt(clu_mean_exon(), value.name = "gene")
    })
    

    output$exon_enri =renderPlotly({
      if (length(input$celltype_enri) != 0){
        long=long_exon_clu()
        line=line_exon_clu()
        long[,2]= line[,2]
        plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                color = ~as.factor(Var1),showlegend = FALSE,
                marker = list(opacity=0.5,width = 2))  %>%
          add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',
                    color = ~as.factor(long$Var1))  %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_predicted()[which(cell_predicted()%in%input$celltype_enri)],
                                                      '</br> Cluster: ', long$Var1,
                                                      '</br> Sample: ',long_exon_clu()$Var2))%>%
          layout(title = paste0( 'Clustered exon expression curve of ',input$celltype_enri),
                 xaxis = list(
                   title = "Pseudotime",
                   showticklabels = FALSE,
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')),
                 yaxis = list(
                   title = 'Clustered exon expression',
                   titlefont = list(
                     size = 16,
                     color = 'rgb(107, 107, 107)'),
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')))
      }
    })    
    
    
    ######## find corresponding cluster exon
    ########## keep only same genes 
    gene_pos_intron = reactive({
      na.omit(match(rownames(rna),rownames(true_intron)))
    })
    
    del_enri_intron = reactive({
      which(is.na(match(rownames(rna),rownames(true_intron))))
    })
    
    rna_intron_gene = reactive({
      rna[-del_enri_intron(),]
    })
    
    clu_gene_intron = reactive({
      clu()[-del_enri_intron()]
    })
    
    enri_intron = reactive({
      if (input$enri_type == 'true'){
        intron_true_dat()
      }else if (input$enri_type == 'pred'){
        intron_pred_dat()
      }else if (input$enri_type == 'both'){
        cbind(intron_pred_dat()[,-rep_intron()],intron_true_dat())
      }
    })
    
    intron_chose_gene= reactive({
      enri_intron()[gene_pos_intron(),pos_dnase()]
    })
    
    intron_chose_raw = reactive({
      dat_enri()[,pos_dnase()]
    })
    
    
    ## construct peusdo time  
    intron_order = reactive({
      set.seed(10)
      intron_mc = exprmclust(intron_chose_raw())
      order_d <- TSCANorder(intron_mc,listbranch = F)
      pos_order_d = match(order_d,colnames(intron_chose_gene()))
      intron_chose_gene()[,pos_order_d]
    })
    
    clu_mean_intron =  reactive({
      t(sapply(c(1:max(clu_gene_intron())), function(i) colMeans(intron_order()[which(clu_gene_intron()==i),])))
    })
    
    
    line_intron_clu = reactive({
      line = t(sapply(c(1:20),function(i) 
        fitted(loess(clu_mean_intron()[i,] ~ c(1:ncol(clu_mean_intron()))))))
      melt(line)
    })
    
    long_intron_clu = reactive({
      melt(clu_mean_intron(), value.name = "gene")
    })

    output$intron_enri =renderPlotly({
      if (length(input$celltype_enri) != 0){
        long=long_intron_clu()
        line=line_intron_clu()
        long[,2]= line[,2]
        plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                color = ~as.factor(Var1),showlegend = FALSE,
                marker = list(opacity=0.5,width = 2))  %>%
          add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',
                    color = ~as.factor(long$Var1))  %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_predicted()[which(cell_predicted()%in%input$celltype_enri)],
                                                      '</br> Cluster: ', long$Var1,
                                                      '</br> Sample: ',long_intron_clu()$Var2))%>%
          layout(title = paste0( 'Clustered intron expression curve of ',input$celltype_enri),
                 xaxis = list(
                   title = "Pseudotime",
                   showticklabels = FALSE,
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')),
                 yaxis = list(
                   title = 'Clustered intron expression',
                   titlefont = list(
                     size = 16,
                     color = 'rgb(107, 107, 107)'),
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')))
      }
    })    
    
    
}

shinyApp(ui = ui, server = server)
