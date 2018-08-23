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
                                
                                h3("Gene expression plots"),
                                checkboxInput("geneexpre", 'Gene expression heatmap'),
                                conditionalPanel(
                                  condition = "input.geneexpre == true",
                                  uiOutput('gene')
                                ),
                                checkboxInput("genecurve", 'Gene expression curve'),
                                conditionalPanel(
                                  condition = "input.genecurve == true",
                                  radioButtons(inputId="data_type_gene", "Data type:",
                                               choices = list("RNA-seq" = "gene_rna"),inline=TRUE),
                                  
                                  uiOutput('gene2'),
                                  uiOutput('celltype')
                                  ),
                                
                                h3("Promoter plots"),
                                checkboxInput("promoter", 'Promoter heatmap'),
                                conditionalPanel(
                                  condition = "input.promoter == true",
                                  uiOutput('gene_pro')
                                ),

                                checkboxInput("procurve", 'Promoter curve'),
                                conditionalPanel(
                                  condition = "input.procurve == true",
                                  radioButtons(inputId="data_type_pro", "Data type:",
                                               choices = list("True DNase" = "pro_true","Predicted DNase" = "pro_pred",
                                                              "Predicted DNase + True DNase-seq" = "pro_both"),inline=TRUE),
                                  uiOutput('gene_pro2'),
                                  uiOutput('celltype_pro')
                                  ),

                                h3("Exon plots"),
                                checkboxInput("exon", 'Exon heatmap'),
                                conditionalPanel(
                                  condition = "input.exon == true",
                                  uiOutput('gene_exon')
                                ),

                                checkboxInput("exoncurve", 'Exon curve'),
                                conditionalPanel(
                                  condition = "input.exoncurve == true",
                                  radioButtons(inputId="data_type_exon", "Data type:",
                                               choices = list("True DNase" = "exon_true","Predicted DNase" = "exon_pred",
                                                              "Predicted DNase + True DNase-seq" = "exon_both"),inline=TRUE),
                                  uiOutput('gene_exon2'),
                                  uiOutput('celltype_exon')
                                  ),

                                h3("Intron plots"),
                                checkboxInput("intron", 'Intron heatmap'),
                                conditionalPanel(
                                  condition = "input.intron == true",
                                  uiOutput('gene_intron')
                                ),

                                checkboxInput("introncurve", 'Intron curve'),
                                conditionalPanel(
                                  condition = "input.introncurve == true",
                                  radioButtons(inputId="data_type_intron", "Data type:",
                                               choices = list("True DNase" = "intron_true","Predicted DNase" = "intron_pred",
                                                              "Predicted DNase + True DNase-seq" = "intron_both"),inline=TRUE),
                                  uiOutput('gene_intron2'),
                                  uiOutput('celltype_intron')
                                  ),

                                h3("Enhancer plots"),
                                checkboxInput("enhancer", 'Enhancer heatmap'),
                                conditionalPanel(
                                  condition = "input.enhancer == true",
                                  uiOutput('gene_enhancer')
                                ),
                                
                                checkboxInput("enhancercurve", 'Enhancer curve'),
                                conditionalPanel(
                                  condition = "input.enhancercurve == true",
                                  radioButtons(inputId="data_type_enhancer", "Data type:",
                                               choices = list("True DNase" = "enhancer_true","Predicted DNase" = "enhancer_pred",
                                                              "Predicted DNase + True DNase-seq" = "enhancer_both"),inline=TRUE),
                                  uiOutput('gene_enhancer2'),
                                  uiOutput('celltype_enhancer')
                                  ),
                                
                    

                                h3("Transcription factor plots"),
                                checkboxInput("tf", 'Transcription factor heatmap'),
                                conditionalPanel(
                                  condition = "input.tf == true",
                                  uiOutput('gene_tf')
                                ),

                                checkboxInput("tfcurve", 'Transcription factor curve'),
                                conditionalPanel(
                                  condition = "input.tfcurve == true",
                                  radioButtons(inputId="data_type_tf", "Data type:",
                                               choices = list("True DNase" = "tf_true","Predicted DNase" = "tf_pred",
                                                              "Predicted DNase + True DNase-seq" = "tf_both"),inline=TRUE),
                                  uiOutput('gene_tf2'),
                                  uiOutput('celltype_tf')
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
                                  
                                  tabPanel("Gene expression", 
                                           h4('Gene expression heatmap'),
                                           plotlyOutput('heat')%>% withSpinner(),
                                           plotlyOutput('heat_pca')%>% withSpinner()
                                  ),
                                  tabPanel("Gene expression curve", 
                                           h4('Gene expression curve'),
                                           verbatimTextOutput("impo"),
                                           
                                           verbatimTextOutput('curve_text'),
                                           verbatimTextOutput('exon_curve_text'),
                                           verbatimTextOutput('intron_curve_text'),
                                           
                                           plotlyOutput('curve')%>% withSpinner(),
                                           plotOutput('pca_curve')%>% withSpinner(),
                                           plotlyOutput('pca_gene')%>% withSpinner()
                                           
                                  ),
                                  
                                  tabPanel("Promoter",
                                           h4('Promoter heatmap'),
                                           plotlyOutput('heat_promoter')%>% withSpinner(),
                                           plotlyOutput('heat_promoter_pca')%>% withSpinner()
                                  ),
                                  tabPanel("Promoter curve",
                                           h4('Promoter curve'),
                                           plotlyOutput('promoter_curve')%>% withSpinner(),
                                           plotOutput('promoter_pca_curve')%>% withSpinner()

                                  ),
                                  tabPanel("Exon",
                                           h4('Exon heatmap'),
                                           plotlyOutput('heat_exon')%>% withSpinner(),
                                           plotlyOutput('heat_exon_pca')%>% withSpinner()
                                  ),
                                  tabPanel("Exon curve",
                                           h4('Exon curve'),

                                           plotlyOutput('exon_curve')%>% withSpinner(),
                                           plotOutput('exon_pca_curve')%>% withSpinner()

                                  ),

                                  tabPanel("Intron",
                                           h4('Intron heatmap'),
                                           plotlyOutput('heat_intron')%>% withSpinner(),
                                           plotlyOutput('heat_intron_pca')%>% withSpinner()
                                  ),
                                  tabPanel("Intron curve",
                                           h4('Intron curve'),
                                           
                                           plotlyOutput('intron_curve')%>% withSpinner(),
                                           plotOutput('intron_pca_curve')%>% withSpinner()

                                  ),
                                  tabPanel("Enhancer",
                                           h4('Enhancer heatmap'),
                                           plotlyOutput('heat_enhancer')%>% withSpinner(),
                                           plotlyOutput('heat_enhancer_pca')%>% withSpinner()
                                  ),
                                  tabPanel("Enhancer curve",
                                           h4('Enhancer curve'),
                                           plotlyOutput('enhancer_curve')%>% withSpinner(),
                                           plotOutput('enhancer_pca_curve')%>% withSpinner()
                                           
                                  ),
                                  
                                  tabPanel("Transcription factor",
                                           h4('Transcription factor heatmap'),
                                           plotlyOutput('heat_tf')%>% withSpinner(),
                                           plotlyOutput('heat_tf_pca')%>% withSpinner()
                                  ),
                                  tabPanel("Transcription factor curve",
                                           h4('Transcription factor curve'),
                                           plotlyOutput('tf_curve')%>% withSpinner(),
                                           plotOutput('tf_pca_curve')%>% withSpinner()

                                  ),
                                  tabPanel("Correction methods evaluation",
                                           # 
                
                                          # plotOutput('dist'),
                                           plotOutput('silhouette'),
                                          # plotOutput('entropies'),
                                           plotOutput('k_unc'),
                                           plotOutput('k_comb'),
                                           plotOutput('k_mnn'),
                                           plotOutput('k_iter')
                                  )
                                  
                                )# end of tabset panel
                              )# end of main panel
                              
                            )
                          )
                          ))
)




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
    gene2 = reactive({
      input$gene2
    })
    output$celltype <- renderUI({
      selectizeInput("celltype",
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
                       options = list(placeholder = 'Select promoters')
        )
      })
      output$gene_pro2 <- renderUI({
        selectizeInput("gene_pro2",
                       label = "Promoter of Interest",
                       choices = rownames(true_promoter),
                       multiple = T,
                       options = list(maxItems = nrow(true_promoter),placeholder = 'Select promoters')
        )
      })

      output$celltype_pro <- renderUI({
        selectizeInput("celltype_pro",
                       label = "Cell types of Interest",
                       choices = cell_pro(),
                       multiple = T,
                       options = list(placeholder = 'Select cell types')
        )
      })


    
    # 
      #######exon
      output$gene_exon <- renderUI({
        selectizeInput("gene_exon",
                       label = "Exon of Interest",
                       choices = rownames(true_exon),
                       multiple = F,
                       options = list(placeholder = 'Select exons')
        )
      })
      output$gene_exon2 <- renderUI({
        selectizeInput("gene_exon2",
                       label = "Exon of Interest",
                       choices = rownames(true_exon),
                       multiple = T,
                       options = list(maxItems = nrow(true_exon),placeholder = 'Select exons')
        )
      })

      output$celltype_exon <- renderUI({
        selectizeInput("celltype_exon",
                       label = "Cell types of Interest",
                       choices = cell_exon(),#change
                       multiple = T,
                       options = list(placeholder = 'Select cell types')
        )
      })
    # 
      #######Intron
      output$gene_intron <- renderUI({
        selectizeInput("gene_intron",
                       label = "Intron of Interest",
                       choices = rownames(true_intron),
                       multiple = F,
                       options = list(placeholder = 'Select introns')
        )
      })
      output$gene_intron2 <- renderUI({
        selectizeInput("gene_intron2",
                       label = "Intron of Interest",
                       choices = rownames(true_intron),
                       multiple = T,
                       options = list(maxItems = nrow(true_intron),placeholder = 'Select introns')
        )
      })

      output$celltype_intron <- renderUI({
        selectizeInput("celltype_intron",
                       label = "Cell types of Interest",
                       choices = cell_intron(),#change
                       multiple = T,
                       options = list(placeholder = 'Select cell types')
        )
      })

      #######Transcription factor
      output$gene_tf <- renderUI({
        selectizeInput("gene_tf",
                       label = "Transcription factor of Interest",
                       choices = rownames(true_tf),
                       multiple = F,
                       options = list(placeholder = 'Select transcription factors')
        )
      })

      output$gene_tf2 <- renderUI({
        selectizeInput("gene_tf2",
                       label = "Transcription factor of Interest",
                       choices = rownames(true_tf),
                       multiple = T,
                       options = list(maxItems = nrow(true_tf),placeholder = 'Select transcription factors')
        )
      })
      output$celltype_tf <- renderUI({
        selectizeInput("celltype_tf",
                       label = "Cell types of Interest",
                       choices = cell_tf(),
                       multiple = T,
                       options = list(placeholder = 'Select cell types')
        )
      })

      ##enhancer
      output$gene_enhancer <- renderUI({
        selectizeInput("gene_enhancer",
                       label = "Enhancer of Interest",
                       choices = rownames(true_enhancer),
                       multiple = F,
                       options = list(placeholder = 'Select Enhancers')
        )
      })
      
      output$gene_enhancer2 <- renderUI({
        selectizeInput("gene_enhancer2",
                       label = "Enhancer of Interest",
                       choices = rownames(true_enhancer),
                       multiple = T,
                       options = list(maxItems = nrow(true_enhancer),placeholder = 'Select Enhancers')
        )
      })
      
      output$celltype_enhancer <- renderUI({
        selectizeInput("celltype_enhancer",
                       label = "Cell types of Interest",
                       choices = cell_enhancer(),
                       multiple = T,
                       options = list(placeholder = 'Select cell types')
        )
      })
      

    
    celltype=reactive({
      input$celltype
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
      match(gene2(),rownames(rna))
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
    
    
    #################################
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
      dd = data_mc[,which(cell_mc%in%celltype())]
      colnames(dd) = colnames(data_mc)[which(cell_mc%in%celltype())]
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
      line = t(sapply(1:length(gene2()),function(i) fitted(loess(time()[i,] ~ c(1:ncol(time()))))))
      melt(line)
    })
    # 
    # output$impo <- renderPrint({
    #   cell_predicted()[which(cell_predicted()%in%celltype())]
    # })
    # output$curve_text <- renderPrint({
    #   cell_pro()[which(cell_pro()%in%celltype_pro())]
    # })
    # 
    # output$exon_curve_text <- renderPrint({
    #   cell_exon()[which(cell_exon()%in%celltype_pro())]
    # })
    # output$intron_curve_text <- renderPrint({
    #  mc_pro()
    # })
    
    output$curve =renderPlotly({
      if (length(input$gene2) == 1){
        line1 = fitted(loess(time() ~ c(1:length(time()))))
        line = data.frame(cbind(Var1 = 1,Var2 = c(1:length(line1)) , value = line1))
        long= data.frame(cbind(Var1 = input$gene2,Var2 = names(time()) , gene = time()))
        rownames(long)=NULL
        long[,2]= line[,2]
        plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='', name=  'WASH7P',
                marker = list(opacity=0.5,width = 2)) %>% 
          add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers'
          ) %>%          
          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_predicted()[which(cell_predicted()%in%celltype())],
                                                                     '</br> Gene: ', long$Var1,
                                                                     '</br> Sample: ',long$Var2))%>%
          layout(title = paste0( 'Gene expression curve of ',celltype()),
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
                color = ~ Var1,
                marker = list(opacity=0.5,width = 2)) %>%
          add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',
                    color = ~long$Var1,
                    showlegend = F)  %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_predicted()[which(cell_predicted()%in%celltype())],
                                                      '</br> Gene: ', long$Var1,
                                                      '</br> Sample: ',long$Var2))%>%
          layout(title = paste0( 'Gene expression curve of ',celltype()),
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
      dd = dat_predicted()[,which(cell_predicted()%in%celltype())]
      colnames(dd) = colnames(data_mc)[which(cell_mc%in%celltype())]
      cell = cell_predicted()[which(cell_predicted()%in%celltype())]
      pca = as.data.frame(prcomp(t(dd),scale=T)$x)
      batch = (batch_var()[[which(batch_var()=='Predicted DNase')]])[which(cell_predicted()%in%celltype())]
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
        cell_time = cell_predicted()[which(cell_predicted()%in%celltype())]
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
    
    output$heat =renderPlotly({
      if (!is.null(tsne())){
        tsne=tsne()[which(batch_var()=='Predicted DNase'),]
        pos = pos()[pos() <= ncol(rna)]
        z = rna[pos_gene(),pos]
        cell = cell_var()[which(batch_var()=='Predicted DNase')]
        x = tsne$Tsne1
        y =tsne$Tsne2
        s <- interp(x,y,z)
        d <- melt(s$z, na.rm = TRUE)
        names(d) <- c("x", "y", "z")
        d$Tsne1 <- s$x[d$x]
        d$Tsne2 <- s$y[d$y]
        mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
        g = ggplot(data = d, aes(x = Tsne1, y = Tsne2, fill = z,text = paste("Gene expression level:",z)))+
          geom_raster()+
          scale_fill_gradientn(colours = mycol,limits=range(z),'Gene expression level')+
          theme_classic()+
          ggtitle(paste0('Gene expression Tsne heatmap of ',input$gene))
        ggplotly(g,tooltip = "text")
      }
    })
    
    output$heat_pca =renderPlotly({
      if (!is.null(pca())){
        
        pca=pca()[which(batch_var()=='Predicted DNase'),]
        pos = pos()[pos() <= ncol(rna)]
        z = rna[pos_gene(),pos]
        cell = cell_var()[which(batch_var()=='Predicted DNase')]
        x = pca$PC1
        y =pca$PC2
        s <- interp(x,y,z)
        d <- melt(s$z, na.rm = TRUE)
        names(d) <- c("x", "y", "z")
        d$PC1 <- s$x[d$x]
        d$PC2 <- s$y[d$y]
        mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
        g = ggplot(data = d, aes(x = PC1, y = PC2, fill = z,text = paste("Gene expression level:",z)))+
          geom_raster()+
          scale_fill_gradientn(colours = mycol,limits=range(z),'Gene expression level')+
          theme_classic()+
          ggtitle(paste0('Gene expression PCA heatmap of ',input$gene))
        ggplotly(g,tooltip = "text")
      }
    })
    ############promoter
      gene_pro2 = reactive({
        input$gene_pro2
      })
      celltype_pro = reactive({
        input$celltype_pro
      })

      pos_gene_pro <- reactive ({
        match(input$gene_pro,rownames(true_promoter))
      })

      gene_curve_pro = reactive({
        match(gene_pro2(),rownames(true_promoter))
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
        dd = data_mc[,which(cell_mc%in%celltype_pro())]
        colnames(dd) = colnames(data_mc)[which(cell_mc%in%celltype_pro())]
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
        line = t(sapply(1:length(gene_pro2()),function(i) fitted(loess(time_pro()[i,] ~ c(1:ncol(time_pro()))))))
        #  rownames(line) = colnames(time())
        melt(line)
      })


      output$promoter_curve =renderPlotly({
        if (length(input$gene_pro2) == 1){
          line1 = fitted(loess(time_pro() ~ c(1:length(time_pro()))))
          line = data.frame(cbind(Var1 = 1,Var2 = c(1:length(line1)) , value = line1))
          long= data.frame(cbind(Var1 = input$gene_pro2,Var2 = names(time_pro()) , gene = time_pro()))
          rownames(long)=NULL
          long[,2]= line[,2]
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='', name=  'WASH7P',
                  marker = list(opacity=0.5,width = 2)) %>% 
            add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers'
            ) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_pro()[which(cell_pro()%in%celltype_pro())],
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = paste0( 'Promoter expression curve of ',celltype_pro()),
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
          long=long_pro()
          line=line_pro()
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                  color = ~Var1,
                  marker = list(opacity=0.5,width = 2)) %>%
            add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',
                      color = ~long$Var1,
                      showlegend = F)  %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_pro()[which(cell_pro()%in%celltype_pro())],
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = paste0( 'Promoter expression curve of ',celltype_pro()),
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

      output$promoter_pca_curve <- renderPlot({
        if (!is.null(input$celltype_pro)){
          cell_time = cell_pro()[which(cell_pro()%in%celltype_pro())]
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
          tsne=tsne()
          z = cbind(pro_pred_dat(),pro_true_dat())[pos_gene_pro(),]
          cell = cell_var()
          x = tsne$Tsne1
          y = tsne$Tsne2
          s <- interp(x,y,z)
          d <- melt(s$z, na.rm = TRUE)
          names(d) <- c("x", "y", "z")
          d$Tsne1 <- s$x[d$x]
          d$Tsne2 <- s$y[d$y]
          mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
          g = ggplot(data = d, aes(x = Tsne1, y = Tsne2, fill = z,text = paste("Promoter level:",z)))+
            geom_raster()+
            scale_fill_gradientn(colours = mycol,limits=range(z),'Promoter level')+
            theme_classic()+
            ggtitle(paste0('Promoter Tsne heatmap of ',input$gene_pro))
          ggplotly(g,tooltip = "text")
        }
      })

      output$heat_promoter_pca =renderPlotly({
        if (!is.null(pca())){
          pca=pca()
          z =cbind(pro_pred_dat(),pro_true_dat())[pos_gene_pro(),]
          cell = cell_var()
          x = pca$PC1
          y =pca$PC2
          s <- interp(x,y,z)
          d <- melt(s$z, na.rm = TRUE)
          names(d) <- c("x", "y", "z")
          d$PC1 <- s$x[d$x]
          d$PC2 <- s$y[d$y]
          mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
          g = ggplot(data = d, aes(x = PC1, y = PC2, fill = z,text = paste("Promoter level:",z)))+
            geom_raster()+
            scale_fill_gradientn(colours = mycol,limits=range(z),'Promoter level')+
            theme_classic()+
            ggtitle(paste0('Promoter PCA heatmap of ',input$gene_pro))
          ggplotly(g,tooltip = "text")
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
        match(gene_exon2(),rownames(true_exon))
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
        dd = data_mc[,which(cell_mc%in%celltype_exon())]
        colnames(dd) = colnames(data_mc)[which(cell_mc%in%celltype_exon())]
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
        line = t(sapply(1:length(gene_exon2()),function(i) fitted(loess(time_exon()[i,] ~ c(1:ncol(time_exon()))))))
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
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='', name=  'WASH7P',
                  marker = list(opacity=0.5,width = 2)) %>% 
            add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers'
            ) %>%
            
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_exon()[which(cell_exon()%in%celltype_exon())],
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = paste0( 'Exon expression curve of ',celltype_exon()),
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
            add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',
                      color = ~long$Var1,
                      showlegend = F)  %>%         
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_exon()[which(cell_exon()%in%celltype_exon())],
                                                                                                  '</br> Gene: ', long$Var1,
                                                                                                  '</br> Sample: ',long$Var2))%>%
            layout(title = paste0( 'Exon expression curve of ',celltype_exon()),
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
          cell_time = cell_exon()[which(cell_exon()%in%celltype_exon())]
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
          tsne=tsne()
          z = cbind(exon_pred_dat(),exon_true_dat())[pos_gene_exon(),]
          cell = cell_var()
          x = tsne$Tsne1
          y = tsne$Tsne2
          s <- interp(x,y,z)
          d <- melt(s$z, na.rm = TRUE)
          names(d) <- c("x", "y", "z")
          d$Tsne1 <- s$x[d$x]
          d$Tsne2 <- s$y[d$y]
          mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
          g = ggplot(data = d, aes(x = Tsne1, y = Tsne2, fill = z,text = paste("Exon level:",z)))+
            geom_raster()+
            scale_fill_gradientn(colours = mycol,limits=range(z),'Exon level')+
            theme_classic()+
            ggtitle(paste0('Exon Tsne heatmap of ',input$gene_exon))
          ggplotly(g,tooltip = "text")
        }
      })

      output$heat_exon_pca =renderPlotly({
        if (!is.null(pca())){
          pca=pca()
          z =cbind(exon_pred_dat(),exon_true_dat())[pos_gene_exon(),]
          cell = cell_var()
          x = pca$PC1
          y =pca$PC2
          s <- interp(x,y,z)
          d <- melt(s$z, na.rm = TRUE)
          names(d) <- c("x", "y", "z")
          d$PC1 <- s$x[d$x]
          d$PC2 <- s$y[d$y]
          mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
          g = ggplot(data = d, aes(x = PC1, y = PC2, fill = z,text = paste("Exon level:",z)))+
            geom_raster()+
            scale_fill_gradientn(colours = mycol,limits=range(z),'Exon level')+
            theme_classic()+
            ggtitle(paste0('Exon PCA heatmap of ',input$gene_exon))
          ggplotly(g,tooltip = "text")
        }
      })


      #################intron
      ###########intron
      gene_intron2 = reactive({
        input$gene_intron2
      })
      celltype_intron = reactive({
        input$celltype_intron
      })

      pos_gene_intron <- reactive ({
        match(input$gene_intron,rownames(true_intron))
      })

      gene_curve_intron = reactive({
        match(gene_intron2(),rownames(true_intron))
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
        dd = data_mc[,which(cell_mc%in%celltype_intron())]
        colnames(dd) = colnames(data_mc)[which(cell_mc%in%celltype_intron())]
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
        line = t(sapply(1:length(gene_intron2()),function(i) fitted(loess(time_intron()[i,] ~ c(1:ncol(time_intron()))))))
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
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='', name=  'WASH7P',
                  marker = list(opacity=0.5,width = 2)) %>% 
            add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers'
            ) %>%
            
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_intron()[which(cell_intron()%in%celltype_intron())],
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = paste0( 'Intron curve of ',celltype_intron()),
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
            add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',
                      color = ~long$Var1,
                      showlegend = F)  %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_intron()[which(cell_intron()%in%celltype_intron())],
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
              layout(title = paste0( 'Intron expression curve of ',celltype_intron()),
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
          cell_time = cell_intron()[which(cell_intron()%in%celltype_intron())]
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
          z = cbind(intron_pred_dat(),intron_true_dat())[pos_gene_intron(),]
          cell = cell_var()
          x = tsne$Tsne1
          y = tsne$Tsne2
          s <- interp(x,y,z)
          d <- melt(s$z, na.rm = TRUE)
          names(d) <- c("x", "y", "z")
          d$Tsne1 <- s$x[d$x]
          d$Tsne2 <- s$y[d$y]
          mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
          g = ggplot(data = d, aes(x = Tsne1, y = Tsne2, fill = z,text = paste("Intron level:",z)))+
            geom_raster()+
            scale_fill_gradientn(colours = mycol,limits=range(z),'Intron level')+
            theme_classic()+
            ggtitle(paste0('Intron Tsne heatmap of ',input$gene_intron))
          ggplotly(g,tooltip = "text")
        }
      })

      output$heat_intron_pca =renderPlotly({
        if (!is.null(pca())){
          pca=pca()
          z =cbind(intron_pred_dat(),intron_true_dat())[pos_gene_intron(),]
          cell = cell_var()
          x = pca$PC1
          y =pca$PC2
          s <- interp(x,y,z)
          d <- melt(s$z, na.rm = TRUE)
          names(d) <- c("x", "y", "z")
          d$PC1 <- s$x[d$x]
          d$PC2 <- s$y[d$y]
          mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
          g = ggplot(data = d, aes(x = PC1, y = PC2, fill = z,text = paste("Intron level:",z)))+
            geom_raster()+
            scale_fill_gradientn(colours = mycol,limits=range(z),'Intron level')+
            theme_classic()+
            ggtitle(paste0('Intron PCA heatmap of ',input$gene_intron))
          ggplotly(g,tooltip = "text")
        }
      })

      ###########enhancer
      gene_enhancer2 = reactive({
        input$gene_enhancer2
      })
      celltype_enhancer = reactive({
        input$celltype_enhancer
      })

      pos_gene_enhancer <- reactive ({
        match(input$gene_enhancer,rownames(true_enhancer))
      })

      gene_curve_enhancer = reactive({
        match(gene_enhancer2(),rownames(true_enhancer))
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
        dd = data_mc[,which(cell_mc%in%celltype_enhancer())]
        colnames(dd) = colnames(data_mc)[which(cell_mc%in%celltype_enhancer())]
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
        line = t(sapply(1:length(gene_enhancer2()),function(i) fitted(loess(time_enhancer()[i,] ~ c(1:ncol(time_enhancer()))))))
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
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='', name=  'WASH7P',
                  marker = list(opacity=0.5,width = 2)) %>% 
            add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers'
            )  %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_enhancer()[which(cell_enhancer()%in%celltype_enhancer())],
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = paste0( 'Enhancer curve of ',celltype_enhancer()),
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
            add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',
                      color = ~long$Var1,
                      showlegend = F)  %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_enhancer()[which(cell_enhancer()%in%celltype_enhancer())],
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = paste0( 'Enhancer expression curve of ',celltype_enhancer()),
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
          cell_time = cell_enhancer()[which(cell_enhancer()%in%celltype_enhancer())]
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
          z = cbind(enhancer_pred_dat(),enhancer_true_dat())[pos_gene_enhancer(),]
          cell = cell_var()
          x = tsne$Tsne1
          y = tsne$Tsne2
          s <- interp(x,y,z)
          d <- melt(s$z, na.rm = TRUE)
          names(d) <- c("x", "y", "z")
          d$Tsne1 <- s$x[d$x]
          d$Tsne2 <- s$y[d$y]
          mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
          g = ggplot(data = d, aes(x = Tsne1, y = Tsne2, fill = z,text = paste("Enhancer level:",z)))+
            geom_raster()+
            scale_fill_gradientn(colours = mycol,limits=range(z),'Enhancer level')+
            theme_classic()+
            ggtitle(paste0('Enhancer Tsne heatmap of ',input$gene_enhancer))
          ggplotly(g,tooltip = "text")
        }
      })

      output$heat_enhancer_pca =renderPlotly({
        if (!is.null(pca())){
          pca=pca()
          z =cbind(enhancer_pred_dat(),enhancer_true_dat())[pos_gene_enhancer(),]
          cell = cell_var()
          x = pca$PC1
          y =pca$PC2
          s <- interp(x,y,z)
          d <- melt(s$z, na.rm = TRUE)
          names(d) <- c("x", "y", "z")
          d$PC1 <- s$x[d$x]
          d$PC2 <- s$y[d$y]
          mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
          g = ggplot(data = d, aes(x = PC1, y = PC2, fill = z,text = paste("enhancer level:",z)))+
            geom_raster()+
            scale_fill_gradientn(colours = mycol,limits=range(z),'enhancer level')+
            theme_classic()+
            ggtitle(paste0('enhancer PCA heatmap of ',input$gene_enhancer))
          ggplotly(g,tooltip = "text")
        }
      })

    ###########Transcription factor
    gene_tf2 = reactive({
      input$gene_tf2
    })
    celltype_tf = reactive({
      input$celltype_tf
    })

    pos_gene_tf <- reactive ({
      match(input$gene_tf,rownames(true_tf))
    })

    gene_curve_tf = reactive({
      match(gene_tf2(),rownames(true_tf))
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
      dd = data_mc[,which(cell_mc%in%celltype_tf())]
      colnames(dd) = colnames(data_mc)[which(cell_mc%in%celltype_tf())]
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
      line = t(sapply(1:length(gene_tf2()),function(i) fitted(loess(time_tf()[i,] ~ c(1:ncol(time_tf()))))))
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
        plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='', name=  'WASH7P',
                marker = list(opacity=0.5,width = 2)) %>% 
          add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers'
          ) %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_tf()[which(cell_tf()%in%celltype_tf())],
                                                      '</br> Gene: ', long$Var1,
                                                      '</br> Sample: ',long$Var2))%>%
          layout(title = paste0( 'Transcription factor curve of ',celltype_tf()),
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
          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', cell_tf()[which(cell_tf()%in%celltype_tf())],
                                                      '</br> Gene: ', long$Var1,
                                                      '</br> Sample: ',long$Var2))%>%
          layout(title = paste0( 'Transcription factor expression curve of ',celltype_tf()),
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
        cell_time = cell_tf()[which(cell_tf()%in%celltype_tf())]
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
        z = cbind(tf_pred_dat(),tf_true_dat())[pos_gene_tf(),]
        cell = cell_var()
        x = tsne$Tsne1
        y = tsne$Tsne2
        s <- interp(x,y,z)
        d <- melt(s$z, na.rm = TRUE)
        names(d) <- c("x", "y", "z")
        d$Tsne1 <- s$x[d$x]
        d$Tsne2 <- s$y[d$y]
        mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
        g = ggplot(data = d, aes(x = Tsne1, y = Tsne2, fill = z,text = paste("tf level:",z)))+
          geom_raster()+
          scale_fill_gradientn(colours = mycol,limits=range(z),'tf level')+
          theme_classic()+
          ggtitle(paste0('Transcription factor Tsne heatmap of ',input$gene_tf))
        ggplotly(g,tooltip = "text")
      }
    })

    output$heat_tf_pca =renderPlotly({
      if (!is.null(pca())){
        pca=pca()
        z =cbind(tf_pred_dat(),tf_true_dat())[pos_gene_tf(),]
        cell = cell_var()
        x = pca$PC1
        y =pca$PC2
        s <- interp(x,y,z)
        d <- melt(s$z, na.rm = TRUE)
        names(d) <- c("x", "y", "z")
        d$PC1 <- s$x[d$x]
        d$PC2 <- s$y[d$y]
        mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
        g = ggplot(data = d, aes(x = PC1, y = PC2, fill = z,text = paste("tf level:",z)))+
          geom_raster()+
          scale_fill_gradientn(colours = mycol,limits=range(z),'tf level')+
          theme_classic()+
          ggtitle(paste0('Transcription factor PCA heatmap of ',input$gene_tf))
        ggplotly(g,tooltip = "text")
      }
    })
    
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
    
    dat = reactive({
      dat_sel()[,keep()]
    })
    
   
    
    #for uncorrected data
    dd_unc = reactive({
      set.seed(10)
      dat_ts = (raw[,pos_sel()])[,keep()]
      rtsne=Rtsne(prcomp(t((dat_ts)),scale=T)$x[,1:50])$Y
      rownames(rtsne)=colnames(dat_ts)
      tsne=as.data.frame(rtsne)
      colnames(tsne)=c('Tsne1','Tsne2')
      as.matrix(dist(tsne)) #all.dists2_unc#
    })
    
    sil_unc = reactive({
      score_sil <- (cluster::silhouette(as.numeric(factor(cell_var())), dd_unc()))
      score_sil[,3] 
    })
    
    dd_mnn= reactive({
      set.seed(10)
      dat_ts = (mnn[,pos_sel()])[,keep()]
      rtsne=Rtsne(prcomp(t((dat_ts)),scale=T)$x[,1:50])$Y
      rownames(rtsne)=colnames(dat_ts)
      tsne=as.data.frame(rtsne)
      colnames(tsne)=c('Tsne1','Tsne2')
       as.matrix(dist(tsne))
    })
    
    sil_mnn = reactive({
      score_sil <- (cluster::silhouette(as.numeric(factor(cell_var())), dd_mnn()))
      score_sil[,3] 
    })
    
    dd_comb= reactive({
      set.seed(10)
      dat_ts = (comb[,pos_sel()])[,keep()]
      rtsne=Rtsne(prcomp(t((dat_ts)),scale=T)$x[,1:50])$Y
      rownames(rtsne)=colnames(dat_ts)
      tsne=as.data.frame(rtsne)
      colnames(tsne)=c('Tsne1','Tsne2')
      as.matrix(dist(tsne))
    })
    sil_comb = reactive({
      score_sil <- (cluster::silhouette(as.numeric(factor(cell_var())), dd_comb()))
      score_sil[,3] 
    })

    dd_iter= reactive({
      set.seed(10)
      dat_ts = (iter[,pos_sel()])[,keep()]
      rtsne=Rtsne(prcomp(t((dat_ts)),scale=T)$x[,1:50])$Y
      rownames(rtsne)=colnames(dat_ts)
      tsne=as.data.frame(rtsne)
      colnames(tsne)=c('Tsne1','Tsne2')
      as.matrix(dist(tsne))
    })
    sil_iter = reactive({
      score_sil <- (cluster::silhouette(as.numeric(factor(cell_var())), dd_iter()))
      score_sil[,3] 
    })


    output$silhouette  <- renderPlot({
      sils<-cbind(sil_unc(),sil_comb(),sil_mnn(),sil_iter())
      boxplot(sils,main="Silhouette coefficient",names=c("Uncorrected","MNN",'iter','combat'),lwd=4,ylab="Silhouette coefficient")
    })
    
    # ################## Batch mixing entropy on PCAs
    # BatchEntropy <- function(dataset, batch0, L=100, M=100, k=500) {
    #   require(RANN) 
    #   nbatches<-length(unique(batch0))
    #   entropy<-matrix(0,L,1)
    #   set.seed(0) 
    #   for (boot in 1:L) {
    #     bootsamples<-sample(1:nrow(dataset),M)
    #     W21<-nn2(dataset,query=dataset[bootsamples,],k)
    #     for (i in 1:length(bootsamples)){
    #       
    #       for (j in 1:nbatches) {
    #         xi<-max(1,sum(batch0[W21$nn.idx[i,]]==j))
    #         entropy[boot]<-entropy[boot]+xi*log(xi)
    #       }
    #     }
    #   }
    #   return( (-1)*entropy/length(bootsamples) )
    # }
    # entropy_unc<-reactive({
    #   dat_en = (raw[,pos_sel()])[,keep()]
    #   pca_en=prcomp(t((dat_en)),scale=T)$x
    #   pca_en = as.data.frame(pca_en)
    #   BatchEntropy(pca_en[,1:2],batch_var())
    # })
    # 
    # entropy_mnn<-reactive({
    #   dat_en = (mnn[,pos_sel()])[,keep()]
    #   pca_en=prcomp(t((dat_en)),scale=T)$x
    #   pca_en = as.data.frame(pca_en)
    #   BatchEntropy(pca_en[,1:2],batch_var())
    # })
    # 
    # entropy_iter<-reactive({
    #   dat_en = (iter[,pos_sel()])[,keep()]
    #   pca_en=prcomp(t((dat_en)),scale=T)$x
    #   pca_en = as.data.frame(pca_en)
    #   BatchEntropy(pca_en[,1:2],batch_var())
    # })
    # 
    # entropy_comb<-reactive({
    #   dat_en = (comb[,pos_sel()])[,keep()]
    #   pca_en=prcomp(t((dat_en)),scale=T)$x
    #   pca_en = as.data.frame(pca_en)
    #   BatchEntropy(pca_en[,1:2],batch_var())
    # })
    #   
    # output$entropies  <- renderPlot({
    #   entropies<-cbind(entropy_unc(),entropy_comb(),entropy_mnn(),entropy_iter())
    #   boxplot(entropies,main="",names=c("Uncorrected","MNN",'iter','combat'),lwd=4,ylab="Batch mixing entropy")#,col="Yellow",ylab="Alpha dists")
    # })
    
    
    
###################
#### kbet 
    kbet_unc <- reactive({
      dat = (raw[,pos_sel()])[,keep()]
      kBET(dat, batch_var(),plot=FALSE)
    })
    
    output$k_unc  <- renderPlot({
      kbet_unc = kbet_unc()
    plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                      each=length(kbet_unc$stats$kBET.observed)), 
                            data =  c(kbet_unc$stats$kBET.observed,
                                      kbet_unc$stats$kBET.expected))
   ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
      labs(x='Test', y='Rejection rate',title='uncorrected kBET test results') +
      theme_bw() +  
      scale_y_continuous(limits=c(0,1))
    })
    kbet_mnn <- reactive({
      dat = (mnn[,pos_sel()])[,keep()]
      kBET(dat, batch_var(),plot=FALSE)
    })
    
    output$k_mnn  <- renderPlot({
      kbet_mnn = kbet_mnn()
      plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                        each=length(kbet_mnn$stats$kBET.observed)), 
                              data =  c(kbet_mnn$stats$kBET.observed,
                                        kbet_mnn$stats$kBET.expected))
      ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
        labs(x='Test', y='Rejection rate',title='MNN kBET test results') +
        theme_bw() +  
        scale_y_continuous(limits=c(0,1))
    })
    
    kbet_comb <- reactive({
      dat = (comb[,pos_sel()])[,keep()]
      kBET(dat, batch_var(),plot=FALSE)
    })
    
    output$k_comb  <- renderPlot({
      kbet_comb = kbet_comb()
      plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                        each=length(kbet_comb$stats$kBET.observed)), 
                              data =  c(kbet_comb$stats$kBET.observed,
                                        kbet_comb$stats$kBET.expected))
      ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
        labs(x='Test', y='Rejection rate',title='Combat kBET test results') +
        theme_bw() +  
        scale_y_continuous(limits=c(0,1))
    })
    
    kbet_iter <- reactive({
      dat = (iter[,pos_sel()])[,keep()]
      kBET(dat, batch_var(),plot=FALSE)
    })
    
    output$k_iter  <- renderPlot({
      kbet_iter = kbet_iter()
      plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                        each=length(kbet_iter$stats$kBET.observed)), 
                              data =  c(kbet_iter$stats$kBET.observed,
                                        kbet_iter$stats$kBET.expected))
      ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
        labs(x='Test', y='Rejection rate',title='iterated MNN kBET test results') +
        theme_bw() +  
        scale_y_continuous(limits=c(0,1))
    })
    ## neighbor 
 #    min.n <- function(x,n){ 
 #      s <- sort(x, index.return=TRUE) 
 #      s$ix[c(1:n)]
 #    }
 #    
 #    closedist_unc=reactive({
 #      dd_unc = dd_unc()
 #      diag(dd_unc)=10000
 #      min_unc = lapply(c(1:nrow(dd_unc)), function(i) FUN=min.n(dd_unc[,i],10))
 #      close_unc=data.frame(lapply(c(1:nrow(dd_unc)), function(i) ifelse(cell_k[min_unc[[i]]]==cell_var[i],10,-1)))
 #      ave_unc=colMeans(close_unc)
 #      colMeans(data.frame(lapply(c(1:nrow(dd_unc)), function(i) close_unc[,i]/dd_unc[min_unc[[i]],i])))
 #    })
 #    
 #    closedist_mnn= reactive({
 #      dd_mnn = dd_mnn()
 #      diag(dd_mnn)=10000
 #      min_mnn = lapply(c(1:nrow(dd_mnn)), function(i) FUN=min.n(dd_mnn[,i],10))
 #      close_mnn=data.frame(lapply(c(1:nrow(dd_mnn)), function(i) ifelse(cell_k[min_mnn[[i]]]==cell_var[i],10,-1)))
 #      ave_mnn=colMeans(close_mnn)
 #      colMeans(data.frame(lapply(c(1:nrow(dd_mnn)), function(i) close_mnn[,i]/dd_mnn[min_mnn[[i]],i])))
 #    })
 #    
 #    closedist_iter=reactive({
 #      dd_iter = dd_iter()
 #      diag(dd_iter)=10000
 #      min_iter = lapply(c(1:nrow(dd_iter)), function(i) FUN=min.n(dd_iter[,i],10))
 #      close_iter=data.frame(lapply(c(1:nrow(dd_iter)), function(i) ifelse(cell_k[min_iter[[i]]]==cell_var[i],10,-1)))
 #      ave_iter=colMeans(close_iter)
 #      colMeans(data.frame(lapply(c(1:nrow(dd_iter)), function(i) close_iter[,i]/dd_iter[min_iter[[i]],i])))
 #    })
 #    
 #    closedist_comb=reactive({
 #      dd_comb = dd_comb()
 #      diag(dd_comb)=10000
 #      min_comb = lapply(c(1:nrow(dd_comb)), function(i) FUN=min.n(dd_comb[,i],10))
 #      close_comb=data.frame(lapply(c(1:nrow(dd_comb)), function(i) ifelse(cell_k[min_comb[[i]]]==cell_var[i],10,-1)))
 #      ave_comb=colMeans(close_comb)
 #      colMeans(data.frame(lapply(c(1:nrow(dd_comb)), function(i) close_comb[,i]/dd_comb[min_comb[[i]],i])))
 #    })
 # 
 # output$dist = renderPlot({
 #    ave_dist<-cbind(closedist_unc(),closedist_comb(),closedist_mnn(),closedist_iter())
 #    boxplot(ave_dist,main="",names=c("Uncorrected",'Combat',"MNN",'iter'),lwd=4,ylim=c(-5,5),
 #            ylab="10 nearest neighbor distance (10,-1),k=15",
 #            xlab = paste0('Mean: unc =',round(mean(closedist_unc()),4), 
 #                          ' combat =',round(mean(closedist_comb()),4), 
 #                          ' mnn=',round(mean(closedist_mnn()),4), 
 #                          ' iter=',round(mean(closedist_iter()),4)
 #            ))
 # })
 
 #
    # output$compare  <- renderPlot({
    #   closedist_unc = ave_dist[,1]
    #   closedist_mnn = ave_dist[,2]
    #   closedist_iter = ave_dist[,3]
    #   boxplot(ave_dist,main="",names=c("Uncorrected","MNN",'iter'),lwd=4,ylim=c(-10,10),ylab="10 nearest neighbor distance (1,-1),k=15",
    #           xlab = paste0('Mean: unc =',round(mean(closedist_unc),4),
    #                         ' mnn=',round(mean(closedist_mnn),4),
    #                         ' iter=',round(mean(closedist_iter),4)
    #           ))
    # })
    #

}

shinyApp(ui = ui, server = server)
