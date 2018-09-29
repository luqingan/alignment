##all
options(repos = BiocInstaller::biocinstallRepos())
getOption("repos")
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

load("rnaseq.rda") ## clustered, 2000 row
load("info.rda")
load("optimized.rda") # k=25,sigma = 0.297540831779859
# load("mnn.rda") # k=25,sigma = 0.297540831779859
# load("X_mnn.rda") # k=25,sigma = 0.297540831779859
# load("X_mnn2.rda") # k=25,sigma = 0.297540831779859
# load("X_mnn3.rda") # k=25,sigma = 0.297540831779859

load('true_promoter.rda')
load('pred_promoter.rda')
load('true_enhancer.rda')
load('pred_enhancer.rda')
load('true_tf_qn.rda')
load('pred_tf_qn.rda')

rna <- readRDS("RNA_data_norm_hg19_all_name.rds") # 

ui <-  shinyUI(navbarPage("APP",
                          tabPanel("APP", fluidPage(
                            #headerPanel("GrowthAnalyst"),
                            sidebarLayout(
                              sidebarPanel(
                                h3('Data processing'),
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
                                # h3('Method selection'),
                                # checkboxInput("method", 'Correction methods'),
                                # conditionalPanel(
                                #   condition = "input.method== true",
                                #   selectInput("corr", "Batch correction method", 
                                #               list( 'Optimized MNN' = 'opt',"uncorrected"="uncorrected", "combat"="comb","mnn"="mnn","iterated mnn"="iter")
                                #   )),
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
                                    uiOutput('gene'),
                                    uiOutput('gene_pro'),
                                    uiOutput('gene_enhancer'),
                                    uiOutput('gene_tf')
                                ), 
                                
                                h3("Curves"),
                                checkboxInput("curve", 'Curves'),
                                conditionalPanel(
                                  condition = "input.curve == true",
                                  radioButtons(inputId="pseudo", "Construct pseudo time with:",
                                               choices = list("RNA-seq" = "rna_pseudo",
                                                              'Predicted DNase-seq' = 'pred_pseudo',
                                                              'True DNase-seq'= 'true_pseudo',
                                                              'Predicted DNase-seq + True DNase-seq'='comb_pseudo'),
                                               inline=TRUE),
                                  uiOutput('celltype'),
                                  
                                  h5('Gene expression curve'),
                                  uiOutput('gene2'),
                                  
                                  h5('Promoter curve'),
                                  
                                  uiOutput('gene_pro2'),

                                  h5('Transcription factor curve'),
                              
                                  uiOutput('gene_tf2'),
                                  h5('Enhancer curve'),
                                  uiOutput('gene_enhancer2')
  
                                
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
                                           plotlyOutput('promoter_curve_pred')%>% withSpinner(),
                                           h2(),
                                           plotlyOutput('promoter_curve_true')%>% withSpinner(),
                                           h2(),
                                           plotlyOutput('promoter_curve_comb')%>% withSpinner(),
                                           
                                           
                                           h4('Transcription factor curve'),
                                           plotlyOutput('tf_curve_pred')%>% withSpinner(),
                                           h2(),
                                           plotlyOutput('tf_curve_true')%>% withSpinner(),
                                           h2(),
                                           plotlyOutput('tf_curve_comb')%>% withSpinner(),
                                           
                                           h4('Enhancer curve'),
                                           plotlyOutput('enhancer_curve_pred')%>% withSpinner(),
                                           h2(),
                                           plotlyOutput('enhancer_curve_true')%>% withSpinner(),
                                           h2(),
                                           plotlyOutput('enhancer_curve_comb')%>% withSpinner()
                                       
                                  )
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
                     choices = colnames(optimized),
                     multiple = T,
                     options = list(maxItems = nrow(optimized), placeholder = 'Select a sample')
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


      #######Transcription factor
      output$gene_tf <- renderUI({
        selectizeInput("gene_tf",
                       label = "Transcription factor of Interest",
                       choices = rownames(x_q),
                       multiple = F,
                       #selected = input$gene,
                       options = list(placeholder = 'Select transcription factors')
        )
      })

      output$gene_tf2 <- renderUI({
        selectizeInput("gene_tf2",
                       label = "Transcription factor of Interest",
                       choices = rownames(true_tf_qn),
                       multiple = T,
                       #selected = input$gene2,
                       options = list(maxItems = nrow(true_tf_qn),placeholder = 'Select transcription factors')
        )
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
        optimized[,pos_sel()]
      })
    
    
    # dat_sel = reactive({
    #   if (input$corr=='opt'){
    #     optimized[,pos_sel()]
    #   }else if (input$corr=='uncorrected'){
    #     X_mnn[,pos_sel()]
    #   }else if (input$corr=='comb'){
    #     X_mnn2[,pos_sel()]
    #   }else if (input$corr=='mnn'){
    #     mnn[,pos_sel()]
    #   }else{
    #     X_mnn3[,pos_sel()]
    #   }
    # })
    
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
        d = true_promoter[,na.omit(match(pos(),which(batch=='True DNase')))]
        colnames(d) = colnames(true_promoter)[na.omit(match(pos(),which(batch=='True DNase')))]
        d
      })

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

    ##########tf
    tf_pred_dat = reactive({
      d =  pred_tf_qn[,na.omit(match(pos(),which(batch=='Predicted DNase')))]
      colnames(d) = colnames(pred_tf_qn)[na.omit(match(pos(),which(batch=='Predicted DNase')))]
      d
    })

    tf_true_dat = reactive({
      d =true_tf_qn[,na.omit(match(pos(),which(batch=='True DNase')))]
      colnames(d) = colnames(true_tf_qn)[na.omit(match(pos(),which(batch=='True DNase')))]
      d
    })

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

    rna_k = reactive({
      rna_dat()[,which(cell_predicted()%in%input$celltype)]
    })
    pred_k =  reactive({
      dat_predicted()[,which(cell_predicted()%in%input$celltype)]
    })
    true_k =  reactive({
      dat_true()[,which(cell_true()%in%input$celltype)]
    })
    
    mc = reactive({
      if (input$pseudo == 'rna_pseudo'){
        data_mc = rna_k()
        cell_mc = cell_predicted()
      }else if (input$pseudo == 'true_pseudo'){
        data_mc = true_k()
        cell_mc = cell_true()
      }else if (input$pseudo == 'pred_pseudo'){
        data_mc = pred_k()
        cell_mc = cell_predicted()
      }else if (input$pseudo == 'comb_pseudo'){
        ### replicate sample, use true data 
        data_mc = cbind(pred_k(),true_k())
        cell_mc = cell_var()
      }
      set.seed(10)
      exprmclust(data_mc)
    })
    
   
    order <- reactive({
      TSCANorder(mc(),listbranch = F)
    })

    
    order_rna = reactive({
       if (input$pseudo == 'rna_pseudo'){
         order()
      }else if (input$pseudo == 'true_pseudo'){
        if (is.na(match(colnames(true_k()),order()))){
          true_k = true_k()[,-which(is.na(match(colnames(true_k()),order())))]
        }else{
          true_k=true_k()
        }
        pred_k = pred_k()
        corr = cor(pred_k,true_k)
        pos_max =  apply(corr, 1, which.max) 
        insert = rep(0,ncol(rna_k()))
        for (i in c(1:length(pos_max))){
          if (1<pos_max[i] & pos_max[i]<length(order())){
            cor_in = corr[i,pos_max[i]]
            cor_prev = corr[i,pos_max[i]-1]
            cor_next = corr[i,pos_max[i]+1]
            if (cor_prev>cor_next){
              insert[i] = pos_max[i]-1 + cor_in
            }else{
              insert[i] = pos_max[i]+1- cor_in
            }
          }else if (pos_max[i]==1){
            insert[i] = pos_max[i]+1- cor_in
          }else{
            insert[i] = pos_max[i]-1 + cor_in
          }
        }
        names(insert) = names(pos_max)
        order = sort(insert)
        names(order)
      }else if (input$pseudo == 'pred_pseudo'){
        order()
      }else if (input$pseudo == 'comb_pseudo'){
        order()[na.omit(match(colnames(rna_k()),order()))]
      }
    })
    
    
    order_pred = reactive({
      order_rna()
    })
    
    order_true = reactive({
      if (input$pseudo == 'true_pseudo'){
        order()
      }else if (input$pseudo == 'comb_pseudo'){
        order()[na.omit(match(colnames(true_k()),order()))]
      }else{
        true_k = true_k()
        if (is.na(match(colnames(pred_k()),order()))){
          pred_k = pred_k()[,-which(is.na(match(colnames(pred_k()),order())))]
        }else{
          pred_k=pred_k()
        }
        corr = cor(true_k,pred_k)
        pos_max =  apply(corr, 1, which.max) 
        insert = rep(0,ncol(true_k()))
        for (i in c(1:length(pos_max))){
          if (1<pos_max[i] & pos_max[i]<length(order())){
            cor_in = corr[i,pos_max[i]]
            cor_prev = corr[i,pos_max[i]-1]
            cor_next = corr[i,pos_max[i]+1]
            if (cor_prev>cor_next){
              insert[i] = pos_max[i]-1 + cor_in
            }else{
              insert[i] = pos_max[i]+1- cor_in
            }
          }else if (pos_max[i]==1){
            insert[i] = pos_max[i]+1- cor_in
          }else{
            insert[i] = pos_max[i]-1 + cor_in
          }
        }
        names(insert) = names(pos_max)
        order_true = sort(insert)
        names(order_true)
      }
    })
    
    order_comb = reactive({
      order = order()
      if (input$pseudo == 'true_pseudo'){
        if (is.na(match(colnames(true_k()),order()))){
          true_k = true_k()[,-which(is.na(match(colnames(true_k()),order())))]
        }else{
          true_k=true_k()
        }
        pred_k = pred_k()
        corr = cor(pred_k,true_k)
        pos_max =  apply(corr, 1, which.max) 
        insert = rep(0,ncol(pred_k()))
        for (i in c(1:length(pos_max))){
          if (1<pos_max[i] & pos_max[i]<length(order())){
            cor_in = corr[i,pos_max[i]]
            cor_prev = corr[i,pos_max[i]-1]
            cor_next = corr[i,pos_max[i]+1]
            if (cor_prev>cor_next){
              insert[i] = pos_max[i]-1 + cor_in
            }else{
              insert[i] = pos_max[i]+1- cor_in
            }
          }else if (pos_max[i]==1){
            insert[i] = pos_max[i]+1- cor_in
          }else{
            insert[i] = pos_max[i]-1 + cor_in
          }
        }
        names(insert) = names(pos_max)
        base = c(1:length(order()))
        names(base) = order()
        order_comb = sort(c(insert,base))
        names(order_comb)
    }else if (input$pseudo == 'comb_pseudo'){
      order()
    }else{
      true_k = true_k()
      if (is.na(match(colnames(pred_k()),order()))){
        pred_k = pred_k()[,-which(is.na(match(colnames(pred_k()),order())))]
      }else{
        pred_k=pred_k()
      }

        corr = cor(true_k,pred_k)
        pos_max =  apply(corr, 1, which.max) 
        insert = rep(0,ncol(true_k()))
        
        for (i in 1:length(pos_max)){
          if (1<pos_max[i] & pos_max[i]<length(order())){
            cor_in = corr[i,pos_max[i]]
            cor_prev = corr[i,pos_max[i]-1]
            cor_next = corr[i,pos_max[i]+1]
            if (cor_prev>cor_next){
              insert[i] = pos_max[i]-1 + cor_in
            }else{
              insert[i] = pos_max[i]+1- cor_in
            }
          }else if (pos_max[i]==1){
            insert[i] = pos_max[i]+1- cor_in
          }else{
            insert[i] = pos_max[i]-1 + cor_in
          }
        }
        names(insert) = names(pos_max)
        base = c(1:length(order()))
        names(base) = order()
        order_comb = sort(c(insert,base))
        names(order_comb)
      }
    })
    
    # one curve for gene expression use order_rna
 
    time = reactive({
      if (length(input$gene2) == 1){
        type = rna[gene_curve(),which(colnames(rna_dat())%in%order_rna())]
        order = na.omit(match(order_rna(),names(type)))
        type[order]
      }else{
        type = rna[gene_curve(),which(colnames(rna_dat())%in%order_rna())]
        order = na.omit(match(order_rna(),colnames(type)))
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

          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_predicted()[match(order(), colnames(optimized))],each=length(input$gene2)),
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

      }else{
        long=long()
        line=line()

        plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                color = ~Var1,
                marker = list(opacity=0.5,width = 2)) %>%
          add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',showlegend = F,
                    color = ~long$Var1)  %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_predicted()[match(order_rna(), colnames(optimized))],each=length(input$gene2)),
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
        if (input$pseudo == 'rna_pseudo'){
          cell_time = cell_predicted()[which(cell_predicted()%in%input$celltype)]
        }else if (input$pseudo == 'true_pseudo'){
          cell_time = cell_true()[which(cell_true()%in%input$celltype)]
        }else if (input$pseudo == 'pred_pseudo'){
          cell_time = cell_predicted()[which(cell_predicted()%in%input$celltype)]
        }else if (input$pseudo == 'comb_pseudo'){
          cell_time = cell_var()[which(cell_var()%in%input$celltype)]
        }
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


    pro_comb_dat = reactive({
        cbind(pro_pred_dat(),pro_true_dat())
    })
    
      time_pro_pred = reactive({
       if (length(input$gene_pro2) == 1){
          type = pro_pred_dat()[gene_curve_pro(),which(colnames(pro_pred_dat())%in%order_pred())]
          order = na.omit(match(order_pred(),names(type)))
          type[order]
        }else{
          type = pro_pred_dat()[gene_curve_pro(),which(colnames(pro_pred_dat())%in%order_pred())]
          order = na.omit(match(order_pred(),colnames(type)))
          type[,order]
        }
      })
      time_pro_true = reactive({
        if (length(input$gene_pro2) == 1){
          type = pro_true_dat()[gene_curve_pro(),which(colnames(pro_true_dat())%in%order_true())]
          order = na.omit(match(order_true(),names(type)))
          type[order]
        }else{
          type = pro_true_dat()[gene_curve_pro(),which(colnames(pro_true_dat())%in%order_true())]
          order = na.omit(match(order_true(),colnames(type)))
          type[,order]
        }
      })
 
      time_pro_comb = reactive({
        if (length(input$gene_pro2) == 1){
          type = pro_comb_dat()[gene_curve_pro(),which(colnames(pro_comb_dat())%in%order_comb())]
          order = na.omit(match(order_comb(),names(type)))
          type[order]
        }else{
          type = pro_comb_dat()[gene_curve_pro(),which(colnames(pro_comb_dat())%in%order_comb())]
          order = na.omit(match(order_comb(),colnames(type)))
          type[,order]
        }
      })
      
      long_pro_pred = reactive({
          long = melt(time_pro_pred(), value.name = "gene")
          long
        })
      
      long_pro_true = reactive({
        long = melt(time_pro_true(), value.name = "gene")
        long
      })
      long_pro_comb = reactive({
        long = melt(time_pro_comb(), value.name = "gene")
        long
      })
      

      
      line_pro_pred = reactive({
        line = t(sapply(1:length(input$gene_pro2),function(i) fitted(loess(time_pro_pred()[i,] ~ c(1:ncol(time_pro_pred()))))))
        melt(line)
      })
      line_pro_true = reactive({
        line = t(sapply(1:length(input$gene_pro2),function(i) fitted(loess(time_pro_true()[i,] ~ c(1:ncol(time_pro_true()))))))
        melt(line)
      })
      
      line_pro_comb = reactive({
        line = t(sapply(1:length(input$gene_pro2),function(i) fitted(loess(time_pro_comb()[i,] ~ c(1:ncol(time_pro_comb()))))))
        melt(line)
      })
      
      output$promoter_curve_pred =renderPlotly({
        if (length(input$gene_pro2) == 1){
          line1 = fitted(loess(time_pro_pred() ~ c(1:length(time_pro_pred()))))
          line = data.frame(cbind(Var1 = 1,Var2 = c(1:length(line1)) , value = line1))
          long= data.frame(cbind(Var1 = input$gene_pro2,Var2 = names(time_pro_pred()) , gene = time_pro_pred()))
          rownames(long)=NULL
          long[,2]= line[,2]


          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                  marker = list(opacity=0.5,width = 2)) %>%
            add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers',showlegend = F) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_predicted()[match(order_pred(), colnames(optimized))],each=length(input$gene_pro2)),
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = paste0( 'Promoter expression curve of ',input$celltype),
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
          long=long_pro_pred()
          line=line_pro_pred()

          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                  color = ~Var1,
                  marker = list(opacity=0.5,width = 2)) %>%
            add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',showlegend = F,
                      color = ~long$Var1)  %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_predicted()[match(order_pred(), colnames(optimized))],each=length(input$gene_pro2)),
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = 'Promoter Expression Curve of Predicted DNase',
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
      
      # output$impo = renderPrint({
      # })
      # 
      # output$curve_text = renderPrint({
      # })
      # 
      # output$exon_curve_text = renderPrint({
      # })
      # 
      # output$intron_curve_text = renderPrint({
      # })

      output$promoter_curve_true =renderPlotly({
        if (length(input$gene_pro2) == 1){
          line1 = fitted(loess(time_pro_true() ~ c(1:length(time_pro_true()))))
          line = data.frame(cbind(Var1 = 1,Var2 = c(1:length(line1)) , value = line1))
          long= data.frame(cbind(Var1 = input$gene_pro2,Var2 = names(time_pro_true()) , gene = time_pro_true()))
          rownames(long)=NULL
          long[,2]= line[,2]
                plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                  marker = list(opacity=0.5,width = 2)) %>%
            add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers',showlegend = F) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_true()[match(order_true(), colnames(optimized)[445:852])],each=length(input$gene2)),
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title =  'Promoter Expression Curve of True DNase-seq ',
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
          long=long_pro_true()
          line=line_pro_true()
          
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                  color = ~Var1,
                  marker = list(opacity=0.5,width = 2)) %>%
            add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',showlegend = F,
                      color = ~long$Var1)  %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_true()[match(order_true(), colnames(optimized)[445:852])],each=length(input$gene2)),
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = 'Promoter expression curve of True DNase-seq',
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
      
      output$promoter_curve_comb =renderPlotly({
        if (length(input$gene_pro2) == 1){
          line1 = fitted(loess(time_pro_comb() ~ c(1:length(time_pro_comb()))))
          line = data.frame(cbind(Var1 = 1,Var2 = c(1:length(line1)) , value = line1))
          long= data.frame(cbind(Var1 = input$gene_pro2,Var2 = names(time_pro_comb()) , gene = time_pro_comb()))
          rownames(long)=NULL
          long[,2]= line[,2]
          
          
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                  marker = list(opacity=0.5,width = 2)) %>%
            add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers',showlegend = F) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_var()[match(order_comb(), colnames(optimized))],each=length(input$gene_pro2)),
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title =  'Promoter expression curve of aligned predicted and true DNase-seq',
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
          long=long_pro_comb()
          line=line_pro_comb()
          
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                  color = ~Var1,
                  marker = list(opacity=0.5,width = 2)) %>%
            add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',showlegend = F,
                      color = ~long$Var1)  %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_var()[match(order_comb(), colnames(optimized))],each=length(input$gene_pro2)),
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = 'Promoter Expression Curve of aligned predicted and true DNase-seq',
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


      ###########enhancer

      pos_gene_enhancer <- reactive ({
        match(input$gene_enhancer,rownames(true_enhancer))
      })

      gene_curve_enhancer = reactive({
        match(input$gene_enhancer2,rownames(true_enhancer))
      })
      
      
      enhancer_comb_dat = reactive({
        cbind(enhancer_pred_dat(),enhancer_true_dat())
      })
      
      time_enhancer_pred = reactive({
        if (length(input$gene_enhancer2) == 1){
          type = enhancer_pred_dat()[gene_curve_enhancer(),which(colnames(enhancer_pred_dat())%in%order_pred())]
          order = na.omit(match(order_pred(),names(type)))
          type[order]
        }else{
          type = enhancer_pred_dat()[gene_curve_enhancer(),which(colnames(enhancer_pred_dat())%in%order_pred())]
          order = na.omit(match(order_pred(),colnames(type)))
          type[,order]
        }
      })
      
      time_enhancer_true = reactive({
        if (length(input$gene_enhancer2) == 1){
          type = enhancer_true_dat()[gene_curve_enhancer(),which(colnames(enhancer_true_dat())%in%order_true())]
          order = na.omit(match(order_true(),names(type)))
          type[order]
        }else{
          type = enhancer_true_dat()[gene_curve_enhancer(),which(colnames(enhancer_true_dat())%in%order_true())]
          order = na.omit(match(order_true(),colnames(type)))
          type[,order]
        }
      })
      
      time_enhancer_comb = reactive({
        if (length(input$gene_enhancer2) == 1){
          type = enhancer_comb_dat()[gene_curve_enhancer(),which(colnames(enhancer_comb_dat())%in%order_comb())]
          order = na.omit(match(order_comb(),names(type)))
          type[order]
        }else{
          type = enhancer_comb_dat()[gene_curve_enhancer(),which(colnames(enhancer_comb_dat())%in%order_comb())]
          order = na.omit(match(order_comb(),colnames(type)))
          type[,order]
        }
      })
      
      long_enhancer_pred = reactive({
        long = melt(time_enhancer_pred(), value.name = "gene")
        long
      })
      long_enhancer_true = reactive({
        long = melt(time_enhancer_true(), value.name = "gene")
        long
      })
      long_enhancer_comb = reactive({
        long = melt(time_enhancer_comb(), value.name = "gene")
        long
      })
      
      
      
      line_enhancer_pred = reactive({
        line = t(sapply(1:length(input$gene_enhancer2),function(i) fitted(loess(time_enhancer_pred()[i,] ~ c(1:ncol(time_enhancer_pred()))))))
        melt(line)
      })
      line_enhancer_true = reactive({
        line = t(sapply(1:length(input$gene_enhancer2),function(i) fitted(loess(time_enhancer_true()[i,] ~ c(1:ncol(time_enhancer_true()))))))
        melt(line)
      })
      
      line_enhancer_comb = reactive({
        line = t(sapply(1:length(input$gene_enhancer2),function(i) fitted(loess(time_enhancer_comb()[i,] ~ c(1:ncol(time_enhancer_comb()))))))
        melt(line)
      })
      
      output$enhancer_curve_pred =renderPlotly({
        if (length(input$gene_enhancer2) == 1){
          line1 = fitted(loess(time_enhancer_pred() ~ c(1:length(time_enhancer_pred()))))
          line = data.frame(cbind(Var1 = 1,Var2 = c(1:length(line1)) , value = line1))
          long= data.frame(cbind(Var1 = input$gene_enhancer2,Var2 = names(time_enhancer_pred()) , gene = time_enhancer_pred()))
          rownames(long)=NULL
          long[,2]= line[,2]
          
          
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                  marker = list(opacity=0.5,width = 2)) %>%
            add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers',showlegend = F) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_predicted()[match(order_pred(), colnames(optimized))],each=length(input$gene_enhancer2)),
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = 'Enhancer expression curve of Predicted DNase',
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
          long=long_enhancer_pred()
          line=line_enhancer_pred()
          
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                  color = ~Var1,
                  marker = list(opacity=0.5,width = 2)) %>%
            add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',showlegend = F,
                      color = ~long$Var1)  %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_predicted()[match(order_pred(), colnames(optimized))],each=length(input$gene_enhancer2)),
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = 'Enhancer Expression Curve of Predicted DNase',
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
        }
      })
      
      output$enhancer_curve_true =renderPlotly({
        if (length(input$gene_enhancer2) == 1){
          line1 = fitted(loess(time_enhancer_true() ~ c(1:length(time_enhancer_true()))))
          line = data.frame(cbind(Var1 = 1,Var2 = c(1:length(line1)) , value = line1))
          long= data.frame(cbind(Var1 = input$gene_enhancer2,Var2 = names(time_enhancer_true()) , gene = time_enhancer_true()))
          rownames(long)=NULL
          long[,2]= line[,2]
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                  marker = list(opacity=0.5,width = 2)) %>%
            add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers',showlegend = F) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_true()[match(order_true(), colnames(optimized)[445:852])],each=length(input$gene_enhancer2)),
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title =  'Enhancer Expression Curve of True DNase-seq ',
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
          long=long_enhancer_true()
          line=line_enhancer_true()
          
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                  color = ~Var1,
                  marker = list(opacity=0.5,width = 2)) %>%
            add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',showlegend = F,
                      color = ~long$Var1)  %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_true()[match(order_true(), colnames(optimized)[445:852])],each=length(input$gene_enhancer2)),
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = 'Enhancer expression curve of True DNase-seq',
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
        }
      })
      
      output$enhancer_curve_comb =renderPlotly({
        if (length(input$gene_enhancer2) == 1){
          line1 = fitted(loess(time_enhancer_comb() ~ c(1:length(time_enhancer_comb()))))
          line = data.frame(cbind(Var1 = 1,Var2 = c(1:length(line1)) , value = line1))
          long= data.frame(cbind(Var1 = input$gene_enhancer2,Var2 = names(time_enhancer_comb()) , gene = time_enhancer_comb()))
          rownames(long)=NULL
          long[,2]= line[,2]  
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                  marker = list(opacity=0.5,width = 2)) %>%
            add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers',showlegend = F) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_var()[match(order_comb(), colnames(optimized))],each=length(input$gene_enhancer2)),
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = 'Enhancer expression curve of aligned predicted and true DNase-seq',
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
          long=long_enhancer_comb()
          line=line_enhancer_comb()
          
          plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                  color = ~Var1,
                  marker = list(opacity=0.5,width = 2)) %>%
            add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',showlegend = F,
                      color = ~long$Var1)  %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_var()[match(order_comb(), colnames(optimized))],each=length(input$gene_enhancer2)),
                                                        '</br> Gene: ', long$Var1,
                                                        '</br> Sample: ',long$Var2))%>%
            layout(title = 'Enhancer expression curve of aligned predicted and true DNase-seq',
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
    # 
    ###########Transcription factor
    pos_gene_tf <- reactive ({
      match(input$gene_tf,rownames(true_tf_qn))
    })

    gene_curve_tf = reactive({
      match(input$gene_tf2,rownames(true_tf_qn))
    })

    
    
    tf_comb_dat = reactive({
      cbind(tf_pred_dat(),tf_true_dat())
    })
    
    time_tf_pred = reactive({
      if (length(input$gene_tf2) == 1){
        type = tf_pred_dat()[gene_curve_tf(),which(colnames(tf_pred_dat())%in%order_pred())]
        order = na.omit(match(order_pred(),names(type)))
        type[order]
      }else{
        type = tf_pred_dat()[gene_curve_tf(),which(colnames(tf_pred_dat())%in%order_pred())]
        order = na.omit(match(order_pred(),colnames(type)))
        type[,order]
      }
    })
    
    time_tf_true = reactive({
      if (length(input$gene_tf2) == 1){
        type = tf_true_dat()[gene_curve_tf(),which(colnames(tf_true_dat())%in%order_true())]
        order = na.omit(match(order_true(),names(type)))
        type[order]
      }else{
        type = tf_true_dat()[gene_curve_tf(),which(colnames(tf_true_dat())%in%order_true())]
        order = na.omit(match(order_true(),colnames(type)))
        type[,order]
      }
    })
    
    time_tf_comb = reactive({
      if (length(input$gene_tf2) == 1){
        type = tf_comb_dat()[gene_curve_tf(),which(colnames(tf_comb_dat())%in%order_comb())]
        order = na.omit(match(order_comb(),names(type)))
        type[order]
      }else{
        type = tf_comb_dat()[gene_curve_tf(),which(colnames(tf_comb_dat())%in%order_comb())]
        order = na.omit(match(order_comb(),colnames(type)))
        type[,order]
      }
    })
    
    long_tf_pred = reactive({
      long = melt(time_tf_pred(), value.name = "gene")
      long
    })
    long_tf_true = reactive({
      long = melt(time_tf_true(), value.name = "gene")
      long
    })
    long_tf_comb = reactive({
      long = melt(time_tf_comb(), value.name = "gene")
      long
    })
    
    
    
    line_tf_pred = reactive({
      line = t(sapply(1:length(input$gene_tf2),function(i) fitted(loess(time_tf_pred()[i,] ~ c(1:ncol(time_tf_pred()))))))
      melt(line)
    })
    line_tf_true = reactive({
      line = t(sapply(1:length(input$gene_tf2),function(i) fitted(loess(time_tf_true()[i,] ~ c(1:ncol(time_tf_true()))))))
      melt(line)
    })
    
    line_tf_comb = reactive({
      line = t(sapply(1:length(input$gene_tf2),function(i) fitted(loess(time_tf_comb()[i,] ~ c(1:ncol(time_tf_comb()))))))
      melt(line)
    })
    
    output$tf_curve_pred =renderPlotly({
      if (length(input$gene_tf2) == 1){
        line1 = fitted(loess(time_tf_pred() ~ c(1:length(time_tf_pred()))))
        line = data.frame(cbind(Var1 = 1,Var2 = c(1:length(line1)) , value = line1))
        long= data.frame(cbind(Var1 = input$gene_tf2,Var2 = names(time_tf_pred()) , gene = time_tf_pred()))
        rownames(long)=NULL
        long[,2]= line[,2]
        
        
        plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                marker = list(opacity=0.5,width = 2)) %>%
          add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers',showlegend = F) %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_predicted()[match(order_pred(), colnames(optimized))],each=length(input$gene_tf2)),
                                                      '</br> Gene: ', long$Var1,
                                                      '</br> Sample: ',long$Var2))%>%
          layout(title = 'Transcription Factor expression curve of Predicted DNase-seq',
                 xaxis = list(
                   title = "Pseudotime",
                   showticklabels = FALSE,
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')),
                 yaxis = list(
                   title = 'Transcription Factor expression',
                   titlefont = list(
                     size = 16,
                     color = 'rgb(107, 107, 107)'),
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')))
        
      }else{
        long=long_tf_pred()
        line=line_tf_pred()
        
        plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                color = ~Var1,
                marker = list(opacity=0.5,width = 2)) %>%
          add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',showlegend = F,
                    color = ~long$Var1)  %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_predicted()[match(order_pred(), colnames(optimized))],each=length(input$gene_tf2)),
                                                      '</br> Gene: ', long$Var1,
                                                      '</br> Sample: ',long$Var2))%>%
          layout(title = 'Transcription Factor Expression Curve of Predicted DNase-seq',
                 xaxis = list(
                   title = "Pseudotime",
                   showticklabels = FALSE,
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')),
                 yaxis = list(
                   title = 'Transcription Factor expression',
                   titlefont = list(
                     size = 16,
                     color = 'rgb(107, 107, 107)'),
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')))
      }
    })
    
    output$tf_curve_true =renderPlotly({
      if (length(input$gene_tf2) == 1){
        line1 = fitted(loess(time_tf_true() ~ c(1:length(time_tf_true()))))
        line = data.frame(cbind(Var1 = 1,Var2 = c(1:length(line1)) , value = line1))
        long= data.frame(cbind(Var1 = input$gene_tf2,Var2 = names(time_tf_true()) , gene = time_tf_true()))
        rownames(long)=NULL
        long[,2]= line[,2]
        plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                marker = list(opacity=0.5,width = 2)) %>%
          add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers',showlegend = F) %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_true()[match(order_true(), colnames(optimized)[445:852])],each=length(input$gene_tf2)),
                                                      '</br> Gene: ', long$Var1,
                                                      '</br> Sample: ',long$Var2))%>%
          layout(title =  'Transcription Factor Expression Curve of True DNase-seq ',
                 xaxis = list(
                   title = "Pseudotime",
                   showticklabels = FALSE,
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')),
                 yaxis = list(
                   title = 'Transcription Factor expression',
                   titlefont = list(
                     size = 16,
                     color = 'rgb(107, 107, 107)'),
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')))
        
      }else{
        long=long_tf_true()
        line=line_tf_true()
        
        plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                color = ~Var1,
                marker = list(opacity=0.5,width = 2)) %>%
          add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',showlegend = F,
                    color = ~long$Var1)  %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_true()[match(order_true(), colnames(optimized)[445:852])],each=length(input$gene_tf2)),
                                                      '</br> Gene: ', long$Var1,
                                                      '</br> Sample: ',long$Var2))%>%
          layout(title = 'Transcription Factor expression curve of True DNase-seq',
                 xaxis = list(
                   title = "Pseudotime",
                   showticklabels = FALSE,
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')),
                 yaxis = list(
                   title = 'Transcription Factor expression',
                   titlefont = list(
                     size = 16,
                     color = 'rgb(107, 107, 107)'),
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')))
      }
    })
    
    output$tf_curve_comb =renderPlotly({
      if (length(input$gene_tf2) == 1){
        line1 = fitted(loess(time_tf_comb() ~ c(1:length(time_tf_comb()))))
        line = data.frame(cbind(Var1 = 1,Var2 = c(1:length(line1)) , value = line1))
        long= data.frame(cbind(Var1 = input$gene_tf2,Var2 = names(time_tf_comb()) , gene = time_tf_comb()))
        rownames(long)=NULL
        long[,2]= line[,2]  
        plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                marker = list(opacity=0.5,width = 2)) %>%
          add_trace(y = ~ line$value, type = 'scatter', mode = 'lines+markers',showlegend = F) %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_var()[match(order_comb(), colnames(optimized))],each=length(input$gene_tf2)),
                                                      '</br> Gene: ', long$Var1,
                                                      '</br> Sample: ',long$Var2))%>%
          layout(title = 'Transcription Factor expression curve of aligned predicted and true DNase-seq',
                 xaxis = list(
                   title = "Pseudotime",
                   showticklabels = FALSE,
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')),
                 yaxis = list(
                   title = 'Transcription Factor expression',
                   titlefont = list(
                     size = 16,
                     color = 'rgb(107, 107, 107)'),
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')))
        
      }else{
        long=long_tf_comb()
        line=line_tf_comb()
        
        plot_ly(long, x = ~Var2, y = ~gene, type = 'scatter',xlab='',
                color = ~Var1,
                marker = list(opacity=0.5,width = 2)) %>%
          add_trace(x =~long$Var2,y = ~ line$value, type = 'scatter', mode = 'lines',showlegend = F,
                    color = ~long$Var1)  %>%
          add_markers(hoverinfo="text" ,text = ~paste('</br> Cell: ', rep(cell_var()[match(order_comb(), colnames(optimized))],each=length(input$gene_tf2)),
                                                      '</br> Gene: ', long$Var1,
                                                      '</br> Sample: ',long$Var2))%>%
          layout(title = 'Transcription Factor expression curve of aligned predicted and true DNase-seq',
                 xaxis = list(
                   title = "Pseudotime",
                   showticklabels = FALSE,
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')),
                 yaxis = list(
                   title = 'Transcription Factor expression',
                   titlefont = list(
                     size = 16,
                     color = 'rgb(107, 107, 107)'),
                   tickfont = list(
                     size = 14,
                     color = 'rgb(107, 107, 107)')))
      }
    })
    
    
    # 
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


    
    
}

shinyApp(ui = ui, server = server)
