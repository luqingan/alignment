library(shiny)
library(ggplot2)
library(reshape2)
library(matrixStats)
library(markdown)
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
# source("https://bioconductor.org/biocLite.R")
# biocLite("preprocessCore")
library(preprocessCore)
library(htmlwidgets)
load("raw.rda")
load("info.rda")
load("mnn.rda")
load("iter.rda")
load('ave_dist.rda')
rna <- readRDS("RNA_data_norm_hg19_all.rds")

ui <-  shinyUI(navbarPage("APP",
                          tabPanel("APP", fluidPage(
                            #headerPanel("GrowthAnalyst"),
                            sidebarLayout(
                              sidebarPanel(
                                #selectInput('search','Search',multiple = T),
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
                                              list("uncorrected"="uncorrected", "mnn"="mnn","iterated mnn"="iter")
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
                                 # checkboxInput('compare','Compare methods')
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
                                 uiOutput('celltype'),
                                 uiOutput('gene2'))
                              ),
                              # Show a plot of the generated distribution
                              mainPanel(
                                  tabsetPanel(
                                    tabPanel("PCA", 
                                             h4('PCA plot'),
                                             plotlyOutput('pca1')
                                    ), 
                                    tabPanel("Tsne", 
                                             verbatimTextOutput("impo"),
                                             h4('Tsne plot'),
                                             plotlyOutput('plot1')
                                             ), 
                                    
                                    tabPanel("Gene expression", 
                                             h4('Gene expression heatmap'),
                                             plotlyOutput('heat'),
                                             plotlyOutput('heat_pca')
                                             ),
                                    tabPanel("Gene expression curve", 
                                             h4('Gene expression curve'),
                                             verbatimTextOutput('curve_text'),
                                             plotlyOutput('curve'),
                                             plotlyOutput('pca_curve')
                                             
                                    ),
                                    tabPanel("Correction methods evaluation", 
                                             plotOutput('kbet'),
                                             plotOutput('compare'))
                                )# end of tabset panel 
                              )# end of main panel 
                             
                                )
                              )
                            ))
)
                          
                            
                          

server <- function(input, output) {
  ##data 
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
                   options = list(maxItems = nrow(rna),placeholder = 'Select cell types')
    )
  })
  
  output$celltype <- renderUI({
    selectizeInput("celltype",
                   label = "Cell types of Interest",
                   choices = cell_var(),
                   multiple = T,
                   options = list(placeholder = 'Select a sample')
    )
  })
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

  
  dat_sel = reactive({
    if (input$corr=='uncorrected'){
     d= normalize.quantiles(as.matrix(raw[,pos_sel()]),copy=TRUE)
     colnames(d) = colnames(raw)
     d
    }else if (input$corr=='mnn'){
      d=normalize.quantiles(as.matrix(mnn[,pos_sel()]),copy=TRUE)
      colnames(d) = colnames(raw)
      d
    }else{
      normalize.quantiles(as.matrix(iter[,pos_sel()]),copy=TRUE)
      colnames(d) = colnames(raw)
      d
    }
  })
  
  cell_sel = reactive({
    cell[pos_sel()]
  })
  
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



 tsne = reactive({
       set.seed(10)
       rtsne=Rtsne(prcomp(t((dat())),scale=T)$x[,1:50])$Y
       rownames(rtsne)=colnames(dat())
       tsne=as.data.frame(rtsne)
       colnames(tsne)=c('Tsne1','Tsne2')
       tsne
 })

dat_cv = reactive({
 raw = dat()
 cv = apply(raw,1,sd)/rowMeans(raw)
 raw[which(cv>0),]
})

pca = reactive({
  pca_unc=prcomp(t((dat_cv())),scale=T)$x
  pca_unc = as.data.frame(pca_unc)
  pca_unc
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
  cell_rna = reactive({
    cell_var()[which(batch_var()=='Predicted DNase')]
  })

  pca_curve = reactive({
    pca()[which(cell_rna()%in%input$celltype),]
  })
  rna_curve = reactive({
    rna[,which(cell_rna()%in%input$celltype)]
  })
  
  dd = reactive({
    if (length(input$gene2)==1){
      dd = data.frame(cbind(pca_curve()$PC1,rna_curve()[gene_curve(),]))
    }else{
    dd = data.frame(cbind(pca_curve()$PC1,t(rna_curve()[gene_curve(),])))
    #cell_d = cell[which(cell_rna()%in%input$celltype)]
    colnames(dd)[1]= 'PC1'
    }
    dd
  })

  long = reactive({
    long = melt(dd(),id='PC1', value.name = "gene")
    long
  })
  line = reactive({
    dd=dd()
    fitted = sapply(1:length(input$gene2),function(i) fitted(loess(dd[,i+1] ~ dd$PC1)))
    line = data.frame(cbind(dd$PC1,fitted))
    colnames(line) = colnames(dd)
    line = melt(line,id='PC1', value.name = "fit")
    line
  })

  # output$impo <- renderPrint({
  #   cate_color()[which(cell_rna()%in%input$celltype)]
  #   })


  # output$curve_text = renderPrint({
  #   if (length(input$gene2==0)&length(input$celltype==0)){
  #     cat('Please select gene and cell type of interest')
  #   }else if (length(input$celltype==0)&length(input$gene2!=0)){
  #     cat('Please select cell type of interest')
  #   }else if (length(input$celltype!=0)&length(input$gene2==0)){
  #     cat('Please select gene of interest')
  #   }else{
  #     cat(paste0('Plot gene ',input$gene2, 'expression curve of ', input$celltype, 'samples'))
  #   }
  # })
  
  output$curve =renderPlotly({
    if (!is.null(long())&!is.null(line())&!is.null(input$celltype)&!is.null(input$gene2)){
    long=long()
    line=line()
    plot_ly(long, x = ~PC1, y = ~gene, type = 'bar',
            color = ~variable,
            marker = list(opacity=0.5,width = 2)) %>%
      add_lines(x =~line$PC1,y = ~ line$fit, type = 'line',
                color = ~line$variable,
                showlegend = F)  %>%
      layout(title = paste0( 'Gene expression curve of',input$celltype),
             xaxis = list(
               title = "PC1",
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

  output$pca_curve <- renderPlotly({
    if (!is.null(pca_curve())){
      
      plot_ly(pca_curve(), x = ~PC1, y = ~PC2,
              color = cell_var()[which(cell_rna()%in%input$celltype)],colors='Set1',
              symbol= ~cell_var()[which(cell_rna()%in%input$celltype)] ,symbols =shape_var()[which(cell_rna()%in%input$celltype)],  marker=list(size=10, opacity=0.7)) %>%
        add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat())[which(cell_rna()%in%input$celltype)],
                                                    '</br> Cell: ', cell_var()[which(cell_rna()%in%input$celltype)],
                                                    '</br> Batch: ', batch_var()[which(cell_rna()%in%input$celltype)]),
                    hovermode="closest"
        )
    }
  })
  
  output$plot1 <- renderPlotly({
    if (!is.null(tsne())){
      
    if (length(input$search) == 0) {
      plot_ly(tsne(), x = ~Tsne1, y = ~Tsne2,
            color = cate_color(),colors='Set1',
            symbol=~cate_shape(),symbols =shape_var(),  marker=list(size=10, opacity=0.7)) %>%
            add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                        '</br> Cell: ', cell_var(),
                                                        '</br> Batch: ', batch_var()),
                        hovermode="closest"
            )
                 }else{
                   group = rep('unselected',ncol(dat()))
                   group[pos_sample()] = cate_color()[pos_sample()]
                   plot_ly(tsne(), x = ~Tsne1, y = ~Tsne2,
                           color = ~group,colors='Set1',
                           symbol= ~cate_shape() ,symbols =shape_var(),  marker=list(size=10, opacity=0.7)) %>%
                     add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                                               '</br> Cell: ', cell_var(),
                                                                               '</br> Batch: ', batch_var()),
                                 hovermode="closest"
                     )

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
                                                                  '</br> Batch: ', batch_var()),
                    hovermode="closest"
        )

           }else{
                   group = rep('unselected',ncol(dat()))
                   group[pos_sample()] = cate_color()[pos_sample()]
                   plot_ly(pca(), x = ~PC1, y = ~PC2,
                           color = ~group,colors='Set1',
                           symbol= ~cate_shape() ,symbols =shape_var(),  marker=list(size=10, opacity=0.7)) %>%
                     add_markers(hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                                               '</br> Cell: ', cell_var(),
                                                                               '</br> Batch: ', batch_var()),
                                 hovermode="closest"
                     )
           }
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
