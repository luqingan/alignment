library(shiny)
library(ggplot2)
library(reshape2)
library(matrixStats)
library(markdown)
library(plotly)
library(kBET)
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
load("raw_var.rda")
load("info_var.rda")
load("mnn_var.rda")
load("iter_var.rda")
load("raw_all.rda")
load("info_all.rda")
load("mnn_all.rda")
load("iter_all.rda")

ui <-  shinyUI(navbarPage("APP",
                          tabPanel("APP", fluidPage(
                            #headerPanel("GrowthAnalyst"),
                            sidebarLayout(
                              sidebarPanel(
                              
                                #selectInput('search','Search',multiple = T),
                                checkboxInput('part','Low variability data',value = TRUE),
                                conditionalPanel(
                                  condition = "input.part == true",
                                  uiOutput('search'),
                                  selectInput("corr_part", "Batch correction", 
                                              list("uncorrected"="uncorrected", "mnn"="mnn","iterated mnn"="iter")
                                  ),
                                  checkboxInput('tsne_var','Tsne plot',value = TRUE),
                                  checkboxInput('pca_var','PCA plot')
                                 # checkboxInput('compare','Cmpare methods')
                                  
                                  
                                ),
                             
                                tags$hr(),
                                tags$hr(),
                                
                                checkboxInput('all','All data'),
                                conditionalPanel(
                                  condition = "input.all == true",
                                  uiOutput('search2'),
                                  selectInput("corr_all", "Batch correction", 
                                              list("uncorrected"="uncorrected", "mnn"="mnn","iterated mnn"="iter")
                                  ),
                                  checkboxInput('tsne_all','Tsne plot',value = TRUE),
                                  checkboxInput('pca_all','PCA plot')
                                ) 
                                ),
                               
                              # Show a plot of the generated distribution
                              mainPanel(
                                verbatimTextOutput("impo"),
                                  conditionalPanel(
                                  condition = 'input.part==true && input.tsne_var==true',
                                  h4('Tsne plot of low variability data'),
                                  plotlyOutput('plot1'),
                                  h4('Tsne plot of low variability data (by type)'),
                                  plotlyOutput('plot2'),
                                  h4('Tsne plot of low variability data (by layer)'),
                                  plotlyOutput('layer_var'),
                                  plotOutput('kbet'),
                                  plotOutput('compare')
                                  ),

                                # conditionalPanel(
                                #   condition = 'input.part==true && input.tsne_var==true',
                                #   h4('Tsne plot of low variability data (by type)'),
                                #   plotlyOutput('plot2')),
                                
                                conditionalPanel(
                                  condition = 'input.part==true && input.pca_var==true',
                                  h4('PCA plot of low variability data'),
                                  plotlyOutput('pca1'),
                                  h4('PCA plot of low variability data (by type)'),
                                  plotlyOutput('pca2'),
                                  h4('PCA plot of low variability data (by layer)'),
                                  plotlyOutput('layer_var_pca')),
                                
                                conditionalPanel(
                                  condition = 'input.all==true && input.tsne_all==true',
                                  h4('Tsne plot of all data'),
                                  plotlyOutput('plot3'),
                                  h4('Tsne plot of all data (by type)'),
                                  plotlyOutput('plot4'),
                                  h4('Tsne plot of all data (by layer)'),
                                  plotlyOutput('layer_all')),
                                conditionalPanel(
                                  condition = 'input.all==true && input.pca_all==true',
                                  h4('PCA plot of all data'),
                                  plotlyOutput('pca3'),
                                  h4('PCA plot of all data (by type)'),
                                  plotlyOutput('pca4'),
                                  h4('PCA plot of all data (by layer)'),
                                  plotlyOutput('layer_all_pca'))
                                )
                              )
                            )
                          )
                            )
                            )
                          

server <- function(input, output) {
  shape_var=as.numeric(unique(info_var[,2]))
  shape_var=ifelse(shape_var<=18,shape_var,shape_var-18)
  pred_var=ifelse(info_var[,3]=='match','3',info_var[,4])
  link_var= gsub("@.*","",info_var[,5])

  urls=link_var
  
  dat = reactive({
    if (input$corr_part=='uncorrected'){
      raw_var
   } else if (input$corr_part=='mnn'){
       mnn_var
    }else{
      iter_var
      }
  })
  tsne=reactive({
    set.seed(10)
    rtsne=Rtsne(prcomp(t((dat())),scale=T)$x[,1:50])$Y
    rownames(rtsne)=colnames(dat())
    tsne_unc=as.data.frame(rtsne)
    tsne_unc
  })
pca = reactive({
  pca_unc=prcomp(t((dat())),scale=T)$x
  pca_unc = as.data.frame(pca_unc)
  pca_unc
})

  shape_all=as.numeric(unique(info_all[,2]))
  shape_all=ifelse(shape_all<=18,shape_all,shape_all-18)
  pred_all=ifelse(info_all[,3]=='match','3',info_all[,4])
  link_all= gsub("@.*","",info_all[,5])

  #urls=link_all
  dat_all = reactive({
    if (input$corr_all=='uncorrected'){
      raw_all
    } else if (input$corr_all=='mnn'){
      mnn_all
    }else{
      iter_all
    }
  })

  
  tsne_all=reactive({
    set.seed(10)
    rtsne=Rtsne(prcomp(t((dat_all())),scale=T)$x[,1:50])$Y
    rownames(rtsne)=colnames(dat_all())
    tsne_unc=as.data.frame(rtsne)
    tsne_unc
  })
  
  pca_all = reactive({
    pca_unc=prcomp(t((dat_all())),scale=T)$x
    pca_unc = as.data.frame(pca_unc)
    pca_unc
  })
  

  output$search <- renderUI({
    selectizeInput("search",
                   label = "Sample of Interest",
                   choices = colnames(dat()),
                   multiple = T,
                   options = list(maxItems = nrow(dat()), placeholder = 'Select a sample')
                   )
  })
  
  output$search2 <- renderUI({
    selectizeInput("search2",
                   label = "Sample of Interest",
                   choices = colnames(dat_all()),
                   multiple = T,
                   options = list(maxItems = 5, placeholder = 'Select a sample')
    )
  })
  # output$impo <- renderPrint({
  #   group = as.character(rep(1,ncol(dat())))
  #   group[pos()] = info_var[,1][pos()]
  #   group
  # })
    
  pos <- reactive ({
    match(input$search,colnames(dat()))
    })
  
  pos_all <- reactive ({
    match(input$search2,colnames(dat_all()))
  })
  

  output$plot1 <- renderPlotly({
    if (length(input$search) == 0) {
     plot_ly(tsne(), x = ~V1, y = ~V2,
              color = info_var[,1],colors='Set1',
              symbol=~info_var[,1] ,symbols =shape_var,  marker=list(size=10, opacity=0.7)
            ) %>%
        add_markers(info = ~link_var, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                                      '</br> Cell: ', info_var[,1],
                                                                      '</br> Batch: ', info_var[,3]),
                    hovermode="closest"
            ) %>%
    onRender("
             function(el, x) {
             el.on('plotly_click', function(d) {
             // d.points is an array of objects which, in this case,
             // is length 1 since the click is tied to 1 point.
             var pt = d.points[0];
             var url = pt.data.info[pt.pointNumber];
             // DISCLAIMER: this won't work from RStudio
             window.open(url);
             });
             }
             ")
      }else{
        group = rep('unselected',ncol(dat()))
        group[pos()] = info_var[,1][pos()]
        plot_ly(tsne(), x = ~V1, y = ~V2,
                color = ~group,colors='Set1',
                symbol= ~info_var[,1] ,symbols =shape_var,  marker=list(size=10, opacity=0.7)) %>%
          add_markers(info = ~link_var, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                                        '</br> Cell: ', info_var[,1],
                                                                        '</br> Batch: ', info_var[,3]),
                      hovermode="closest"
          ) %>%
          onRender("
                   function(el, x) {
                   el.on('plotly_click', function(d) {
                   // d.points is an array of objects which, in this case,
                   // is length 1 since the click is tied to 1 point.
                   var pt = d.points[0];
                   var url = pt.data.info[pt.pointNumber];
                   // DISCLAIMER: this won't work from RStudio
                   window.open(url);
                   });
                   }
                   ")

      }
  })
  
  
  output$plot2 <- renderPlotly({
    if (length(input$search) == 0) {
      plot_ly(tsne(), x = ~V1, y = ~V2,
              color = info_var[,1],colors='Set1',
              symbol=info_var[,3],symbols =unique(pred_var),  marker=list(size=10, opacity=0.7)) %>%
        add_markers(info = ~link_var, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                                      '</br> Cell: ', info_var[,1],
                                                                      '</br> Batch: ', info_var[,3]),
                    hovermode="closest"
        ) %>%
        onRender("
                 function(el, x) {
                 el.on('plotly_click', function(d) {
                 // d.points is an array of objects which, in this case,
                 // is length 1 since the click is tied to 1 point.
                 var pt = d.points[0];
                 var url = pt.data.info[pt.pointNumber];
                 // DISCLAIMER: this won't work from RStudio
                 window.open(url);
                 });
                 }
                 ")
                 }else{
                   group = rep('unselected',ncol(dat()))
                   group[pos()] = info_var[,1][pos()]
                   plot_ly(tsne(), x = ~V1, y = ~V2,
                           color = ~as.factor(group),colors='Set1',
                           symbol=info_var[,3],symbols =unique(pred_var),  marker=list(size=10, opacity=0.7)) %>%
                     add_markers(info = ~link_var, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                                                   '</br> Cell: ', info_var[,1],
                                                                                   '</br> Batch: ', info_var[,3]),
                                 hovermode="closest"
                     ) %>%
                     onRender("
                              function(el, x) {
                              el.on('plotly_click', function(d) {
                              // d.points is an array of objects which, in this case,
                              // is length 1 since the click is tied to 1 point.
                              var pt = d.points[0];
                              var url = pt.data.info[pt.pointNumber];
                              // DISCLAIMER: this won't work from RStudio
                              window.open(url);
                              });
                              }
                              ")
                   
                              }
                              })
  output$layer_var= renderPlotly({
    if (length(input$search) == 0) {
      plot_ly(tsne(), x = ~V1, y = ~V2,
              color = info_var[,1],colors='Set1',
              symbol=info_var[,6],symbols =as.factor(unique(info_var[,6])),  marker=list(size=10, opacity=0.7)) %>%
        add_markers(info = ~link_var, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                                      '</br> Cell: ', info_var[,1],
                                                                      '</br> Batch: ', info_var[,3]),
                    hovermode="closest"
        ) %>%
        onRender("
                 function(el, x) {
                 el.on('plotly_click', function(d) {
                 // d.points is an array of objects which, in this case,
                 // is length 1 since the click is tied to 1 point.
                 var pt = d.points[0];
                 var url = pt.data.info[pt.pointNumber];
                 // DISCLAIMER: this won't work from RStudio
                 window.open(url);
                 });
                 }
                 ")
                 }else{
                   group = rep('unselected',ncol(dat()))
                   group[pos()] = info_var[,1][pos()]
                   plot_ly(tsne(), x = ~V1, y = ~V2,
                           color = ~as.factor(group),colors='Set1',
                           symbol=info_var[,6],symbols =as.factor(unique(info_var[,6])),  marker=list(size=10, opacity=0.7)) %>%
                     add_markers(info = ~link_var, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                                                   '</br> Cell: ', info_var[,1],
                                                                                   '</br> Batch: ', info_var[,3]),
                                 hovermode="closest"
                     ) %>%
                     onRender("
                              function(el, x) {
                              el.on('plotly_click', function(d) {
                              // d.points is an array of objects which, in this case,
                              // is length 1 since the click is tied to 1 point.
                              var pt = d.points[0];
                              var url = pt.data.info[pt.pointNumber];
                              // DISCLAIMER: this won't work from RStudio
                              window.open(url);
                              });
                              }
                              ")
                   
                              }
                              })

  output$plot3 <- renderPlotly({
    if (length(input$search2) == 0) {
      plot_ly(tsne_all(), x = ~V1, y = ~V2,
              color = info_all[,1],colors='Set1',
              symbol=~info_all[,1] ,symbols =shape_all,  marker=list(size=10, opacity=0.7)) %>%
        add_markers(info = ~link_all, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat_all()),
                                                                      '</br> Cell: ', info_all[,1],
                                                                      '</br> Batch: ', info_all[,3]),
                    hovermode="closest"
        ) %>%
        onRender("
                 function(el, x) {
                 el.on('plotly_click', function(d) {
                 // d.points is an array of objects which, in this case,
                 // is length 1 since the click is tied to 1 point.
                 var pt = d.points[0];
                 var url = pt.data.info[pt.pointNumber];
                 // DISCLAIMER: this won't work from RStudio
                 window.open(url);
                 });
                 }
                 ")
                 }else{
                   group = rep('unselected',ncol(dat_all()))
                   group[pos_all()] = info_all[,1][pos_all()]
                   plot_ly(tsne_all(), x = ~V1, y = ~V2,
                           color = ~as.factor(group),colors='Set1',
                           symbol=~info_all[,1] ,symbols =shape_all,  marker=list(size=10, opacity=0.7)) %>%
                     add_markers(info = ~link_all, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat_all()),
                                                                                   '</br> Cell: ', info_all[,1],
                                                                                   '</br> Batch: ', info_all[,3]),
                                 hovermode="closest"
                     ) %>%
                  
                     onRender("
                              function(el, x) {
                              el.on('plotly_click', function(d) {
                              // d.points is an array of objects which, in this case,
                              // is length 1 since the click is tied to 1 point.
                              var pt = d.points[0];
                              var url = pt.data.info[pt.pointNumber];
                              // DISCLAIMER: this won't work from RStudio
                              window.open(url);
                              });
                              }
                              ")
                   
                              }
                              })
  
  
  output$plot4 <- renderPlotly({
    if (length(input$search2) == 0) {
      plot_ly(tsne_all(), x = ~V1, y = ~V2,
              color = info_all[,1],colors='Set1',
              symbol=info_all[,3],symbols =unique(pred_all),  marker=list(size=10, opacity=0.7)) %>%
        add_markers(info = ~link_all, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat_all()),
                                                                      '</br> Cell: ', info_all[,1],
                                                                      '</br> Batch: ', info_all[,3]),
                    hovermode="closest"
        ) %>%
        onRender("
                 function(el, x) {
                 el.on('plotly_click', function(d) {
                 // d.points is an array of objects which, in this case,
                 // is length 1 since the click is tied to 1 point.
                 var pt = d.points[0];
                 var url = pt.data.info[pt.pointNumber];
                 // DISCLAIMER: this won't work from RStudio
                 window.open(url);
                 });
                 }
                 ")
                 }else{
                   group = rep('unselected',ncol(dat_all()))
                   group[pos_all()] = info_all[,1][pos_all()]
                   plot_ly(tsne_all(), x = ~V1, y = ~V2,
                           color = ~as.factor(group),colors='Set1',
                           symbol=pred_all,symbols =unique(pred_all),  marker=list(size=10, opacity=0.7)) %>%
                     add_markers(info = ~link_all, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat_all()),
                                                                                   '</br> Cell: ', info_all[,1],
                                                                                   '</br> Batch: ', info_all[,3]),
                                 hovermode="closest"
                     ) %>%

                     onRender("
                              function(el, x) {
                              el.on('plotly_click', function(d) {
                              // d.points is an array of objects which, in this case,
                              // is length 1 since the click is tied to 1 point.
                              var pt = d.points[0];
                              var url = pt.data.info[pt.pointNumber];
                              // DISCLAIMER: this won't work from RStudio
                              window.open(url);
                              });
                              }
                              ")
                   
                              }
                              })
  
  output$layer_all <- renderPlotly({
    if (length(input$search2) == 0) {
      plot_ly(tsne_all(), x = ~V1, y = ~V2,
              color = info_all[,1],colors='Set1',
              symbol=info_all[,6],symbols =as.factor(unique(info_all[,6])),  marker=list(size=10, opacity=0.7)) %>%
        add_markers(info = ~link_all, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat_all()),
                                                                      '</br> Cell: ', info_all[,1],
                                                                      '</br> Batch: ', info_all[,3]),
                    hovermode="closest"
        ) %>%
        onRender("
                 function(el, x) {
                 el.on('plotly_click', function(d) {
                 // d.points is an array of objects which, in this case,
                 // is length 1 since the click is tied to 1 point.
                 var pt = d.points[0];
                 var url = pt.data.info[pt.pointNumber];
                 // DISCLAIMER: this won't work from RStudio
                 window.open(url);
                 });
                 }
                 ")
                 }else{
                   group = rep('unselected',ncol(dat_all()))
                   group[pos_all()] = info_all[,1][pos_all()]
                   plot_ly(tsne_all(), x = ~V1, y = ~V2,
                           color = ~as.factor(group),colors='Set1',
                           symbol=pred_all,symbols =unique(pred_all),  marker=list(size=10, opacity=0.7)) %>%
                     add_markers(info = ~link_all, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat_all()),
                                                                                   '</br> Cell: ', info_all[,1],
                                                                                   '</br> Batch: ', info_all[,3]),
                                 hovermode="closest"
                     ) %>%
                     
                     onRender("
                              function(el, x) {
                              el.on('plotly_click', function(d) {
                              // d.points is an array of objects which, in this case,
                              // is length 1 since the click is tied to 1 point.
                              var pt = d.points[0];
                              var url = pt.data.info[pt.pointNumber];
                              // DISCLAIMER: this won't work from RStudio
                              window.open(url);
                              });
                              }
                              ")
                   
                              }
         })
  
  ##pca
  
  
  output$pca1 <- renderPlotly({
    if (length(input$search) == 0) {
      plot_ly(pca(), x = ~PC1, y = ~PC2,
              color = info_var[,1],colors='Set1',
              symbol=~info_var[,1] ,symbols =shape_var,  marker=list(size=10, opacity=0.7)) %>%
        add_markers(info = ~link_var, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                                      '</br> Cell: ', info_var[,1],
                                                                      '</br> Batch: ', info_var[,3]),
                    hovermode="closest"
        ) %>%
        onRender("
                 function(el, x) {
                 el.on('plotly_click', function(d) {
                 // d.points is an array of objects which, in this case,
                 // is length 1 since the click is tied to 1 point.
                 var pt = d.points[0];
                 var url = pt.data.info[pt.pointNumber];
                 // DISCLAIMER: this won't work from RStudio
                 window.open(url);
                 });
                 }
                 ")
                 }else{
                   group = rep('unselected',ncol(dat()))
                   group[pos()] = info_var[,1][pos()]
                   plot_ly(pca(), x = ~PC1, y = ~PC2,
                           color = ~group,colors='Set1',
                           symbol=~info_var[,1] ,symbols =shape_var,  marker=list(size=10, opacity=0.7)) %>%
                     add_markers(info = ~link_var, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                                                   '</br> Cell: ', info_var[,1],
                                                                                   '</br> Batch: ', info_var[,3]),
                                 hovermode="closest"
                     ) %>%
              
                     onRender("
                              function(el, x) {
                              el.on('plotly_click', function(d) {
                              // d.points is an array of objects which, in this case,
                              // is length 1 since the click is tied to 1 point.
                              var pt = d.points[0];
                              var url = pt.data.info[pt.pointNumber];
                              // DISCLAIMER: this won't work from RStudio
                              window.open(url);
                              });
                              }
                              ")
                   
                              }
                              })
  
  
  output$pca2 <- renderPlotly({
    if (length(input$search) == 0) {
      plot_ly(pca(), x = ~PC1, y = ~PC2,
              color = info_var[,1],colors='Set1',
              symbol=info_var[,3],symbols =unique(pred_var),  marker=list(size=10, opacity=0.7)) %>%
        add_markers(info = ~link_var, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                                      '</br> Cell: ', info_var[,1],
                                                                      '</br> Batch: ', info_var[,3]),
                    hovermode="closest"
        ) %>%

        onRender("
                 function(el, x) {
                 el.on('plotly_click', function(d) {
                 // d.points is an array of objects which, in this case,
                 // is length 1 since the click is tied to 1 point.
                 var pt = d.points[0];
                 var url = pt.data.info[pt.pointNumber];
                 // DISCLAIMER: this won't work from RStudio
                 window.open(url);
                 });
                 }
                 ")
                 }else{
                   group = rep('unselected',ncol(dat()))
                   group[pos()] = info_var[,1][pos()]
                   plot_ly(pca(), x = ~PC1, y = ~PC2,
                           color = ~as.factor(group),colors='Set1',
                           symbol=info_var[,3],symbols =unique(pred_var),  marker=list(size=10, opacity=0.7)) %>%
                     add_markers(info = ~link_var, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                                                   '</br> Cell: ', info_var[,1],
                                                                                   '</br> Batch: ', info_var[,3]),
                                 hovermode="closest"
                     ) %>%

                     onRender("
                              function(el, x) {
                              el.on('plotly_click', function(d) {
                              // d.points is an array of objects which, in this case,
                              // is length 1 since the click is tied to 1 point.
                              var pt = d.points[0];
                              var url = pt.data.info[pt.pointNumber];
                              // DISCLAIMER: this won't work from RStudio
                              window.open(url);
                              });
                              }
                              ")
                   
                              }
                              })
  
  output$layer_var_pca <- renderPlotly({
    if (length(input$search) == 0) {
      plot_ly(pca(), x = ~PC1, y = ~PC2,
              color = info_var[,1],colors='Set1',
              symbol=info_var[,6],symbols =as.factor(unique(info_var[,6])),  marker=list(size=10, opacity=0.7)) %>%
        add_markers(info = ~link_var, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                                      '</br> Cell: ', info_var[,1],
                                                                      '</br> Batch: ', info_var[,3]),
                    hovermode="closest"
        ) %>%
        
        onRender("
                 function(el, x) {
                 el.on('plotly_click', function(d) {
                 // d.points is an array of objects which, in this case,
                 // is length 1 since the click is tied to 1 point.
                 var pt = d.points[0];
                 var url = pt.data.info[pt.pointNumber];
                 // DISCLAIMER: this won't work from RStudio
                 window.open(url);
                 });
                 }
                 ")
                 }else{
                   group = rep('unselected',ncol(dat()))
                   group[pos()] = info_var[,1][pos()]
                   plot_ly(pca(), x = ~PC1, y = ~PC2,
                           color = ~as.factor(group),colors='Set1',
                           symbol=info_var[,6],symbols =as.factor(unique(info_var[,6])),  marker=list(size=10, opacity=0.7)) %>%
                           add_markers(info = ~link_var, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat()),
                                                                                   '</br> Cell: ', info_var[,1],
                                                                                   '</br> Batch: ', info_var[,3]),
                                 hovermode="closest"
                     ) %>%
                     
                     onRender("
                              function(el, x) {
                              el.on('plotly_click', function(d) {
                              // d.points is an array of objects which, in this case,
                              // is length 1 since the click is tied to 1 point.
                              var pt = d.points[0];
                              var url = pt.data.info[pt.pointNumber];
                              // DISCLAIMER: this won't work from RStudio
                              window.open(url);
                              });
                              }
                              ")
                   
                              }
                              })
  output$pca3 <- renderPlotly({
    if (length(input$search2) == 0) {
      
      plot_ly(pca_all(), x = ~PC1, y = ~PC2,
              color = info_all[,1],colors='Set1',
              symbol=~info_all[,1] ,symbols =shape_all,  marker=list(size=10, opacity=0.7)) %>%
        add_markers(info = ~link_all, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat_all()),
                                                                      '</br> Cell: ', info_all[,1],
                                                                      '</br> Batch: ', info_all[,3]),
                    hovermode="closest"
        ) %>%

        onRender("
                 function(el, x) {
                 el.on('plotly_click', function(d) {
                 // d.points is an array of objects which, in this case,
                 // is length 1 since the click is tied to 1 point.
                 var pt = d.points[0];
                 var url = pt.data.info[pt.pointNumber];
                 // DISCLAIMER: this won't work from RStudio
                 window.open(url);
                 });
                 }
                 ")
                 }else{
                   group = rep('unselected',ncol(dat_all()))
                   group[pos_all()] = info_all[,1][pos_all()]
                   plot_ly(pca_all(), x = ~PC1, y = ~PC2,
                           color = ~as.factor(group),colors='Set1',
                           symbol=~info_all[,1] ,symbols =shape_all,  marker=list(size=10, opacity=0.7)) %>%
                     add_markers(info = ~link_all, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat_all()),
                                                                                   '</br> Cell: ', info_all[,1],
                                                                                   '</br> Batch: ', info_all[,3]),
                                 hovermode="closest"
                     ) %>%
            
                     
                     onRender("
                              function(el, x) {
                              el.on('plotly_click', function(d) {
                              // d.points is an array of objects which, in this case,
                              // is length 1 since the click is tied to 1 point.
                              var pt = d.points[0];
                              var url = pt.data.info[pt.pointNumber];
                              // DISCLAIMER: this won't work from RStudio
                              window.open(url);
                              });
                              }
                              ")
                   
                              }
                              })
  
  
  output$pca4 <- renderPlotly({
    if (length(input$search2) == 0) {
      plot_ly(pca_all(), x = ~PC1, y = ~PC2,
              color = info_all[,1],colors='Set1',
              symbol=info_all[,3],symbols =unique(pred_all),  marker=list(size=10, opacity=0.7)) %>%
        add_markers(info = ~link_all, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat_all()),
                                                                      '</br> Cell: ', info_all[,1],
                                                                      '</br> Batch: ', info_all[,3]),
                    hovermode="closest"
        ) %>%

        onRender("
                 function(el, x) {
                 el.on('plotly_click', function(d) {
                 // d.points is an array of objects which, in this case,
                 // is length 1 since the click is tied to 1 point.
                 var pt = d.points[0];
                 var url = pt.data.info[pt.pointNumber];
                 // DISCLAIMER: this won't work from RStudio
                 window.open(url);
                 });
                 }
                 ")
                 }else{
                   group = rep('unselected',ncol(dat_all()))
                   group[pos_all()] = info_all[,1][pos_all()]
                   plot_ly(pca_all(), x = ~PC1, y = ~PC2,
                           color = ~as.factor(group),colors='Set1',
                           symbol=info_all[,3],symbols =unique(pred_all),  marker=list(size=10, opacity=0.7)) %>%
                     add_markers(info = ~link_all, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat_all()),
                                                                                   '</br> Cell: ', info_all[,1],
                                                                                   '</br> Batch: ', info_all[,3]),
                                 hovermode="closest"
                     ) %>%
                   
                     onRender("
                              function(el, x) {
                              el.on('plotly_click', function(d) {
                              // d.points is an array of objects which, in this case,
                              // is length 1 since the click is tied to 1 point.
                              var pt = d.points[0];
                              var url = pt.data.info[pt.pointNumber];
                              // DISCLAIMER: this won't work from RStudio
                              window.open(url);
                              });
                              }
                              ")
                   
                 }
                 })
  
  output$layer_all_pca <- renderPlotly({
    if (length(input$search2) == 0) {
      plot_ly(pca_all(), x = ~PC1, y = ~PC2,
              color = info_all[,1],colors='Set1',
              symbol=info_all[,6],symbols =as.factor(unique(info_all[,6])),  marker=list(size=10, opacity=0.7)) %>%
        add_markers(info = ~link_all, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat_all()),
                                                                      '</br> Cell: ', info_all[,1],
                                                                      '</br> Batch: ', info_all[,3]),
                    hovermode="closest"
        ) %>%
        
        onRender("
                 function(el, x) {
                 el.on('plotly_click', function(d) {
                 // d.points is an array of objects which, in this case,
                 // is length 1 since the click is tied to 1 point.
                 var pt = d.points[0];
                 var url = pt.data.info[pt.pointNumber];
                 // DISCLAIMER: this won't work from RStudio
                 window.open(url);
                 });
                 }
                 ")
                 }else{
                   group = rep('unselected',ncol(dat_all()))
                   group[pos_all()] = info_all[,1][pos_all()]
                   plot_ly(pca_all(), x = ~PC1, y = ~PC2,
                           color = ~as.factor(group),colors='Set1',
                           symbol=info_all[,6],symbols =as.factor(unique(info_all[,6])),  marker=list(size=10, opacity=0.7)) %>%
                     
                     add_markers(info = ~link_all, hoverinfo="text" ,text = ~paste('</br> Sample: ',colnames(dat_all()),
                                                                                   '</br> Cell: ', info_all[,1],
                                                                                   '</br> Batch: ', info_all[,3]),
                                 hovermode="closest"
                     ) %>%
                     
                     onRender("
                              function(el, x) {
                              el.on('plotly_click', function(d) {
                              // d.points is an array of objects which, in this case,
                              // is length 1 since the click is tied to 1 point.
                              var pt = d.points[0];
                              var url = pt.data.info[pt.pointNumber];
                              // DISCLAIMER: this won't work from RStudio
                              window.open(url);
                              });
                              }
                              ")
                   
                              }
                              })
  
  
  output$kbet  <- renderPlot({
    kbet<- kBET(dat(), info_var[,4])
    plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                      each=length(kbet$stats$kBET.observed)), 
                            data =  c(kbet$stats$kBET.observed,
                                      kbet$stats$kBET.expected))
    g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
      labs(x='Test', y='Rejection rate',title='kBET test results') +
      theme_bw() +  
      scale_y_continuous(limits=c(0,1))
    g
  })
  
  
  output$compare  <- renderPlot({
    set.seed(10)
    rtsne=Rtsne(prcomp(t((raw_var)),scale=T)$x[,1:50])$Y  
    rownames(rtsne)=colnames(raw_var)
    tsne_unc=as.data.frame(rtsne)
    
    set.seed(10)
    rtsne=Rtsne(prcomp(t((mnn_var)),scale=T)$x[,1:50])$Y  
    rownames(rtsne)=colnames(mnn_var)
    tsne_mnn=as.data.frame(rtsne)
    
    set.seed(10)
    rtsne=Rtsne(prcomp(t((iter_var)),scale=T)$x[,1:50])$Y  
    rownames(rtsne)=colnames(iter_var)
    tsne_iter=as.data.frame(rtsne)
    
    dd_unc <- as.matrix(dist(tsne_unc)) #all.dists2_unc#
    dd_mnn <- as.matrix(dist(tsne_mnn))#all.dists2.c#
    dd_iter <- as.matrix(dist(tsne_iter))#all.dists2.c#
    diag(dd_unc)=10000
    min_unc=apply(dd_unc, 1, FUN=which.min)
    n=c(1:length(min_unc))
    close=unlist(lapply(n, function(i) ifelse(info_var[,1][min_unc[i]]==info_var[,1][i],1,0)))
    
    
    diag(dd_mnn)=10000
    min_mnn=apply(dd_mnn, 1, FUN=which.min)
    
    close=unlist(lapply(n, function(i) ifelse(info_var[,1][min_mnn[i]]==info_var[,1][i],1,0)))
    mean(close)
    
    diag(dd_iter)=10000
    min_iter=apply(dd_iter, 1, FUN=which.min)
    close=unlist(lapply(n, function(i) ifelse(info_var[,1][min_iter[i]]==info_var[,1][i],1,0)))
    mean(close)
    
    ## neighbor 
    min.n <- function(x,n){ 
      s <- sort(x, index.return=TRUE) 
      s$ix[c(1:n)]
    }
    min_unc=lapply(c(1:nrow(dd_unc)), function(i) FUN=min.n(dd_unc[,i],10))
    close_unc=data.frame(lapply(c(1:nrow(dd_unc)), function(i) ifelse(info_var[,1][min_unc[[i]]]==info_var[,1][i],1,-1)))
    ave_unc=colMeans(close_unc)
    closedist_unc = colMeans(data.frame(lapply(c(1:nrow(dd_unc)), function(i) close_unc[,i]/dd_unc[min_unc[[i]],i])))
    
    min_mnn=lapply(c(1:nrow(dd_mnn)), function(i) FUN=min.n(dd_mnn[,i],10))
    close_mnn=data.frame(lapply(c(1:nrow(dd_mnn)), function(i) ifelse(info_var[,1][min_mnn[[i]]]==info_var[,1][i],1,-1)))
    ave_mnn=colMeans(close_mnn)
    closedist_mnn = colMeans(data.frame(lapply(c(1:nrow(dd_mnn)), function(i) close_mnn[,i]/dd_mnn[min_mnn[[i]],i])))
    
    min_iter=lapply(c(1:nrow(dd_iter)), function(i) FUN=min.n(dd_iter[,i],10))
    close_iter=data.frame(lapply(c(1:nrow(dd_iter)), function(i) ifelse(info_var[,1][min_iter[[i]]]==info_var[,1][i],1,-1)))
    ave_iter=colMeans(close_iter)
    closedist_iter = colMeans(data.frame(lapply(c(1:nrow(dd_iter)), function(i) close_iter[,i]/dd_iter[min_iter[[i]],i])))
    

    ave_dist<-cbind(closedist_unc,closedist_mnn,closedist_iter)
    
    boxplot(ave_dist,main="",names=c("Uncorrected","MNN",'iter'),lwd=4,,ylim=c(-10,10),ylab="10 nearest neighbor distance (1,-1),k=15",
            xlab = paste0('Mean: unc =',round(mean(closedist_unc),4), 
                          ' mnn=',round(mean(closedist_mnn),4), 
                          ' iter=',round(mean(closedist_iter),4)
            ))
  })
    
    
}

shinyApp(ui = ui, server = server)
