library(shiny)
library(htmltools)
library(plotly)
library(DT)
library(shinythemes)
library(readr)
library(readxl)
library(ggpubr)
library(tidyr)
library(dplyr)
library(Seurat)
library(svgPanZoom)
library(showtext)
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
showtext_auto()
font_add("Arial", "www/Arial.ttf")
colors <- c(
  "#726BAE", "#C7AED5", "#F7DBF0", "#EA945A", "#60A897", "#F5BC6E", "#86A667", 
  "#276D9F", "#6488B9", "#ACD48A", "#A5AA99", "#E5948E", "#E05F48", "#5D88BF", 
  "#78C4D4", "#8AB6D6", "#E7C5DB", "#FCEDDC", "#FCE8E3", "#CCD7DD", "#D5E7AC",
  "#67C9F2", "#F4D4D1", "#F9E7F0", "#F6EDF3", "#F3F8E5"  
)
ui <- fluidPage(
  tags$head(
    tags$title("scTExplorer"),
    tags$link(
      rel = "icon",
      type = "image/png",
      href = "logo.png"
    ),
    tags$style(HTML("
      body {
        background-color: white;
        color: black;
        overflow-x: auto !important; 
      }

      /*-----------NAV BAR STYLE-----------*/
      .nav-bar {
        background: linear-gradient(90deg, #1e1f29, #3a3d5c, #6d6fa3);
        display: flex;
        padding: 10px 0;
      }

      .navbar-left {
        font-size: 26px;
        font-weight: bold;
        color: #f3eaff;
        margin-left: 20px;
      }

      .nav-bar button {
        background-color: transparent;
        color: #f8f5ff;
        border: none;
        font-size: 18px;
        padding: 10px;
      }
      .nav-bar button:hover {
        color: #D0A9FF;
        text-shadow: 0px 0px 6px #ffbfff;
      }
      .container-fluid {
          padding-left: 0 !important;
          padding-right: 0 !important;
      }

      /*-----------HOME HEADER-----------*/
      .home-header {
        width: 100%;
        text-align: center;
        padding: 35px 0 25px 0;
        background: linear-gradient(90deg, #1e1f29, #3a3d5c, #6d6fa3);
        border-bottom: 1px solid #ffffff22;
      }

      .header-title {
        font-size: 40px;
        font-weight: bold;
        background: linear-gradient(90deg, #e3b3ff, #cfa9ff, #9d8cff);
        -webkit-background-clip: text;
        color: transparent;
      }

      .header-subtitle {
        font-size: 18px;
        margin-top: 5px;
        color: #e8e8e8;
        opacity: 0.85;
      }

      /*-----------PAGE CONTENT-----------*/
      .content {
        text-align: center;
        margin-top: 50px;
        color: black;
      }
      
      .main-content-container {
          min-width: 1300px; 
          margin: 0 auto;
          display: block;
      }
      .plot-card-small {
          display: flex;
          flex-direction: column;
          align-items: center;
          justify-content: center;
      }
      .plot-card-bar {
          height: 380px;
      }
      /*-----------CARD CONTENT-----------*/
      .plot-card-small {  /*卡片风格设计*/
          background: #ffffff !important;
          border: 1px solid #eef0f2;
          border-radius: 16px !important;
          box-shadow: 0 4px 20px rgba(0, 0, 0, 0.05) !important;
          transition: all 0.3s ease-in-out;
          position: relative;
          display: block;
          overflow: hidden !important;
      }
      .plot-card-small:hover {  /*悬停效果*/
          transform: translateY(-5px);
          box-shadow: 0 12px 30px rgba(109, 111, 163, 0.15) !important;
          border-color: #6d6fa3;
      }
      .plot-card-small h4 { /*标题*/
          font-size: 2rem;
          font-weight: 700 !important;
          color: #344767 !important;
          margin-bottom: 20px;
          padding-left: 12px;
          border-left: 4px solid #6d6fa3;
      }
      .svg-container {  /*图片*/
          background: #ffffff;
          overflow-x: auto !important;
          overflow-y: hidden;
          width: 100%;           
          border-radius: 8px;
          display: flex !important;
          justify-content: center;
          align-items: center;
      }
      .svg-container1 {  /*图片*/
          background: #ffffff;
          overflow-x: auto !important;
          overflow-y: hidden;
          width: 100%;           
          border-radius: 8px;
          display: block !important;
          align-items: center;
      }
      .svg-container img {
          margin: 0 !important;
          max-width: none !important; 
          height: auto;
          vertical-align: middle;
          display: inline-block !important;
          align-items: center !important;
      }
      .card-header-custom { /*标题和按钮位置设计*/
          display: flex;
          justify-content: space-between;
          align-items: center;
          margin-bottom: 20px;
          gap: 15px;
      }
      .card-header-custom h4 {
          margin-bottom: 0 !important; 
      }
      .btn-download-custom { /*下载按钮*/
          background-color: transparent;
          color: #6d6fa3;
          border: 1px solid #6d6fa3;
          border-radius: 15px;
          padding: 4px 12px;
          font-size: 12px;
          transition: all 0.2s;
      }
      .btn-download-custom:hover { /*下载按钮悬停效果*/
          background-color: #6d6fa3;
          color: white;
          padding-left: 5px;
          box-shadow: 0 2px 8px rgba(109, 111, 163, 0.4);
      }
      
    "))
  ),
  
  #NAVIGATION BAR------------------------
  div(class = "nav-bar",
      style = "justify-content: space-between; align-items: center;width:1900px;",
      div(class = "navbar-left", "scTExplorer"),
      div(class = "navbar-right", style = "margin-right:20px;",
          actionButton("home", "Home"),
          actionButton("oscc_progression", "OSCC Progression"),
          actionButton("oral_diseases", "Oral Diseases"),
          actionButton("pancancer", "Pan-Cancer")
          #actionButton("download", "Download"),
          #actionButton("document", "Document")
      )
  ),
  
  #HOME HEADER AREA----------------------
  div(
    class = "home-header",style = "width:1900px;",
    uiOutput("mainHeader")
  ),
  
  #HOME INTRO CONTENT--------------------
  div(
    uiOutput("HomeContent"),
    style = "width:1900px;"
  ),
  
  #DYNAMIC PAGE CONTENT--------------------
  div(
    uiOutput("pageContent"),
    style = "width:1900px;"
  )
  
)
  #Server-------------------
server <- function(input, output, session) {
  
  activePage <- reactiveVal("home")
  
  # Main Header
  output$mainHeader <- renderUI({
    if (activePage() == "home") {
      tagList(
        div(class = "header-title", "scTExplorer"),
        div(class = "header-subtitle", "A comprehensive database for T cell subtypes and states")
      )
    }
  })
  
  #Home page description
  output$HomeContent <- renderUI({
    if (activePage() == "home") {
      div(
        style = "padding:30px; font-size: 12pt; line-height:1.6; color:black;",
        h3(style = "text-align: center; color: #2c3e50; border-bottom: 2px solid #ff007f; padding-bottom: 10px; font-weight: bold;",
           "Introduction"
        ),
        HTML("scTExplorer is a comprehensive resource focusing on T cell single-cell transcriptomic data, covering OSCC progression, oral diseases, and pan-cancer datasets."),
        br(), br(),
        HTML("With scTExplorer, you can:"),
        tags$ul(
          tags$li("Investigate T cell dynamics and functional changes during OSCC progression."),
          tags$li("Compare T cell features across oral diseases."),
          tags$li("Explore T cell states across pan-cancer datasets.")
        ),
        br(),
        #三大模块-------------------
        div(
          style = "
              width:80%;
              margin:auto;
              margin-top:30px;
              display:flex;
              flex-direction:column;
              gap:30px;
            ",
          #第一行：单个OSCC卡片--------------
          div(
            style = "display:flex; justify-content:center; width:100%;",
            div(
              style = "
                width:100%;
                padding:20px;
                border-radius:18px;
                background: linear-gradient(135deg, #ffffff, #f6f3ff);
                box-shadow:0 6px 20px rgba(0,0,0,0.15);
                transition: all 0.35s ease-in-out;
                cursor:pointer;
              ",
              onmouseover = "this.style.transform='scale(1.03)'; this.style.boxShadow='0 12px 30px rgba(0,0,0,0.25)';",
              onmouseout  = "this.style.transform='scale(1.0)';  this.style.boxShadow='0 6px 20px rgba(0,0,0,0.15)';",
              actionLink(
                "OSCC_link",
                div(
                  div(
                    style='font-size:24px; font-weight:bold; color:#2c3e50; margin-bottom:15px;',
                    "OSCC Progression"
                  ),
                  #固定高度的图片容器：高度更低，不再铺满
                  div(
                    style = "
                    overflow:hidden;
                    border-radius:12px;
                    height:230px;
                    width:100%;
                    display:flex;
                    align-items:center;
                    justify-content:center;
                    background:#ffffff;
                  ",
                    tags$img(
                      src = "OSCC progression .png",
                      style="
                    max-height:100%;
                    max-width:100%;
                    object-fit:contain;
                  "
                    )
                  )
                )
              ),
              p(
                style='color:#333; margin-top:12px; font-size:14px;',
                strong("OSCC Progression: "),
                "Explore T cell dynamics and functional changes from normal (NL), leukoplakia (LP), carcinoma (CA), to lymph node metastasis (LN)."
              )
            )
          ),
          #第二行：Oral+Pan-Cancer两个卡片-----------------
          div(
            style = "display:flex; justify-content:space-between; width:100%;",
            #Oral Diseases--------
            div(
              style = "
                width:48%;
                padding:20px;
                border-radius:18px;
                background: linear-gradient(135deg, #ffffff, #f2f6ff);
                box-shadow:0 5px 18px rgba(0,0,0,0.12);
                transition: all 0.35s ease-in-out;
                cursor:pointer;
              ",
              onmouseover = "this.style.transform='scale(1.03)'; this.style.boxShadow='0 12px 30px rgba(0,0,0,0.25)';",
              onmouseout  = "this.style.transform='scale(1.0)';  this.style.boxShadow='0 5px 18px rgba(0,0,0,0.12)';",
              actionLink(
                "Oral_link",
                div(
                  div(
                    style='font-size:22px; font-weight:bold; color:#2c3e50; margin-bottom:12px;',
                    "Oral Diseases"
                  ),
                  #和OSCC一样高度的图片容器
                  div(
                    style = "
                    overflow:hidden;
                    border-radius:12px;
                    height:230px;
                    width:100%;
                    display:flex;
                    align-items:center;
                    justify-content:center;
                    background:#ffffff;
                     ",
                    tags$img(
                      src = "Oral diseases.png",
                      style="
                    max-height:100%;
                    max-width:100%;
                    object-fit:contain;
                      "
                    ))
                )
              ),
              p(
                style='color:#333; margin-top:12px; font-size:14px;',
                strong("Oral Diseases: "),
                "Compare T cell features across various oral diseases, including periodontitis, pulpitis, and oral lichen planus."
              )
            ),
            #Pan-Cancer--------
            div(
              style = "
              width:48%;
              padding:20px;
              border-radius:18px;
              background: linear-gradient(135deg, #ffffff, #f2fbff);
              box-shadow:0 5px 18px rgba(0,0,0,0.12);
              transition: all 0.35s ease-in-out;
              cursor:pointer;
            ",
              onmouseover = "this.style.transform='scale(1.03)'; this.style.boxShadow='0 12px 30px rgba(0,0,0,0.25)';",
              onmouseout  = "this.style.transform='scale(1.0)';  this.style.boxShadow='0 5px 18px rgba(0,0,0,0.12)';",
              
              actionLink(
                "Pancancer_link",
                div(
                  div(
                    style='font-size:22px; font-weight:bold; color:#2c3e50; margin-bottom:12px;',
                    "Pan-Cancer"
                  ),
                  div(
                    style = "
                  overflow:hidden;
                  border-radius:12px;
                  height:230px;
                  display:flex;
                  align-items:center;
                  justify-content:center;
                  background:#ffffff;
                    ",
                    tags$img(
                      src = "pancancer.png",  
                      style = "
                      max-height:100%;
                      max-width:100%;
                      object-fit:contain;
                    "
                    )
                    #imageOutput("Pancancer_figure", height = "220px", width = "auto")
                  )
                )
              ),
              p(
                style='color:#333; margin-top:12px; font-size:14px;',
                strong("Pan-Cancer: "),
                "Analyze T cell subset distribution and functional states across multiple cancer types."
              )
            )
          )
        ),
        br(),
        h3(style = "text-align: center; color: #2c3e50; border-bottom: 2px solid #ff007f; padding-bottom: 10px; font-weight: bold;",
           "Contact"
        ),
        #em(strong("Contact:")),
        p("Taiwen Li", a("litaiwen@scu.edu.cn", href = "mailto:litaiwen@scu.edu.cn"))
      )
    }
  })
  #Dynamic navigation pages
  observeEvent(input$home, {
    activePage("home")
    output$pageContent <- renderUI({
      div(class = "content",
          h4("Welcome to scTExplorer"),
          p("Explore the database for T cell subtypes and their states.")
      )
    })
  })
  ##OSCC_progression UI---------------------
  oscc_progression <- function(){
    output$pageContent <- renderUI({
      tagList(
        div(class = "main-content-container",style = "width:1900px;",
            div(
              style="width:100%; text-align:center; padding:25px 0; background:linear-gradient(90deg,#1e1f29,#3a3d5c,#6d6fa3); color:white; border-radius:2px; margin-bottom:25px;",
                h2("OSCC Progression", style="font-weight:bold; margin:0;"),
                p("Explore T cell dynamic changes across OSCC progression stages", style="opacity:0.8; font-size:15px; margin-top:5px;")
            ),
            div(style="display:flex; gap:25px; align-items: flex-start;",
                div(
                  style="width:300px;height:320px; flex-shrink:0; background:white; border-radius:16px; margin-left:20px;padding:20px; box-shadow:0 4px 15px rgba(0,0,0,0.1); color:black;",
                    h4("Data Selection", style="color:#2c3e50; font-weight:bold;"),
                    radioButtons("OSCC_dataset", tags$b("Choose Dataset:"), choices=c("In_house","GSE181919","HRA001006"), selected="In_house"),
                    radioButtons("OSCC_Celltype", tags$b("Choose T cell type:"), choices=c("CD4T","CD8T"), selected="CD4T"),
                    actionButton("confirm_button","Confirm", class="btn btn-primary", style="width:100%; background:#6d6fa3; border:none;")
                ),
                div(style="flex-grow:1; min-width:700px;",
                    conditionalPanel(
                      condition = "output.isDataLoaded_tab1 == false",
                      div(style="height:320px;background:white; border-radius:16px; padding:100px; text-align:center; box-shadow:0 4px 15px rgba(0,0,0,0.1); color:#999;",
                          icon("hand-pointer", class="fa-3x"),
                          h3("Please select parameters and click 'Confirm' to load data."),
                          p("The analysis results will appear here.")
                      )
                    ),
                    conditionalPanel(
                      condition = "output.isDataLoaded_tab1 == true",
                      div(style="background:white; border-radius:16px; padding:15px; box-shadow:0 4px 15px rgba(0,0,0,0.1);flex-shrink: 0 !important;width:100%;",
                          tabsetPanel(
                            id = "oscc_tabs",
                            #Overview Panel------
                            tabPanel(
                              title = tags$span(icon("chart-area"), "Overview"),
                              #第一行：UMAP+Pie Chart
                              div(style="display:flex; flex-shrink: 0 !important;gap:25px; margin-top:20px; justify-content: center;",
                                  div(class="plot-card-small", style="width:620px; padding:25px;flex-shrink: 0 !important;",
                                      div(
                                        class="card-header-custom",
                                        h4("Cluster UMAP Projection"),
                                        downloadButton("download_UMAP_OSCC", "Download PDF", class="btn-download-custom")
                                        ),
                                      div(class="svg-container", uiOutput("UMAP_plot"))
                                  ),
                                  div(class="plot-card-small", style="width:620px; padding:25px;flex-shrink: 0 !important;",
                                      div(
                                        class="card-header-custom",
                                        h4("Cell Composition Ratio"),
                                        downloadButton("download_Pie_chart_OSCC", "Download PDF", class="btn-download-custom")
                                        ),
                                      div(class="svg-container", uiOutput("Pie_chart"))
                                  )
                              ),
                              #第二行：Dotplot
                              div(style="display:flex; justify-content: center; margin-top:25px;flex-shrink: 0;",
                                  div(class="plot-card-small", style="width:1265px; padding:25px;",
                                      div(
                                        class="card-header-custom",
                                        h4("Marker Gene Expression Dotplot"),
                                        downloadButton("download_Dot_plot_OSCC", "Download PDF", class="btn-download-custom")
                                        ),
                                      div(class="svg-container", uiOutput("Dot_plot"))
                                  )
                              ),
                              #第三行：Conditions+Patient Barplots
                              div(
                                style="display:flex; gap:25px; margin-top:25px; justify-content: center;",
                                  div(
                                    class="plot-card-small", style="width:430px; padding:25px;",
                                      div(
                                        class="card-header-custom",
                                        h4("Distribution by Condition"),
                                        downloadButton("download_Barplot_conditions_OSCC", "Download PDF", class="btn-download-custom")
                                        ),
                                      div(class="svg-container1", uiOutput("Barplot_conditions"))
                                  ),
                                  div(
                                    class="plot-card-small", style="width:810px; padding:25px;",
                                      div(
                                        class="card-header-custom",
                                        h4("Individual Patient Profiles"),
                                        downloadButton("download_Barplot_patients_OSCC", "Download PDF", class="btn-download-custom")
                                        ),
                                      div(class="svg-container1", uiOutput("Barplot_patients"))
                                  )
                              ),
                              div(style="height: 50px;")
                            ),
                            ##Gene Panel----------
                            tabPanel(
                              title = tags$span(icon("dna"), "Gene"),
                              div(
                                style = "background: #ffffff; border-radius: 16px; padding: 20px; margin: 20px 0; box-shadow: 0 4px 15px rgba(0,0,0,0.05);",
                                div(
                                  style = "display: flex; gap: 20px; align-items: flex-end; justify-content: center;",
                                  div(style = "width: 400px;font-size:15pt;",
                                      selectizeInput("gene_input_OSCC", tags$b("Select Gene:"), choices = NULL, width = "100%")
                                  ),
                                  div(style = "width: 150px;",
                                      actionButton(
                                        "Submit_Gene_OSCC", "Submit",
                                        icon = icon("paper-plane"),
                                        style = "background: #6d6fa3; color: white; width: 100%; border: none; height: 38px; margin-bottom: 15px;"
                                      )
                                  )
                                )
                              ),
                              div(
                                style = "display: flex; gap: 25px; align-items: flex-start; justify-content: center;",
                                #FeaturePlot
                                div(
                                  class = "plot-card-small",
                                  style = "width: 600px;height: 550px; padding: 25px; flex-shrink: 0;",
                                  div(class = "card-header-custom",
                                      h4("Expression Distribution (Feature Plot)"),
                                      uiOutput("download_btn_Gene_FeaturePlot_OSCC")
                                  ),
                                  div(
                                    style = "background: #ffffff; min-height: 420px; display: flex; align-items: center; justify-content: center;",
                                    plotOutput("gene_expression_plot_OSCC", height = "420px", width = "460px")
                                  )
                                ),
                                #VioPlot
                                div(
                                  class = "plot-card-small",
                                  style = "width: 900px;height: 550px; padding: 25px; flex-shrink: 0;",
                                  div(class = "card-header-custom",
                                      div(class = "card-header-custom",
                                        h4("Comparison by Group"),
                                        uiOutput("download_btn_Gene_VlnPlot_OSCC")
                                      ),
                                      div(style = "width: 150px;",
                                          selectInput("gene_compare_by", NULL, 
                                                      choices = c("Celltype", "Tissue"), 
                                                      selected = "Celltype", width = "100%")
                                      )
                                  ),
                                  div(
                                    style = "background: #ffffff; min-height: 420px; width: 100%; margin-top: 15px;display: flex; flex-direction: column; align-items: center; justify-content: center;overflow-x: auto !important;",
                                    div(
                                      style = "display: inline-block; margin: 0 auto;", 
                                      plotOutput("violin_plot_OSCC", height = "420px", width = "100%")
                                    )
                                  )
                                )
                              ),
                              div(style = "height: 50px;")
                            ),
                            ##DEG Panel--------
                            tabPanel(
                              title = tags$span(icon("table"), "DEGs"),
                              div(
                                style = "background: #ffffff; border-radius: 16px; padding: 20px; margin: 20px 0; box-shadow: 0 4px 15px rgba(0,0,0,0.05);",
                                div(
                                  style = "display: flex; gap: 20px; align-items: flex-end; justify-content: center;",
                                  div(style = "width: 300px;",
                                      selectizeInput("DEG_Celltype_OSCC", tags$b("Select Celltype:"),
                                                     choices = NULL, width = "100%",
                                                     options = list(placeholder = 'Select a cell type...'))
                                  ),
                                  tags$div(style = "margin-bottom: 25px; color: #ccc;", icon("chevron-right")),
                                  div(style = "width: 300px;",
                                      selectizeInput("Condition1_OSCC", tags$b("Condition 1:"), choices = NULL, 
                                                     options = list(placeholder = 'Select a Condition 1...'),
                                                     width = "100%")
                                  ),
                                  tags$div(style = "margin-bottom: 25px; font-weight: bold; color: #6d6fa3;", "VS"),
                                  div(style = "width: 300px;",
                                      selectizeInput("Condition2_OSCC", tags$b("Condition 2:"), choices = NULL,
                                                     options = list(placeholder = 'Select a Condition 2...'),
                                                     width = "100%")
                                  ),
                                  div(style = "width: 150px;",
                                      actionButton(
                                        "Compare_DEGs", "Compare",
                                        icon = icon("sync-alt"),
                                        style = "background: #6d6fa3; color: white; width: 100%; border: none; height: 38px; margin-bottom: 15px;"
                                      )
                                  )
                                )
                              ),
                              div(
                                style = "display: flex; gap: 25px; align-items: flex-start; justify-content: center;",
                                div(
                                  class = "plot-card-small",
                                  style = "width: 700px; padding: 25px; flex-shrink: 0;",
                                  div(class = "card-header-custom",
                                      h4("Differential Expression Gene List"),
                                      downloadButton("download_DEG_table", "Export CSV", class = "btn-download-custom")
                                  ),
                                  div(
                                    style = "background: #ffffff; min-height: 500px; margin-top: 15px;",
                                    DTOutput("DEG_table")
                                  )
                                ),
                                div(
                                  class = "plot-card-small",
                                  style = "width: 800px; padding: 25px; flex-shrink: 0;",
                                  div(class = "card-header-custom",
                                      h4("Volcano Visualization"),
                                      downloadButton("download_Volcano_OSCC", "PDF", class = "btn-download-custom")
                                  ),
                                  div(
                                    style = "background: #ffffff; height: 500px; align-items: center; justify-content: center; margin-top: 15px;display: flex;overflow: hidden;",
                                    plotOutput("DEG_Vol_OSCC", height = "480px", width = "560px")
                                  )
                                )
                              ),
                              div(style = "height: 50px;")
                            ),
                            ##Signature Score Panel--------
                            tabPanel(
                              title = tags$span(icon("sliders-h"), "Signature Score"),
                              
                              fluidRow(
                                style = "margin: 20px 0;",
                                # 左侧控制面板：上传与预览
                                column(4,
                                       div(class = "plot-card-small", style = "padding: 25px; min-height: 650px;",
                                           h4(icon("upload"), "1. Upload Gene Set File"),
                                           tags$p("Upload a file containing a list of genes(e.g., .tsv, .txt, .xlsx, .csv etc.).", 
                                                  style = "color: #888; font-size: 13px;"),
                                           
                                           fileInput("gene_file", NULL, accept = c(".tsv", ".txt", ".xlsx", ".csv"), width = "100%"),
                                           
                                           div(style = "background: #f8f9fa; border-radius: 8px; padding: 15px; margin-top: 10px;",
                                               tags$b("Detected Genes:"),
                                               textOutput("gene_list_text")
                                           ),
                                           
                                           tags$hr(),
                                           
                                           h4(icon("cogs"), "2. Execution"),
                                           div(style = "display: flex; gap: 10px; margin-top: 15px;",
                                               actionButton("Upload_Signature", "Analyze Signature", 
                                                            icon = icon("play"),
                                                            style = "background: #6d6fa3; color: white; border: none; flex-grow: 1; height: 45px;"),
                                               # 重置按钮
                                               actionButton("Reset_Signature", "Reset", 
                                                            icon = icon("undo"),
                                                            style = "background: #f1f1f1; color: #333; border: none; width: 100px;")
                                           )
                                       )
                                ),
                                
                                # 右侧结果面板：动态显示分析结果
                                column(8,
                                       uiOutput("Signature_Score")
                                )
                              )
                            )
                                            # #Signature Score Panel--------
                                            # tabPanel(
                                            #   title = tags$span(icon("sliders-h"), "Signature Score"),
                                            # 
                                            #   fileInput("gene_file","上传基因文件（tsv/txt/xlsx/csv）"),
                                            #   verbatimTextOutput("gene_list_text"),
                                            #   actionButton("Upload_Signature","Upload",
                                            #                style="background:#6d6fa3;color:white;border:none;"),
                                            # 
                                            #   tags$hr(),
                                            # 
                                            #   uiOutput("Signature_Score")
                                            # )
                          )
                      )
                    )
                )
            )
        )
      )
    })
  }
  observeEvent(input$oscc_progression, {
    activePage("oscc_progression")
    oscc_progression()
  })
  ##Oral_diseases UI---------------------
  oral_diseases <- function(){
    output$pageContent <- renderUI({
      tagList(
        
        div(class = "main-content-container",style = "width:1900px;",
            div(
              style="width:100%; text-align:center; padding:25px 0; background:linear-gradient(90deg,#1e1f29,#3a3d5c,#6d6fa3); color:white; border-radius:2px; margin-bottom:25px;",
              h2("Oral Disease T Cell Exploration", style="font-weight:bold; margin:0;"),
              p("Explore T cell dynamic changes across different oral diseases.", style="opacity:0.8; font-size:15px; margin-top:5px;")
            ),
            div(style="display:flex; gap:25px; align-items: flex-start;",
                div(
                  style="width:300px;height:650px; flex-shrink:0; background:white; border-radius:16px; margin-left:20px;padding:20px; box-shadow:0 4px 15px rgba(0,0,0,0.1); color:black;",
                  h4("Data Selection", style="color:#2c3e50; font-weight:bold;"),
                  radioButtons("Oral_Diseases_Names", tags$b("Choose Dataset:"), choices=c(
                    "CAP_GSE181688","CAP_GSE197680","HNSCC_GSE139324","HNSCC_GSE172577","HNSCC_GSE173468","HNSCC_GSE185965",
                    "HNSCC_GSE188737","HNSCC_GSE200996","HNSCC_GSE215403","HNSCC_GSE234933","HNSCC_GSE243359",
                    "OLP_GSE211630", "Periodontitis_GSE152042", "Periodontitis_GSE164241","Periodontitis_GSE171213", 
                    "Periodontitis_GSE207502", "Periodontitis_GSE266897"
                  ), selected="CAP_GSE181688"),
                  radioButtons("Oral_Diseases_Celltype", tags$b("Choose T cell type:"), choices=c("CD4T","CD8T"), selected="CD4T"),
                  actionButton("confirm_button_Oral_Disease","Confirm", class="btn btn-primary", style="width:100%; background:#6d6fa3; border:none;")
                ),
                div(style="flex-grow:1; min-width:700px;",
                    conditionalPanel(
                      condition = "output.isDataLoaded_tab2 == false",
                      div(style="height:650px;background:white; border-radius:16px; text-align:center; box-shadow:0 4px 15px rgba(0,0,0,0.1); color:#999;display: flex; flex-direction: column; justify-content: center; align-items: center; padding: 20px;",
                          icon("hand-pointer", class="fa-3x"),
                          h3("Please select parameters and click 'Confirm' to load data."),
                          p("The analysis results will appear here.")
                      )
                    ),
                    conditionalPanel(
                      condition = "output.isDataLoaded_tab2 == true",
                      div(style="background:white; border-radius:16px; padding:15px; box-shadow:0 4px 15px rgba(0,0,0,0.1);flex-shrink: 0 !important;width:100%;",
                          tabsetPanel(
                            id = "Oral_Diseases_tabs",
                            #Overview Panel------
                            tabPanel(
                              title = tags$span(icon("chart-area"), "Overview"),
                              #第一行：UMAP+Pie Chart
                              div(style="display:flex; flex-shrink: 0 !important;gap:25px; margin-top:20px; justify-content: center;",
                                  div(class="plot-card-small", style="width:620px; padding:25px;flex-shrink: 0 !important;",
                                      div(
                                        class="card-header-custom",
                                        h4("Cluster UMAP Projection"),
                                        downloadButton("download_UMAP_Oral_Diseases", "Download PDF", class="btn-download-custom")
                                      ),
                                      div(class="svg-container", uiOutput("UMAP_plot_Oral_Diseases"))
                                  ),
                                  div(class="plot-card-small", style="width:620px; padding:25px;flex-shrink: 0 !important;",
                                      div(
                                        class="card-header-custom",
                                        h4("Cell Composition Ratio"),
                                        downloadButton("download_Pie_chart_Oral_Diseases", "Download PDF", class="btn-download-custom")
                                      ),
                                      div(class="svg-container", uiOutput("Pie_chart_Oral_Diseases"))
                                  )
                              ),
                              #第二行：Dotplot
                              div(style="display:flex; justify-content: center; margin-top:25px;flex-shrink: 0;",
                                  div(class="plot-card-small", style="width:1265px; padding:25px;",
                                      div(
                                        class="card-header-custom",
                                        h4("Marker Gene Expression Dotplot"),
                                        downloadButton("download_Dot_plot_Oral_Diseases", "Download PDF", class="btn-download-custom")
                                      ),
                                      div(class="svg-container", uiOutput("Dot_plot_Oral_Diseases"))
                                  )
                              ),
                              #第三行：Conditions+Patient Barplots
                              div(
                                style="display:flex; gap:25px; margin-top:25px; justify-content: center;",
                                div(
                                  class="plot-card-small", style="width:430px; padding:25px;",
                                  div(
                                    class="card-header-custom",
                                    h4("Distribution by Condition"),
                                    downloadButton("download_Barplot_conditions_Oral_Diseases", "Download PDF", class="btn-download-custom")
                                  ),
                                  div(class="svg-container1", uiOutput("Barplot_conditions_Oral_Diseases"))
                                ),
                                div(
                                  class="plot-card-small", style="width:810px; padding:25px;",
                                  div(
                                    class="card-header-custom",
                                    h4("Individual Patient Profiles"),
                                    downloadButton("download_Barplot_patients_Oral_Diseases", "Download PDF", class="btn-download-custom")
                                  ),
                                  div(class="svg-container1", uiOutput("Barplot_patients_Oral_Diseases"))
                                )
                              ),
                              div(style="height: 50px;")
                            ),
                            ##Gene Panel----------
                            tabPanel(
                              title = tags$span(icon("dna"), "Gene"),
                              div(
                                style = "background: #ffffff; border-radius: 16px; padding: 20px; margin: 20px 0; box-shadow: 0 4px 15px rgba(0,0,0,0.05);",
                                div(
                                  style = "display: flex; gap: 20px; align-items: flex-end; justify-content: center;",
                                  div(style = "width: 400px;font-size:15pt;",
                                      selectizeInput("gene_input_Oral_Diseases", tags$b("Select Gene:"), choices = NULL, width = "100%")
                                  ),
                                  div(style = "width: 150px;",
                                      actionButton(
                                        "Submit_Gene_Oral_Diseases", "Submit",
                                        icon = icon("paper-plane"),
                                        style = "background: #6d6fa3; color: white; width: 100%; border: none; height: 38px; margin-bottom: 15px;"
                                      )
                                  )
                                )
                              ),
                              div(
                                style = "display: flex; gap: 25px; align-items: flex-start; justify-content: center;",
                                #FeaturePlot
                                div(
                                  class = "plot-card-small",
                                  style = "width: 600px;height: 550px; padding: 25px; flex-shrink: 0;",
                                  div(class = "card-header-custom",
                                      h4("Expression Distribution (Feature Plot)"),
                                      uiOutput("download_btn_Gene_FeaturePlot_Oral_Diseases")
                                  ),
                                  div(
                                    style = "background: #ffffff; min-height: 420px; display: flex; align-items: center; justify-content: center;",
                                    plotOutput("gene_expression_plot_Oral_Diseases", height = "420px", width = "460px")
                                  )
                                ),
                                #VioPlot
                                div(
                                  class = "plot-card-small",
                                  style = "width: 900px;height: 550px; padding: 25px; flex-shrink: 0;",
                                  div(class = "card-header-custom",
                                      div(class = "card-header-custom",
                                          h4("Comparison by Group"),
                                          uiOutput("download_btn_Gene_VlnPlot_Oral_Diseases")
                                      ),
                                      div(style = "width: 150px;",
                                          selectInput("gene_compare_Oral_Diseases", NULL, 
                                                      choices = c("Celltype", "Tissue"), 
                                                      selected = "Celltype", width = "100%")
                                      )
                                  ),
                                  div(
                                    style = "background: #ffffff; min-height: 420px; width: 100%; margin-top: 15px;display: flex; flex-direction: column; align-items: center; justify-content: center;overflow-x: auto !important;",
                                    div(
                                      style = "display: inline-block; margin: 0 auto;", 
                                      plotOutput("violin_plot_Oral_Diseases", height = "420px", width = "100%")
                                    )
                                  )
                                )
                              ),
                              div(style = "height: 50px;")
                            ),
                            ##DEG Panel--------
                            tabPanel(
                              title = tags$span(icon("table"), "DEGs"),
                              div(
                                style = "background: #ffffff; border-radius: 16px; padding: 20px; margin: 20px 0; box-shadow: 0 4px 15px rgba(0,0,0,0.05);",
                                div(
                                  style = "display: flex; gap: 20px; align-items: flex-end; justify-content: center;",
                                  div(style = "width: 300px;",
                                      selectizeInput("DEG_Celltype_Oral_Diseases", tags$b("Select Celltype:"),
                                                     choices = NULL, width = "100%",
                                                     options = list(placeholder = 'Select a cell type...'))
                                  ),
                                  tags$div(style = "margin-bottom: 25px; color: #ccc;", icon("chevron-right")),
                                  div(style = "width: 300px;",
                                      selectizeInput("Condition1_Oral_Diseases", tags$b("Condition 1:"), choices = NULL, 
                                                     options = list(placeholder = 'Select a Condition 1...'),
                                                     width = "100%")
                                  ),
                                  tags$div(style = "margin-bottom: 25px; font-weight: bold; color: #6d6fa3;", "VS"),
                                  div(style = "width: 300px;",
                                      selectizeInput("Condition2_Oral_Diseases", tags$b("Condition 2:"), choices = NULL,
                                                     options = list(placeholder = 'Select a Condition 2...'),
                                                     width = "100%")
                                  ),
                                  div(style = "width: 150px;",
                                      actionButton(
                                        "Compare_DEGs_Oral_Diseases", "Compare",
                                        icon = icon("sync-alt"),
                                        style = "background: #6d6fa3; color: white; width: 100%; border: none; height: 38px; margin-bottom: 15px;"
                                      )
                                  )
                                )
                              ),
                              div(
                                style = "display: flex; gap: 25px; align-items: flex-start; justify-content: center;",
                                div(
                                  class = "plot-card-small",
                                  style = "width: 700px; padding: 25px; flex-shrink: 0;",
                                  div(class = "card-header-custom",
                                      h4("Differential Expression Gene List"),
                                      downloadButton("download_DEG_table_Oral_Diseases", "Export CSV", class = "btn-download-custom")
                                  ),
                                  div(
                                    style = "background: #ffffff; min-height: 500px; margin-top: 15px;",
                                    DTOutput("DEG_table_Oral_Diseases")
                                  )
                                ),
                                div(
                                  class = "plot-card-small",
                                  style = "width: 800px; padding: 25px; flex-shrink: 0;",
                                  div(class = "card-header-custom",
                                      h4("Volcano Visualization"),
                                      downloadButton("download_Volcano_Oral_Diseases", "PDF", class = "btn-download-custom")
                                  ),
                                  div(
                                    style = "background: #ffffff; height: 500px; align-items: center; justify-content: center; margin-top: 15px;display: flex;overflow: hidden;",
                                    plotOutput("DEG_Vol_Oral_Diseases", height = "480px", width = "560px")
                                  )
                                )
                              ),
                              div(style = "height: 50px;")
                            ),
                            ##Signature Score Panel--------
                            tabPanel(
                              title = tags$span(icon("sliders-h"), "Signature Score"),
                              fluidRow(
                                style = "margin: 20px 0;",
                                column(4,
                                       div(class = "plot-card-small", style = "padding: 25px; min-height: 650px;",
                                           h4(icon("upload"), "1. Upload Gene Set File"),
                                           tags$p("Upload a file containing a list of genes(e.g., .tsv, .txt, .xlsx, .csv etc.).", 
                                                  style = "color: #888; font-size: 13px;"),
                                           fileInput("gene_file_Oral_Diseases", NULL, accept = c(".tsv", ".txt", ".xlsx", ".csv"), width = "100%"),
                                           div(style = "background: #f8f9fa; border-radius: 8px; padding: 15px; margin-top: 10px;",
                                               tags$b("Detected Genes:"),
                                               textOutput("gene_list_text_Oral_Diseases")
                                           ),
                                           tags$hr(),
                                           h4(icon("cogs"), "2. Execution"),
                                           div(style = "display: flex; gap: 10px; margin-top: 15px;",
                                               actionButton("Upload_Signature_Oral_Diseases", "Analyze Signature", 
                                                            icon = icon("play"),
                                                            style = "background: #6d6fa3; color: white; border: none; flex-grow: 1; height: 45px;"),
                                               # 重置按钮
                                               actionButton("Reset_Signature_Oral_Diseases", "Reset", 
                                                            icon = icon("undo"),
                                                            style = "background: #f1f1f1; color: #333; border: none; width: 100px;")
                                           )
                                       )
                                ),
                                column(8,
                                       uiOutput("Signature_Score_Oral_Diseases")
                                )
                              )
                            )

                          )
                      )
                    )
                )
            )
        ),
        tags$hr(),
        # fluidRow(
        #   column(
        #     width = 3, style = "width:300px;", 
        #     wellPanel(
        #       style = "font-size: 12pt; line-height:1.6;color:black;",
        #       h4("Oral Disease T Cell Exploration"),
        #       radioButtons(
        #         inputId = "Oral_Diseases_Names",
        #         label = "Choose Dataset:",
        #         choices = c("CAP_GSE181688","CAP_GSE197680","HNSCC_GSE139324","HNSCC_GSE172577","HNSCC_GSE173468",
        #                     "HNSCC_GSE188737","HNSCC_GSE200996","HNSCC_GSE215403","HNSCC_GSE234933","HNSCC_GSE243359",
        #                     "OLP_GSE211630", "Periodontitis_GSE152042", "Periodontitis_GSE164241","Periodontitis_GSE171213", 
        #                     "Periodontitis_GSE207502", "Periodontitis_GSE266897"),
        #         selected = "CAP_GSE181688"
        #       ),
        #       radioButtons(
        #         inputId = "Oral_Diseases_Celltype",
        #         label = "Choose T cell type:",
        #         choices = c("CD4T", "CD8T"),
        #         selected = "CD4T"
        #       ),
        #       actionButton("confirm_button_Oral_Disease", "Confirm", class = "btn-primary")
        #     )
        #   ),
        #   # 右侧动态显示分面框
        #   column(
        #     width = 9,
        #     conditionalPanel(
        #       condition = "output.isDataLoaded_tab2 == true",
        #       div(
        #         tabsetPanel(
        #           id = "Oral_Diseases_tabs",
        #           type = "tabs",
        #           #Overview panel----------
        #             tabPanel(
        #               title = tags$span(icon("chart-area"), "Overview"),
        #               #第一行：UMAP+Pie Chart
        #               div(style="display:flex; gap:25px; margin-top:20px; justify-content: center;",
        #                   div(class="plot-card-small", style="width:620px; padding:25px;",
        #                       div(
        #                         class="card-header-custom",
        #                         h4("Cluster UMAP Projection"),
        #                         downloadButton("download_UMAP_OSCC", "Download PDF", class="btn-download-custom")
        #                       ),
        #                       div(class="svg-container", uiOutput("UMAP_plot_Oral_Diseases"))
        #                   ),
        #                   div(class="plot-card-small", style="width:620px; padding:25px;",
        #                       div(
        #                         class="card-header-custom",
        #                         h4("Cell Composition Ratio"),
        #                         downloadButton("download_Pie_chart_OSCC", "Download PDF", class="btn-download-custom")
        #                       ),
        #                       div(class="svg-container", uiOutput("Pie_chart_Oral_Diseases"))
        #                   )
        #               ),
        #               #第二行：Dotplot
        #               div(style="display:flex; justify-content: center; margin-top:25px;",
        #                   div(class="plot-card-small", style="width:1265px; padding:25px;",
        #                       div(
        #                         class="card-header-custom",
        #                         h4("Marker Gene Expression Dotplot"),
        #                         downloadButton("download_Dot_plot_OSCC", "Download PDF", class="btn-download-custom")
        #                       ),
        #                       div(class="svg-container", uiOutput("Dot_plot_Oral_Diseases"))
        #                   )
        #               ),
        #               #第三行：Conditions+Patient Barplots
        #               div(style="display:flex; gap:25px; margin-top:25px; justify-content: center;",
        #                   div(class="plot-card-small", style="width:430px; padding:25px;",
        #                       div(
        #                         class="card-header-custom",
        #                         h4("Distribution by Condition"),
        #                         downloadButton("download_Barplot_conditions_OSCC", "Download PDF", class="btn-download-custom")
        #                       ),
        #                       div(class="svg-container", uiOutput("Barplot_conditions_Oral_Diseases"))
        #                   ),
        #                   div(class="plot-card-small", style="width:810px; padding:25px;",
        #                       div(
        #                         class="card-header-custom",
        #                         h4("Individual Patient Profiles"),
        #                         downloadButton("download_Barplot_patients_OSCC", "Download PDF", class="btn-download-custom")
        #                       ),
        #                       div(class="svg-container", uiOutput("Barplot_patients_Oral_Diseases"))
        #                   )
        #               ),
        #               #底部留白
        #               div(style="height: 50px;")
        #             ),
        #           #Gene panel-----------
        #           tabPanel(
        #             "Gene",
        #             fluidRow(
        #               column(
        #                 width = 4,
        #                 selectizeInput("gene_input_Oral_Diseases",
        #                                "Gene (Please input a gene):", choices = NULL)
        #               ),
        #               column(
        #                 width = 4,
        #                 actionButton(
        #                   inputId = "Submit_Gene_Oral_Diseases",
        #                   label = "Submit",
        #                   icon = icon("paper-plane"),
        #                   style = "background-color: #87CEEB; color: white !important;
        #                          border: 2px solid #87CEEB; font-size:15px"
        #                 )
        #               )
        #             ),
        #             
        #             fluidRow(
        #               column(
        #                 width = 12,
        #                 plotlyOutput("gene_expression_plot_Oral_Diseases",
        #                              height = "500px", width = "500px")
        #               )
        #             ),
        #             
        #             tags$hr(),
        #             
        #             fluidRow(
        #               column(
        #                 width = 12,
        #                 selectInput(
        #                   "gene_compare_Oral_Diseases",
        #                   "The Violin plot of gene expression, compared by:",
        #                   choices = c("celltype", "Tissue")
        #                 ),
        #                 plotOutput("violin_plot_Oral_Diseases",
        #                            height = "300px", width = "auto")
        #               )
        #             )
        #           )
        #         )
        #       )
        #     )
        #   )
        # )
      )
    })
  }
  observeEvent(input$oral_diseases, {
    activePage("oral_diseases")
    oral_diseases()
  })
  ##Pancancer---------------------
  pancancer <- function(){
    output$pageContent <- renderUI({
      tagList(
        div(class = "main-content-container",style = "width:1900px;",
            div(
              style="width:100%; text-align:center; padding:25px 0; background:linear-gradient(90deg,#1e1f29,#3a3d5c,#6d6fa3); color:white; border-radius:2px; margin-bottom:25px;",
              h2("Pan-Cancer", style="font-weight:bold; margin:0;"),
              p("Explore pan-cancer insights related to T cell states.", style="opacity:0.8; font-size:15px; margin-top:5px;")
            ),
            div(style="display:flex; gap:25px; align-items: flex-start;",
                div(
                  style="width:300px;height:280px; flex-shrink:0; background:white; border-radius:16px; margin-left:20px;padding:20px; box-shadow:0 4px 15px rgba(0,0,0,0.1); color:black;",
                  h4("Data Selection", style="color:#2c3e50; font-weight:bold;"),
                  radioButtons("Pancancer_Names", tags$b("Choose Dataset:"), choices=c(
                    "Pancancer"
                  ), selected="Pancancer"),
                  radioButtons("Pancancer_Celltype", tags$b("Choose T cell type:"), choices=c("CD4T","CD8T"), selected="CD4T"),
                  actionButton("confirm_button_Pancancer","Confirm", class="btn btn-primary", style="width:100%; background:#6d6fa3; border:none;")
                ),
                div(style="flex-grow:1; min-width:700px;",
                    conditionalPanel(
                      condition = "output.isDataLoaded_tab3 == false",
                      div(style="height:650px;background:white; border-radius:16px; text-align:center; box-shadow:0 4px 15px rgba(0,0,0,0.1); color:#999;display: flex; flex-direction: column; justify-content: center; align-items: center; padding: 20px;",
                          icon("hand-pointer", class="fa-3x"),
                          h3("Please select parameters and click 'Confirm' to load data."),
                          p("The analysis results will appear here.")
                      )
                    ),
                    conditionalPanel(
                      condition = "output.isDataLoaded_tab3 == true",
                      div(style="background:white; border-radius:16px; padding:15px; box-shadow:0 4px 15px rgba(0,0,0,0.1);flex-shrink: 0 !important;width:100%;",
                          tabsetPanel(
                            id = "Pancancer_tabs",
                            #Overview Panel------
                            tabPanel(
                              title = tags$span(icon("chart-area"), "Overview"),
                              #第一行：UMAP
                              div(style="display:flex; flex-shrink: 0 !important;gap:25px; margin-top:20px; justify-content: center;",
                                  div(class="plot-card-small", style="width:1265px; padding:25px;flex-shrink: 0 !important;",
                                      div(
                                        class="card-header-custom",
                                        h4("Cluster UMAP Projection"),
                                        downloadButton("download_UMAP_Pancancer", "Download PNG", class="btn-download-custom")
                                      ),
                                      div(class="svg-container", uiOutput("UMAP_plot_Pancancer"))
                                  )
                              ),
                              #第二行：Dotplot
                              div(style="display:flex; justify-content: center; margin-top:25px;flex-shrink: 0;",
                                  div(class="plot-card-small", style="width:1265px; padding:25px;",
                                      div(
                                        class="card-header-custom",
                                        h4("Marker Gene Expression Dotplot"),
                                        downloadButton("download_Dot_plot_Pancancer", "Download PNG", class="btn-download-custom")
                                      ),
                                      div(class="svg-container", uiOutput("Dot_plot_Pancancer"))
                                  )
                              ),
                              #第三行：Conditions+Patient Barplots
                              div(
                                style="display:flex; gap:25px; margin-top:25px; justify-content: center;",
                                div(
                                  class="plot-card-small", style="width:490px; padding:25px;",
                                  div(
                                    class="card-header-custom",
                                    h4("Distribution by Condition"),
                                    downloadButton("download_Barplot_conditions_Pancancer", "Download PNG", class="btn-download-custom")
                                  ),
                                  div(class="svg-container1", uiOutput("Barplot_conditions_Pancancer"))
                                ),
                                div(
                                  class="plot-card-small", style="width:750px; padding:25px;",
                                  div(
                                    class="card-header-custom",
                                    h4("Individual Patient Profiles"),
                                    downloadButton("download_Barplot_patients_Pancancer", "Download PNG", class="btn-download-custom")
                                  ),
                                  div(class="svg-container1", uiOutput("Barplot_patients_Pancancer"))
                                )
                              ),
                              div(style="height: 50px;")
                            ),
                            ##Gene Panel----------
                            tabPanel(
                              title = tags$span(icon("dna"), "Gene"),
                              div(
                                style = "background: #ffffff; border-radius: 16px; padding: 20px; margin: 20px 0; box-shadow: 0 4px 15px rgba(0,0,0,0.05);",
                                div(
                                  style = "display: flex; gap: 20px; align-items: flex-end; justify-content: center;",
                                  div(style = "width: 400px;font-size:15pt;",
                                      selectizeInput("gene_input_Pancancer", tags$b("Select Gene:"), choices = NULL, width = "100%")
                                  ),
                                  div(style = "width: 150px;",
                                      actionButton(
                                        "Submit_Gene_Pancancer", "Submit",
                                        icon = icon("paper-plane"),
                                        style = "background: #6d6fa3; color: white; width: 100%; border: none; height: 38px; margin-bottom: 15px;"
                                      )
                                  )
                                )
                              ),
                              div(
                                style = "display: flex; gap: 25px; align-items: flex-start; justify-content: center;",
                                #FeaturePlot
                                div(
                                  class = "plot-card-small",
                                  style = "width: 600px;height: 580px; padding: 25px; flex-shrink: 0;",
                                  div(class = "card-header-custom",
                                      h4("Expression Distribution (Feature Plot)"),
                                      uiOutput("download_btn_Gene_FeaturePlot_Pancancer")
                                  ),
                                  div(
                                    style = "background: #ffffff; min-height: 420px; display: flex; align-items: center; justify-content: center;",
                                    uiOutput("gene_expression_plot_Pancancer", height = "420px", width = "460px")
                                  )
                                ),
                                #VioPlot
                                div(
                                  class = "plot-card-small",
                                  style = "width: 900px;height: 580px; padding: 25px; flex-shrink: 0;",
                                  div(class = "card-header-custom",
                                      div(class = "card-header-custom",
                                          h4("Comparison by Group"),
                                          uiOutput("download_btn_Gene_VlnPlot_Pancancer")
                                      ),
                                      div(style = "width: 150px;",
                                          selectInput("gene_compare_Pancancer", NULL,
                                                      choices = c("Celltype", "Cancer_type"),
                                                      selected = "Celltype", width = "100%")
                                      )
                                  ),
                                  div(
                                    style = "background: #ffffff; min-height: 420px; width: 100%; margin-top: 15px;display: flex; flex-direction: column; align-items: center; justify-content: center;overflow-x: auto !important;",
                                    div(
                                      style = "display: inline-block; margin: 0 auto;",
                                      uiOutput("violin_plot_Pancancer", height = "420px", width = "100%")
                                    )
                                  )
                                )
                              ),
                              div(style = "height: 50px;")
                            )
                            
                          )
                      )
                    )
                )
            )
        )
      )
    })
  }
  observeEvent(input$pancancer, {
    activePage("pancancer")
    pancancer()
  })
  #主页三大板块图片---------
  output$OSCC_figure <- renderImage({
    list(src = "www/OSCC progression .png", contentType = "image/png", width = "100%")
  }, deleteFile = FALSE)
  
  output$Oral_figure <- renderImage({
    list(src = "www/Oral diseases.png", contentType = "image/png", width = "100%")
  }, deleteFile = FALSE)
  
  output$Pancancer_figure <- renderImage({
    list(src = "www/pancancer.png", contentType = "image/png", width = "100%")
  }, deleteFile = FALSE)
  
  #图片点击事件
  observeEvent(input$OSCC_link, {
    activePage("oscc_progression")
    output$pageContent <- renderUI({
      oscc_progression()
    })
  })
  
  observeEvent(input$Oral_link, {
    activePage("oral_diseases")
    oral_diseases()
  })
  
  observeEvent(input$Pancancer_link, {
    activePage("pancancer")
    pancancer()
  })
  #初始化数据加载状态
  isDataLoaded_tab1 <- reactiveVal(FALSE)
  isDataLoaded_tab2 <- reactiveVal(FALSE)
  isDataLoaded_tab3 <- reactiveVal(FALSE)
  ##OSCC Progression 后端--------------------
  #数据加载状态供前端使用
  output$isDataLoaded_tab1 <- reactive({
    isDataLoaded_tab1()
  })
  outputOptions(output, "isDataLoaded_tab1", suspendWhenHidden = FALSE)
  #点击确认按钮后加载数据
  observeEvent(input$confirm_button, {
    #更新数据加载状态为TRUE
    isDataLoaded_tab1(TRUE)
    
    #gene name选项信息
    load(paste0("Dataset/OSCC/",input$OSCC_dataset,"/",input$OSCC_Celltype,"/gene_names.Rdata"))
    load(paste0("Dataset/OSCC/",input$OSCC_dataset,"/",input$OSCC_Celltype,"/metadata.Rdata"))
    updateSelectizeInput(session = session,inputId = "gene_input_OSCC",choices=c("",sort(gene_names)),server=TRUE)
    updateSelectizeInput(session = session,inputId = "DEG_Celltype_OSCC",choices=c("",unique(as.character(metadata$Celltype))),server=TRUE)
    updateSelectizeInput(session = session,inputId = "Condition1_OSCC",choices=c("",unique(as.character(metadata$Tissue))),server=TRUE)
    updateSelectizeInput(session = session,inputId = "Condition2_OSCC",choices=c("",unique(as.character(metadata$Tissue))),server=TRUE)
    #UMAP图
    output$UMAP_plot <- renderUI({
      tags$img(src = paste0("Dataset/OSCC/",input$OSCC_dataset,"/",input$OSCC_Celltype,"/UMAP.svg"),
               style = "width: auto; height: 450px; object-fit: contain;"
               )
    })
    output$download_UMAP_OSCC <- downloadHandler(
      filename = function() {
        paste0("UMAP_", input$OSCC_dataset, "_", input$OSCC_Celltype, ".pdf")
      },
      content = function(file) {
        path <- paste0("www/Dataset/OSCC/", input$OSCC_dataset, "/", input$OSCC_Celltype, "/UMAP.pdf")
        file.copy(path, file)
      }
    )
    #饼图
    output$Pie_chart <- renderUI({
      #svgPanZoomOutput("mysvg", width = "100%", height = "500px")
      tags$img(
        src = paste0("Dataset/OSCC/",input$OSCC_dataset,"/",input$OSCC_Celltype,"/Pieplot.svg"),
        style = "width: auto; height: 450px; object-fit: contain;"
      )
    })
    # output$mysvg <- renderSvgPanZoom({
    #   svgPanZoom(paste0("www/Dataset/OSCC/",input$OSCC_dataset,"/",input$OSCC_Celltype,"/Pieplot.svg"), controlIconsEnabled = TRUE)
    # })
    output$download_Pie_chart_OSCC <- downloadHandler(
      filename = function() {
        paste0("Pie_chart_", input$OSCC_dataset, "_", input$OSCC_Celltype, ".pdf")
      },
      content = function(file) {
        path <- paste0("www/Dataset/OSCC/", input$OSCC_dataset, "/", input$OSCC_Celltype, "/Pieplot.pdf")
        file.copy(path, file)
      }
    )
    #Dotplot
    output$Dot_plot <- renderUI({
      tags$img(
        src = paste0("Dataset/OSCC/",input$OSCC_dataset,"/",input$OSCC_Celltype,"/Dotplot.svg"),
        style = "width: 1000px; height: auto; object-fit: contain;"
      )
    })
    output$download_Dot_plot_OSCC <- downloadHandler(
      filename = function() {
        paste0("Dot_plot_", input$OSCC_dataset, "_", input$OSCC_Celltype, ".pdf")
      },
      content = function(file) {
        path <- paste0("www/Dataset/OSCC/", input$OSCC_dataset, "/", input$OSCC_Celltype, "/Dotplot.pdf")
        file.copy(path, file)
      }
    )
    #柱状图：按条件分组
    output$Barplot_conditions <- renderUI({
      tags$img(src = paste0("Dataset/OSCC/",input$OSCC_dataset,"/",input$OSCC_Celltype,"/Barplot_Condition.svg"), 
               #type = "image/svg+xml",
               style = "height: 400px;flex-shrink: 0;"
      )
    })
    output$download_Barplot_conditions_OSCC <- downloadHandler(
      filename = function() {
        paste0("Barplot_condition_", input$OSCC_dataset, "_", input$OSCC_Celltype, ".pdf")
      },
      content = function(file) {
        path <- paste0("www/Dataset/OSCC/", input$OSCC_dataset, "/", input$OSCC_Celltype, "/Barplot_Condition.pdf")
        file.copy(path, file)
      }
    )
    #柱状图：按患者分组
    output$Barplot_patients <- renderUI({
      tags$img(src = paste0("Dataset/OSCC/",input$OSCC_dataset,"/",input$OSCC_Celltype,"/Barplot_Patient.svg"), 
               #type = "image/svg+xml",
               style = "height: 400px;flex-shrink: 0;"
      )
    })
    output$download_Barplot_patients_OSCC <- downloadHandler(
      filename = function() {
        paste0("Barplot_patient_", input$OSCC_dataset, "_", input$OSCC_Celltype, ".pdf")
      },
      content = function(file) {
        path <- paste0("www/Dataset/OSCC/", input$OSCC_dataset, "/", input$OSCC_Celltype, "/Barplot_Patient.pdf")
        file.copy(path, file)
      }
    )
    #Gene panel--------------------------
    show_download <- reactiveVal(FALSE)
    observeEvent(input$Submit_Gene_OSCC,{
      if(input$gene_input_OSCC != ""){
        req(input$gene_input_OSCC != "")
        file_path = paste0("Dataset/OSCC/",input$OSCC_dataset,"/",input$OSCC_Celltype,"/seurat_plot.Rdata")
        load(file_path)
        dataset <- seurat_plot
      }
      output$gene_expression_plot_OSCC <- renderPlot({
        req(input$gene_input_OSCC != "")
        if((input$gene_input_OSCC != "")&(input$gene_input_OSCC %in% rownames(dataset))){
          gene_data <- FetchData(dataset, vars = input$gene_input_OSCC, slot = "data")[,1]
          min_val <- 0
          mid_val <- quantile(gene_data, 0.8) 
          max_val <- max(gene_data)
          rescale_mid <- (mid_val - min_val) / (max_val - min_val)
          p <- FeaturePlot(
            dataset, features = input$gene_input_OSCC, order = TRUE, 
            cols = c("lightgrey", "#ff0000"), pt.size = 0.5, combine = TRUE) + 
            theme_bw() + 
            theme(
              panel.border = element_rect(fill = NA, color = "black", size = 1.2), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_blank(), 
              plot.title = element_text(hjust = 0.5, size = 20, family = "Arial", face = "italic"),
              axis.title = element_text(size = 17, family = "Arial", color = "black"),
              axis.text = element_text(size = 15, family = "Arial", color = "black"),
              legend.position = "right",
              legend.text = element_text(size = 15, family = "Arial"),
              legend.title = element_text(size = 17, family = "Arial")
            ) + 
            xlab("UMAP_1") + ylab("UMAP_2") +
            scale_color_gradientn(
              colors = BlueAndRed(),
              values = c(0, rescale_mid * 0.99, rescale_mid, rescale_mid * 1.01, 1),
              limits = c(min_val, max_val),
              na.value = "lightgray",
              guide = guide_colorbar(
                #title = "Expression",
                frame.colour = "black", 
                ticks.colour = "black"
              )
            )
            #scale_color_gradientn(colors = BlueAndRed())
          p
        }else{}
        })
      output$violin_plot_OSCC <- renderPlot({
        req(input$gene_input_OSCC != "", input$gene_compare_by != "")
        if((input$gene_input_OSCC!="") & (input$gene_input_OSCC %in% rownames(dataset))){
          theme_clean_pub <- function() { 
            theme_classic(base_size = 18) +  
              theme(   
                panel.background = element_rect(fill = "white", color = NA),   
                plot.background = element_rect(fill = "white", color = NA),   
                panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
                axis.line = element_blank(),
                axis.ticks = element_line(color = "black", linewidth = 0.8),
                axis.ticks.length = unit(0.15, "cm"),   
                axis.text = element_text(color = "black", size = 14, family = "Arial"),
                axis.title = element_text(color = "black", size = 16, family = "Arial"), 
                plot.title = element_text(size = 18, face = "bold.italic", hjust = 0.5, family = "Arial"), 
                plot.subtitle = element_text(size = 12, hjust = 0.5, family = "Arial"),   
                strip.background = element_blank(), 
                strip.text = element_text(size = 14, face = "bold", family = "Arial"),   
                panel.grid = element_blank(),
                legend.position = "right",
                legend.title = element_text(size = 12, face = "bold", family = "Arial"),   
                legend.text = element_text(size = 11, family = "Arial"),   
                legend.background = element_blank(),
                legend.key = element_blank()
              )
          }
          if(input$gene_compare_by == "Celltype"){
            genes_to_extract <- input$gene_input_OSCC
            group_col <- input$gene_compare_by
            plot_data <- FetchData(
              object = dataset, 
              vars = c(genes_to_extract, group_col), 
              slot = "data"                         
            )
            p <- ggplot(plot_data, aes(x = .data[[group_col]], y = .data[[genes_to_extract]], fill = .data[[group_col]])) +
              geom_violin(trim = FALSE, 
                          color = NA, 
                          alpha = 0.8, 
                          width = 0.9) +
              geom_boxplot(width = 0.1, 
                           linewidth = 0.8, 
                           color = "black", 
                           outlier.shape = NA, 
                           alpha = 1) +
              stat_summary(fun = mean, 
                           geom = "point", 
                           shape = 23,  size =4,  
                           fill ="white",  color ="black",  stroke =1.4 ) +
              scale_fill_manual(values = colors) + 
              theme_bw()+
              labs(  title ="",  x ="",  y =paste0("The expression levels of ",genes_to_extract)) +
              theme_clean_pub() + theme(legend.position ="none")
            p
          }else if(input$gene_compare_by == "Tissue"){
            genes_to_extract <- input$gene_input_OSCC
            group_x <- "Celltype"  
            group_col <- input$gene_compare_by
            plot_data <- FetchData(
              object = dataset, 
              vars = c(genes_to_extract, group_x, group_col), 
              slot = "data"
            )
            p <- ggplot(plot_data, aes(x = .data[[group_x]], 
                                       y = .data[[genes_to_extract]], 
                                       fill = .data[[group_col]])) +
              geom_violin(trim = FALSE, 
                          color = "black", 
                          linewidth = 0.3,
                          alpha = 0.7, 
                          width = 0.8,
                          position = position_dodge(width = 0.9)) + 
              geom_boxplot(
                aes(group = interaction(.data[[group_x]], .data[[group_col]])),
                width = 0.15, 
                linewidth = 0.5, 
                color = "black", 
                outlier.shape = NA, 
                fill = "white",
                position = position_dodge(width = 0.9)) + 
              stat_summary(aes(group = .data[[group_col]]),
                           fun = mean, 
                           geom = "point", 
                           shape = 23, size = 1.5,  
                           fill = "white", color = "black",
                           position = position_dodge(width = 0.9)) +
              
              scale_fill_manual(values = colors) + 
              labs(x = "", y = paste0("Expression levels of ", genes_to_extract), fill = "Tissue") +
              theme_clean_pub() + 
              theme(
                legend.position = "right",
                axis.text.x = element_text(angle = 45, hjust = 1)
              )
            p
          }
        }
      }, width = function() {
        if(input$gene_compare_by == "Celltype"){
          n_groups <- length(unique(dataset@meta.data[[input$gene_compare_by]]))
          calculated_width <- n_groups * 80 + 200
        }else if(input$gene_compare_by == "Tissue"){
          n_groups <- length(unique(dataset@meta.data[[input$gene_compare_by]]))
          m_groups <- length(unique(dataset@meta.data[['Celltype']]))
          calculated_width <- n_groups * m_groups * 80 + 200
        }
        return(max(600, calculated_width))
      }
      )
      show_download(TRUE)
    })
    output$download_btn_Gene_FeaturePlot_OSCC <- renderUI({
      if (show_download()) {
        downloadButton(
          "download_Gene_FeaturePlot_OSCC", 
          "Download PDF", 
          class = "btn-download-custom"
        )
      } else {
        return(NULL) 
      }
    })
    output$download_Gene_FeaturePlot_OSCC <- downloadHandler(
      filename = function() {
        paste0("FeaturePlot_", input$gene_input_OSCC, "_", input$OSCC_dataset, "_", input$OSCC_Celltype, ".pdf")
      },
      content = function(file) {
        req(input$gene_input_OSCC != "")
        file_path = paste0("Dataset/OSCC/",input$OSCC_dataset,"/",input$OSCC_Celltype,"/seurat.Rdata")
        load(file_path)
        dataset <- seurat
        gene_data <- FetchData(dataset, vars = input$gene_input_OSCC, slot = "data")[,1]
        min_value <- quantile(gene_data, 0.8) 
        max_value <- max(gene_data)
        p <- FeaturePlot(
          dataset, features = input$gene_input_OSCC, order = TRUE, 
          cols = c("lightgrey", "#ff0000"), pt.size = 0.5, combine = TRUE) + 
          theme_bw() + 
          theme(
            panel.border = element_rect(fill = NA, color = "black", size = 1.2), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_blank(), 
            plot.title = element_text(hjust = 0.5, size = 20,  face = "italic", family = "Arial"),
            axis.title = element_text(size = 17,  color = "black", family = "Arial"),
            axis.text = element_text(size = 15,  color = "black", family = "Arial"),
            legend.position = "right",
            legend.text = element_text(size = 15, family = "Arial"),
            legend.title = element_text(size = 17, family = "Arial")
          ) + 
          xlab("UMAP_1") + ylab("UMAP_2") +
          scale_color_gradientn(colors = BlueAndRed(),
                                limits = c(min_value, max_value),       
                                na.value = "lightgray")
        ggsave(file, plot = p, device = "pdf", width = 5.5, height = 5, units = "in")
      }
    )
    output$download_btn_Gene_VlnPlot_OSCC <- renderUI({
      if (show_download()) {
        downloadButton(
          "download_Gene_VlnPlot_OSCC", 
          "Download PDF", 
          class = "btn-download-custom"
        )
      } else {
        return(NULL) 
      }
    })
    output$download_Gene_VlnPlot_OSCC <- downloadHandler(
      filename = function() {
        paste0("VlnPlot_", input$gene_input_OSCC, "_", input$OSCC_dataset, "_", input$OSCC_Celltype, ".pdf")
      },
      content = function(file) {
        req(input$gene_input_OSCC != "")
        file_path = paste0("Dataset/OSCC/",input$OSCC_dataset,"/",input$OSCC_Celltype,"/seurat.Rdata")
        load(file_path)
        dataset <- seurat
        theme_clean_pub <- function() { 
          theme_classic(base_size = 18) +  
            theme(   
              panel.background = element_rect(fill = "white", color = NA),   
              plot.background = element_rect(fill = "white", color = NA),   
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
              axis.line = element_blank(),
              axis.ticks = element_line(color = "black", linewidth = 0.8),
              axis.ticks.length = unit(0.15, "cm"),   
              axis.text = element_text(color = "black", size = 14, family = "Arial"),
              axis.title = element_text(color = "black", size = 16, family = "Arial"), 
              plot.title = element_text(size = 18, face = "bold.italic", hjust = 0.5, family = "Arial"), 
              plot.subtitle = element_text(size = 12, hjust = 0.5, family = "Arial"),   
              strip.background = element_blank(), 
              strip.text = element_text(size = 14, face = "bold", family = "Arial"),   
              panel.grid = element_blank(),
              legend.position = "right",
              legend.title = element_text(size = 12, face = "bold", family = "Arial"),   
              legend.text = element_text(size = 11, family = "Arial"),   
              legend.background = element_blank(),
              legend.key = element_blank()
            )
        }
        if(input$gene_compare_by == "Celltype"){
          genes_to_extract <- input$gene_input_OSCC
          group_col <- input$gene_compare_by
          plot_data <- FetchData(
            object = dataset, 
            vars = c(genes_to_extract, group_col), 
            slot = "data"                         
          )
          p <- ggplot(plot_data, aes(x = .data[[group_col]], y = .data[[genes_to_extract]], fill = .data[[group_col]])) +
            geom_violin(trim = FALSE, 
                        color = NA, 
                        alpha = 0.8, 
                        width = 0.9) +
            geom_boxplot(width = 0.1, 
                         linewidth = 0.8, 
                         color = "black", 
                         outlier.shape = NA, 
                         alpha = 1) +
            stat_summary(fun = mean, 
                         geom = "point", 
                         shape = 23,  size =4,  
                         fill ="white",  color ="black",  stroke =1.4 ) +
            scale_fill_manual(values = colors) + 
            theme_bw()+
            labs(  title ="",  x ="",  y =paste0("The expression levels of ",genes_to_extract)) +
            theme_clean_pub() + theme(legend.position ="none")
        }else if(input$gene_compare_by == "Tissue"){
          genes_to_extract <- input$gene_input_OSCC
          group_x <- "Celltype"  
          group_col <- input$gene_compare_by
          plot_data <- FetchData(
            object = dataset, 
            vars = c(genes_to_extract, group_x, group_col), 
            slot = "data"
          )
          p <- ggplot(plot_data, aes(x = .data[[group_x]], 
                                     y = .data[[genes_to_extract]], 
                                     fill = .data[[group_col]])) +
            geom_violin(trim = FALSE, 
                        color = "black", 
                        linewidth = 0.3,
                        alpha = 0.7, 
                        width = 0.8,
                        position = position_dodge(width = 0.9)) + 
            geom_boxplot(
              aes(group = interaction(.data[[group_x]], .data[[group_col]])),
              width = 0.15, 
              linewidth = 0.5, 
              color = "black", 
              outlier.shape = NA, 
              fill = "white",
              position = position_dodge(width = 0.9)) + 
            stat_summary(aes(group = .data[[group_col]]),
                         fun = mean, 
                         geom = "point", 
                         shape = 23, size = 1.5,  
                         fill = "white", color = "black",
                         position = position_dodge(width = 0.9)) +
            
            scale_fill_manual(values = colors) + 
            labs(x = "", y = paste0("Expression levels of ", genes_to_extract), fill = "Tissue") +
            theme_clean_pub() + 
            theme(
              legend.position = "right",
              axis.text.x = element_text(angle = 45, hjust = 1)
            )
        }
        if(input$gene_compare_by == "Celltype"){
          n_groups <- length(unique(dataset@meta.data[[input$gene_compare_by]]))
          calculated_width <- n_groups * 80 + 200
        }else if(input$gene_compare_by == "Tissue"){
          n_groups <- length(unique(dataset@meta.data[[input$gene_compare_by]]))
          m_groups <- length(unique(dataset@meta.data[['Celltype']]))
          calculated_width <- n_groups * m_groups * 80 + 200
        }
        ggsave(file, plot = p, device = "pdf", width = calculated_width/72, height = 450/72, units = "in")
      }
    )
    ##DEGs Panel---------
    observeEvent(input$Compare_DEGs,{
      req(input$Condition1_OSCC != "", input$Condition2_OSCC != "")
      if(is.null(input$DEG_Celltype_OSCC)) {
        showNotification(
          "Please select Celltype！", 
          type = "warning",
          duration = 3
        )
        return()
      }
      
      if(input$Condition1_OSCC == input$Condition2_OSCC) {
        showNotification(
          "Please select two different conditions for comparison！", 
          type = "warning",
          duration = 3
        )
        return()
      }
      ##DEGs table------
      output$DEG_table <- renderDT({
        req()
        filename = paste0("Dataset/OSCC/",input$OSCC_dataset,"/",input$OSCC_Celltype,"/DEGs/",input$DEG_Celltype_OSCC,"/",input$DEG_Celltype_OSCC,"_",input$Condition1_OSCC,"_vs_",input$Condition2_OSCC,".Rdata")
        load(filename)
        #DEGs <- FindMarkers(dataset,group.by = "Tissue",ident.1 = input$Condition1,ident.2 = input$Condition2,min.pct=0.25,logfc.threshold = 0.25)
        #DEGs <- DEGs[which(abs(DEGs$avg_log2FC) >= 0.5 & DEGs$p_val<0.05),]
        #DEGs$gene <- rownames(DEGs)
        #DEGs <- DEGs[,c("gene","avg_log2FC","p_val","p_val_adj","pct.1","pct.2")]
        #DEGs[, sapply(DEGs, is.numeric)] <- round(DEGs[, sapply(DEGs, is.numeric)], 2)
        DEGs
      },escape = FALSE,
      rownames = FALSE,
      options = list(
        paging = TRUE,
        searching = TRUE,
        info = TRUE,
        lengthMenu = NULL,
        columnDefs = list(
          list(width = '80px',targets = c(0,1,2,3,4,5))
        )
      )
      )
      ##Volcano---------
      output$DEG_Vol_OSCC <- renderPlot({
        req()
        filename = paste0("Dataset/OSCC/",input$OSCC_dataset,"/",input$OSCC_Celltype,"/DEGs/",input$DEG_Celltype_OSCC,"/",input$DEG_Celltype_OSCC,"_",input$Condition1_OSCC,"_vs_",input$Condition2_OSCC,".Rdata")
        load(filename)
        DEGs$p_val_adj <- as.numeric(DEGs$p_val_adj)
        DEGs$avg_log2FC <- as.numeric(DEGs$avg_log2FC)
        DEGs <- DEGs[order(DEGs$p_val_adj,DEGs$avg_log2FC,decreasing = c(FALSE,TRUE)),]
        DEGs_up <- DEGs[which(DEGs$avg_log2FC >= 1 & DEGs$p_val_adj<0.05),]
        DEGs_down <- DEGs[which(DEGs$avg_log2FC <= -1 & DEGs$p_val_adj<0.05),]
        DEGs_up <- DEGs_up[order(DEGs_up$avg_log2FC,decreasing = TRUE),]
        DEGs_down <- DEGs_down[order(DEGs_down$avg_log2FC,decreasing = FALSE),]
        Upvals <- rownames(DEGs_up)[1:5]
        Downvals <- rownames(DEGs_down)[1:5]
        vals <- c(Upvals,Downvals)
        group<-ifelse(
          DEGs$avg_log2FC<(-0.5)&DEGs$p_val_adj<0.05,'#4D4398',
          ifelse(DEGs$avg_log2FC>(0.5)&DEGs$p_val_adj<0.05,'indianred1',
                 '#b5b5b5'))
        max <- max(abs(na.omit(DEGs[,c("avg_log2FC")])))
        x_max <- max*(11/10)
        x_min <- -x_max
        y_max <- -log10(min(DEGs$p_val_adj[DEGs$p_val_adj > 0]))+1
        group[is.na(group)]<-'#b5b5b5'
        names(group)[group=='indianred1']<-'Up'
        names(group)[group=='#b5b5b5']<-'No_sig'
        names(group)[group=='#4D4398']<-'Down'
        if(nrow(DEGs)>0){
          g <- EnhancedVolcano(DEGs,
                               x = "avg_log2FC",
                               y = "p_val_adj",
                               lab = rownames(DEGs),
                               pCutoff = 0.05,
                               FCcutoff = 0.1,
                               pointSize=c(ifelse(rownames(DEGs) %in% vals,1.5,0.5)),
                               #pointSize = 3,
                               labSize = 5,
                               xlim = c(x_min,x_max),
                               ylim = c(0,y_max),
                               colCustom = group,
                               title = NULL,
                               subtitle = NULL,
                               caption = NULL,
                               legendPosition = "right",
                               selectLab = c(Upvals,Downvals),
                               xlab = bquote(~Log[2]~'fold change'),
                               legendLabSize = 12,
                               legendIconSize = 6,
                               labCol = 'black',
                               drawConnectors = TRUE,typeConnectors = 'closed',lengthConnectors = unit(0, "mm"),
                               widthConnectors = 0.5,
                               max.overlaps = 30

          )+
            theme(
              panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
              plot.title = element_text(
                size = 15,
                face = "bold",
                hjust = 0.5,
                vjust = 2,
                margin = margin(b = 10)
              ),
              axis.title = element_text(size = 14),
              legend.title = element_blank(),
              text = element_text()
            )
          g
        }else{}
      },width = 560,height = 400)
      
      
    })
    
    ##Gene_Signature--------
    parsed_genes <- reactive({
      req(input$gene_file)
      file <- input$gene_file
      ext <- tools::file_ext(file$name)
      df <- switch(ext,
                   "csv"  = read.csv(file$datapath, header = TRUE),
                   "tsv"  = read.table(file$datapath, header = TRUE, sep = "\t"),
                   "txt"  = read.table(file$datapath, header = TRUE, sep = "\t"),
                   "xlsx" = readxl::read_excel(file$datapath),
                   { showNotification("不支持格式", type = "error"); return(NULL) }
      )
      genes <- as.character(df[[1]])
      genes <- unique(trimws(genes))
      genes <- genes[genes != "" & !is.na(genes)]
      return(genes)
    })
    
    output$gene_list_text <- renderText({
      g <- parsed_genes()
      req(g)
      #将基因向量合并成用逗号分隔的字符串
      paste(g, collapse = ", ")
    })
    
    observeEvent(input$Upload_Signature, {
      genes <- parsed_genes()
      valid_genes <- intersect(genes, rownames(dataset))
      if(length(valid_genes) == 0) {
        showNotification("上传的基因均不在当前数据集中，请检查符号。", type = "error")
        return(NULL)
      }
      
      #计算分数(AddModuleScore)
      withProgress(message = 'Calculating Signature Score...', {
        sce_score <- AddModuleScore(
          object = dataset,
          features = list(valid_genes),
          name = "Signature_score"
        )
        colnames(sce_score@meta.data)[colnames(sce_score@meta.data) == "Signature_score1"] <- "Signature_score"
      })
      
      output$Signature_Score <- renderUI({
        tagList(
          div(class = "plot-card-small", style = "padding: 20px; margin-bottom: 20px;",
              h4("Feature Plot"),
              plotOutput("Signature_FeaturePlot", height = "500px")),
          
          div(class = "plot-card-small", style = "padding: 20px; margin-bottom: 20px;",
              h4("Violin Plot"),
              plotOutput("Signature_Vlnplot", height = "400px")),
          
          div(class = "plot-card-small", style = "padding: 20px;",
              h4("Bar Plot"),
              plotOutput("Signature_Barplot", height = "400px"))
        )
      })
      
      output$Signature_FeaturePlot <- renderPlot({
        FeaturePlot(sce_score,features = "Signature_score")
      })
      output$Signature_Barplot <- renderPlot({
        avg_scores <- sce_score@meta.data %>%
          group_by(celltype) %>%
          summarise(mean_score = mean(Signature_score, na.rm = TRUE),
                    sd_score = sd(Signature_score, na.rm = TRUE))
        
        ggplot(avg_scores, aes(x = celltype, y = mean_score, fill = celltype)) +
          geom_bar(stat = "identity", width = 0.7, color = "black") +
          geom_errorbar(aes(ymin = mean_score - sd_score, ymax = mean_score + sd_score),
                        width = 0.2, position = position_dodge(0.9)) +
          
          labs(x = "Cell Type", y = "Mean Signature Score") +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "none")  
        
      })
      output$Signature_Vlnplot <- renderPlot({
        ggviolin(sce_score@meta.data,
                 x="celltype",y="Signature_score",
                 width=0.8,color="black",
                 fill="celltype",
                 xlab=F,
                 add='mean_sd',
                 bxp.errorbar=T,
                 bxp.errorbar.width=0.05,
                 size=0.5,
                 palette="npg",
                 legend="right")
        
      })
      # output$signature_score_barplot <- renderPlotly({
      #   plot_ly(
      #     data = dataset,
      #     x = ~Condition,
      #     y = ~SignatureScore,
      #     type = "bar",
      #     color = ~Condition
      #   ) %>% layout(title = "Signature Score by Condition")
      # })
      # output$signature_score_violin <- renderPlotly({
      #   req(input$signature_compare_by)
      #   plot_ly(
      #     data = dataset,
      #     x = as.formula(paste0("~", input$signature_compare_by)),
      #     y = ~SignatureScore,
      #     type = "violin",
      #     color = as.formula(paste0("~", input$signature_compare_by))
      #   ) %>% layout(title = paste("Violin plot of Signature Score by", input$signature_compare_by))
      # })
    })
    
    # gene_list <- reactive({
    #   req(input$gene_file)
    #   ext <- tools::file_ext(input$gene_file$name)
    #   if (ext == "csv") {
    #     df <- read_csv(input$gene_file$datapath)
    #   } else if (ext == "tsv") {
    #     df <- read_tsv(input$gene_file$datapath)
    #   } else if (ext == "txt") {
    #     df <- read.table(input$gene_file$datapath, header = TRUE, sep = "\t")
    #   } else if (ext == "xlsx") {
    #     df <- read_excel(input$gene_file$datapath)
    #   } else {
    #     showNotification("不支持的文件格式", type = "error")
    #     return(NULL)
    #   }
    #   genes <- df[[1]]
    #   genes <- genes[!is.na(genes) & genes != ""]
    #   return(genes)
    # })
    # observeEvent(input$Upload_Signature,{
    #   req(input$gene_file)
    #   output$Signature_Score <- renderUI({
    #     tagList(
    #       fluidRow(
    #         column(12,
    #                plotOutput("Signature_UMAP")
    #         )
    #       ),
    #       tags$hr(),
    #       fluidRow(
    #         column(12,
    #                plotOutput("Signature_Barplot")
    #         )
    #       ),
    #       tags$hg(),
    #       fluidRow(
    #         column(12,
    #                plotOutput("Signature_Vlnplot")
    #         )
    #       )
    #     )
    #   }
    #   )
    #   output$gene_list_text <- renderPrint({
    #     gene_list()
    #   })
    #   genes <- gene_list()
    #   sce_score <- tryCatch({
    #     AddModuleScore(dataset,
    #                    features = list(c(genes)),
    #                    ctrl = min(100),  # 控制基因数不超过可用基因数
    #                    name = "Signature_score")
    #   }, error = function(e) {
    #     showNotification(paste("计算特征分数出错:", e$message), type = "error")
    #     return(NULL)
    #   })
    #   colnames(sce_score@meta.data)[colnames(sce_score@meta.data) == "Signature_score1"] <- "Signature_score"
    # })
    
  })
  ##Oral Diseases 后端--------------------
  output$isDataLoaded_tab2 <- reactive({
    isDataLoaded_tab2()
  })
  outputOptions(output, "isDataLoaded_tab2", suspendWhenHidden = FALSE)
  observeEvent(input$confirm_button_Oral_Disease, {
    #更新数据加载状态为TRUE
    isDataLoaded_tab2(TRUE)
    #gene name选项信息
    load(paste0("Dataset/Oral/",input$Oral_Diseases_Names,"/",input$Oral_Diseases_Celltype,"/gene_names.Rdata"))
    load(paste0("Dataset/Oral/",input$Oral_Diseases_Names,"/",input$Oral_Diseases_Celltype,"/metadata.Rdata"))
    updateSelectizeInput(session = session,inputId = "gene_input_Oral_Diseases",choices=c("",gene_names),server=TRUE)
    updateSelectizeInput(session = session,inputId = "DEG_Celltype_Oral_Diseases",choices=c("",unique(as.character(metadata$Celltype))),server=TRUE)
    updateSelectizeInput(session = session,inputId = "Condition1_Oral_Diseases",choices=c("",unique(metadata$patient)),server=TRUE)
    updateSelectizeInput(session = session,inputId = "Condition2_Oral_Diseases",choices=c("",unique(metadata$patient)),server=TRUE)
    ##Overview_Oral_Diseases---------
    #UMAP图
    output$UMAP_plot_Oral_Diseases <- renderUI({
      tags$img(src = paste0("Dataset/Oral/",input$Oral_Diseases_Names,"/",input$Oral_Diseases_Celltype,"/UMAP.svg"), 
               style = "width: auto; height: 450px; object-fit: contain;")
    })
    output$download_UMAP_Oral_Diseases <- downloadHandler(
      filename = function() {
        paste0("UMAP_", input$Oral_Diseases_Names, "_", input$Oral_Diseases_Celltype, ".pdf")
      },
      content = function(file) {
        path <- paste0("www/Dataset/Oral/", input$Oral_Diseases_Names, "/", input$Oral_Diseases_Celltype, "/UMAP.pdf")
        file.copy(path, file)
      }
    )
    #饼图
    output$Pie_chart_Oral_Diseases <- renderUI({
      div(
        tags$img(
          src = paste0("Dataset/Oral/",input$Oral_Diseases_Names,"/",input$Oral_Diseases_Celltype,"/Pieplot.svg"),
          style = "width: auto; height: 450px; object-fit: contain;"
        )
      )
    })
    output$download_Pie_chart_Oral_Diseases <- downloadHandler(
      filename = function() {
        paste0("Pie_chart_", input$Oral_Diseases_Names, "_", input$Oral_Diseases_Celltype, ".pdf")
      },
      content = function(file) {
        path <- paste0("www/Dataset/Oral/", input$Oral_Diseases_Names, "/", input$Oral_Diseases_Celltype, "/Pieplot.pdf")
        file.copy(path, file)
      }
    )
    #Dotplot
    output$Dot_plot_Oral_Diseases <- renderUI({
      tags$img(src = paste0("Dataset/Oral/",input$Oral_Diseases_Names,"/",input$Oral_Diseases_Celltype,"/Dotplot.svg"), 
               style = "width: 1000px; height: auto; object-fit: contain;")
    })
    output$download_Dot_plot_Oral_Diseases <- downloadHandler(
      filename = function() {
        paste0("Dot_plot_", input$Oral_Diseases_Names, "_", input$Oral_Diseases_Celltype, ".pdf")
      },
      content = function(file) {
        path <- paste0("www/Dataset/Oral/", input$Oral_Diseases_Names, "/", input$Oral_Diseases_Celltype, "/Dotplot.pdf")
        file.copy(path, file)
      }
    )
    #柱状图：按条件分组
    output$Barplot_conditions_Oral_Diseases <- renderUI({
      tags$img(src = paste0("Dataset/Oral/",input$Oral_Diseases_Names,"/",input$Oral_Diseases_Celltype,"/Barplot_Condition.svg"), 
               style = "height: 400px;flex-shrink: 0;")
    })
    output$download_Barplot_conditions_Oral_Diseases <- downloadHandler(
      filename = function() {
        paste0("Barplot_conditions_", input$Oral_Diseases_Names, "_", input$Oral_Diseases_Celltype, ".pdf")
      },
      content = function(file) {
        path <- paste0("www/Dataset/Oral/", input$Oral_Diseases_Names, "/", input$Oral_Diseases_Celltype, "/Barplot_Condition.pdf")
        file.copy(path, file)
      }
    )
    #柱状图：按患者分组
    output$Barplot_patients_Oral_Diseases <- renderUI({
      tags$img(src = paste0("Dataset/Oral/",input$Oral_Diseases_Names,"/",input$Oral_Diseases_Celltype,"/Barplot_Patient.svg"), 
               style = "height: 400px;flex-shrink: 0;")
    })
    output$download_Barplot_patients_Oral_Diseases <- downloadHandler(
      filename = function() {
        paste0("Barplot_patients_", input$Oral_Diseases_Names, "_", input$Oral_Diseases_Celltype, ".pdf")
      },
      content = function(file) {
        path <- paste0("www/Dataset/Oral/", input$Oral_Diseases_Names, "/", input$Oral_Diseases_Celltype, "/Barplot_Patient.pdf")
        file.copy(path, file)
      }
    )
    ##Gene_Oral_Diseases---------
    show_download_Oral_Diseases <- reactiveVal(FALSE)
    observeEvent(input$Submit_Gene_Oral_Diseases,{
      if(input$gene_input_Oral_Diseases != ""){
        req(input$gene_input_Oral_Diseases != "")
        file_path = paste0("Dataset/Oral/",input$Oral_Diseases_Names,"/",input$Oral_Diseases_Celltype,"/seurat_plot.Rdata")
        load(file_path)
        dataset <- seurat_plot
      }
      output$gene_expression_plot_Oral_Diseases <- renderPlot({
        req(input$gene_input_Oral_Diseases != "")
        if((input$gene_input_Oral_Diseases != "")&(input$gene_input_Oral_Diseases %in% rownames(dataset))){
          p <- FeaturePlot(
            dataset, features = input$gene_input_Oral_Diseases, order = TRUE, 
            cols = c("lightgrey", "#ff0000"), pt.size = 0.5, combine = TRUE) + 
            theme_bw() + 
            theme(
              panel.border = element_rect(fill = NA, color = "black", size = 1.2), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_blank(), 
              plot.title = element_text(hjust = 0.5, size = 20, family = "Arial", face = "italic"),
              axis.title = element_text(size = 17, family = "Arial", color = "black"),
              axis.text = element_text(size = 15, family = "Arial", color = "black"),
              legend.position = "right",
              legend.text = element_text(size = 15, family = "Arial"),
              legend.title = element_text(size = 17, family = "Arial")
            ) + 
            xlab("UMAP_1") + ylab("UMAP_2") +
            scale_color_gradientn(colors = BlueAndRed())
          p
        }else{}
      })
      output$violin_plot_Oral_Diseases <- renderPlot({
        req(input$gene_input_Oral_Diseases != "", input$gene_compare_Oral_Diseases != "")
        if((input$gene_input_Oral_Diseases!="") & (input$gene_input_Oral_Diseases %in% rownames(dataset))){
          theme_clean_pub <- function() { 
            theme_classic(base_size = 18) +  
              theme(   
                panel.background = element_rect(fill = "white", color = NA),   
                plot.background = element_rect(fill = "white", color = NA),   
                panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
                axis.line = element_blank(),
                axis.ticks = element_line(color = "black", linewidth = 0.8),
                axis.ticks.length = unit(0.15, "cm"),   
                axis.text = element_text(color = "black", size = 14, family = "Arial"),
                axis.title = element_text(color = "black", size = 16, family = "Arial"), 
                plot.title = element_text(size = 18, face = "bold.italic", hjust = 0.5, family = "Arial"), 
                plot.subtitle = element_text(size = 12, hjust = 0.5, family = "Arial"),   
                strip.background = element_blank(), 
                strip.text = element_text(size = 14, face = "bold", family = "Arial"),   
                panel.grid = element_blank(),
                legend.position = "right",
                legend.title = element_text(size = 12, face = "bold", family = "Arial"),   
                legend.text = element_text(size = 11, family = "Arial"),   
                legend.background = element_blank(),
                legend.key = element_blank()
              )
          }
          if(input$gene_compare_Oral_Diseases == "Celltype"){
            genes_to_extract <- input$gene_input_Oral_Diseases
            group_col <- input$gene_compare_Oral_Diseases
            plot_data <- FetchData(
              object = dataset, 
              vars = c(genes_to_extract, group_col), 
              slot = "data"                         
            )
            p <- ggplot(plot_data, aes(x = .data[[group_col]], y = .data[[genes_to_extract]], fill = .data[[group_col]])) +
              geom_violin(trim = FALSE, 
                          color = NA, 
                          alpha = 0.8, 
                          width = 0.9) +
              geom_boxplot(width = 0.1, 
                           linewidth = 0.8, 
                           color = "black", 
                           outlier.shape = NA, 
                           alpha = 1) +
              stat_summary(fun = mean, 
                           geom = "point", 
                           shape = 23,  size =4,  
                           fill ="white",  color ="black",  stroke =1.4 ) +
              scale_fill_manual(values = colors) + 
              theme_bw()+
              labs(  title ="",  x ="",  y =paste0("The expression levels of ",genes_to_extract)) +
              theme_clean_pub() + theme(legend.position ="none")
            p
          }else if(input$gene_compare_Oral_Diseases == "Tissue"){
            genes_to_extract <- input$gene_input_Oral_Diseases
            group_x <- "Celltype"  
            group_col <- input$gene_compare_Oral_Diseases
            plot_data <- FetchData(
              object = dataset, 
              vars = c(genes_to_extract, group_x, group_col), 
              slot = "data"
            )
            p <- ggplot(plot_data, aes(x = .data[[group_x]], 
                                       y = .data[[genes_to_extract]], 
                                       fill = .data[[group_col]])) +
              geom_violin(trim = FALSE, 
                          color = "black", 
                          linewidth = 0.3,
                          alpha = 0.7, 
                          width = 0.8,
                          position = position_dodge(width = 0.9)) + 
              geom_boxplot(
                aes(group = interaction(.data[[group_x]], .data[[group_col]])),
                width = 0.15, 
                linewidth = 0.5, 
                color = "black", 
                outlier.shape = NA, 
                fill = "white",
                position = position_dodge(width = 0.9)) + 
              stat_summary(aes(group = .data[[group_col]]),
                           fun = mean, 
                           geom = "point", 
                           shape = 23, size = 1.5,  
                           fill = "white", color = "black",
                           position = position_dodge(width = 0.9)) +
              
              scale_fill_manual(values = colors) + 
              labs(x = "", y = paste0("Expression levels of ", genes_to_extract), fill = "Tissue") +
              theme_clean_pub() + 
              theme(
                legend.position = "right",
                axis.text.x = element_text(angle = 45, hjust = 1)
              )
            p
          }
        }
      }, width = function() {
        if(input$gene_compare_Oral_Diseases == "Celltype"){
          n_groups <- length(unique(dataset@meta.data[[input$gene_compare_Oral_Diseases]]))
          calculated_width <- n_groups * 80 + 200
        }else if(input$gene_compare_Oral_Diseases == "Tissue"){
          n_groups <- length(unique(dataset@meta.data[[input$gene_compare_Oral_Diseases]]))
          m_groups <- length(unique(dataset@meta.data[['Celltype']]))
          calculated_width <- n_groups * m_groups * 80 + 200
        }
        return(max(600, calculated_width))
      }
      )
      show_download_Oral_Diseases(TRUE)
    })
    output$download_btn_Gene_FeaturePlot_Oral_Diseases <- renderUI({
      if (show_download_Oral_Diseases()) {
        downloadButton(
          "download_Gene_FeaturePlot_Oral_Diseases", 
          "Download PDF", 
          class = "btn-download-custom"
        )
      } else {
        return(NULL) 
      }
    })
    output$download_Gene_FeaturePlot_Oral_Diseases <- downloadHandler(
      filename = function() {
        paste0("FeaturePlot_", input$gene_input_Oral_Diseases, "_", input$Oral_Diseases_Names, "_", input$Oral_Diseases_Celltype, ".pdf")
      },
      content = function(file) {
        req(input$gene_input_Oral_Diseases != "")
        file_path = paste0("Dataset/Oral/",input$Oral_Diseases_Names,"/",input$Oral_Diseases_Celltype,"/seurat.Rdata")
        load(file_path)
        dataset <- seurat
        gene_data <- FetchData(dataset, vars = input$gene_input_Oral_Diseases, slot = "data")[,1]
        min_value <- quantile(gene_data, 0.8) 
        max_value <- max(gene_data)
        p <- FeaturePlot(
          dataset, features = input$gene_input_Oral_Diseases, order = TRUE, 
          cols = c("lightgrey", "#ff0000"), pt.size = 0.5, combine = TRUE) + 
          theme_bw() + 
          theme(
            panel.border = element_rect(fill = NA, color = "black", size = 1.2), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_blank(), 
            plot.title = element_text(hjust = 0.5, size = 20,  face = "italic", family = "Arial"),
            axis.title = element_text(size = 17,  color = "black", family = "Arial"),
            axis.text = element_text(size = 15,  color = "black", family = "Arial"),
            legend.position = "right",
            legend.text = element_text(size = 15, family = "Arial"),
            legend.title = element_text(size = 17, family = "Arial")
          ) + 
          xlab("UMAP_1") + ylab("UMAP_2") +
          scale_color_gradientn(colors = BlueAndRed(),
                                limits = c(min_value, max_value),       
                                na.value = "lightgray")
        ggsave(file, plot = p, device = "pdf", width = 5.5, height = 5, units = "in")
      }
    )
    output$download_btn_Gene_VlnPlot_Oral_Diseases <- renderUI({
      if (show_download_Oral_Diseases()) {
        downloadButton(
          "download_Gene_VlnPlot_Oral_Diseases", 
          "Download PDF", 
          class = "btn-download-custom"
        )
      } else {
        return(NULL) 
      }
    })
    output$download_Gene_VlnPlot_Oral_Diseases <- downloadHandler(
      filename = function() {
        paste0("VlnPlot_", input$gene_input_Oral_Diseases, "_", input$Oral_Diseases_Names, "_", input$Oral_Diseases_Celltype, ".pdf")
      },
      content = function(file) {
        req(input$gene_input_Oral_Diseases != "")
        file_path = paste0("Dataset/Oral/",input$Oral_Diseases_Names,"/",input$Oral_Diseases_Celltype,"/seurat.Rdata")
        load(file_path)
        dataset <- seurat
        theme_clean_pub <- function() { 
          theme_classic(base_size = 18) +  
            theme(   
              panel.background = element_rect(fill = "white", color = NA),   
              plot.background = element_rect(fill = "white", color = NA),   
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
              axis.line = element_blank(),
              axis.ticks = element_line(color = "black", linewidth = 0.8),
              axis.ticks.length = unit(0.15, "cm"),   
              axis.text = element_text(color = "black", size = 14, family = "Arial"),
              axis.title = element_text(color = "black", size = 16, family = "Arial"), 
              plot.title = element_text(size = 18, face = "bold.italic", hjust = 0.5, family = "Arial"), 
              plot.subtitle = element_text(size = 12, hjust = 0.5, family = "Arial"),   
              strip.background = element_blank(), 
              strip.text = element_text(size = 14, face = "bold", family = "Arial"),   
              panel.grid = element_blank(),
              legend.position = "right",
              legend.title = element_text(size = 12, face = "bold", family = "Arial"),   
              legend.text = element_text(size = 11, family = "Arial"),   
              legend.background = element_blank(),
              legend.key = element_blank()
            )
        }
        if(input$gene_compare_Oral_Diseases == "Celltype"){
          genes_to_extract <- input$gene_input_Oral_Diseases
          group_col <- input$gene_compare_Oral_Diseases
          plot_data <- FetchData(
            object = dataset, 
            vars = c(genes_to_extract, group_col), 
            slot = "data"                         
          )
          p <- ggplot(plot_data, aes(x = .data[[group_col]], y = .data[[genes_to_extract]], fill = .data[[group_col]])) +
            geom_violin(trim = FALSE, 
                        color = NA, 
                        alpha = 0.8, 
                        width = 0.9) +
            geom_boxplot(width = 0.1, 
                         linewidth = 0.8, 
                         color = "black", 
                         outlier.shape = NA, 
                         alpha = 1) +
            stat_summary(fun = mean, 
                         geom = "point", 
                         shape = 23,  size =4,  
                         fill ="white",  color ="black",  stroke =1.4 ) +
            scale_fill_manual(values = colors) + 
            theme_bw()+
            labs(  title ="",  x ="",  y =paste0("The expression levels of ",genes_to_extract)) +
            theme_clean_pub() + theme(legend.position ="none")
        }else if(input$gene_compare_Oral_Diseases == "Tissue"){
          genes_to_extract <- input$gene_input_Oral_Diseases
          group_x <- "Celltype"  
          group_col <- input$gene_compare_Oral_Diseases
          plot_data <- FetchData(
            object = dataset, 
            vars = c(genes_to_extract, group_x, group_col), 
            slot = "data"
          )
          p <- ggplot(plot_data, aes(x = .data[[group_x]], 
                                     y = .data[[genes_to_extract]], 
                                     fill = .data[[group_col]])) +
            geom_violin(trim = FALSE, 
                        color = "black", 
                        linewidth = 0.3,
                        alpha = 0.7, 
                        width = 0.8,
                        position = position_dodge(width = 0.9)) + 
            geom_boxplot(
              aes(group = interaction(.data[[group_x]], .data[[group_col]])),
              width = 0.15, 
              linewidth = 0.5, 
              color = "black", 
              outlier.shape = NA, 
              fill = "white",
              position = position_dodge(width = 0.9)) + 
            stat_summary(aes(group = .data[[group_col]]),
                         fun = mean, 
                         geom = "point", 
                         shape = 23, size = 1.5,  
                         fill = "white", color = "black",
                         position = position_dodge(width = 0.9)) +
            
            scale_fill_manual(values = colors) + 
            labs(x = "", y = paste0("Expression levels of ", genes_to_extract), fill = "Tissue") +
            theme_clean_pub() + 
            theme(
              legend.position = "right",
              axis.text.x = element_text(angle = 45, hjust = 1)
            )
        }
        if(input$gene_compare_Oral_Diseases == "Celltype"){
          n_groups <- length(unique(dataset@meta.data[[input$gene_compare_Oral_Diseases]]))
          calculated_width <- n_groups * 80 + 200
        }else if(input$gene_compare_Oral_Diseases == "Tissue"){
          n_groups <- length(unique(dataset@meta.data[[input$gene_compare_Oral_Diseases]]))
          m_groups <- length(unique(dataset@meta.data[['Celltype']]))
          calculated_width <- n_groups * m_groups * 80 + 200
        }
        ggsave(file, plot = p, device = "pdf", width = calculated_width/72, height = 450/72, units = "in")
      }
    )
    
    ##DEGs Panel_Oral_Diseases---------
    observeEvent(input$Compare_DEGs_Oral_Diseases,{
      req(input$Condition1_Oral_Diseases != "", input$Condition2_Oral_Diseases != "")
      if (input$Condition1_Oral_Diseases == input$Condition2_Oral_Diseases) {
        showNotification(
          "请选择两个不同的条件进行比较！", 
          type = "warning",
          duration = 3
        )
        return()
      }
      ##DEGs table------
      output$DEG_table_Oral_Diseases <- renderDT({
        req()
        filename = paste0("Dataset/Oral/",input$Oral_Diseases_Names,"/",input$Oral_Diseases_Celltype,"/DEGs/",input$DEG_Celltype_Oral_Diseases,"/",input$DEG_Celltype_Oral_Diseases,"_",input$Condition1_Oral_Diseases,"_vs_",input$Condition2_Oral_Diseases,".Rdata")
        load(filename)
        #DEGs <- FindMarkers(dataset,group.by = "Tissue",ident.1 = input$Condition1,ident.2 = input$Condition2,min.pct=0.25,logfc.threshold = 0.25)
        #DEGs <- DEGs[which(abs(DEGs$avg_log2FC) >= 0.5 & DEGs$p_val<0.05),]
        #DEGs$gene <- rownames(DEGs)
        #DEGs <- DEGs[,c("gene","avg_log2FC","p_val","p_val_adj","pct.1","pct.2")]
        #DEGs[, sapply(DEGs, is.numeric)] <- round(DEGs[, sapply(DEGs, is.numeric)], 2)
        DEGs
      },escape = FALSE,
      rownames = FALSE,
      options = list(
        paging = TRUE,
        searching = TRUE,
        info = TRUE,
        lengthMenu = NULL,
        columnDefs = list(
          list(width = '80px',targets = c(0,1,2,3,4,5))
        )
      )
      )
      ##Volcano---------
      output$DEG_Vol_Oral_Diseases <- renderPlot({
        req()
        filename = paste0("Dataset/Oral/",input$Oral_Diseases_Names,"/",input$Oral_Diseases_Celltype,"/DEGs/",input$DEG_Celltype_Oral_Diseases,"/",input$DEG_Celltype_Oral_Diseases,"_",input$Condition1_Oral_Diseases,"_vs_",input$Condition2_Oral_Diseases,".Rdata")
        load(filename)
        DEGs$p_val_adj <- as.numeric(DEGs$p_val_adj)
        DEGs$avg_log2FC <- as.numeric(DEGs$avg_log2FC)
        DEGs <- DEGs[order(DEGs$p_val_adj,DEGs$avg_log2FC,decreasing = c(FALSE,TRUE)),]
        DEGs_up <- DEGs[which(DEGs$avg_log2FC >= 1 & DEGs$p_val_adj<0.05),]
        DEGs_down <- DEGs[which(DEGs$avg_log2FC <= -1 & DEGs$p_val_adj<0.05),]
        DEGs_up <- DEGs_up[order(DEGs_up$avg_log2FC,decreasing = TRUE),]
        DEGs_down <- DEGs_down[order(DEGs_down$avg_log2FC,decreasing = FALSE),]
        Upvals <- rownames(DEGs_up)[1:5]
        Downvals <- rownames(DEGs_down)[1:5]
        vals <- c(Upvals,Downvals)
        group<-ifelse(
          DEGs$avg_log2FC<(-0.5)&DEGs$p_val_adj<0.05,'#4D4398',
          ifelse(DEGs$avg_log2FC>(0.5)&DEGs$p_val_adj<0.05,'indianred1',
                 '#b5b5b5'))
        max <- max(abs(na.omit(DEGs[,c("avg_log2FC")])))
        x_max <- max*(11/10)
        x_min <- -x_max
        y_max <- -log10(min(DEGs$p_val_adj[DEGs$p_val_adj > 0]))+1
        group[is.na(group)]<-'#b5b5b5'
        names(group)[group=='indianred1']<-'Up'
        names(group)[group=='#b5b5b5']<-'No_sig'
        names(group)[group=='#4D4398']<-'Down'
        if(nrow(DEGs)>0){
          g <- EnhancedVolcano(DEGs,
                               x = "avg_log2FC",
                               y = "p_val_adj",
                               lab = rownames(DEGs),
                               pCutoff = 0.05,
                               FCcutoff = 0.1,
                               pointSize=c(ifelse(rownames(DEGs) %in% vals,1.5,0.5)),
                               #pointSize = 3,
                               labSize = 5,
                               xlim = c(x_min,x_max),
                               ylim = c(0,y_max),
                               colCustom = group,
                               title = NULL,
                               subtitle = NULL,
                               caption = NULL,
                               legendPosition = "right",
                               selectLab = c(Upvals,Downvals),
                               xlab = bquote(~Log[2]~'fold change'),
                               legendLabSize = 12,
                               legendIconSize = 6,
                               labCol = 'black',
                               drawConnectors = TRUE,typeConnectors = 'closed',lengthConnectors = unit(0, "mm"),
                               widthConnectors = 0.5,
                               max.overlaps = 30
                               
          )+
            theme(
              panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
              plot.title = element_text(
                size = 15,
                face = "bold",
                hjust = 0.5,    
                vjust = 2,    
                margin = margin(b = 10) 
              ),
              axis.title = element_text(size = 14),                   
              legend.title = element_blank(),                                       
              text = element_text() 
            )
          g
        }else{}
      },width = 560,height = 400)
    })
    ##Gene_Signature_Oral_Diseases--------
    parsed_genes_Oral_Diseases <- reactive({
      req(input$gene_file_Oral_Diseases)
      file <- input$gene_file_Oral_Diseases
      ext <- tools::file_ext(file$name)
      df <- switch(ext,
                   "csv"  = read.csv(file$datapath, header = TRUE),
                   "tsv"  = read.table(file$datapath, header = TRUE, sep = "\t"),
                   "txt"  = read.table(file$datapath, header = TRUE, sep = "\t"),
                   "xlsx" = readxl::read_excel(file$datapath),
                   { showNotification("不支持格式", type = "error"); return(NULL) }
      )
      genes <- as.character(df[[1]])
      genes <- unique(trimws(genes))
      genes <- genes[genes != "" & !is.na(genes)]
      return(genes)
    })
    
    output$gene_list_text_Oral_Diseases <- renderText({
      g <- parsed_genes_Oral_Diseases()
      req(g)
      #将基因向量合并成用逗号分隔的字符串
      paste(g, collapse = ", ")
    })
    
    observeEvent(input$Upload_Signature_Oral_Diseases, {
      genes <- parsed_genes_Oral_Diseases()
      valid_genes <- intersect(genes, rownames(dataset))
      if(length(valid_genes) == 0) {
        showNotification("上传的基因均不在当前数据集中，请检查符号。", type = "error")
        return(NULL)
      }
      
      #计算分数(AddModuleScore)
      withProgress(message = 'Calculating Signature Score...', {
        sce_score <- AddModuleScore(
          object = dataset,
          features = list(valid_genes),
          name = "Signature_score"
        )
        colnames(sce_score@meta.data)[colnames(sce_score@meta.data) == "Signature_score1"] <- "Signature_score"
      })
      
      output$Signature_Score_Oral_Diseases <- renderUI({
        tagList(
          div(class = "plot-card-small", style = "padding: 20px; margin-bottom: 20px;",
              h4("Feature Plot"),
              plotOutput("Signature_FeaturePlot_Oral_Diseases", height = "500px")),
          
          div(class = "plot-card-small", style = "padding: 20px; margin-bottom: 20px;",
              h4("Violin Plot"),
              plotOutput("Signature_Vlnplot_Oral_Diseases", height = "400px")),
          
          div(class = "plot-card-small", style = "padding: 20px;",
              h4("Bar Plot"),
              plotOutput("Signature_Barplot_Oral_Diseases", height = "400px"))
        )
      })
      
      output$Signature_FeaturePlot_Oral_Diseases <- renderPlot({
        FeaturePlot(sce_score,features = "Signature_score")
      })
      output$Signature_Barplot_Oral_Diseases <- renderPlot({
        avg_scores <- sce_score@meta.data %>%
          group_by(celltype) %>%
          summarise(mean_score = mean(Signature_score, na.rm = TRUE),
                    sd_score = sd(Signature_score, na.rm = TRUE))
        
        ggplot(avg_scores, aes(x = celltype, y = mean_score, fill = celltype)) +
          geom_bar(stat = "identity", width = 0.7, color = "black") +
          geom_errorbar(aes(ymin = mean_score - sd_score, ymax = mean_score + sd_score),
                        width = 0.2, position = position_dodge(0.9)) +
          
          labs(x = "Cell Type", y = "Mean Signature Score") +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "none")  
        
      })
      output$Signature_Vlnplot_Oral_Diseases <- renderPlot({
        ggviolin(sce_score@meta.data,
                 x="celltype",y="Signature_score",
                 width=0.8,color="black",
                 fill="celltype",
                 xlab=F,
                 add='mean_sd',
                 bxp.errorbar=T,
                 bxp.errorbar.width=0.05,
                 size=0.5,
                 palette="npg",
                 legend="right")
        
      })
    })
    
    
    
  })
  
  ##Pancancer 后端-----------
  output$isDataLoaded_tab3 <- reactive({
    isDataLoaded_tab3()
  })
  outputOptions(output, "isDataLoaded_tab3", suspendWhenHidden = FALSE)
  observeEvent(input$confirm_button_Pancancer, {
    #更新数据加载状态为TRUE
    if(input$Pancancer_Celltype == "CD4T"){
      gene_names <- read.table(file = "/home/cuihao/scTExplorer/Dataset/Pan_cancer/CD4_T_genes.txt")
      gene_names <- gene_names$V1
    }else if(input$Pancancer_Celltype == "CD8T"){
      gene_names <- read.table(file = "/home/cuihao/scTExplorer/Dataset/Pan_cancer/CD8_T_genes.txt")
      gene_names <- gene_names$V1
    }
    isDataLoaded_tab3(TRUE)
    #gene name选项信息
    updateSelectizeInput(session = session,inputId = "gene_input_Pancancer",choices=c("",gene_names),server=TRUE)
    ##Overview_Oral_Diseases---------
    #UMAP图
    output$UMAP_plot_Pancancer <- renderUI({
      tags$img(src = paste0("Dataset/Pan_cancer/","UMAP.png"), 
               style = "width: auto; height: 450px; object-fit: contain;")
    })
    output$download_UMAP_Pancancer <- downloadHandler(
      filename = function() {
        paste0("UMAP_", input$Pancancer_Names, "_", input$Pancancer_Celltype, ".png")
      },
      content = function(file) {
        path <- paste0("www/Dataset/Pan_cancer/", "/UMAP.png")
        file.copy(path, file)
      }
    )
    #Dotplot
    output$Dot_plot_Pancancer <- renderUI({
      tags$img(src = paste0("Dataset/Pan_cancer/","/Dotplot.png"), 
               style = "width: 1000px; height: auto; object-fit: contain;")
    })
    output$download_Dot_plot_Pancancer <- downloadHandler(
      filename = function() {
        paste0("Dot_plot_", input$Pancancer_Names, "_", input$Pancancer_Celltype, ".png")
      },
      content = function(file) {
        path <- paste0("www/Dataset/Pan_cancer/", "/Dotplot.png")
        file.copy(path, file)
      }
    )
    #柱状图：按条件分组
    output$Barplot_conditions_Pancancer <- renderUI({
      tags$img(src = paste0("Dataset/Pan_cancer/","/Barplot_Condition.png"), 
               style = "height: 800px;flex-shrink: 0;")
    })
    output$download_Barplot_conditions_Pancancer <- downloadHandler(
      filename = function() {
        paste0("Barplot_conditions_", input$Pancancer_Names, "_", input$Pancancer_Celltype, ".png")
      },
      content = function(file) {
        path <- paste0("www/Dataset/Pan_cancer/", "/Barplot_Condition.png")
        file.copy(path, file)
      }
    )
    #柱状图：按患者分组
    output$Barplot_patients_Pancancer <- renderUI({
      tags$img(src = paste0("Dataset/Pan_cancer/","/Barplot_Patient.png"), 
               style = "height: 800px;flex-shrink: 0;")
    })
    output$download_Barplot_patients_Pancancer <- downloadHandler(
      filename = function() {
        paste0("Barplot_patients_", input$Pancancer_Names, "_", input$Pancancer_Celltype, ".png")
      },
      content = function(file) {
        path <- paste0("www/Dataset/Pan_cancer/", "/Barplot_Patient.png")
        file.copy(path, file)
      }
    )
    ##Gene_Oral_Diseases---------
    show_download_Pancancer <- reactiveVal(FALSE)
    observeEvent(input$Submit_Gene_Pancancer,{
      output$gene_expression_plot_Pancancer <- renderUI({
        tags$img(src = paste0("Dataset/Pan_cancer/",input$Pancancer_Celltype,"/Gene/Featureplot/PNG/",input$gene_input_Pancancer,"_Feature_Plot.png"), 
                 style = "width: auto; height: 450px; object-fit: contain;")
      })
      output$violin_plot_Pancancer <- renderUI({
        if(input$gene_compare_Pancancer == "Celltype"){
          type = "cell_type"
        }else if(input$gene_compare_Pancancer == "Cancer_type"){
          type = "cancer_type"
        }
        tags$img(src = paste0("Dataset/Pan_cancer/",input$Pancancer_Celltype,"/Gene/Vlnplot/",type,"/PNG/",input$gene_input_Pancancer,"_Vln_Plot.png"), 
                 style = "width: auto; height: 450px; object-fit: contain;")
      })

      show_download_Pancancer(TRUE)
    })
    output$download_btn_Gene_FeaturePlot_Pancancer <- renderUI({
      if (show_download_Pancancer()) {
        downloadButton(
          "download_Gene_FeaturePlot_Pancancer",
          "Download PDF",
          class = "btn-download-custom"
        )
      } else {
        return(NULL)
      }
    })
    output$download_Gene_FeaturePlot_Pancancer <- downloadHandler(
      filename = function() {
        paste0("FeaturePlot_", input$gene_input_Pancancer, "_", input$Pancancer_Names, "_", input$Pancancer_Celltype, ".pdf")
      },
      content = function(file) {
        req(input$gene_input_Pancancer != "")
        src = paste0("www/Dataset/Pan_cancer/",input$Pancancer_Celltype,"/Gene/Featureplot/PDF/",input$gene_input_Pancancer,"_Feature_Plot.pdf")
        if (!file.exists(src)) {
          stop("PDF file not found: ", src)
        }
        file.copy(src, file, overwrite = TRUE)
        }
    )
    output$download_btn_Gene_VlnPlot_Pancancer <- renderUI({
      if (show_download_Pancancer()) {
        downloadButton(
          "download_Gene_VlnPlot_Pancancer",
          "Download PDF",
          class = "btn-download-custom"
        )
      } else {
        return(NULL)
      }
    })
    output$download_Gene_VlnPlot_Pancancer <- downloadHandler(
      filename = function() {
        paste0("VlnPlot_", input$gene_input_Pancancer, "_", input$Pancancer_Names, "_", input$Pancancer_Celltype, ".pdf")
      },
      content = function(file) {
        req(input$gene_input_Pancancer != "")
        if(input$gene_compare_Pancancer == "Celltype"){
          type = "cell_type"
        }else if(input$gene_compare_Pancancer == "Cancer_type"){
          type = "cancer_type"
        }
        src = paste0("www/Dataset/Pan_cancer/",input$Pancancer_Celltype,"/Gene/Vlnplot/",type,"/PDF/",input$gene_input_Pancancer,"_Vln_Plot.pdf")
        if (!file.exists(src)) {
          stop("PDF file not found: ", src)
        }
        file.copy(src, file, overwrite = TRUE)
      }
    )
    
  })
  
}

shinyApp(ui = ui, server = server)
