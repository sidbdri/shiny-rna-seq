library(shiny)
library(tidyverse)
library(pheatmap)
library(plotly)
library(readxl)
library(shinycssloaders)
library(DT)
library(conflicted)
conflict_prefer("dataTableOutput", "DT")

shinyUI(
  fluidPage(
    tags$style(HTML("
      #first {
          border: 16px double white;
      }
      #second {
          border: 2px dashed blue;
      }
    ")),
    fluidRow(id = 'first',
      uiOutput('data_select_ui'),
      uiOutput('comp_select_ui'),
      height = '150%'
    ),
    fillRow(
      tabsetPanel(type = 'tabs', id = 'resTabs',
                  tabPanel('DEA Results', br(), 
                           fluidRow(id = 'first',
                             downloadButton('downloadDEResults', 'Download .xlsx')  %>% withSpinner(color="grey"), br(), 
                             dataTableOutput('de_tab') %>% withSpinner(color="grey")
                           )
                          ),
                  tabPanel('Plots',
                           tabsetPanel(type = 'tabs', id = 'plotTabs',
                                       tabPanel('MA plot', br(),
                                                fluidRow(
                                                  column(width = 3, textInput(inputId = 'ma_cutoff', label = 'P-value cutoff:', value = 0.05)),
                                                  column(width = 6, sliderInput(inputId = 'ma_limits', label = 'Log-fold change limits:',
                                                                     value = c(-1, 1), min = -3, max = 3, step = 0.1))
                                                ),
                                                fluidRow(id = 'first',
                                                  plotOutput('ma_plot') %>% withSpinner(color="grey")
                                                  )
                                                ),
                                       tabPanel('Volcano plot', br(),
                                                fluidRow(
                                                  column(width = 2, textInput(inputId = 'vulcano_cutoff', label = 'P-value cutoff:', value = 0.05)),
                                                  column(width = 4, uiOutput(outputId = 'vplot_xlim')),
                                                  column(width = 2, uiOutput(outputId = 'vplot_ylim'))
                                                ),
                                                fluidRow(id = 'first',
                                                  plotlyOutput('vplot') %>% withSpinner(color="grey")
                                                )
                                              ),
                                       tabPanel('Expression plot', br(),
                                                sidebarLayout(
                                                  sidebarPanel(
                                                    textInput(inputId = 'exp_cutoff', label = 'P-value cutoff:', value = 0.05),
                                                    uiOutput('exp_gene_selector') %>% withSpinner(color="grey")
                                                  ),
                                                  mainPanel(
                                                    #plotOutput('expplot') %>% withSpinner(color="grey"),
                                                    plotlyOutput('expplotly') %>% withSpinner(color="grey")
                                                  )
                                                )
                                                
                                              ),
                                       tabPanel('Gene plot', br(),
                                                fluidRow(id = 'first',
                                                  column(width = 4, uiOutput(outputId = 'gene_plot_selector') %>% withSpinner(color="grey")),
                                                  column(width = 2, radioButtons(inputId = 'gene_plot_value', label = 'Select value:', choices = c('Counts' = 'counts', 'FPKM' = 'fpkm'))),
                                                  column(width = 2, radioButtons(inputId = 'gene_plot_type', label = 'Select plot type:', choices = c('Scatter plot' = 'scatter', 'Box plot' = 'box', 'Violin plot' = 'violin')))
                                                ),
                                                fluidRow(id = 'first',
                                                         #tableOutput(outputId = 'gene_plot')
                                                         plotlyOutput(outputId = 'gene_plot')
                                                )
                                       )
                                       )
                           ),
                  tabPanel('Sample PCA', br(),
                           sidebarLayout(
                             sidebarPanel(
                               textInput(inputId = 'clust_ngenes', label = 'Number of genes:', value = 500),
                               uiOutput('clust_samples')
                             ),
                             mainPanel(
                               plotOutput(outputId = 'pca_plot') %>% withSpinner(color="grey")
                             )
                           )),
                  tabPanel('GO', br(),
                           fluidRow(id = 'first',
                             column(width = 4,
                               radioButtons(inputId = 'go_category', label = 'GO Category:',
                                            choices = c('BP' = 'bp', 'CC' = 'cc', 'MF' = 'mf'),
                                            selected = 'bp', inline = T),
                               radioButtons(inputId = 'go_regulation', label = 'Gene regulation:',
                                            choices = c('All DE' = '', 'Upregulated' = '.up', 'Downregulated' = '.down'),
                                            selected = '', inline = T)
                             ),
                             column(width = 6,
                                    plotOutput('go_plot')) %>% withSpinner(color="grey")
                           ),
                           fluidRow(id = 'first',
                             downloadButton('downloadGO', 'Download .xlsx') %>% withSpinner(color="grey"), br(),
                             dataTableOutput(outputId = 'go_table') %>% withSpinner(color="grey")
                           )
                          ),
                  tabPanel('GSA', br(),
                           uiOutput('gsa_cat'),
                           downloadButton('downloadGSA', 'Download .xlsx') %>% withSpinner(color="grey"), br(),
                           dataTableOutput('gsa_table'))
                  )
    )
  )
)