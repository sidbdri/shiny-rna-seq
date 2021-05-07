library(shiny)
library(tidyverse)
library(pheatmap)
library(plotly)
library(readxl)
library(shinycssloaders)
library(DT)

shinyUI(
  fluidPage(
    fluidRow(
      uiOutput('data_select_ui'),
      uiOutput('comp_select_ui'),
      height = '150%'
    ),
    fillRow(
      tabsetPanel(type = 'tabs', id = 'resTabs',
                  tabPanel('DEA Results', br(), 
                           fluidRow(
                             dataTableOutput('de_tab') %>% withSpinner(color="grey")
                           )
                          ),
                  tabPanel('Plots',
                           tabsetPanel(type = 'tabs', id = 'plotTabs',
                                       tabPanel('MA plot', br(),
                                                textInput(inputId = 'ma_cutoff', label = 'P-value cutoff:', value = 0.05),
                                                plotOutput('ma_plot') %>% withSpinner(color="grey")),
                                       tabPanel('Volcano plot', br(),
                                                textInput(inputId = 'vulcano_cutoff', label = 'P-value cutoff:', value = 0.05),
                                                plotOutput('vplot') %>% withSpinner(color="grey")
                                                ),
                                       tabPanel('Expression plot', br(),
                                                sidebarLayout(
                                                  sidebarPanel(
                                                    textInput(inputId = 'exp_cutoff', label = 'P-value cutoff:', value = 0.05),
                                                    uiOutput('exp_gene_selector') %>% withSpinner(color="grey")
                                                  ),
                                                  mainPanel(
                                                    plotOutput('expplot') %>% withSpinner(color="grey")
                                                  )
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
                           fluidRow(
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
                           fluidRow(
                             dataTableOutput(outputId = 'go_table')
                           )
                          ),
                  tabPanel('GSA', br(),
                           radioButtons(inputId = 'gsa_cat', label = 'Category:',
                                        choices = c('Cell type' = 'CELL_TYPE',
                                                    'Curated' = 'CURATED',
                                                    'GO' = 'GO',
                                                    'Motif' = 'MOTIF'),
                                        inline = T),
                           dataTableOutput('gsa_table'))
                  )
    )
  )
)