library(shiny)
library(tidyverse)
library(pheatmap)
library(plotly)
library(readxl)
library(shinyWidgets)
library(ggrepel)
library(shinycssloaders)
library(matrixStats)
library(DT)
library(conflicted)

conflict_prefer("renderDataTable", "DT")
conflict_prefer("filter", "dplyr")
conflict_prefer("arrange", "dplyr")

shinyServer(function(input, output, session){
  
  datasets <- list.dirs('./Data', full.names = F, recursive = F)
  output$data_select_ui <- renderUI(selectInput('data_select', 'Select dataset:', choices = datasets, selected = datasets[1]))
  
  de_data <- reactive({
    req(input$data_select)
    fname <- list.files(str_c('./Data/', input$data_select), pattern = 'fpkm', full.names = T)
    d <- read.csv(fname)
    cnames <- colnames(d)
    cnames <- sapply(cnames, function(x){
      if(str_sub(x, 1, 1) == 'X'){
        x <- str_sub(x, 2, -1)
      }
      return(x)
    })
    colnames(d) <- cnames
    d
  })
  
  de_summary <- reactive({
    req(input$data_select)
    fname <- list.files(str_c('./Data/', input$data_select), pattern = 'summary', full.names = T)
    read.csv(fname)
  })
  
  comps <- reactive({
    comps <- de_summary()$Comparison
    comps
  })
  output$comp_select_ui <- renderUI(selectInput(inputId = 'comp_select', label = 'Select comparison:', choices = comps(), selected = comps()[1]))
  
  #DE data table
  de_tab_data <- reactive({
    req(input$data_select, input$comp_select)
    dtb <- de_data()
    dtb <- cbind(dtb[, c('gene', 'gene_name', 'description', 'gene_type')],
                 dtb[, str_detect(colnames(dtb), str_c(input$comp_select, '\\.'))])
    dtb$description <- str_sub(dtb$description, 1, (str_locate(dtb$description, '\\[') -2)[, 1])
    colnames(dtb) <- c('Gene', 'Gene name', 'Description', 'Gene type', 'Log2FC', 'P-value', 'P-adj', 'Raw Log2FC')
    dtb
  })
  rtd <- reactive({
    req(de_tab_data())
    de_tab_data() %>% na.omit %>% arrange(`P-adj`)
  })
  output$de_tab <- renderDataTable(rtd())
  
  output$downloadDEResults <- downloadHandler(
    filename = function() {
      str_c(input$comp_select, '_DE_results.xlsx')
    },
    content = function(file) {
      WriteXLS::WriteXLS(x = rtd(), ExcelFileName = file, row.names = F, col.names = T, AdjWidth = T, FreezeRow = 1)
    }
  )
  
  #MA plot
  ma_plot <- reactive({
    req(input$ma_cutoff, input$ma_limits)
    samples <- c(
      de_summary() %>% filter(Comparison == input$comp_select) %>% pull(Sample_names_in_base_level_condition) %>% str_split(',') %>% unlist %>% str_c('_fpkm'),
      de_summary() %>% filter(Comparison == input$comp_select) %>% pull(Sample_names_in_comparison_level_condition) %>% str_split(',') %>% unlist %>% str_c('_fpkm')
    )
    avgs <- de_data()[, samples] %>% rowMeans()
    ma_data <- de_tab_data()
    ma_data$avg_fpkm <- avgs
    ma_cutoff <- input$ma_cutoff %>% as.numeric
    ma_data <- ma_data %>% na.omit
    plot_ma_data <- data.frame(`M` = ma_data$avg_fpkm, `A` = ma_data$Log2FC, `isde` = ma_data$`P-adj` <= ma_cutoff)
    geneplotter::plotMA(plot_ma_data, ylim = c(input$ma_limits[1], input$ma_limits[2]))
  })
  output$ma_plot <- renderPlot(ma_plot())
                                
  #Volcano plot
  vplot_xlim <- reactive({
    req(input$vulcano_cutoff)
    sliderInput(inputId = 'vulcano_xlim', label = 'Log-fold change limits:', step = 0.5,
                value = c(min(rtd()$Log2FC) %>% floor, max(rtd()$Log2FC) %>% ceiling), min = -8, max = 8)
                #value = c(min(de_tab_data()$Log2FC), max(de_tab_data()$Log2FC)),
                #min = 1.1*min(de_tab_data()$Log2FC), max = 1.1*max(de_tab_data()$Log2FC))
  })
  output$vplot_xlim <- renderUI({vplot_xlim()})
  
  vplot_ylim <- reactive({
    max_value <- rtd() %>% filter(`P-adj` > 0) %>% pull(`P-adj`) %>% min
    max_value <- -log10(max_value)/100
    max_value <- 100*ceiling(max_value)
    #max_value <- 100*ceiling(max(-log10(de_tab_data() %>% filter(`P-ajd` > 0) %>% pull(`P-adj`)))/100)
    textInput(inputId = 'vulcano_ylim', value = max_value, label = 'y-limit:')
  })
  output$vplot_ylim <- renderUI({vplot_ylim()})
  
  vplot <- reactive({
    req(input$vulcano_cutoff, input$vulcano_xlim, input$vulcano_ylim)
    p_cutoff <- as.numeric(input$vulcano_cutoff)
    vplot_reg <- rep('Not sig.', nrow(de_tab_data()))
    vplot_reg[de_tab_data()$Log2FC > 0 & de_tab_data()$`P-adj` <= p_cutoff] <- 'Up reg.'
    vplot_reg[de_tab_data()$Log2FC < 0 & de_tab_data()$`P-adj` <= p_cutoff] <- 'Down reg.'
    vdata <- de_tab_data()
    vdata$regulation <- vplot_reg
    p <- ggplot(vdata, aes(x = Log2FC, y = -log10(`P-adj`), color = regulation, Gene = `Gene name`)) + 
      scale_color_manual(values = c('Not sig.' = 'gray', 'Up reg.' = 'darkred', 'Down reg.' = 'navy')) +
      geom_point() + theme_bw() + ylim(0, 1.5*max(-log10(vdata$`P-adj`))) +
      xlim(input$vulcano_xlim[1], input$vulcano_xlim[2]) + ylim(0, (input$vulcano_ylim %>% as.numeric))
    ggplotly(p)
  })
  output$vplot <- renderPlotly(vplot())
  
  #Expression plot
  genelist <- reactive({
    req(input$comp_select)
    de_tab_data() %>% na.omit %>% pull(`Gene name`) %>% sort
  })
  output$exp_gene_selector <- renderUI({
    multiInput(inputId = 'exp_genes', label = 'Select genes:', choices = genelist(), options = list(enable_search = T))
  })
  
  expplot <- reactive({
    #req(input$exp_genes)
    #dtb <- de_tab_data() %>% na.omit
    exp_cutoff <- input$exp_cutoff %>% as.numeric
    base_samp <- de_summary() %>% filter(Comparison == input$comp_select) %>% pull(Sample_names_in_base_level_condition) %>% str_split(',') %>% unlist %>% str_c('_fpkm')
    comp_samp <- de_summary() %>% filter(Comparison == input$comp_select) %>% pull(Sample_names_in_comparison_level_condition) %>% str_split(',') %>% unlist %>% str_c('_fpkm')
    base_avgs <- de_data()[, base_samp] %>% rowMeans()
    comp_avgs <- de_data()[, comp_samp] %>% rowMeans()
    base_cond <- de_summary() %>% filter(Comparison == input$comp_select) %>% pull(Base_level_condition)
    base_cond <- str_c('Average ', base_cond, ' fpkm')
    comp_cond <- de_summary() %>% filter(Comparison == input$comp_select) %>% pull(Comparison_level_condition)
    comp_cond <- str_c('Average ', comp_cond, ' fpkm')
    exp_data <- data.frame(Gene = de_tab_data()$`Gene name`,
                           base_avgs = base_avgs,
                           comp_avgs = comp_avgs,
                           signif = de_tab_data()$`P-adj` <= exp_cutoff
                           ) %>% na.omit()
    #exp_data <- exp_data[!is.na(exp_data$signif), ]
    #exp_data$signif <- exp_data$signif <= exp_cutoff
    
    label_data <- exp_data[exp_data$Gene %in% input$exp_genes, ]
    #exp_data$Gene[!(exp_data$Gene %in% input$exp_genes)]
    #exp_data <- exp_data[!(is.na(exp_data$signif)), ]
    
    p <- ggplot(exp_data, aes(x = base_avgs, y = comp_avgs, colour = signif, Gene = `Gene`)) + geom_point() +
      scale_color_manual(values = c('TRUE' = 'darkred', 'FALSE' = 'gray')) + theme_bw() + 
      xlab(base_cond) + ylab(comp_cond) + scale_x_log10() + scale_y_log10() + 
      theme(legend.position = "none", text = element_text(size = 16)) + geom_label_repel(data = label_data, aes(x = base_avgs, y = comp_avgs, label = Gene))+
      geom_point(data = label_data, aes(x = base_avgs, y = comp_avgs, fill = signif), colour = 'black', shape = 21)
    p
  })
  output$expplotly <- renderPlotly(ggplotly(expplot()))
  output$expplot <- renderPlot(expplot())
  
  
  ### PCA plot
  clust_samples <- reactive({
    de_summary() %>% filter(Comparison == input$comp_select) %>% select_at(c(7, 10)) %>% as.character() %>% str_split(',') %>% unlist
  })
  output$clust_samples <- renderUI(multiInput(inputId = 'clust_samples', label = 'Select samples:',
                                              choices = clust_samples(), selected = clust_samples()))
  
  clust_plot <- reactive({
    req(clust_samples(), input$clust_samples)
    df <- de_summary() %>% filter(Comparison == input$comp_select)
    samples <- clust_samples()
    conds <- factor(c(rep(df %>% pull(5), df %>% pull(6)), rep(df %>% pull(8), df %>% pull(9))),
                    levels = c(df %>% pull(5), df %>% pull(8)))
    names(conds) <- samples
    
    # filtering just selected samples and their conditions
    samples <- input$clust_samples
    conds <- conds[input$clust_samples]
    
    clust_df <- de_data()[, c(str_c(samples, '_fpkm'))] %>% as.matrix # initial df with comparison samples
    colnames(clust_df) <- samples
    # n most variant genes
    clust_df <- clust_df[order(rowVars(clust_df[, -1]), decreasing = T), ][1:input$clust_ngenes, ]
    clust_pca <- prcomp(t(clust_df), scale = T) # do PCA on samples
    pca_df <- cbind(samples, Condition = conds, clust_pca$x[, 1:2] %>% as.data.frame)
    x_var <- str_c('PC1: ', round(100*summary(clust_pca)$importance[2, 1], 2), '% variance explained')
    y_var <- str_c('PC2: ', round(100*summary(clust_pca)$importance[2, 2], 2), '% variance explained')
    
    ggplot(pca_df, aes(x = PC1, y = PC2, colour = Condition)) +
      geom_point() + geom_text_repel(aes(label = samples)) + theme_bw() + xlab(x_var) + ylab(y_var) +
      theme(text = element_text(size = 16))
  })
  
  output$pca_plot <- renderPlot(clust_plot())
  
  #GO
  
  go_table <- reactive({
    species <- list.dirs(str_c('Data/', input$data_select, '/enrichment_tests'), recursive = F, full.names = F)
    go_file = str_c('Data/', input$data_select, '/enrichment_tests/', species, '/', input$comp_select, '/',
                    input$comp_select, input$go_regulation, '_', 'go', '_',
                    input$go_category, '.csv'
                    )
    if(file.exists(go_file)){
      go_tab <- read.csv(go_file)
    }else{
      go_tab <- data.frame()
    }
    go_tab
  })
  
  
  
  output$go_table <- renderDataTable(datatable(go_table()[, -7], selection = 'single'),
                                     server = T) # Renders the table without gene list
  
  output$downloadGO <- downloadHandler(
    filename = function() {
      if(input$go_regulation == ''){
        go.reg <- 'all'
      }else{
        go.reg <- c('.up' = 'up', '.down' = 'down')[input$go_regulation]
      }
      
      str_c(input$comp_select, '_', input$go_category, '_', go.reg, '.xlsx')
    },
    content = function(file) {
      WriteXLS::WriteXLS(x = go_table(), ExcelFileName = file, row.names = F, col.names = T, AdjWidth = T, FreezeRow = 1)
    }
  )
  
  go_plot <- reactive({
    go_genes <- go_table()[input$go_table_rows_selected, 8]
    go_genes <- str_split(go_genes, ', ') %>% unlist
    #req(input$exp_genes)
    #dtb <- de_tab_data() %>% na.omit
    exp_cutoff <- input$exp_cutoff %>% as.numeric
    base_samp <- de_summary() %>% filter(Comparison == input$comp_select) %>% pull(Sample_names_in_base_level_condition) %>% str_split(',') %>% unlist %>% str_c('_fpkm')
    comp_samp <- de_summary() %>% filter(Comparison == input$comp_select) %>% pull(Sample_names_in_comparison_level_condition) %>% str_split(',') %>% unlist %>% str_c('_fpkm')
    base_avgs <- de_data()[, base_samp] %>% rowMeans()
    comp_avgs <- de_data()[, comp_samp] %>% rowMeans()
    base_cond <- de_summary() %>% filter(Comparison == input$comp_select) %>% pull(Base_level_condition)
    base_cond <- str_c('Average ', base_cond, ' fpkm')
    comp_cond <- de_summary() %>% filter(Comparison == input$comp_select) %>% pull(Comparison_level_condition)
    comp_cond <- str_c('Average ', comp_cond, ' fpkm')
    exp_data <- data.frame(Gene = de_tab_data()$`Gene name`,
                           base_avgs = base_avgs,
                           comp_avgs = comp_avgs,
                           signif = de_tab_data()$`P-adj` <= exp_cutoff
    ) %>% na.omit()

    label_data <- exp_data[exp_data$Gene %in% go_genes, ]
    
    ggplot(exp_data, aes(x = base_avgs, y = comp_avgs, colour = signif)) + geom_point() +
      scale_color_manual(values = c('TRUE' = 'darkred', 'FALSE' = 'gray')) + theme_bw() + 
      xlab(base_cond) + ylab(comp_cond) + scale_x_log10() + scale_y_log10() +
      theme(legend.position = "none", text = element_text(size = 16)) + geom_label_repel(data = label_data, aes(x = base_avgs, y = comp_avgs, label = Gene))+
      geom_point(data = label_data, aes(x = base_avgs, y = comp_avgs, fill = signif), colour = 'black', shape = 21) +
      ggtitle(go_table()[input$go_table_rows_selected, 2])
  })
  
  output$go_plot <- renderPlot(go_plot())
  
  # GSA
  
  gsa_table <- reactive({
    species <- list.dirs(str_c('Data/', input$data_select, '/gene_set_tests'), recursive = F, full.names = F)
    gsa_file <- str_c('Data/', input$data_select, '/gene_set_tests/', species, '/', input$comp_select, '/',
                      input$comp_select, '-', input$gsa_cat, '_sets.csv')
    gsa_tab <- read.csv(gsa_file)
    gsa_tab
  })
  
  output$gsa_table <- renderDataTable(gsa_table())
  
  output$downloadGSA <- downloadHandler(
    filename = function() {
      str_c(input$comp_select, '_', input$gsa_cat, '.xlsx')
    },
    content = function(file) {
      WriteXLS::WriteXLS(x = gsa_table(), ExcelFileName = file, row.names = F, col.names = T, AdjWidth = T, FreezeRow = 1)
    }
  )
  
})


