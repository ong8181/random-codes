####
#### Shiny app to convert Quant-it data
#### 2020.8.28 v1 Ushio
#### 2020.9.1 v2 Ushio
####

# Load library
library(shiny); packageVersion("shiny") # 1.5.0, 2020.8.28
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.8.29

# Load settings
well_position <- read.csv("data/WellPosition.csv")
well_conc <- read.csv("data/WellPosi_Convert.txt",stringsAsFactors = F, header = T)


# Define UI for app that draws a histogram ----
ui <- fluidPage(
  # Sidebar settings
  sidebarLayout(
    sidebarPanel(
      # Input Quant-it data
      fileInput("file", "1. Quant-it データを選択", accept = c(".csv")),
      
      # Select STD column
      checkboxGroupInput("std_columns", 
                         "2. 標準 DNA/RNA の列を選択", 
                         choices = list("Column 1" = 2, 
                                        "Column 2" = 3, 
                                        "Column 3" = 4, 
                                        "Column 4" = 5, 
                                        "Column 5" = 6, 
                                        "Column 6" = 7, 
                                        "Column 7" = 8, 
                                        "Column 8" = 9, 
                                        "Column 9" = 10, 
                                        "Column 10" = 11, 
                                        "Column 11" = 12, 
                                        "Column 12" = 13),
                         selected = NULL),
      
      # Specify STD concentrations
      textInput("std_dna_conc", "3. 標準 DNA/RNA の濃度を入力 (ng/µl)", 
                value = "0,50,100,200,400,600,800,1000"),
      
      # Specify volume of standard
      numericInput("std_vol", "4. 標準 DNA/RNA 溶液の量を入力", 10),
      
      # Show standard curve
      actionButton("show_standard_curve", "5. 検量線をチェック"),
      
      # Specify volume of samples
      numericInput("sample_vol", "6. サンプル溶液の量を入力", 5),
      
      # Show standard curve
      actionButton("convert_data1", "7. ウェルあたりの DNA 量を計算 (ng/well)"),
      
      # Show standard curve
      actionButton("convert_data2", "8. DNA 濃度を計算 (ng/µl in sample)"),
      
      # Specify NC wells
      textInput("nc_wells", "[オプション] ネガコンの濃度を入力 (Ex. 1.0, 2.1, ...)",  value = "0"),
      
      # Download results
      downloadButton('data1', label = "DNA/RNA 濃度データ"),
      downloadButton('data2', label = "[オプション] ネガコンを差し引いたデータ"),
      
      # Input data for bioTEC
      fileInput("input_for_biotec", "9. DNA/RNA 濃度データを選択", accept = c(".csv")),
      numericInput("target_conc", "10. 揃えたい DNA/RNA 濃度", 0.5),
      
      # Show standard curve
      actionButton("generate_biotec", "11. BioTEC 用シートを生成"),
      downloadButton('data3', label = "BioTEC シート"),
    ),
    
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Quant-it データ", tableOutput('table')),
                  tabPanel("標準 DNA/RNA データ", tableOutput('std_data')),
                  tabPanel("検量線", plotOutput('std_plot')),
                  tabPanel("ウェルあたりの DNA/RNA 量 (ng/well)", tableOutput('amount_data')),
                  tabPanel("DNA/RNA 濃度", tableOutput('not_corrected_data')),
                  tabPanel("DNA/RNA 濃度 (ネガコン補正後)", tableOutput('corrected_data')),
                  tabPanel("BioTEC 用シート", tableOutput('biotec_sheet'))
      )
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  csv_file <- reactive(read.csv(input$file$datapath))
  std_vals <- reactive(as.numeric(unlist(strsplit(input$std_dna_conc,","))))
  std_concs <- reactive(data.frame(std = std_vals()))
  nc_wells <- reactive(as.numeric(unlist(strsplit(input$nc_wells,","))))
  
  # Show original data
  observeEvent(input$file, {
    output$table <- renderTable(csv_file())
  })
  
  # Show standard data
  observeEvent(input$file, {
    output$std_data <- renderTable(cbind(std_concs(),
                                         csv_file()[,as.numeric(input$std_columns)]))
  })
  
  # Show standard curve
  observeEvent(input$show_standard_curve, {
    std_df <- cbind(std_concs(), csv_file()[,as.numeric(input$std_columns)])
    d_long <- pivot_longer(std_df, cols = -std)
    
    output$std_plot = renderPlot({
      ggplot(d_long, aes(x = value, y = std)) +
        geom_point() + geom_smooth(method = lm, formula = y ~ x + I(x^2) + I(x^3) + 0) +
        xlab("Quant-it values") + ylab("DNA concentration")
    })
  })
  
  observeEvent(input$convert_data1, {
    d <- read.csv(input$file$datapath)
    std_col <- as.numeric(input$std_columns)
    std_df <- cbind(as.numeric(input$std_vol)*std_concs(), d[,std_col])
    d_long <- pivot_longer(std_df, cols = -std)
    
    # Perform regression
    lm_model <- summary(lm(std ~ value + I(value^2) + I(value^3) + 0, data = d_long))
    coeffs <- coefficients(lm_model)[,1]
    d_samples <- d[,-c(1,std_col)]
    d_amount <- coeffs[1] * d_samples + coeffs[2] * d_samples^2 + coeffs[3] * d_samples^3
    
    # Substitude negative controls
    d_amount[d_amount < 0] <- 0
    
    d_amount <- cbind(d[,1], d_amount)
    colnames(d_amount)[1] <- "X"
    
    output$amount_data = renderTable({
      d_amount
    })
  })
  
  
  observeEvent(input$convert_data2, {
    d <- read.csv(input$file$datapath)
    std_col <- as.numeric(input$std_columns)
    std_df <- cbind(as.numeric(input$std_vol)*std_concs(), d[,std_col])
    d_long <- pivot_longer(std_df, cols = -std)
    
    # Perform regression
    lm_model <- summary(lm(std ~ value + I(value^2) + I(value^3) + 0, data = d_long))
    coeffs <- coefficients(lm_model)[,1]
    d_samples <- d[,-c(1,std_col)]
    d_amount <- coeffs[1] * d_samples + coeffs[2] * d_samples^2 + coeffs[3] * d_samples^3
    
    # Substitude negative controls
    d_amount[d_amount < 0] <- 0
    d_concs <- d_amount/as.numeric(input$sample_vol)
    
    # Substitude negative controls
    d_concs_corrected <- d_concs - mean(nc_wells())
    d_concs_corrected[d_concs_corrected < 0] <- 0
    #d_concs_corrected <- d_amount_corrected/as.numeric(input$sample_vol)
    
    d_concs <- cbind(d[,1], d_concs)
    d_concs_corrected <- cbind(d[,1], d_concs_corrected)
    colnames(d_concs)[1] <- colnames(d_concs_corrected)[1] <- "X"
    
    output$not_corrected_data = renderTable({
      d_concs
    })
    
    output$corrected_data = renderTable({
      d_concs_corrected
    })
    
    # Downloadable csv of selected dataset ----
    output$data1 <- downloadHandler(
      filename = function() {
        "not_corrected_dna_conc.csv"
      },
      content = function(file) {
        write.csv(d_concs, file, row.names = FALSE)
      }
    )
    
    output$data2 <- downloadHandler(
      filename = function() {
        "corrected_dna_conc.csv"
      },
      content = function(file) {
        write.csv(d_concs_corrected, file, row.names = FALSE)
      }
    )
    
  })
  
  # Generate biotec data sheet
  observeEvent(input$generate_biotec, {
    dna <- read.csv(input$input_for_biotec$datapath, row.names = 1)
    target_conc <- input$target_conc
    total_vol <- 24 # Volume of diluted sample
    dispense_vol <- 0
    transfer_manner <- 1 # 1,2,3,4から選択
    dilution_H2O_vol <- 450
    
    # Convert "dna" file as an acceptable input
    well_position$row2 <- substr(well_position$well_name, 1, 1)
    biotec_input <- well_position[,c("row2", "col")]
    
    ## Assign values
    biotec_input$dna <- NA
    for(i in 1:nrow(biotec_input)){
      if(!is.null(dna[biotec_input$row2[i],biotec_input$col[i]])){
        biotec_input$dna[i] <- dna[biotec_input$row2[i],biotec_input$col[i]]
      }
    }
    
    biotec_input <- na.omit(biotec_input)
    colnames(biotec_input) <- c("RowPos", "ColPos", "TECAN_conc")
    
    # Convert biotec_input for TECAN
    tmp <- matrix(0, nrow = nrow(biotec_input), ncol = 3)
    colnames(tmp) <- c("pos","conc","total_vol")
    tmp[,"pos"] <- sprintf("%s%s",biotec_input$RowPos,biotec_input$ColPos)
    tmp[,"conc"] <- biotec_input$TECAN_conc
    tmp[,"total_vol"] = rep(total_vol, nrow(biotec_input))
    
    # >>>> Ushio added, 2020.8.24 >>>>
    #tmp[is.na(tmp[,"conc"]),"conc"] <- 0
    #tmp[,"conc"] <- as.numeric(tmp[,"conc"])
    #tmp[tmp[,"conc"] < target_conc,"conc"] <- 0.5
    # <<<< Ushio added <<<<
    
    biotec_input <- tmp
    biotec_output <- matrix("",ncol = 12, nrow = nrow(biotec_input))
    
    # Format biotec_output file
    biotec_output[,1] <- biotec_input[,"pos"]
    biotec_output[,2] <- total_vol * target_conc/as.numeric(biotec_input[,"conc"])
    biotec_output[,2][as.numeric(biotec_output[,2]) > total_vol] <- total_vol
    biotec_output[,3] <- biotec_output[,1] # Copy column 1 to column 3
    biotec_output[,4] <- biotec_output[,2] # Copy column 2 to column 4
    
    biotec_output[,12] <- rep(total_vol, nrow(biotec_output)) - as.numeric(biotec_output[,4])
    # Pipetting volume
    biotec_output[,5] = (as.numeric(biotec_output[,2]) + as.numeric(biotec_output[,12]))*0.8
    # 
    biotec_output[,6] <- rep(3, nrow(biotec_output))
    # 
    biotec_output[,7] <- biotec_output[,1]
    # 
    biotec_output[,8] <- rep(dispense_vol, nrow(biotec_output))
    biotec_output[,9] <- as.vector(well_conc[sprintf(biotec_output[,1]),transfer_manner])
    biotec_output[,10] <- biotec_output[,8]
    used_water <- 0
    dilution_tube_pos <- c("A1","A2","A3","A4","A5","A6","A7","A8")
    dilution_tube_pos_num <- 1
    
    # Format biotec_output file 2
    for(i in 1:nrow(biotec_output)){
      used_water = used_water + as.numeric(biotec_output[i,12])
      biotec_output[i,11] = dilution_tube_pos[dilution_tube_pos_num]
      if(used_water > dilution_H2O_vol - 80){
        used_water = 0
        dilution_tube_pos_num = dilution_tube_pos_num + 1        
      }
    }
    for(i in 1:nrow(biotec_output)){
      if(biotec_output[i,8]==0){
        biotec_output[i,6] = ""
        biotec_output[i,7] = ""
        biotec_output[i,8] = ""
      }
    }
    for(i in 1:nrow(biotec_output)){
      if(biotec_output[i,10]==0){
        biotec_output[i,9] = ""
        biotec_output[i,10] = ""
      }
    }
    for(i in 1:nrow(biotec_output)){
      if(biotec_output[i,12]==0){
        biotec_output[i,11] = ""
        biotec_output[i,12] = ""
      }
    }
    
    # Generate biotec_output file
    biotec_output_colnames <- c("Well","Vol(ul)","Well","Vol(ul)","Mix Vol(ul)","Stage","Well","Vol(ul)","Well","Vol(ul)","Well","Vol(ul)")
    #biotec_output2 <- data.frame(biotec_output)
    #colnames(biotec_output2) <- biotec_output_colnames
    biotec_output2 <- as.data.frame(rbind(biotec_output_colnames, biotec_output))
    
    output$biotec_sheet = renderTable({
      biotec_output2
    })
    
    # Downloadable csv of selected dataset ----
    output$data3 <- downloadHandler(
      filename = function() {
        "boitec_sheet.csv"
      },
      content = function(file) {
        write.csv(biotec_output2, file, row.names = FALSE)
      }
    )
  })
  
  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)

## Publish App
#library(rsconnect)
#rsconnect::deployApp('/Users/ushio/Work/Company/Clockmics/20200825_TECAN_Shiny/Quantit_Calculator_v2')
