library(shiny)
library(ggplot2)
library(Seurat)
library(ComplexHeatmap)
library(RColorBrewer)
library(plotly)
library(FNN)
library(InteractiveComplexHeatmap)

load("/media/hd/pimenta2/R/tirosh_gbm_spatial_Data/metadata_tirosh_allGBM_score_cell_type.rda")
load("/media/hd/pimenta2/R/tirosh_gbm_spatial_Data/score_all_pathways_tirosh_all_gbm.rda")
pathwaylist <- colnames(score_all_pathways_tirosh_all_gbm)
pathwaylist <- pathwaylist[order(pathwaylist)]

color <- c(
  "AC" = "#E41A1C",
  "Chromatin.Reg" = "#984EA3",
  "Inflammatory.Mac" = "#4DAF4A",
  "low.quality" = "#999999",
  "Mac" = "#A65628",
  "MES" = "#FF7F00",
  "MES.Ast" = "#FFD700",
  "MES.Hyp" = "#D95F02",
  "Neuron" = "#006666",
  "NPC" = "#999900",
  "OC.NPC1" = "#1F78B4",
  "Oligo" = "#A6CEE3",
  "OPC" = "#FFFF99",
  "Prolif.Metab" = "#FF69B4",
  "Reactive.Ast" = "#008080",
  "Vasc" = "#4B0082"
)

my_shiny_function_all_samples_color_pval <- function(metadata, patient_column, list_of_scores_columns,
                                                     cell_type_column, proximity_split = "", k=5){
  metadata <- metadata[!is.na(metadata[,cell_type_column]),]
  metadata[,cell_type_column] <- as.character(metadata[,cell_type_column])
  if(proximity_split != ""){
    metadata[,proximity_split] <- as.character(metadata[,proximity_split])
  }

  library(grid)
  plot_heatmap_pvalue_spatial_score <- function(output_function){
    aux <- output_function[[1]][["adj_p_value"]]
    aux <- as.matrix(aux)
    p <- Heatmap(as.matrix(output_function[[1]][["normalized_score"]]), cluster_rows = FALSE, cluster_columns = FALSE, column_title = names(output_function)[[1]],
                 cell_fun = function(j, i, x, y, width, height, fill, p_val_mat = aux) {
                   # Check the value from text_matrix to decide the color
                   text_color <- if (p_val_mat[i, j] < 0.05) {
                     "green"
                   } else {
                     "black"
                   }
                   
                   # Use text_matrix[i, j] for the label and the determined color for the text
                   grid.text(
                     sprintf("%.3f", p_val_mat[i, j]), 
                     x, y, 
                     gp = gpar(fontsize = 10, col = text_color)
                   )
                 })
    p
  }
  
  proximity_metabolic_score_pvalue_perm_metab <- function(meta_from_seurat, k = 5,
                                                          sep = "GSM", cell_type_column,
                                                          list_metabolic, n_permutations = 1000,
                                                          seed = 42) {
    if(k<3){
      print("select k > 3")
    }
    # --- Initial Setup ---
    # Set a seed for random number generation to ensure reproducibility
    if (!is.null(seed)) {
      set.seed(seed)
    }
    list_permutations <- list()
    for(i in 1:length(list_metabolic)){
      list_permutations[[i]] <- list()
    }
    names(list_permutations) <- list_metabolic
    meta_from_seurat <- as.data.frame(meta_from_seurat)
    meta_from_seurat <- meta_from_seurat[!(is.na(meta_from_seurat[,cell_type_column])),]
    
    meta_from_seurat$cells_row_names <- rownames(meta_from_seurat)
    meta_from_seurat <- meta_from_seurat[,c(sep, cell_type_column, list_metabolic, "x", "y", "cells_row_names")]
    cell_types <- unique(meta_from_seurat[,cell_type_column])
    base_matrix <- matrix(0,nrow = length(cell_types), ncol = length(cell_types))
    rownames(base_matrix) <- cell_types
    colnames(base_matrix) <- cell_types
    
    
    # ============================== KNN in metadata =========================
    
    #first puttin the knn information back on metadata, to calculate just one time
    list_metas <- split(meta_from_seurat, meta_from_seurat[,sep])
    for(i in 1:length(list_metas)){
      # Extract spatial data and cell type information from the Seurat object
      spatial_data <- as.data.frame(list_metas[[i]][, c('x', 'y')])
      
      # Find k-nearest neighbors for each cell
      nearest_neighbors <- get.knnx(
        data = spatial_data[, c("x", "y")],
        query = spatial_data[, c("x", "y")],
        k = k + 1  # Include the cell itself
      )
      
      # Extract neighbor indices (excluding the cell itself)
      neighbor_indices <- nearest_neighbors$nn.index
      #butting back the knn information
      list_metas[[i]] <- cbind.data.frame(neighbor_indices, list_metas[[i]])
      for(rows in 1:nrow(list_metas[[i]])){
        for(cols in 1:(k+1)){
          list_metas[[i]][rows,cols] <- list_metas[[i]][as.numeric(list_metas[[i]][rows,cols]),cell_type_column]
        }
      }
    }
    metadata_knn_info <- list_metas[[1]]
    if(length(list_metas)>1){
      for(i in 2:length(list_metas)){
        metadata_knn_info <- rbind.data.frame(metadata_knn_info, list_metas[[i]])
      }
    }
    #now i have a metadata with all the info i need to proceed with the spatial score
    aux <- metadata_knn_info
    metadata_knn_info <- aux[,-c(3:(k+1))]
    colnames(metadata_knn_info) <- c("reference_cell_type", "interaction_cell_type", colnames(aux)[-c(1:(k+1))])
    for(new_df in 3:k){
      aux2 <- aux[,-(c(2:(k+1))[-which(c(2:(k+1)) == new_df)])]
      colnames(aux2) <- c("reference_cell_type", "interaction_cell_type", colnames(aux)[-c(1:(k+1))])
      metadata_knn_info <- rbind.data.frame(metadata_knn_info, aux2)
    }
    
    
    #===================== original value part ===================
    list_original_score <- split(metadata_knn_info, metadata_knn_info$reference_cell_type)
    for(cell_types in 1:length(list_original_score)){
      list_original_score[[cell_types]] <- split(list_original_score[[cell_types]], list_original_score[[cell_types]][,"interaction_cell_type"])
    }
    
    #metabolic list
    list_metabolic_original <- list()
    for(i in 1:length(list_metabolic)){
      list_metabolic_original[[i]] <- list()
    }
    names(list_metabolic_original) <- list_metabolic
    for(metab in 1:length(list_metabolic)){
      original_score <- base_matrix
      for(reference in 1:length(list_original_score)){
        for(interaction in 1:length(list_original_score[[reference]])){
          original_score[names(list_original_score)[reference],
                         names(list_original_score[[reference]])[interaction]] <- 
            sum((list_original_score[[reference]][[interaction]][,which(colnames(list_original_score[[reference]][[interaction]]) == list_metabolic[metab])]))/nrow(list_original_score[[reference]][[interaction]])
          list_metabolic_original[[metab]] <- original_score
        }
      }
    }
    
    
    # ========================= permutation score part ==============================
    #getting matrix each simulation
    list_score_permutation <- split(metadata_knn_info, metadata_knn_info$reference_cell_type)
    for(p in 1:n_permutations){
      list_score_permutation_p <- list_score_permutation
      for(cell_types in 1:length(list_score_permutation_p)){
        for(metabs in 1:length(list_metabolic)){
          list_score_permutation_p[[cell_types]][, list_metabolic[metabs]] <- sample(list_score_permutation_p[[cell_types]][, list_metabolic[metabs]])
        }
        list_score_permutation_p[[cell_types]] <- split(list_score_permutation_p[[cell_types]], list_score_permutation_p[[cell_types]][,"interaction_cell_type"])
      }
      
      list_metabolic_permutation_p <- list()
      for(i in 1:length(list_metabolic)){
        list_metabolic_permutation_p[[i]] <- list()
      }
      names(list_metabolic_permutation_p) <- list_metabolic
      for(metab in 1:length(list_metabolic)){
        original_score <- base_matrix
        for(reference in 1:length(list_score_permutation_p)){
          for(interaction in 1:length(list_score_permutation_p[[reference]])){
            original_score[names(list_score_permutation_p)[reference],
                           names(list_score_permutation_p[[reference]])[interaction]] <- 
              sum((list_score_permutation_p[[reference]][[interaction]][,which(colnames(list_score_permutation_p[[reference]][[interaction]]) == list_metabolic[metab])]))/nrow(list_score_permutation_p[[reference]][[interaction]])
            list_metabolic_permutation_p[[metab]] <- original_score
          }
        }
      }
      
      for(j in 1:length(list_metabolic)){
        list_permutations[[j]][[p]] <- list_metabolic_permutation_p[[j]]
      }
    }
    
    #=================== pvalue =======================================
    pvalue_list <- list()
    for(i in 1:length(list_metabolic)){
      pvalue_list[[i]] <- list_metabolic_original[[1]]*0
    }
    names(pvalue_list) <- list_metabolic
    
    for(metabs in 1:length(list_metabolic)){
      for(row_name in rownames(list_metabolic_original[[metabs]])){
        for(col_names in colnames(list_metabolic_original[[metabs]])){
          real_score1 <- list_metabolic_original[[metabs]][row_name, col_names] 
          permutations_scores <- c()
          for(p in 1:n_permutations){
            permutations_scores[[p]] <- list_permutations[[metabs]][[p]][row_name, col_names]
          }
          count_greater <- sum(permutations_scores >= real_score1)
          count_less <- sum(permutations_scores <= real_score1)
          pvalue <- min((2*min(count_greater, count_less) +1)/(n_permutations+1), 1.0)
          pvalue_list[[metabs]][row_name, col_names] <- pvalue
        }}}
    
    # ============================ pvalue adj ===========================
    pvalue_adj_list <- list()
    for(i in 1:length(list_metabolic)){
      pvalue_adj_list[[i]] <- list()
    }
    names(pvalue_adj_list) <- list_metabolic
    for(metabs in 1:length(list_metabolic)){
      adj_p_value_matrix <- matrix(
        p.adjust(as.vector(as.matrix(pvalue_list[[metabs]])), method = "BH"),
        nrow = nrow(pvalue_list[[metabs]])
      )
      dimnames(adj_p_value_matrix) <- dimnames(pvalue_list[[metabs]])
      
      pvalue_adj_list[[metabs]] <- adj_p_value_matrix
    }
    
    final_output <- list()
    for(i in 1:length(list_metabolic)){
      final_output[[i]] <- list()
    }
    names(final_output) <- list_metabolic
    for(metabs in 1:length(list_metabolic)){
      final_output[[metabs]] <- list(list_metabolic_original[[metabs]], pvalue_list[[metabs]], pvalue_adj_list[[metabs]])
      names(final_output[[metabs]]) <- c("normalized_score", "p_value", "adj_p_value")
    }
    return(final_output)
  }
  
  
  ui <- fluidPage(
    tags$head(
      tags$style(HTML("
      /* changing size of checkboxes */
      #dynamic_checkboxes_ui .checkbox label {
        font-size: 10px;
      }
    "))
    ),
    titlePanel("Analysis of Spatially distributed scores"),
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h4("1. Filter Data by slice:"),
        selectInput("gsm_filter_selector", "Select slice:",
                    choices = c("All", unique(metadata[,patient_column]))),
        hr(),
        h4("2. Select Scores to plot:"),
        uiOutput("dynamic_checkboxes_ui"),
        hr(),
        # MODIFICATION: Added numericInput for dot size
        numericInput("dot_size_selector", "Select Dot Size for Scatter Plots:",
                     value = 2, min = 0.5, max = 10, step = 0.5),
        hr(),
        actionButton("process_button", "Generate Plots & Heatmaps"),
        sliderInput("summary_heatmap_row_slider", "Select Displayed Row Range:", min = 1, max = 10, value = c(1, 10), step = 1),
        actionButton("generate_subset_heatmap_button", "Generate Subset Heatmap from Range"), 
        hr(),
        actionButton("all_samples_interaction_button", "All Samples Score Interaction"),
        hr(),
        h4("Selected Area Info:"),
        verbatimTextOutput("brush_info")
      ),
      mainPanel(
        width = 8,
        h4("Cell Type Spatial Distribution (select interesting area)"),
        plotOutput("fixed_cell_type_plot", height = "450px",
                   brush = brushOpts(id = "fixed_plot_brush", resetOnNew = TRUE)),
        hr(),
        uiOutput("brushed_heatmaps_ui"),
        hr(),
        uiOutput("main_dynamic_plots_ui"),
        hr(), # Added a separator
        h4("All Samples Interaction Heatmaps"), # Title for the new section
        uiOutput("all_samples_heatmaps_ui")
      )
    )
  )
  
  server <- function(input, output, session) {
    
    output$dynamic_checkboxes_ui <- renderUI({ lapply(list_of_scores_columns, function(col_name) {
      checkboxInput(inputId = paste0("chk_", make.names(col_name)), label = col_name, value = FALSE) }) %>% 
        do.call(tagList, .) })
    
    
    processed_data_and_scores <- eventReactive(input$process_button, {
      selected_cols <- character(0)
      for (col_name in list_of_scores_columns) {
        input_id <- paste0("chk_", make.names(col_name))
        if (!is.null(input[[input_id]]) && input[[input_id]] == TRUE) {
          selected_cols <- c(selected_cols, col_name)
        }
      }
      chosen_gsm <- input$gsm_filter_selector
      data_for_plots <- metadata
      if (chosen_gsm != "All") {
        data_for_plots <- data_for_plots[data_for_plots[,patient_column] == chosen_gsm, ]
      }
      proximity_function_output_main <- NULL
      if (nrow(data_for_plots) > 0 && length(selected_cols) > 0 && cell_type_column %in% names(data_for_plots)) {
        tryCatch({
          proximity_function_output_main <- proximity_metabolic_score_pvalue_perm_metab(
            meta_from_seurat = data_for_plots, cell_type_column = cell_type_column, list_metabolic = selected_cols, 
            sep = patient_column, k=k
          )
        }, error = function(e) { warning(paste("Error in main proximity_metabolic_score_all_matix:", e$message)) })
      }
      #   list(
      #     data_for_plotting = data_for_plots,
      #     columns_checked = selected_cols,
      #     active_gsm_filter = chosen_gsm,
      #     proximity_scores_main = proximity_function_output_main
      #   )
      # })
      full_heatmap_matrix <- NULL
      if (nrow(data_for_plots) > 0 && cell_type_column %in% names(data_for_plots)) {
        candidate_cols <- intersect(list_of_scores_columns, names(data_for_plots)); numeric_pathway_cols <- character(0)
        if(length(candidate_cols) > 0){ are_numeric <- sapply(data_for_plots[, candidate_cols, drop = FALSE], is.numeric); numeric_pathway_cols <- candidate_cols[are_numeric] }
        if (length(numeric_pathway_cols) > 0) { full_heatmap_matrix <- tryCatch({ t(data_for_plots[, numeric_pathway_cols, drop = FALSE]) }, error = function(e){warning(paste("Error creating full summary matrix:",e$message));NULL})
        } else { warning("No numeric pathway columns for summary heatmap.") }
      } else { warning("Not enough data for summary heatmap matrix.") }
      list( data_for_plotting = data_for_plots, columns_checked = selected_cols, active_gsm_filter = chosen_gsm, proximity_scores_main = proximity_function_output_main, summary_heatmap_full_matrix = full_heatmap_matrix )
    })
    
    main_heatmap_drawn_info_rv <- reactiveVal(NULL)
    observe({ # Updates slider based on drawn main heatmap
      info <- main_heatmap_drawn_info_rv()
      num_total_displayed_rows <- 0
      if(!is.null(info) && !is.null(info$pathway_names_in_display_order)){ num_total_displayed_rows <- length(info$pathway_names_in_display_order) }
      slider_max_val <- max(1, num_total_displayed_rows)
      current_slider_values <- isolate(input$summary_heatmap_row_slider)
      new_slider_value <- c(min(current_slider_values[1], slider_max_val), min(current_slider_values[2], slider_max_val))
      if(new_slider_value[1] > new_slider_value[2]) new_slider_value[1] <- new_slider_value[2]
      if(new_slider_value[1] < 1 && slider_max_val > 0) new_slider_value[1] <- 1
      if(new_slider_value[2] < 1 && slider_max_val > 0) new_slider_value[2] <- slider_max_val
      
      # Check if an update is truly needed to prevent infinite loops
      if (is.null(isolate(input$summary_heatmap_row_slider)) || 
          isolate(input$summary_heatmap_row_slider[1]) != new_slider_value[1] || 
          isolate(input$summary_heatmap_row_slider[2]) != new_slider_value[2] ||
          (exists("last_slider_max_val", envir = session$userData) && session$userData$last_slider_max_val != slider_max_val) ||
          !exists("last_slider_max_val", envir = session$userData) ) {
        
        updateSliderInput(session, "summary_heatmap_row_slider", max = slider_max_val, value = new_slider_value)
        session$userData$last_slider_max_val <- slider_max_val # Store last max in session data
      }
    })
    
    
    
    brushed_data_proximity_scores <- reactive({
      req(input$fixed_plot_brush, processed_data_and_scores())
      main_plot_data <- processed_data_and_scores()$data_for_plotting
      pathways_checked_for_brushed <- processed_data_and_scores()$columns_checked
      current_brush_object <- input$fixed_plot_brush
      if (nrow(main_plot_data) == 0) {
        return(list(brushed_data_exists = FALSE, brushed_df_for_summary = data.frame(), brushed_points_count = 0, rescaled_brush_coords = NULL))
      }
      actual_x_range <- if(sum(!is.na(main_plot_data$x)) > 0) range(main_plot_data$x, na.rm = TRUE) else c(0,1)
      actual_y_range <- if(sum(!is.na(main_plot_data$y)) > 0) range(main_plot_data$y, na.rm = TRUE) else c(0,1)
      brush_domain <- current_brush_object$domain
      brush_xmin_norm <- current_brush_object$xmin; brush_xmax_norm <- current_brush_object$xmax
      brush_ymin_norm <- current_brush_object$ymin; brush_ymax_norm <- current_brush_object$ymax
      rescaled_xmin <- NA; rescaled_xmax <- NA; rescaled_ymin <- NA; rescaled_ymax <- NA
      if (!is.null(brush_domain$left) && !is.null(brush_domain$right) && (brush_domain$right - brush_domain$left != 0)) {
        domain_x_span <- brush_domain$right - brush_domain$left; data_x_span <- actual_x_range[2] - actual_x_range[1]
        rescaled_xmin <- actual_x_range[1] + ((brush_xmin_norm - brush_domain$left) / domain_x_span) * data_x_span
        rescaled_xmax <- actual_x_range[1] + ((brush_xmax_norm - brush_domain$left) / domain_x_span) * data_x_span
      } else { warning("Brush X domain is invalid for rescaling.") }
      if (!is.null(brush_domain$bottom) && !is.null(brush_domain$top) && (brush_domain$top - brush_domain$bottom != 0)) {
        domain_y_span <- brush_domain$top - brush_domain$bottom; data_y_span <- actual_y_range[2] - actual_y_range[1]
        rescaled_ymin <- actual_y_range[1] + ((brush_ymin_norm - brush_domain$bottom) / domain_y_span) * data_y_span
        rescaled_ymax <- actual_y_range[1] + ((brush_ymax_norm - brush_domain$bottom) / domain_y_span) * data_y_span
      } else { warning("Brush Y domain is invalid for rescaling.") }
      rescaled_coords_list <- list(xmin=rescaled_xmin, xmax=rescaled_xmax, ymin=rescaled_ymin, ymax=rescaled_ymax)
      brushed_df <- data.frame()
      if (!anyNA(c(rescaled_xmin, rescaled_xmax, rescaled_ymin, rescaled_ymax))) {
        brushed_df <- main_plot_data[ main_plot_data$x >= rescaled_xmin & main_plot_data$x <= rescaled_xmax & main_plot_data$y >= rescaled_ymin & main_plot_data$y <= rescaled_ymax & !is.na(main_plot_data$x) & !is.na(main_plot_data$y), , ]
      } else { warning("Could not rescale brush coordinates properly; brushed_df will be empty.") }
      if (cell_type_column %in% names(brushed_df) && is.factor(brushed_df[,cell_type_column])) {
        brushed_df[,cell_type_column] <- droplevels(brushed_df[,cell_type_column])
      }
      data_exists_for_processing <- nrow(brushed_df) > 0 && cell_type_column %in% names(brushed_df)
      proximity_output_brushed <- NULL
      if (data_exists_for_processing && length(pathways_checked_for_brushed) > 0) {
        tryCatch({
          proximity_output_brushed <- proximity_metabolic_score_pvalue_perm_metab( meta_from_seurat = brushed_df, k=k,
                                                                                   cell_type_column = cell_type_column, 
                                                                                   list_metabolic = pathways_checked_for_brushed, sep = patient_column)
        }, error = function(e) { warning(paste("Error in brushed proximity_metabolic_score_all_matix:", e$message)) })
      }
      return(list(
        brushed_data_exists = data_exists_for_processing,
        brushed_df_for_summary = brushed_df,
        proximity_scores_brushed = proximity_output_brushed,
        columns_for_pathway_specific_brushed_heatmap = pathways_checked_for_brushed,
        brushed_points_count = nrow(brushed_df),
        rescaled_brush_coords = rescaled_coords_list
      ))
    })
    
    output$brush_info <- renderPrint({
      current_brush <- input$fixed_plot_brush
      if (is.null(current_brush)) { return("Brush the 'Overall Cell Type Distribution' plot to select points for new heatmaps.") }
      results <- brushed_data_proximity_scores()
      if (!is.null(results)) {
        count <- results$brushed_points_count
        cat("--- Brush Selection ---\n")
        if(!is.null(results$rescaled_brush_coords) && !anyNA(unlist(results$rescaled_brush_coords))) {
          cat(sprintf("Rescaled Brush Area (approx.):\nX: %.2f to %.2f\nY: %.2f to %.2f\n", results$rescaled_brush_coords$xmin, results$rescaled_brush_coords$xmax, results$rescaled_brush_coords$ymin, results$rescaled_brush_coords$ymax))
        } else { cat("Brush coordinates could not be fully rescaled.\nRaw brush (normalized-like):\n"); print(current_brush[c("xmin", "xmax", "ymin", "ymax")]) }
        if (isTRUE(results$brushed_data_exists)) {
          cat(paste("\nSelected area includes", count, "points usable for heatmaps.\n"))
          if(count > 0 && is.null(results$proximity_scores_brushed) && length(results$columns_for_pathway_specific_brushed_heatmap) > 0){ cat("Note: Proximity function (pathway-specific) on brushed data might have returned NULL or failed.\n") }
        } else { cat(paste("\nSelected area includes", count, "points. Conditions not met for generating heatmaps from brushed data (e.g., 0 points actually selected, or critical 'cell_type_column' column missing in selection).\n")) }
      } else { "Processing brush..." }
    })
    
    output$fixed_cell_type_plot <- renderPlot({
      req(processed_data_and_scores()); plot_details <- processed_data_and_scores(); data_to_use_scatter <- plot_details$data_for_plotting; active_filter_name <- plot_details$active_gsm_filter
      if (nrow(data_to_use_scatter) == 0) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, paste("No data for patient filter '", active_filter_name, "'."), cex = 1.2, col="red"); return() }
      if (!cell_type_column %in% names(data_to_use_scatter)) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, "Error: cell_type_column column missing.", cex = 1.2, col="red"); return() }
      plot_title_fixed <- "Scatter Plot by Cell Type"; if (active_filter_name != "All") plot_title_fixed <- paste(plot_title_fixed, "\n(patient Filter:", active_filter_name, ")")
      # MODIFICATION: Use input$dot_size_selector
      p_fixed <- ggplot(data_to_use_scatter, aes(x = .data$x, y = .data$y, color = factor(.data[[cell_type_column]]))) +
        geom_point(size = input$dot_size_selector) +
        labs(title = plot_title_fixed, color = "Cell Type") +
        theme_classic() +
        theme(plot.title = element_text(size = 10)) +
        coord_cartesian(xlim = if(sum(!is.na(data_to_use_scatter$x)) > 0) range(data_to_use_scatter$x, na.rm = TRUE) else NULL, 
                        ylim = if(sum(!is.na(data_to_use_scatter$y)) > 0) range(data_to_use_scatter$y, na.rm = TRUE) else NULL, 
                        expand = TRUE) +
        scale_color_manual(values = color) # THIS LINE IS NEW
      
      print(p_fixed)
    })
    
    output$brushed_heatmaps_ui <- renderUI({
      current_brush_val_ui <- input$fixed_plot_brush; if(is.null(current_brush_val_ui)){ return(NULL) }
      brushed_results <- brushed_data_proximity_scores()
      all_brushed_ui_elements <- list()
      if (isTRUE(brushed_results$brushed_data_exists) && !is.null(brushed_results$brushed_df_for_summary) && nrow(brushed_results$brushed_df_for_summary) > 0) {
        all_brushed_ui_elements[[length(all_brushed_ui_elements) + 1]] <- h4("Heatmap all Scores, selected cells")
        all_brushed_ui_elements[[length(all_brushed_ui_elements) + 1]] <- plotOutput(outputId = "brushed_summary_pathway_heatmap", height = "600px")
        all_brushed_ui_elements[[length(all_brushed_ui_elements) + 1]] <- hr()
      }
      if (isTRUE(brushed_results$brushed_data_exists) && !is.null(brushed_results$proximity_scores_brushed) && !is.null(brushed_results$columns_for_pathway_specific_brushed_heatmap) && length(brushed_results$columns_for_pathway_specific_brushed_heatmap) > 0) {
        all_brushed_ui_elements[[length(all_brushed_ui_elements) + 1]] <- h4("Heatmaps of Cell-cell Scores Spatially Interacting (just selected cells)")
        pathway_specific_brushed_heatmaps_ui <- lapply(brushed_results$columns_for_pathway_specific_brushed_heatmap, function(col_name) {
          specific_brushed_heatmap_data <- brushed_results$proximity_scores_brushed[[col_name]]
          if (!is.null(specific_brushed_heatmap_data) && is.list(specific_brushed_heatmap_data) && length(specific_brushed_heatmap_data) >= 1 && (is.data.frame(specific_brushed_heatmap_data[[1]]) || is.matrix(specific_brushed_heatmap_data[[1]])) ) {
            plotOutput(outputId = paste0("brushed_heatmap_", make.names(col_name)), height = "500px")
          } else { tags$div(style="color:orange;", paste("Could not prepare UI for pathway-specific brushed heatmap:", col_name)) }
        })
        all_brushed_ui_elements <- c(all_brushed_ui_elements, Filter(Negate(is.null), pathway_specific_brushed_heatmaps_ui))
        if(length(pathway_specific_brushed_heatmaps_ui) > 0 && sum(!sapply(pathway_specific_brushed_heatmaps_ui, is.null)) > 0 ) all_brushed_ui_elements[[length(all_brushed_ui_elements) + 1]] <- hr()
      }
      if (length(all_brushed_ui_elements) > 0) {
        return(do.call(tagList, all_brushed_ui_elements))
      } else if (!is.null(current_brush_val_ui)) {
        return(tags$p(style="font-style: italic;", "Brush selection active. If heatmaps from selection are not shown, it may be due to insufficient points, no pathways chosen for specific heatmaps, or a data processing error for the selected area."))
      }
      return(NULL)
    })
    
    output$main_dynamic_plots_ui <- renderUI({
      req(input$process_button > 0); plot_info <- processed_data_and_scores()
      all_main_ui_elements <- list()
      cols_to_plot_scatter <- plot_info$columns_checked
      if (!is.null(cols_to_plot_scatter) && length(cols_to_plot_scatter) > 0) {
        if (nrow(plot_info$data_for_plotting) > 0) { all_main_ui_elements[[length(all_main_ui_elements) + 1]] <- h4("Spatial plots of Selected Scores") }
        scatter_plot_ui <- lapply(cols_to_plot_scatter, function(col_name) { tagList(plotlyOutput(outputId = paste0("plot_", make.names(col_name)), height = "450px"), hr()) })
        all_main_ui_elements <- c(all_main_ui_elements, scatter_plot_ui)
      } else if (input$process_button > 0) { all_main_ui_elements[[length(all_main_ui_elements) + 1]] <- tags$p("No pathways CHECKED for dynamic scatter plots.") }
      cols_for_main_heatmap <- plot_info$columns_checked
      main_proximity_data_available <- !is.null(plot_info$proximity_scores_main)
      if (main_proximity_data_available && !is.null(cols_for_main_heatmap) && length(cols_for_main_heatmap) > 0) {
        all_main_ui_elements[[length(all_main_ui_elements) + 1]] <- h4("Heatmaps of Cell-cell Scores Spatially Interacting (row is active when close to column)")
        main_heatmap_ui_elements <- lapply(cols_for_main_heatmap, function(col_name) {
          specific_main_heatmap_data <- plot_info$proximity_scores_main[[col_name]]
          if (!is.null(specific_main_heatmap_data) && is.list(specific_main_heatmap_data) && length(specific_main_heatmap_data) >= 1 && (is.data.frame(specific_main_heatmap_data[[1]]) || is.matrix(specific_main_heatmap_data[[1]])) ) {
            tagList(plotOutput(outputId = paste0("main_heatmap_", make.names(col_name)), height = "500px"), hr())
          } else { tags$div(tags$h5(paste("Main Heatmap for:", col_name)), tags$p(style="color:orange;", "Data for this main proximity heatmap is not available or in an unexpected format.")) }
        })
        all_main_ui_elements <- c(all_main_ui_elements, Filter(Negate(is.null), main_heatmap_ui_elements))
      } else if (input$process_button > 0 && !is.null(cols_for_main_heatmap) && length(cols_for_main_heatmap) > 0 && !main_proximity_data_available) {
        all_main_ui_elements[[length(all_main_ui_elements) + 1]] <- tags$p(style="color:red;", "Main proximity score data (for checked pathways) could not be generated.")
      }
      
      #changed here ==================================================
      # --- UI for Main Summary Heatmap (static) AND the NEW Sub-Heatmap based on slider ---
      plot_info_check_main <- isolate(processed_data_and_scores()) 
      if (input$process_button > 0 && !is.null(plot_info_check_main$data_for_plotting) && nrow(plot_info_check_main$data_for_plotting) > 0) {
        if (!is.null(plot_info_check_main$summary_heatmap_full_matrix) && nrow(plot_info_check_main$summary_heatmap_full_matrix) > 0) { # Check if matrix exists
          all_main_ui_elements[[length(all_main_ui_elements) + 1]] <- h4("Heatmap all Scores and all Cells")
          all_main_ui_elements[[length(all_main_ui_elements) + 1]] <- plotOutput(outputId = "summary_pathway_heatmap", height = "700px") 
          all_main_ui_elements[[length(all_main_ui_elements) + 1]] <- hr()
          all_main_ui_elements[[length(all_main_ui_elements) + 1]] <- h4("Heatmap selected Scores and all Cells")
          all_main_ui_elements[[length(all_main_ui_elements) + 1]] <- plotOutput(outputId = "subset_summary_heatmap", height = "600px")
          all_main_ui_elements[[length(all_main_ui_elements) + 1]] <- hr()
        } else {
          all_main_ui_elements[[length(all_main_ui_elements) + 1]] <- tags$p("Not enough data or pathways to generate the overall summary heatmap.")
        }
      } else if (input$process_button > 0) {
        all_main_ui_elements[[length(all_main_ui_elements) + 1]] <- tags$p("Not enough data after GSM filter to generate any summary heatmaps.")
      }
      do.call(tagList, all_main_ui_elements)
    })
    
    # This eventReactive will now hold the data needed for the subset heatmap,
    # triggered by the new button.
    subset_heatmap_data_trigger <- eventReactive(input$generate_subset_heatmap_button, {
      info <- main_heatmap_drawn_info_rv()
      slider_vals <- input$summary_heatmap_row_slider
      req(info, info$pathway_names_in_display_order, info$full_matrix_original_order, slider_vals)
      
      ordered_pathways <- info$pathway_names_in_display_order
      full_matrix <- info$full_matrix_original_order
      df_for_annotation <- info$df_for_annotation # Get the df used for main heatmap annotation
      
      num_displayed_rows <- length(ordered_pathways)
      start_idx <- max(1, slider_vals[1])
      end_idx <- min(num_displayed_rows, slider_vals[2])
      if (start_idx > end_idx || num_displayed_rows == 0) return(NULL) # Invalid range or no rows
      
      pathways_for_subset <- ordered_pathways[start_idx:end_idx]
      sub_matrix <- full_matrix[pathways_for_subset, , drop = FALSE]
      
      return(list(
        sub_matrix = sub_matrix, 
        df_for_annotation = df_for_annotation, # Pass the correct df for annotation
        active_gsm_filter = isolate(processed_data_and_scores()$active_gsm_filter), # Get current GSM filter
        slider_display_range = c(start_idx, end_idx)
      ))
    })
    #   # ================================== maybe change in here
    #   if (input$process_button > 0 && nrow(plot_info$data_for_plotting) > 0) {
    #     all_main_ui_elements[[length(all_main_ui_elements) + 1]] <- h4("Heatmap all Scores and all Cells")
    #     all_main_ui_elements[[length(all_main_ui_elements) + 1]] <- plotOutput(outputId = "summary_pathway_heatmap", height = "600px")
    #     all_main_ui_elements[[length(all_main_ui_elements) + 1]] <- hr()
    #   }
    #   do.call(tagList, all_main_ui_elements)
    # })
    
    observe({
      plot_details_main_obs <- processed_data_and_scores() # Triggered by process_button
      if (!is.null(plot_details_main_obs)) {
        data_to_use_scatter <- plot_details_main_obs$data_for_plotting
        cols_for_dynamic_plots <- plot_details_main_obs$columns_checked
        active_filter_name <- plot_details_main_obs$active_gsm_filter
        proximity_scores_main_data <- plot_details_main_obs$proximity_scores_main
        # --- Render Dynamic Scatter Plots (ggplotly) ---
        if (!is.null(cols_for_dynamic_plots) && length(cols_for_dynamic_plots) > 0) {
          for (col_name_iter_scatter in cols_for_dynamic_plots) { local({ current_col_scatter <- col_name_iter_scatter; safe_col_name_scatter <- make.names(current_col_scatter); output[[paste0("plot_", safe_col_name_scatter)]] <- renderPlotly({ if (nrow(data_to_use_scatter) == 0) { return(plotly_empty(type = "scatter", mode = "markers") %>% layout(title = "No data to display")) } ; if (!current_col_scatter %in% names(data_to_use_scatter) || !is.numeric(data_to_use_scatter[[current_col_scatter]])) {  return(plotly_empty(type = "scatter", mode = "markers") %>% layout(title = paste("Cannot plot '", current_col_scatter, "'.\nNot numeric or not found."))) }; if(all(is.na(data_to_use_scatter[[current_col_scatter]]))) { return(plotly_empty(type = "scatter", mode = "markers") %>% layout(title = paste("All values for '", current_col_scatter, "' are NA."))) }; midpoint_val <- median(data_to_use_scatter[[current_col_scatter]], na.rm = TRUE); plot_title_dynamic <- paste("Scatter: x vs y, by", current_col_scatter); if (active_filter_name != "All") plot_title_dynamic <- paste(plot_title_dynamic, "\n(GSM Filter:", active_filter_name, ")"); if (!"cell" %in% names(data_to_use_scatter)) data_to_use_scatter$cell <- rownames(data_to_use_scatter); p_dynamic <- ggplot(data_to_use_scatter, aes(x = x, y = y, color = !!sym(current_col_scatter), text = paste(#"Cell:", cell,
            "<br>Cell Type:", .data[[cell_type_column]], "<br>", current_col_scatter, ":", round(!!sym(current_col_scatter), 2) ))) + geom_point(size = input$dot_size_selector) + scale_color_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = midpoint_val, na.value="grey50") + labs(title = plot_title_dynamic, color = "Score") + theme_classic() + theme(plot.title = element_text(size = 10)); ggplotly(p_dynamic, tooltip = "text") }) }) }
        }
        # --- Render Main Proximity Heatmaps ---
        if (!is.null(proximity_scores_main_data) && !is.null(cols_for_dynamic_plots) && length(cols_for_dynamic_plots) > 0) {
          for (col_name_iter_heatmap in cols_for_dynamic_plots) { local({ current_col_heatmap <- col_name_iter_heatmap; safe_col_name_heatmap <- make.names(current_col_heatmap); output[[paste0("main_heatmap_", safe_col_name_heatmap)]] <- renderPlot({ if (is.null(proximity_scores_main_data[[current_col_heatmap]]) || !is.list(proximity_scores_main_data[[current_col_heatmap]]) || length(proximity_scores_main_data[[current_col_heatmap]]) < 1 || !(is.data.frame(proximity_scores_main_data[[current_col_heatmap]][[1]]) || is.matrix(proximity_scores_main_data[[current_col_heatmap]][[1]])) ) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, paste("Data for main heatmap '", current_col_heatmap, "'\nis missing or not in expected format."), cex = 1.2, col="red"); return() }; 
            matrix_to_plot <- as.matrix(proximity_scores_main_data[[current_col_heatmap]][[1]]); 
            annotation_matrix <- as.matrix(proximity_scores_main_data[[current_col_heatmap]][["adj_p_value"]]);
            if (nrow(matrix_to_plot) == 0 || ncol(matrix_to_plot) == 0) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, paste("Matrix for main heatmap '", current_col_heatmap, "' is empty."), cex = 1.2, col="orange"); return() }; heatmap_title <- paste("Main HM:", current_col_heatmap); if (active_filter_name != "All") { heatmap_title <- paste(heatmap_title, "(GSM:", active_filter_name,")") }; tryCatch({ ht <- Heatmap(matrix_to_plot, cluster_rows = FALSE, cluster_columns = FALSE, name = "Score", column_title = heatmap_title, column_title_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
                                                                                                                                                                                                                                                                                                                                                                                                                                                   cell_fun = function(j, i, x, y, width, height, fill, p_val_mat = annotation_matrix) {
                                                                                                                                                                                                                                                                                                                                                                                                                                                     # Check the value from text_matrix to decide the color
                                                                                                                                                                                                                                                                                                                                                                                                                                                     text_color <- if (p_val_mat[i, j] < 0.05) {
                                                                                                                                                                                                                                                                                                                                                                                                                                                       "green"
                                                                                                                                                                                                                                                                                                                                                                                                                                                     } else {
                                                                                                                                                                                                                                                                                                                                                                                                                                                       "black"
                                                                                                                                                                                                                                                                                                                                                                                                                                                     }
                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                                                                                                                                                                                                                     # Use text_matrix[i, j] for the label and the determined color for the text
                                                                                                                                                                                                                                                                                                                                                                                                                                                     grid.text(
                                                                                                                                                                                                                                                                                                                                                                                                                                                       sprintf("%.3f", p_val_mat[i, j]), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                       x, y, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                       gp = gpar(fontsize = 10, col = text_color)
                                                                                                                                                                                                                                                                                                                                                                                                                                                     )
                                                                                                                                                                                                                                                                                                                                                                                                                                                   }); draw(ht) }, error = function(e_hm) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, paste("Error for main heatmap '", current_col_heatmap, "':\n", e_hm$message), cex = 0.9, col="darkred") }) }) }) }
        }
        
        output$summary_pathway_heatmap <- renderPlot({
          # This req is important because this renderPlot now depends on the button too
          req(plot_details_main_obs$data_for_plotting) 
          
          df_filtred_main_hm <- plot_details_main_obs$data_for_plotting
          full_heatmap_matrix_main <- plot_details_main_obs$summary_heatmap_full_matrix 
          active_filter_name_sum_main <- plot_details_main_obs$active_gsm_filter
          
          if (is.null(full_heatmap_matrix_main) || nrow(full_heatmap_matrix_main) == 0 ) {
            plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, "Not enough data or pathways for the main summary heatmap.", col="red"); 
            main_heatmap_drawn_info_rv(NULL) # Clear/set empty if no data
            return()
          }
          if (!cell_type_column %in% names(df_filtred_main_hm)) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, "'cell_type_column' column missing for summary heatmap.", col="red"); return() }
          
          cell_types_for_annotation_main <- factor(df_filtred_main_hm[[cell_type_column]])
          unique_cts_main <- levels(cell_types_for_annotation_main)
          column_ha_main <- HeatmapAnnotation(df=data.frame(cell_type = cell_types_for_annotation_main), col = list(cell_type = color),annotation_name_gp = gpar(fontsize = 8), simple_anno_size = unit(0.3, "cm"))
          # Color generation for annotation (same as before)
          summary_heatmap_title_main <- "Scores x cells, separated by Cell Type"; if (active_filter_name_sum_main != "All") { summary_heatmap_title_main <- paste(summary_heatmap_title_main, "\n( Patient/slice:", active_filter_name_sum_main, ")") }
          ht_summary_obj_main <- NULL 
          tryCatch({ 
            ht_summary_obj_main <- Heatmap(full_heatmap_matrix_main, name = "Pathways Score", column_title = summary_heatmap_title_main, 
                                           show_column_names = FALSE, 
                                           cluster_rows = TRUE, # Let Heatmap do the clustering
                                           cluster_columns = TRUE, row_names_gp = gpar(fontsize=6), use_raster = FALSE, 
                                           top_annotation = column_ha_main, column_split = cell_types_for_annotation_main, 
                                           column_title_gp = gpar(fontsize = 10))
            drawn_ht_main <- draw(ht_summary_obj_main) 
            
            final_row_order_indices_displayed <- unlist(ComplexHeatmap::row_order(drawn_ht_main)) 
            path_names_in_display <- rownames(full_heatmap_matrix_main)[final_row_order_indices_displayed]
            
            main_heatmap_drawn_info_rv(list(
              pathway_names_in_display_order = path_names_in_display,
              full_matrix_original_order = full_heatmap_matrix_main,
              df_for_annotation = df_filtred_main_hm # Store the df used for annotation
            ))
          }, error = function(e_sum_hm) { 
            plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, paste("Error generating summary heatmap:\n", e_sum_hm$message), cex = 0.9, col="darkred") 
            main_heatmap_drawn_info_rv(list(pathway_names_in_display_order = character(0), full_matrix_original_order = NULL, df_for_annotation = NULL))
          })
        }) # End summary_pathway_heatmap renderPlot
      } # End if (!is.null(plot_details_main_obs))
      current_brush_val_obs <- input$fixed_plot_brush
      if (!is.null(current_brush_val_obs)) {
        brushed_results <- brushed_data_proximity_scores()
        if (isTRUE(brushed_results$brushed_data_exists)) {
          active_filter_name_for_brushed_title <- processed_data_and_scores()$active_gsm_filter
          output$brushed_summary_pathway_heatmap <- renderPlot({
            req(brushed_results$brushed_df_for_summary)
            df_brushed_summary <- brushed_results$brushed_df_for_summary
            if (nrow(df_brushed_summary) == 0) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, "No points selected in brush for summary heatmap.", col="orange"); return() }
            if (!cell_type_column %in% names(df_brushed_summary)) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, "cell_type_column column missing in brushed data.", col="red"); return() }
            candidate_cols_brushed <- intersect(list_of_scores_columns, names(df_brushed_summary))
            numeric_pathway_cols_brushed <- character(0)
            if(length(candidate_cols_brushed) > 0){ numeric_pathway_cols_brushed <- candidate_cols_brushed[sapply(df_brushed_summary[, candidate_cols_brushed, drop = FALSE], is.numeric)] }
            if (length(numeric_pathway_cols_brushed) == 0) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, "No numeric pathway columns in brushed data for summary heatmap.", col="orange"); return() }
            heatmap_matrix_brushed_summary <- tryCatch({ t(df_brushed_summary[, numeric_pathway_cols_brushed, drop = FALSE]) }, error = function(e){ NULL })
            if (is.null(heatmap_matrix_brushed_summary) || nrow(heatmap_matrix_brushed_summary) == 0 || ncol(heatmap_matrix_brushed_summary) == 0) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, "Matrix for brushed summary heatmap is empty or could not be prepared.", col="orange"); return() }
            cell_types_for_anno_brushed <- factor(df_brushed_summary[,cell_type_column])
            unique_cts_brushed <- levels(cell_types_for_anno_brushed)
            column_ha_brushed <- HeatmapAnnotation(cell_type = cell_types_for_anno_brushed, col = list(cell_type = color), annotation_name_gp = gpar(fontsize = 8), simple_anno_size = unit(0.3, "cm"))
            
            title_brushed_summary <- "Scores x selected cells, separated by Cell Type"; if (active_filter_name_for_brushed_title != "All") { title_brushed_summary <- paste(title_brushed_summary, "\n( Patient/slice:", active_filter_name_for_brushed_title, ")") }; title_brushed_summary <- paste(title_brushed_summary, "-", brushed_results$brushed_points_count, "pts selected")
            tryCatch({ ht_brushed_summary <- Heatmap(heatmap_matrix_brushed_summary, name = "Pathways Score", column_title = title_brushed_summary, show_column_names = FALSE, cluster_rows = TRUE, cluster_columns = TRUE, row_names_gp = gpar(fontsize=6), use_raster = FALSE, top_annotation = column_ha_brushed, column_split = cell_types_for_anno_brushed, column_title_gp = gpar(fontsize = 10)); draw(ht_brushed_summary)
            }, error = function(e_brush_sum_hm) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, paste("Error generating brushed summary heatmap:\n", e_brush_sum_hm$message), cex = 0.9, col="darkred") })
          })
          if (!is.null(brushed_results$proximity_scores_brushed) && !is.null(brushed_results$columns_for_pathway_specific_brushed_heatmap) && length(brushed_results$columns_for_pathway_specific_brushed_heatmap) > 0) {
            for (col_name_iter_brushed_hm in brushed_results$columns_for_pathway_specific_brushed_heatmap) {
              local({
                current_col_brushed_hm <- col_name_iter_brushed_hm; safe_col_name_brushed_hm <- make.names(current_col_brushed_hm)
                output[[paste0("brushed_heatmap_", safe_col_name_brushed_hm)]] <- renderPlot({
                  if (is.null(brushed_results$proximity_scores_brushed[[current_col_brushed_hm]]) || !is.list(brushed_results$proximity_scores_brushed[[current_col_brushed_hm]]) || length(brushed_results$proximity_scores_brushed[[current_col_brushed_hm]]) < 1 || !(is.data.frame(brushed_results$proximity_scores_brushed[[current_col_brushed_hm]][[1]]) || is.matrix(brushed_results$proximity_scores_brushed[[current_col_brushed_hm]][[1]])) ) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, paste("Data for pathway-specific brushed heatmap '", current_col_brushed_hm, "'\nis missing or not in expected format."), cex = 1.2, col="red"); return() }
                  matrix_to_plot_brushed <- as.matrix(brushed_results$proximity_scores_brushed[[current_col_brushed_hm]][[1]])
                  annotation_matrix_brushed <- as.matrix(brushed_results$proximity_scores_brushed[[current_col_brushed_hm]][["adj_p_value"]])
                  if (nrow(matrix_to_plot_brushed) == 0 || ncol(matrix_to_plot_brushed) == 0) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, paste("Matrix for pathway-specific brushed heatmap '", current_col_brushed_hm, "' is empty."), cex = 1.2, col="orange"); return() }
                  heatmap_title_brushed_prox <- paste("Score Interacted:", current_col_brushed_hm); if (active_filter_name_for_brushed_title != "All") { heatmap_title_brushed_prox <- paste(heatmap_title_brushed_prox, "( Patient/slice:", active_filter_name_for_brushed_title,")") }; heatmap_title_brushed_prox <- paste(heatmap_title_brushed_prox, "-", brushed_results$brushed_points_count, "pts")
                  tryCatch({ ht_brushed <- Heatmap(matrix_to_plot_brushed, cluster_rows = FALSE, cluster_columns = FALSE, name = "Score", column_title = heatmap_title_brushed_prox, column_title_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
                                                   cell_fun = function(j, i, x, y, width, height, fill, p_val_mat = annotation_matrix_brushed) {
                                                     # Check the value from text_matrix to decide the color
                                                     text_color <- if (p_val_mat[i, j] < 0.05) {
                                                       "green"
                                                     } else {
                                                       "black"
                                                     }
                                                     
                                                     # Use text_matrix[i, j] for the label and the determined color for the text
                                                     grid.text(
                                                       sprintf("%.3f", p_val_mat[i, j]), 
                                                       x, y, 
                                                       gp = gpar(fontsize = 10, col = text_color)
                                                     )
                                                   }); draw(ht_brushed)
                  }, error = function(e_hm_b) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, paste("Error for pathway-specific brushed heatmap '", current_col_brushed_hm, "':\n", e_hm_b$message), cex = 0.9, col="darkred") })
                })
              })
            }
          }
        }
      }
    })
    # --- Render SUBSET of Main Summary Heatmap (based on slider AND new button) ---
    # This is triggered by the new button `generate_subset_heatmap_button`
    observeEvent(input$generate_subset_heatmap_button, {
      output$subset_summary_heatmap <- renderPlot({
        plot_info_for_subset <- main_heatmap_drawn_info_rv() 
        slider_vals <- input$summary_heatmap_row_slider
        
        req(plot_info_for_subset, 
            plot_info_for_subset$pathway_names_in_display_order, 
            plot_info_for_subset$full_matrix_original_order,
            plot_info_for_subset$df_for_annotation, # Make sure this is passed
            slider_vals)
        
        ordered_pathways_main_heatmap <- plot_info_for_subset$pathway_names_in_display_order
        full_matrix_main_heatmap <- plot_info_for_subset$full_matrix_original_order
        df_annotation_for_subset <- plot_info_for_subset$df_for_annotation # Use this
        
        num_displayed_rows_main <- length(ordered_pathways_main_heatmap)
        start_idx_display <- max(1, slider_vals[1])
        end_idx_display <- min(num_displayed_rows_main, slider_vals[2])
        
        if (start_idx_display > end_idx_display || num_displayed_rows_main == 0 || start_idx_display > num_displayed_rows_main ) { # Added check for start_idx > total
          plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, "Invalid row range selected or no rows to display in subset heatmap.", col="orange"); return()
        }
        
        pathways_for_subset <- ordered_pathways_main_heatmap[start_idx_display:end_idx_display]
        if (length(pathways_for_subset) == 0) { # if range results in no pathways
          plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, "Selected range yields no pathways.", col="orange"); return()
        }
        
        sub_matrix <- full_matrix_main_heatmap[pathways_for_subset, , drop = FALSE]
        
        if (nrow(sub_matrix) == 0 || ncol(sub_matrix) == 0) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, "Selected row range results in an empty sub-heatmap.", col="orange"); return() }
        
        req(df_annotation_for_subset) # df_for_annotation passed from main_heatmap_drawn_info_rv
        active_filter_sub <- isolate(processed_data_and_scores()$active_gsm_filter) # Get current GSM filter
        
        if (!cell_type_column %in% names(df_annotation_for_subset)) { plot(1,type="n",axes=F,xlab="",ylab=""); text(1,1,"'cell_type_column' missing for sub-heatmap anno.", col="red"); return()}
        
        # Column annotation logic needs to use df_annotation_for_subset which corresponds to all cells in full_matrix_main_heatmap
        cell_types_for_sub_annotation <- factor(df_annotation_for_subset[[cell_type_column]])
        # Ensure it matches the columns of sub_matrix (which are all cells from full_matrix_main_heatmap)
        if(length(cell_types_for_sub_annotation) != ncol(sub_matrix)){
          plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, "Cell annotation length mismatch for sub-heatmap.", col="red"); return()
        }
        
        unique_cts_sub <- levels(cell_types_for_sub_annotation)
        column_ha_sub <- HeatmapAnnotation(df=data.frame(cell_type = cell_types_for_sub_annotation), col = list(cell_type = color), annotation_name_gp = gpar(fontsize = 8), simple_anno_size = unit(0.3, "cm"))
        # Color generation logic for annotation
        
        title_sub_heatmap <- paste("Subset of Pathways (Displayed Rows", start_idx_display, "-", end_idx_display, " from Main Clustered Heatmap)");
        if (active_filter_sub != "All") { title_sub_heatmap <- paste(title_sub_heatmap, "\n(GSM Filter:", active_filter_sub, ")") }
        
        tryCatch({
          ht_sub <- Heatmap(sub_matrix, name = "Pathways Score", column_title = title_sub_heatmap, 
                            show_column_names = FALSE, 
                            cluster_rows = FALSE, 
                            row_order = 1:nrow(sub_matrix), # Display in the selected order
                            cluster_columns = TRUE, 
                            row_names_gp = gpar(fontsize = 8), use_raster = FALSE, 
                            top_annotation = column_ha_sub, 
                            column_split = cell_types_for_sub_annotation, 
                            column_title_gp = gpar(fontsize = 10))
          draw(ht_sub)
        }, error = function(e_sub_hm) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, paste("Error generating sub-heatmap:\n", e_sub_hm$message), cex = 0.9, col="darkred") })
      }) # End renderPlot for subset_summary_heatmap
    }) # End observeEvent for generate_subset_heatmap_button
    
    
    #new functionality of score for every category in selected column proximity_split
    # --- EventReactive and Observe/renderPlot for "All Samples Score Interaction" ---
    all_samples_processed_cache_rv <- reactiveVal(NULL) # Moved definition to top of server
    
    observeEvent(input$all_samples_interaction_button, {
      message(paste("--- 'All Samples Interaction' Button Clicked at:", Sys.time(), "---"))
      current_selected_cols <- character(0)
      for (col_name in list_of_scores_columns) { input_id <- paste0("chk_", make.names(col_name)); if (!is.null(input[[input_id]]) && input[[input_id]] == TRUE) { current_selected_cols <- c(current_selected_cols, col_name) } }
      if (length(current_selected_cols) == 0) { showNotification("Please select at least one pathway/score for 'All Samples Interaction'.", type = "warning", duration = 5); all_samples_processed_cache_rv(list(type="error", message="No pathways selected")); return() }
      
      original_metadata_as <- metadata 
      tumor_col_name_as <- proximity_split 
      tumor_col_exists_as <- nzchar(tumor_col_name_as) && tumor_col_name_as %in% names(original_metadata_as)
      unique_tumor_values_as <- NULL
      if(tumor_col_exists_as) { unique_tumor_values_as <- unique(original_metadata_as[,tumor_col_name_as][!is.na(original_metadata_as[,tumor_col_name_as])]) }
      
      results_payload <- NULL
      if (tumor_col_exists_as && length(unique_tumor_values_as) > 0) {
        split_by_tumor_as <- split(original_metadata_as, original_metadata_as[,tumor_col_name_as])
        results_by_tumor_as <- list(); tumor_names_processed_as <- character(0)
        for(tumor_name_iter_as in names(split_by_tumor_as)){
          tumor_subset_df_as <- split_by_tumor_as[[tumor_name_iter_as]]
          if (nrow(tumor_subset_df_as) > 0 && cell_type_column %in% names(tumor_subset_df_as) && all(c("x", "y") %in% names(tumor_subset_df_as))) {
            proximity_output_as <- tryCatch({ proximity_metabolic_score_pvalue_perm_metab( k=k, meta_from_seurat = tumor_subset_df_as, cell_type_column = cell_type_column, list_metabolic = current_selected_cols, sep = patient_column ) }, error = function(e){ warning(paste("Error in prox function for group '", tumor_name_iter_as, "': ", e$message)); NULL })
            if(!is.null(proximity_output_as) && length(proximity_output_as) > 0) { results_by_tumor_as[[tumor_name_iter_as]] <- proximity_output_as; tumor_names_processed_as <- c(tumor_names_processed_as, tumor_name_iter_as) }
          } else { warning(paste("Skipping group '", tumor_name_iter_as, "' due to no data or missing columns.")) }
        }
        if(length(results_by_tumor_as) > 0){ results_payload <- list(type = "split_by_group", results = results_by_tumor_as, selected_pathways = current_selected_cols, group_names = tumor_names_processed_as, group_column_name = tumor_col_name_as)
        } else { showNotification(paste0("No groups based on '", tumor_col_name_as, "' could be processed."), type = "warning", duration = 5); results_payload <- list(type="error", message=paste0("No groups based on '", tumor_col_name_as, "' processed."))}
      } else {
        if (nzchar(tumor_col_name_as) && tumor_col_exists_as && length(unique_tumor_values_as) == 0) { message(paste0("Processing full dataset as '", tumor_col_name_as, "' column has no distinct non-NA values."))
        } else if (!nzchar(tumor_col_name_as) || !tumor_col_exists_as) { message(paste0("Processing full dataset as 'proximity_split' column ('", tumor_col_name_as, "') not provided or not found.")) }
        proximity_output_full_as <- NULL
        if (nrow(original_metadata_as) > 0 && cell_type_column %in% names(original_metadata_as) && all(c("x", "y") %in% names(original_metadata_as))) {
          proximity_output_full_as <- tryCatch({ proximity_metabolic_score_pvalue_perm_metab( k=k, meta_from_seurat = original_metadata_as, cell_type_column = cell_type_column, list_metabolic = current_selected_cols, sep = patient_column ) }, error = function(e){ warning(paste("Error in prox function for full dataset: ", e$message)); NULL })
        }
        if(!is.null(proximity_output_full_as) && length(proximity_output_full_as) > 0){ results_payload <- list(type = "full_dataset", results = proximity_output_full_as, selected_pathways = current_selected_cols)
        } else {showNotification("Could not process full dataset for 'All Samples Interaction'.", type = "error", duration=5); results_payload <- list(type="error", message="Full dataset processing failed or empty.")}
      }
      all_samples_processed_cache_rv(results_payload) # Update the reactiveVal
    })
    
    output$all_samples_heatmaps_ui <- renderUI({
      all_samples_info <- all_samples_processed_cache_rv() # Depend on the reactiveVal
      req(all_samples_info) 
      if(identical(all_samples_info$type, "error")) { return(tags$p(style="color:red;", all_samples_info$message)) }
      req(all_samples_info$selected_pathways, length(all_samples_info$selected_pathways) > 0)
      
      ui_elements_all_samples <- list()
      if (all_samples_info$type == "split_by_group") {
        req(all_samples_info$results, all_samples_info$group_names)
        if(length(all_samples_info$results) == 0) return(tags$p(paste0("No '", all_samples_info$group_column_name,"' groups processed.")))
        for (group_name_iter_ui in all_samples_info$group_names) {
          if (!is.null(all_samples_info$results[[group_name_iter_ui]]) && length(all_samples_info$results[[group_name_iter_ui]]) > 0) {
            ui_elements_all_samples[[length(ui_elements_all_samples) + 1]] <- h4(paste("Interaction Heatmaps for ",all_samples_info$group_column_name,":", group_name_iter_ui))
            for (pathway_name_iter_ui in all_samples_info$selected_pathways) {
              if (!is.null(all_samples_info$results[[group_name_iter_ui]][[pathway_name_iter_ui]])) {
                ui_elements_all_samples[[length(ui_elements_all_samples) + 1]] <- plotOutput(outputId = paste0("all_samples_plot_", make.names(group_name_iter_ui), "_", make.names(pathway_name_iter_ui)), height = "500px")
                ui_elements_all_samples[[length(ui_elements_all_samples) + 1]] <- hr()
              } else { ui_elements_all_samples[[length(ui_elements_all_samples) + 1]] <- tags$p(paste("No data for pathway", pathway_name_iter_ui, "in group", group_name_iter_ui)) }
            }
          }
        }
      } else if (all_samples_info$type == "full_dataset") {
        req(all_samples_info$results); if(length(all_samples_info$results) == 0 && length(all_samples_info$selected_pathways) > 0) { return(tags$p("Proximity function returned no results for selected pathways on full dataset.")) }
        ui_elements_all_samples[[length(ui_elements_all_samples) + 1]] <- h4("Interaction Heatmaps for Full Dataset")
        for (pathway_name_iter_ui_full in all_samples_info$selected_pathways) {
          if (!is.null(all_samples_info$results[[pathway_name_iter_ui_full]])) {
            ui_elements_all_samples[[length(ui_elements_all_samples) + 1]] <- plotOutput(outputId = paste0("all_samples_plot_full_", make.names(pathway_name_iter_ui_full)), height = "500px")
            ui_elements_all_samples[[length(ui_elements_all_samples) + 1]] <- hr()
          } else { ui_elements_all_samples[[length(ui_elements_all_samples) + 1]] <- tags$p(paste("No data for pathway", pathway_name_iter_ui_full, "in full dataset.")) }
        }
      }
      if(length(ui_elements_all_samples) > 0) { do.call(tagList, ui_elements_all_samples) } 
      else { tags$p("No heatmaps for 'All Samples Score Interaction'. Check selections and data.") }
    })
    
    observeEvent(all_samples_processed_cache_rv(), { # Changed from all_samples_data_rval
      processed_all_samples_data <- all_samples_processed_cache_rv() # Use the reactiveVal
      req(processed_all_samples_data)
      if(identical(processed_all_samples_data$type, "error")) { return() }
      req(processed_all_samples_data$selected_pathways)
      
      if (processed_all_samples_data$type == "split_by_group") {
        req(processed_all_samples_data$results, processed_all_samples_data$group_names)
        for (group_name_obs in processed_all_samples_data$group_names) {
          tumor_specific_results <- processed_all_samples_data$results[[group_name_obs]]
          if (!is.null(tumor_specific_results)) {
            for (pathway_name_obs in processed_all_samples_data$selected_pathways) {
              local({ 
                current_group_name_local <- group_name_obs; current_pathway_name_local <- pathway_name_obs
                plot_id <- paste0("all_samples_plot_", make.names(current_group_name_local), "_", make.names(current_pathway_name_local))
                output[[plot_id]] <- renderPlot({
                  pathway_data_list <- processed_all_samples_data$results[[current_group_name_local]][[current_pathway_name_local]]
                  req(pathway_data_list, pathway_data_list[[1]]); matrix_to_plot <- as.matrix(pathway_data_list[[1]])
                  annotation_to_plot <- as.matrix(pathway_data_list[["adj_p_value"]])
                  if (nrow(matrix_to_plot) == 0 || ncol(matrix_to_plot) == 0) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, paste("Matrix empty for", current_pathway_name_local, "\nin group", current_group_name_local), col="orange"); return() }
                  heatmap_title <- paste("Score Interacted:", current_pathway_name_local, "\n(",processed_all_samples_data$group_column_name,":", current_group_name_local, ")"); tryCatch({ ht <- Heatmap(matrix_to_plot, cluster_rows = FALSE, cluster_columns = FALSE, name = "Score", column_title = heatmap_title, column_title_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
                                                                                                                                                                                                           cell_fun = function(j, i, x, y, width, height, fill, p_val_mat = annotation_to_plot) {
                                                                                                                                                                                                             # Check the value from text_matrix to decide the color
                                                                                                                                                                                                             text_color <- if (p_val_mat[i, j] < 0.05) {
                                                                                                                                                                                                               "green"
                                                                                                                                                                                                             } else {
                                                                                                                                                                                                               "black"
                                                                                                                                                                                                             }
                                                                                                                                                                                                             
                                                                                                                                                                                                             # Use text_matrix[i, j] for the label and the determined color for the text
                                                                                                                                                                                                             grid.text(
                                                                                                                                                                                                               sprintf("%.3f", p_val_mat[i, j]), 
                                                                                                                                                                                                               x, y, 
                                                                                                                                                                                                               gp = gpar(fontsize = 10, col = text_color)
                                                                                                                                                                                                             )
                                                                                                                                                                                                           }); draw(ht)
                  }, error = function(e) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, paste("Error plotting heatmap for", current_pathway_name_local, "\nin group", current_group_name_local, ":\n", e$message), cex = 0.9, col="darkred") }) }) }) } } }
      } else if (processed_all_samples_data$type == "full_dataset") {
        req(processed_all_samples_data$results); full_dataset_results <- processed_all_samples_data$results
        for (pathway_name_obs_full in processed_all_samples_data$selected_pathways) {
          local({ 
            current_pathway_name_full_local <- pathway_name_obs_full
            plot_id_full <- paste0("all_samples_plot_full_", make.names(current_pathway_name_full_local))
            output[[plot_id_full]] <- renderPlot({
              pathway_data_list_full <- full_dataset_results[[current_pathway_name_full_local]]
              req(pathway_data_list_full, pathway_data_list_full[[1]]); matrix_to_plot_full <- as.matrix(pathway_data_list_full[[1]])
              annotation_to_full_plot <- as.matrix(pathway_data_list_full[["adj_p_value"]])
              if (nrow(matrix_to_plot_full) == 0 || ncol(matrix_to_plot_full) == 0) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, paste("Matrix empty for", current_pathway_name_full_local, "\nin full dataset"), col="orange"); return() }
              heatmap_title_full <- paste("Score Interacted:", current_pathway_name_full_local, "\n(Full Dataset)"); tryCatch({ ht_full <- Heatmap(matrix_to_plot_full, cluster_rows = FALSE, cluster_columns = FALSE, name = "Score", column_title = heatmap_title_full, column_title_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
                                                                                                                                                   cell_fun = function(j, i, x, y, width, height, fill, p_val_mat = annotation_to_full_plot) {
                                                                                                                                                     # Check the value from text_matrix to decide the color
                                                                                                                                                     text_color <- if (p_val_mat[i, j] < 0.05) {
                                                                                                                                                       "green"
                                                                                                                                                     } else {
                                                                                                                                                       "black"
                                                                                                                                                     }
                                                                                                                                                     
                                                                                                                                                     # Use text_matrix[i, j] for the label and the determined color for the text
                                                                                                                                                     grid.text(
                                                                                                                                                       sprintf("%.3f", p_val_mat[i, j]), 
                                                                                                                                                       x, y, 
                                                                                                                                                       gp = gpar(fontsize = 10, col = text_color)
                                                                                                                                                     )
                                                                                                                                                   }); draw(ht_full)
              }, error = function(e) { plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, paste("Error plotting heatmap for", current_pathway_name_full_local, "\n(Full Dataset):\n", e$message), cex = 0.9, col="darkred") }) }) }) } }
    }, ignoreNULL = TRUE, ignoreInit = TRUE) 
    # --- END NEW observeEvent ---
    
  }
  
  shinyApp(ui = ui, server = server)
  
}
my_shiny_function_all_samples_color_pval(metadata = metadata_tirosh_allGBM_score_cell_type, patient_column = "GSM", 
                                              list_of_scores_columns = pathwaylist,
                                              cell_type_column = "mp",
                                              proximity_split = "", k=4
)
