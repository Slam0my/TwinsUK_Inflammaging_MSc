#RShiny Web for Protein-Protein network graphs
#Contains genetically and environmentally significant Protein-Protein pairs in a network
#You can select and filter parameters with the sliders and dropdowns to observe specific modules
#Network graphs are based on Leiden clustering for modules and uses the force directed layout from ForceAtlas2
#Weights are based on |ρG| or |ρE|


#Load libraries
if (!require("shiny")) install.packages("shiny")
if (!require("shinydashboard")) install.packages("shinydashboard")
if (!require("visNetwork")) install.packages("visNetwork")
if (!require("visNetwork")) install.packages("visNetwork")
if (!require("dplyr")) install.packages("dplyr")
if (!require("igraph")) install.packages("igraph")
if (!require("viridis")) install.packages("viridis")

library(shiny)
library(shinydashboard)
library(dplyr)
library(visNetwork) #Network graph
library(igraph) #Network graph
library(viridis) #Colour


#Load data
#Loading in the protein-protein solar analysis results 
#Make sure this RData is within the same file
load("prot_prot_filtering_results.RData")


#Static clusters based on Leiden
#static for stability
#Genetics
g_base_gen <- graph_from_data_frame(prot_gen_0.1_hits %>% filter(abs(rhoG) >= 0.3) %>% mutate(weight = abs(rhoG)) %>% select(trait1, trait2, weight), directed = FALSE)
leiden_gen <- cluster_leiden(g_base_gen, objective_function = "modularity", weights = E(g_base_gen)$weight)
gen_module_lookup <- data.frame(name = V(g_base_gen)$name, module = as.character(leiden_gen$membership))

#Environmental
g_base_env <- graph_from_data_frame(prot_env_hits %>% filter(abs(rhoE) >= 0.3) %>% mutate(weight = abs(rhoE)) %>% select(trait1, trait2, weight), directed = FALSE)
leiden_env <- cluster_leiden(g_base_env, objective_function = "modularity", weights = E(g_base_env)$weight)
env_module_lookup <- data.frame(name = V(g_base_env)$name, module = as.character(leiden_env$membership))

#Module colours the same 
all_possible_modules <- sort(unique(c(gen_module_lookup$module, env_module_lookup$module)))
global_module_colours <- setNames(viridis_pal(option = "viridis")(length(all_possible_modules)), all_possible_modules)


#Key/legend
build_custom_legend <- function(is_genetic = TRUE) {
  link_title <- ifelse(is_genetic, "Links (ρG)", "Links (ρE)")
  active_modules <- if(is_genetic) sort(unique(gen_module_lookup$module)) else sort(unique(env_module_lookup$module))
  active_colors <- global_module_colours[active_modules]
  
  module_items <- paste0(
    '<div style="display: flex; align-items: center; margin-bottom: 8px;">',
    '<span style="display: inline-block; width: 18px; height: 18px; border-radius: 50%; background-color: ', unname(active_colors), '; margin-right: 10px; flex-shrink: 0;"></span>',
    '<span style="font-size: 15px; line-height: 1.2;">Module ', names(active_colors), '</span>',
    '</div>', collapse = ""
  )
  
  HTML(paste0(
    '<div style="padding: 10px; background-color: #f9f9f9; border-radius: 5px; height: 710px; overflow-y: auto; border: 1px solid #ddd;">',
    
    #Node size and centrality 
    '<h4 style="text-align: center; font-weight: bold; margin-top: 0; font-size: 20px;">Nodes & Centrality</h4>',
    '<div style="text-align: center; font-size: 14px; margin-bottom: 15px;">',
    '<p style="margin-bottom: 8px;">Node Size = Total Number of Connections</p>',
    '<div style="display: flex; align-items: flex-end; justify-content: center; gap: 20px; height: 60px;">',
    '<div style="text-align: center;"><div style="width: 10px; height: 10px; background-color: #7d7d7d; border-radius: 50%; margin: 0 auto 5px auto;"></div>1</div>',
    '<div style="text-align: center;"><div style="width: 20px; height: 20px; background-color: #7d7d7d; border-radius: 50%; margin: 0 auto 5px auto;"></div>5</div>',
    '<div style="text-align: center;"><div style="width: 35px; height: 35px; background-color: #7d7d7d; border-radius: 50%; margin: 0 auto 5px auto;"></div>10+</div>',
    '</div>',
    '</div><hr style="border-top: 2px solid #ddd;">',
    
    #Links
    '<h4 style="text-align: center; font-weight: bold; font-size: 20px;">', link_title, '</h4>',
    '<div style="font-size: 14px; text-align: center; margin-bottom: 15px;">',
    '<p style="margin-bottom: 10px;">Line Thickness = Correlation Strength</p>',
    '<div style="display: flex; align-items: center; justify-content: center; gap: 15px; margin-bottom: 15px;">',
    '<div style="text-align: center;"><div style="width: 30px; height: 1px; background-color: #7d7d7d; margin: 0 auto 5px auto;"></div>Weak</div>',
    '<div style="text-align: center;"><div style="width: 30px; height: 4px; background-color: #7d7d7d; margin: 0 auto 5px auto;"></div>Med</div>',
    '<div style="text-align: center;"><div style="width: 30px; height: 8px; background-color: #7d7d7d; margin: 0 auto 5px auto;"></div>Strong</div>',
    '</div>',
    '<div style="margin-bottom: 8px; font-size: 15px;"><b>&mdash;&mdash;&mdash;</b> Solid: Positive (+ρ)</div>',
    '<div style="font-size: 15px;"><b>- - - -</b> Dotted: Negative (-ρ)</div>',
    '</div><hr style="border-top: 2px solid #ddd;">',
    
    
    #Modules
    '<h4 style="text-align: center; font-weight: bold; margin-bottom: 15px; font-size: 20px;">Leiden Modules</h4>',
    '<div style="column-count: 2; column-gap: 20px;">', module_items, '</div>',
    '</div>'
  ))
}




#UI
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "Protein-Protein Networks"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Genetically Significant", tabName = "gen_pairs", icon = icon("dna")),
      menuItem("Environmentally Significant", tabName = "env_pairs", icon = icon("leaf"))
    )
  ),
  
  dashboardBody(
    tabItems(
      #GENETICS
      tabItem(tabName = "gen_pairs",
              h2("Genetically Significant Protein-Protein Pairs", align = "center"),
              p("Analysis based on shared genetic architecture (ρG). This network is based on Leiden clustering with weights using |ρG|. ForceAtlas2 is used as a force directed layout. Please bare with the time for loading graphs. For these genetically signficant pairs, h² > 0.1 & ρG FDR < 0.05 was implemented. The maximum |ρG| between a pair reached 0.9. ", align = "center"),
              #Sliders
              fluidRow(
                box(title = "Filtering:", width = 4, status = "primary", solidHeader = TRUE,
                    sliderInput("gen_sharing", "Min Connections per Protein:", min = 1, max = 15, value = 3, step = 1),
                    sliderInput("gen_thresh", "Min |ρG| Strength:", min = 0.0, max = 0.9, value = 0.4, step = 0.05)),
                box(title = "Selection:", width = 8, status = "primary", solidHeader = TRUE,
                    fluidRow(
                      column(6, selectInput("gen_mod", "Select Leiden Module:", choices = "All")),
                      column(6, selectInput("gen_prot", "Select Specific Protein:", choices = "All"))
                    ), hr(), htmlOutput("gen_module_summary"))
              ),
              #Dropdowns
              fluidRow(
                box(width = 9, status = "primary", visNetworkOutput("network_graph_gen", height = "750px")),
                box(title = "Network Key", width = 3, status = "primary", solidHeader = TRUE, build_custom_legend(TRUE))
              )
      ),
      
      #ENVIRONMENTAL
      tabItem(tabName = "env_pairs",
              h2("Environmentally Significant Protein-Protein Pairs", align = "center"),
              p("Analysis based on environmental architecture (ρE). This network is based on Leiden clustering with weights using |ρE|. ForceAtlas2 is used as a force directed layout. Please bare with the time for loading graphs. For these environmentally signficant pairs, ρE FDR < 0.05 was implemented. The maximum |ρE| between a pair reached 0.5. ", align = "center"),
              #Sliders
              fluidRow(
                box(title = "Filtering:", width = 4, status = "success", solidHeader = TRUE,
                    sliderInput("env_sharing", "Min Connections per Protein:", min = 1, max = 15, value = 3, step = 1),
                    sliderInput("env_thresh", "Min |ρE| Strength:", min = 0.0, max = 0.5, value = 0.4, step = 0.05)),
                box(title = "Selection:", width = 8, status = "success", solidHeader = TRUE,
                    fluidRow(
                      column(6, selectInput("env_mod", "Select Leiden Module:", choices = "All")),
                      column(6, selectInput("env_prot", "Select Specific Protein:", choices = "All"))
                    ), hr(), htmlOutput("env_module_summary"))
              ),
              #Dropdowns
              fluidRow(
                box(width = 9, status = "success", visNetworkOutput("network_graph_env", height = "750px")),
                box(title = "Network Key", width = 3, status = "success", solidHeader = TRUE, build_custom_legend(FALSE))
              )
      )
    )
  )
)

#Server

server <- function(input, output, session) {
  
  #REACTIVES
  #GENETIC NETWORK
  #Reactive calculations based on input slider
  network_data_gen <- reactive({
    req(input$gen_thresh, input$gen_sharing)
    
    hub_data_all <- prot_gen_0.1_hits %>%
      filter(abs(rhoG) >= input$gen_thresh) %>% 
      mutate(abs_rhoG = abs(rhoG), edge_sign = ifelse(rhoG >= 0, "Positive", "Negative"))
    
    #Fail safe returns 0 
    if(nrow(hub_data_all) == 0) return(NULL)
    
    #Centrality 
    node_counts <- data.frame(node = c(hub_data_all$trait1, hub_data_all$trait2)) %>% count(node, name = "original_centrality")
    
    #nodes based off slider
    valid_nodes <- node_counts %>% filter(original_centrality >= input$gen_sharing) %>% pull(node)
    hub_data_filtered <- hub_data_all %>% filter(trait1 %in% valid_nodes & trait2 %in% valid_nodes)
    if(nrow(hub_data_filtered) == 0) return(NULL)
    
    nodes_df <- data.frame(name = unique(c(hub_data_filtered$trait1, hub_data_filtered$trait2))) %>% 
      left_join(gen_module_lookup, by = "name") %>%
      left_join(node_counts, by = c("name" = "node")) %>% 
      mutate(module = ifelse(is.na(module), "Unassigned", module))
    
    g <- graph_from_data_frame(d = hub_data_filtered %>% select(trait1, trait2, rhoG, abs_rhoG, edge_sign), directed = FALSE, vertices = nodes_df)
    
    #the centrality present 
    V(g)$visible_centrality <- degree(g)
    
    vis_nodes <- igraph::as_data_frame(g, what = "vertices") %>% rename(id = name) %>%
      mutate(label = id, 
             value = original_centrality, 
             color = ifelse(module %in% names(global_module_colours), global_module_colours[module], "#cccccc"), 
             original_color = color, 
             title = paste0("<p><b>", id, "</b><br>Module: ", module, 
                            "<br>Total Connections (≥ ", input$gen_thresh, " |ρG|): <b>", original_centrality, "</b>",
                            "<br>Visible Hub Connections: <b>", visible_centrality, "</b></p>"))
    
    
    vis_edges <- igraph::as_data_frame(g, what = "edges") %>% rename(from = from, to = to) %>%
      mutate(id = row_number(), value = abs_rhoG, dashes = rhoG < 0, color = "#a9a9a9", original_color = color,
             title = paste0("<p><b>ρG:</b> ", round(rhoG, 3), "<br><b>Direction:</b> ", edge_sign, "</p>"))
    list(nodes = vis_nodes, edges = vis_edges)
  })
  
  output$network_graph_gen <- renderVisNetwork({
    data <- network_data_gen(); req(data)
    visNetwork(data$nodes, data$edges) %>%
      visPhysics(solver = "forceAtlas2Based", forceAtlas2Based = list(gravitationalConstant = -200, springLength = 150), stabilization = list(enabled = TRUE, iterations = 1000)) %>%
      visEvents(stabilizationIterationsDone = "function () {this.setOptions({physics: false});}") %>%
      visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE), nodesIdSelection = list(enabled = TRUE, main = "Search by Name")) %>%
      visEdges(smooth = list(enabled = TRUE, type = "continuous")) %>% visInteraction(navigationButtons = TRUE)
  })
  
  
  
  
  
  
  
  #ENVIRONMENT NETWORK
  network_data_env <- reactive({
    req(input$env_thresh, input$env_sharing)
    
    hub_data_all <- prot_env_hits %>%
      filter(abs(rhoE) >= input$env_thresh) %>% 
      mutate(abs_rhoE = abs(rhoE), edge_sign = ifelse(rhoE >= 0, "Positive", "Negative"))
    
    if(nrow(hub_data_all) == 0) return(NULL)
    
    #Centrality calculations 
    node_counts <- data.frame(node = c(hub_data_all$trait1, hub_data_all$trait2)) %>% count(node, name = "original_centrality")
    
    #Nodes from slider 
    valid_nodes <- node_counts %>% filter(original_centrality >= input$env_sharing) %>% pull(node)
    hub_data_filtered <- hub_data_all %>% filter(trait1 %in% valid_nodes & trait2 %in% valid_nodes)
    if(nrow(hub_data_filtered) == 0) return(NULL)
    
    nodes_df <- data.frame(name = unique(c(hub_data_filtered$trait1, hub_data_filtered$trait2))) %>% 
      left_join(env_module_lookup, by = "name") %>%
      left_join(node_counts, by = c("name" = "node")) %>% 
      mutate(module = ifelse(is.na(module), "Unassigned", module))
    
    g <- graph_from_data_frame(d = hub_data_filtered %>% select(trait1, trait2, rhoE, abs_rhoE, edge_sign), directed = FALSE, vertices = nodes_df)
    
    
    #centrality based on what is visible 
    V(g)$visible_centrality <- degree(g)
    
    vis_nodes <- igraph::as_data_frame(g, what = "vertices") %>% rename(id = name) %>%
      mutate(label = id, value = original_centrality,
             color = ifelse(module %in% names(global_module_colours), global_module_colours[module], "#cccccc"),
             original_color = color, 
             title = paste0("<p><b>", id, "</b><br>Module: ", module, 
                            "<br>Total Connections (≥ ", input$env_thresh, " |ρE|): <b>", original_centrality, "</b>",
                            "<br>Visible Hub Connections: <b>", visible_centrality, "</b></p>"))
    
    vis_edges <- igraph::as_data_frame(g, what = "edges") %>% rename(from = from, to = to) %>%
      mutate(id = row_number(), value = abs_rhoE, dashes = rhoE < 0, color = "#a9a9a9", original_color = color,
             title = paste0("<p><b>ρE:</b> ", round(rhoE, 3), "<br><b>Direction:</b> ", edge_sign, "</p>"))
    list(nodes = vis_nodes, edges = vis_edges)
  })
  
  
  output$network_graph_env <- renderVisNetwork({
    data <- network_data_env(); req(data)
    visNetwork(data$nodes, data$edges) %>%
      visPhysics(solver = "forceAtlas2Based", forceAtlas2Based = list(gravitationalConstant = -200, springLength = 150), stabilization = list(enabled = TRUE, iterations = 1000)) %>%
      visEvents(stabilizationIterationsDone = "function () {this.setOptions({physics: false});}") %>%
      visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE), nodesIdSelection = list(enabled = TRUE, main = "Search by Name")) %>%
      visEdges(smooth = list(enabled = TRUE, type = "continuous")) %>% visInteraction(navigationButtons = TRUE)
  })
  
  #Shared drpdowns/selections
  observe({
    data_g <- network_data_gen(); req(data_g)
    mod <- input$gen_mod; prot <- input$gen_prot
    v_mod <- data_g$nodes; if(prot != "All") v_mod <- v_mod %>% filter(id %in% c((data_g$edges %>% filter(from == prot | to == prot))$from, (data_g$edges %>% filter(from == prot | to == prot))$to))
    v_prot <- data_g$nodes; if(mod != "All") v_prot <- v_prot %>% filter(module == mod)
    updateSelectInput(session, "gen_mod", choices = c("All", sort(unique(v_mod$module))), selected = ifelse(mod %in% v_mod$module, mod, "All"))
    updateSelectInput(session, "gen_prot", choices = c("All", sort(unique(v_prot$id))), selected = ifelse(prot %in% v_prot$id, prot, "All"))
    
    data_e <- network_data_env(); req(data_e)
    mod_e <- input$env_mod; prot_e <- input$env_prot
    v_mod_e <- data_e$nodes; if(prot_e != "All") v_mod_e <- v_mod_e %>% filter(id %in% c((data_e$edges %>% filter(from == prot_e | to == prot_e))$from, (data_e$edges %>% filter(from == prot_e | to == prot_e))$to))
    v_prot_e <- data_e$nodes; if(mod_e != "All") v_prot_e <- v_prot_e %>% filter(module == mod_e)
    updateSelectInput(session, "env_mod", choices = c("All", sort(unique(v_mod_e$module))), selected = ifelse(mod_e %in% v_mod_e$module, mod_e, "All"))
    updateSelectInput(session, "env_prot", choices = c("All", sort(unique(v_prot_e$id))), selected = ifelse(prot_e %in% v_prot_e$id, prot_e, "All"))
  })
  
  
  
  #grey out gen
  observe({
    data <- network_data_gen(); req(data, input$gen_mod, input$gen_prot)
    a_nodes <- data$nodes; if(input$gen_mod != "All") a_nodes <- a_nodes %>% filter(module == input$gen_mod); if(input$gen_prot != "All") a_nodes <- a_nodes %>% filter(id == input$gen_prot)
    t_ids <- if(input$gen_mod == "All" && input$gen_prot == "All") data$nodes$id else a_nodes$id
    
    if(input$gen_prot != "All") {
      a_edges <- data$edges %>% filter(from == input$gen_prot | to == input$gen_prot)
      t_ids <- unique(c(a_edges$from, a_edges$to))
    }
    
    t_e_ids <- data$edges %>% filter(from %in% t_ids & to %in% t_ids) %>% pull(id)
    if(input$gen_mod == "All" && input$gen_prot == "All") t_e_ids <- data$edges$id
    
    visNetworkProxy("network_graph_gen") %>%
      visUpdateNodes(nodes = data$nodes %>% select(id, original_color) %>% mutate(color = ifelse(id %in% t_ids, original_color, "rgba(200,200,200,0.1)"))) %>%
      visUpdateEdges(edges = data$edges %>% select(id, original_color) %>% mutate(color = ifelse(id %in% t_e_ids, original_color, "rgba(200,200,200,0.05)")))
  })
  
  
  
  
  #grey out env
  observe({
    data <- network_data_env(); req(data, input$env_mod, input$env_prot)
    a_nodes <- data$nodes; if(input$env_mod != "All") a_nodes <- a_nodes %>% filter(module == input$env_mod); if(input$env_prot != "All") a_nodes <- a_nodes %>% filter(id == input$env_prot)
    t_ids <- if(input$env_mod == "All" && input$env_prot == "All") data$nodes$id else a_nodes$id
    
    if(input$env_prot != "All") {
      a_edges <- data$edges %>% filter(from == input$env_prot | to == input$env_prot)
      t_ids <- unique(c(a_edges$from, a_edges$to))
    }
    
    t_e_ids <- data$edges %>% filter(from %in% t_ids & to %in% t_ids) %>% pull(id)
    if(input$env_mod == "All" && input$env_prot == "All") t_e_ids <- data$edges$id
    
    visNetworkProxy("network_graph_env") %>%
      visUpdateNodes(nodes = data$nodes %>% select(id, original_color) %>% mutate(color = ifelse(id %in% t_ids, original_color, "rgba(200,200,200,0.1)"))) %>%
      visUpdateEdges(edges = data$edges %>% select(id, original_color) %>% mutate(color = ifelse(id %in% t_e_ids, original_color, "rgba(200,200,200,0.05)")))
  })
  
  
  #sumamries for modules
  output$gen_module_summary <- renderUI({
    data <- network_data_gen(); req(data); mod <- input$gen_mod
    if(mod == "All") return(HTML("<p style='color: #555; font-style: italic;'>Select a Module to see its proteins.</p>"))
    mod_prots <- data$nodes %>% filter(module == mod) %>% arrange(desc(value)) 
    if(nrow(mod_prots) == 0) return(HTML("<p>No proteins currently in this module.</p>"))
    
    HTML(paste0(
      "<div style='font-size: 15px;'>",
      "<span style='color: #1e88e5;'><b>Total Proteins in Module ", mod, " (", nrow(mod_prots), "):</b></span><br>", 
      paste(mod_prots$id, collapse = ", "),
      "</div>"
    ))
  })
  
  output$env_module_summary <- renderUI({
    data <- network_data_env(); req(data); mod <- input$env_mod
    if(mod == "All") return(HTML("<p style='color: #555; font-style: italic;'>Select a Module to see its proteins.</p>"))
    mod_prots <- data$nodes %>% filter(module == mod) %>% arrange(desc(value))
    if(nrow(mod_prots) == 0) return(HTML("<p>No proteins currently in this module.</p>"))
    
    HTML(paste0(
      "<div style='font-size: 15px;'>",
      "<span style='color: #00897b;'><b>Total Proteins in Module ", mod, " (", nrow(mod_prots), "):</b></span><br>", 
      paste(mod_prots$id, collapse = ", "),
      "</div>"
    ))
  })
}


shinyApp(ui = ui, server = server)