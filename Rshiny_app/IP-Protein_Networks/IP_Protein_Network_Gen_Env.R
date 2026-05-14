#RShiny Web for IP-Protein network graphs
#Contains genetically and environmentally significant IP-Protein pairs in a network
#You can select and filter parameters with the sliders and dropdowns to observe specific modules, lineages and proteins
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
#Loading in the solar results analysis filtering.
#Make sure this RData is within the same file
load("solar_analysis_filtering_results.RData")


#Lineage colour palette for legend/key
all_possible_lineages <- unique(c(genetic_0.1_hits$Sub.Lineage, environmental_hits$Sub.Lineage))
all_possible_lineages <- sort(all_possible_lineages[!is.na(all_possible_lineages)])
global_lineage_colours <- setNames(viridis_pal(option = "viridis")(length(all_possible_lineages)), all_possible_lineages)


#Static clusters based on Leiden
#static for stability 
#Genetic
g_base_gen <- graph_from_data_frame(genetic_0.1_hits %>% mutate(abs_rhoG = abs(rhoG)) %>% select(trait1, trait2, abs_rhoG), directed = FALSE)
leiden_gen <- cluster_leiden(g_base_gen, objective_function = "modularity", weights = E(g_base_gen)$abs_rhoG)
gen_module_lookup <- data.frame(name = V(g_base_gen)$name, module = as.character(leiden_gen$membership))

#Environmental
g_base_env <- graph_from_data_frame(environmental_hits %>% mutate(abs_rhoE = abs(rhoE)) %>% select(trait1, trait2, abs_rhoE), directed = FALSE)
leiden_env <- cluster_leiden(g_base_env, objective_function = "modularity", weights = E(g_base_env)$abs_rhoE)
env_module_lookup <- data.frame(name = V(g_base_env)$name, module = as.character(leiden_env$membership))


#Key/legend
build_custom_legend <- function(is_genetic = TRUE) {
  link_title <- ifelse(is_genetic, "Links (ρG)", "Links (ρE)")
  
  #Legend appears based on gen or env
  if(is_genetic) {
    active_lineages <- sort(unique(genetic_0.1_hits$Sub.Lineage))
  } else {
    active_lineages <- sort(unique(environmental_hits$Sub.Lineage))
  }
  active_lineages <- active_lineages[!is.na(active_lineages)]
  
  #colours for lineages
  active_colors <- global_lineage_colours[active_lineages]
  
  #Listing key/legend
  lineage_items <- paste0(
    '<div style="display: flex; align-items: center; margin-bottom: 8px;">',
    '<span style="display: inline-block; width: 18px; height: 18px; border-radius: 50%; background-color: ', unname(active_colors), '; margin-right: 10px; flex-shrink: 0;"></span>',
    '<span style="font-size: 16px; line-height: 1.2;">', names(active_colors), '</span>',
    '</div>',
    collapse = ""
  )
  
  HTML(paste0(
    '<div style="padding: 8px; background-color: #f9f9f9; border-radius: 2px; height: 710px; overflow-y: auto; border: 1px solid #ddd;">',
    
    #Node 
    '<h4 style="text-align: center; font-weight: bold; margin-top: 0; font-size: 20px;">Nodes</h4>',
    '<div style="display: flex; justify-content: space-around; margin-bottom: 15px; font-size: 16px;">',
    '<div style="text-align: center;"><div style="width: 18px; height: 18px; background-color: black; transform: rotate(45deg); margin: 0 auto 5px auto;"></div>Protein Hub</div>',
    '<div style="text-align: center;"><div style="width: 20px; height: 20px; background-color: #7d7d7d; border-radius: 50%; margin: 0 auto 5px auto;"></div>IP Node</div>',
    '</div>',
    '<hr style="border-top: 2px solid #ddd;">',
    
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
    
    
    
    #Lineages
    '<h4 style="text-align: center; font-weight: bold; margin-bottom: 15px; font-size: 20px;">Lineages</h4>',
    '<div style="column-count: 2; column-gap: 20px;">',
    lineage_items,
    '</div>',
    
    '</div>'
  ))
}


#UI
ui <- dashboardPage(
  skin = "purple",
  
  dashboardHeader(title = "IP-Protein Networks"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Genetic", tabName = "gen_pairs", icon = icon("dna")),
      menuItem("Environmental", tabName = "env_pairs", icon = icon("leaf"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .main-sidebar .sidebar .sidebar-menu a {
          font-size: 18px;   
          padding: 20px 15px; 
        }
        .main-sidebar .sidebar .sidebar-menu a i {
          font-size: 20px;
          margin-right: 10px;
        }
      "))
    ),
    
    tabItems(
      #GENETICS
      tabItem(tabName = "gen_pairs",
              h2("Genetically Significant IP-Protein Pairs", align = "center"),
              p("Analysis based on shared genetic architecture (ρG) between Immunophenotypes (IPs) and Circulating Proteins. This network is based on Leiden clustering with weights using |ρG|. ForceAtlas2 is used as a force directed layout. Please bare with the time for loading graphs. For these genetically signficant pairs, h² > 0.1 & ρG FDR < 0.05 was implemented. The maximum |ρG| between a pair reached 0.83. The maximum number of proteins shared by IPs was 9.", align = "center"),
              
              fluidRow(
                #Sliders
                box(title = "Filtering:", 
                    width = 4, 
                    status = "primary", 
                    solidHeader = TRUE,
                    sliderInput("gen_sharing", 
                                "Proteins per IP:", 
                                min = 1, 
                                max = 9, 
                                value = 2, 
                                step = 1),
                    sliderInput("gen_thresh", 
                                "Min |ρG| Strength:", 
                                min = 0.0, 
                                max = 0.8, 
                                value = 0.3,
                                step = 0.05)
                ),
                #Dropdowns
                box(title = "Selection:", 
                    width = 8, 
                    status = "primary", 
                    solidHeader = TRUE,
                    fluidRow(
                      column(4, selectInput("gen_mod", "Select Leiden Module:", choices = "All")),
                      column(4, selectInput("gen_lin", "Select Immune Lineage:", choices = "All")),
                      column(4, selectInput("gen_prot", "Select Protein Hub:", choices = "All"))
                    ),
                    hr(),
                    htmlOutput("gen_module_summary"))
              ),
              fluidRow(
                box(width = 9, status = "primary", visNetworkOutput("network_graph_gen", height = "800px")),
                box(title = "Key", width = 3, status = "primary", solidHeader = TRUE, build_custom_legend(TRUE))
              )
      ),
      
      #ENVIRONMENTAL
      tabItem(tabName = "env_pairs",
              h2("Environmentally Significant IP-Protein Pairs", align = "center"),
              p("Analysis based on environmental architecture (ρE). This network is based on Leiden clustering with weights using |ρE|. ForceAtlas2 is used as a force directed layout. Please bare with the time for loading graphs. For these environmentally signficant pairs, ρE FDR < 0.05 was implemented. The maximum |ρE| between a pair reached 0.45. The maximum number of proteins shared by IPs was 6.", align = "center"),
              fluidRow(
                #Sliders
                box(title = "Filtering:", 
                    width = 4, 
                    status = "success", 
                    solidHeader = TRUE,
                    sliderInput("env_sharing", 
                                "Proteins per IP:", 
                                min = 1, 
                                max = 6, 
                                value = 2, 
                                step = 1),
                    sliderInput("env_thresh", 
                                "Min |ρE| Strength:", 
                                min = 0.0, 
                                max = 0.45, 
                                value = 0.3, 
                                step = 0.05)
                ),
                #Dropdowns
                box(title = "Selection:", 
                    width = 8, 
                    status = "success", 
                    solidHeader = TRUE,
                    fluidRow(
                      column(4, selectInput("env_mod", "Select Leiden Module:", choices = "All")),
                      column(4, selectInput("env_lin", "Select Immune Lineage:", choices = "All")),
                      column(4, selectInput("env_prot", "Select Protein Hub:", choices = "All"))
                    ),
                    hr(),
                    htmlOutput("env_module_summary"))
              ),
              fluidRow(
                box(width = 9, status = "success", visNetworkOutput("network_graph_env", height = "800px")),
                box(title = "Key", width = 3, status = "success", solidHeader = TRUE, build_custom_legend(FALSE))
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
    hub_data <- genetic_0.1_hits %>%
      filter(abs(rhoG) >= input$gen_thresh) %>% 
      group_by(trait1) %>%
      filter(n() >= input$gen_sharing) %>%
      ungroup() %>%
      mutate(abs_rhoG = abs(rhoG), edge_sign = ifelse(rhoG >= 0, "Positive", "Negative"))
    
    #Fail safe returns 0 
    if(nrow(hub_data) == 0) return(NULL)
    
    nodes_df <- data.frame(name = unique(c(hub_data$trait1, hub_data$trait2))) %>% 
      left_join(genetic_0.1_hits %>% select(trait1, Sub.Lineage) %>% distinct(trait1, .keep_all=TRUE), by = c("name" = "trait1")) %>%
      left_join(gen_module_lookup, by = "name") %>%
      mutate(type = ifelse(!is.na(Sub.Lineage), "IP", "Protein"), shape_map = ifelse(type == "Protein", "Protein", Sub.Lineage), module = ifelse(is.na(module), "Unassigned", module))
    
    g <- graph_from_data_frame(d = hub_data %>% 
                                 select(trait1, trait2, rhoG, abs_rhoG, edge_sign), 
                               directed = FALSE, vertices = nodes_df)
    V(g)$module <- as.character(cluster_leiden(g, objective_function = "modularity", 
                                               weights = E(g)$abs_rhoG)$membership)
    
    vis_nodes <- igraph::as_data_frame(g, what = "vertices") %>% rename(id = name) %>%
      mutate(label = ifelse(type == "Protein", id, NA), font.size = ifelse(type == "Protein", 40, 10), font.color = ifelse(type == "Protein", "black", "transparent"),
             shape = ifelse(type == "Protein", "diamond", "dot"), size = ifelse(type == "Protein", 40, 20), color = ifelse(type == "Protein", "black", global_lineage_colours[shape_map]), 
             original_color = color, original_font_color = font.color, title = paste0("<p><b>", id, "</b><br>Type: ", type, "<br>Module: ", module, "<br>Lineage: ", shape_map, "</p>"))
    
    
    vis_edges <- igraph::as_data_frame(g, what = "edges") %>% 
      rename(from = from, to = to) %>%
      mutate(id = row_number(), value = abs(rhoG), dashes = rhoG < 0, 
             color = "#7d7d7d", original_color = "#7d7d7d", 
             title = paste0("<p><b>rhoG:</b> ", round(rhoG, 3), "<br><b>Direction:</b> ", edge_sign, "</p>"))
    list(nodes = vis_nodes, edges = vis_edges)
  })
  
  
  
  #Initial network
  output$network_graph_gen <- renderVisNetwork({
    data <- network_data_gen(); req(data)
    visNetwork(data$nodes, data$edges) %>%
      #Force atlas based 
      visPhysics(solver = "forceAtlas2Based", 
                 forceAtlas2Based = list(gravitationalConstant = -200, springLength = 150), 
                 stabilization = list(enabled = TRUE, iterations = 1000)) %>%
      visEvents(stabilizationIterationsDone = "function () {this.setOptions({physics: false});}") %>%
      visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE), 
                 nodesIdSelection = list(enabled = TRUE, main = "Search by Name", values = data$nodes %>% 
                                           arrange(desc(type == "Protein"), id) %>% pull(id))) %>%
      
      visEdges(smooth = list(enabled = TRUE, type = "continuous")) %>% 
      visExport(type = "png", name = "IP_Protein_Gen") %>%
      visInteraction(navigationButtons = TRUE, tooltipDelay = 0)
  })
  
  
  
  #Dropdown changes based on observed
  observe({
    data <- network_data_gen(); req(data)
    mod <- input$gen_mod; lin <- input$gen_lin; prot <- input$gen_prot
    
    valid_mod <- data$nodes %>% 
      filter(type == "IP"); if(lin != "All") valid_mod <- valid_mod %>% filter(shape_map == lin); if(prot != "All") valid_mod <- valid_mod %>% filter(id %in% c((data$edges %>% filter(from == prot | to == prot))$from, (data$edges %>% filter(from == prot | to == prot))$to))
    valid_lin <- data$nodes %>% filter(type == "IP"); if(mod != "All") valid_lin <- valid_lin %>% filter(module == mod); if(prot != "All") valid_lin <- valid_lin %>% filter(id %in% c((data$edges %>% filter(from == prot | to == prot))$from, (data$edges %>% filter(from == prot | to == prot))$to))
    valid_prot <- data$nodes %>% filter(type == "IP"); if(mod != "All") valid_prot <- valid_prot %>% filter(module == mod); if(lin != "All") valid_prot <- valid_prot %>% filter(shape_map == lin)
    eprot <- data$edges %>% filter(from %in% valid_prot$id | to %in% valid_prot$id)
    
    updateSelectInput(session, "gen_mod", choices = c("All", as.character(sort(unique(valid_mod$module[!is.na(valid_mod$module)])))), selected = ifelse(mod %in% valid_mod$module, mod, "All"))
    updateSelectInput(session, "gen_lin", choices = c("All", as.character(sort(unique(valid_lin$shape_map)))), selected = ifelse(lin %in% valid_lin$shape_map, lin, "All"))
    updateSelectInput(session, "gen_prot", choices = c("All", as.character(sort(unique(data$nodes$id[data$nodes$type == "Protein" & data$nodes$id %in% c(eprot$from, eprot$to)])))), selected = ifelse(prot %in% c(eprot$from, eprot$to), prot, "All"))
  })
  
  #Grey out based on dropdown observed
  observe({
    data <- network_data_gen(); req(data, input$gen_mod, input$gen_lin, input$gen_prot)
    a_ips <- data$nodes %>% filter(type == "IP"); if(input$gen_mod != "All") a_ips <- a_ips %>% filter(module == input$gen_mod); if(input$gen_lin != "All") a_ips <- a_ips %>% filter(shape_map == input$gen_lin)
    a_prots <- data$nodes %>% filter(type == "Protein"); if(input$gen_prot != "All") a_prots <- a_prots %>% filter(id == input$gen_prot)
    a_edges <- data$edges %>% filter((from %in% a_ips$id & to %in% a_prots$id) | (to %in% a_ips$id & from %in% a_prots$id))
    
    t_ids <- if(input$gen_mod == "All" && input$gen_lin == "All" && input$gen_prot == "All") data$nodes$id else unique(c(a_edges$from, a_edges$to))
    t_e_ids <- if(input$gen_mod == "All" && input$gen_lin == "All" && input$gen_prot == "All") data$edges$id else a_edges$id
    
    visNetworkProxy("network_graph_gen") %>%
      visUpdateNodes(nodes = data$nodes %>% select(id, original_color, original_font_color) %>% mutate(color = ifelse(id %in% t_ids, original_color, "rgba(200,200,200,0.1)"), font.color = ifelse(id %in% t_ids, original_font_color, "rgba(200,200,200,0)"))) %>%
      visUpdateEdges(edges = data$edges %>% select(id, original_color) %>% mutate(color = ifelse(id %in% t_e_ids, original_color, "rgba(200,200,200,0.05)")))
  })
  
  #Module summary
  output$gen_module_summary <- renderUI({
    data <- network_data_gen(); req(data)
    mod <- input$gen_mod
    
    if(mod == "All") {
      return(HTML("<p style='color: #555; font-style: italic;'>Select a Module for information.</p>"))
    }
    
    #All Ips
    mod_ips <- data$nodes %>% filter(module == mod, type == "IP")
    if(nrow(mod_ips) == 0) return(HTML("<p>No IPs in this module.</p>"))
    
    #IP summary
    overall_lin_counts <- mod_ips %>% count(shape_map) %>% arrange(desc(n))
    overall_lin_text <- paste0(overall_lin_counts$shape_map, " (", overall_lin_counts$n, ")", collapse = ", ")
    
    #Find edges connection
    mod_edges <- data$edges %>% filter(from %in% mod_ips$id | to %in% mod_ips$id)
    
    #Ips to proteins
    edge_details <- mod_edges %>%
      left_join(data$nodes %>% select(id, type, shape_map), by = c("from" = "id")) %>% rename(from_type = type, from_lineage = shape_map) %>%
      left_join(data$nodes %>% select(id, type, shape_map), by = c("to" = "id")) %>% rename(to_type = type, to_lineage = shape_map)
    
    connections <- edge_details %>%
      mutate(
        Protein = ifelse(from_type == "Protein", from, to),
        IP = ifelse(from_type == "IP", from, to),
        Lineage = ifelse(from_type == "IP", from_lineage, to_lineage)
      ) %>%
      filter(IP %in% mod_ips$id) 
    
    if(nrow(connections) > 0) {
      prot_summary <- connections %>%
        group_by(Protein, Lineage) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(Protein) %>%
        summarise(details = paste0(Lineage, " (", n, ")", collapse = ", "), .groups = "drop") %>%
        mutate(html_string = paste0("<li><b>", Protein, ":</b> ", details, "</li>")) %>%
        pull(html_string) %>% paste(collapse = "")
      
      prot_text <- paste0("<ul style='margin-top: 5px; margin-bottom: 0; padding-left: 20px;'>", prot_summary, "</ul>")
      unique_prots_count <- length(unique(connections$Protein))
    } else {
      prot_text <- "None"
      unique_prots_count <- 0
    }
    
    HTML(paste0(
      "<div style='font-size: 15px;'>",
      "<span style='color: #1e88e5;'><b>Total IPs in Module ", mod, " (", nrow(mod_ips), "):</b></span><br>", overall_lin_text, "<br><br>",
      "<span style='color: #1e88e5;'><b>Total Proteins (", unique_prots_count, "):</b></span>", prot_text,
      "</div>"
    ))
  })
  
  
  
  
  #ENVIRONMENT NETWORK
  network_data_env <- reactive({
    req(input$env_thresh, input$env_sharing)
    hub_data <- environmental_hits %>%
      filter(abs(rhoE) >= input$env_thresh) %>% 
      group_by(trait1) %>%
      filter(n() >= input$env_sharing) %>% 
      ungroup() %>%
      mutate(abs_rhoE = abs(rhoE), edge_sign = ifelse(rhoE >= 0, "Positive", "Negative"))
    
    if(nrow(hub_data) == 0) return(NULL) 
    
    nodes_df <- data.frame(name = unique(c(hub_data$trait1, hub_data$trait2))) %>% 
      left_join(environmental_hits %>% select(trait1, Sub.Lineage) %>% distinct(trait1, .keep_all=TRUE), by = c("name" = "trait1")) %>%
      left_join(env_module_lookup, by = "name") %>%
      mutate(type = ifelse(!is.na(Sub.Lineage), "IP", "Protein"), shape_map = ifelse(type == "Protein", "Protein", Sub.Lineage), module = ifelse(is.na(module), "Unassigned", module))
    
    g <- graph_from_data_frame(d = hub_data %>% 
                                 select(trait1, trait2, rhoE, abs_rhoE, edge_sign), 
                               directed = FALSE, vertices = nodes_df)
    V(g)$module <- as.character(cluster_leiden(g, objective_function = "modularity", 
                                               weights = E(g)$abs_rhoE)$membership)
    
    vis_nodes <- igraph::as_data_frame(g, what = "vertices") %>% rename(id = name) %>%
      mutate(label = ifelse(type == "Protein", id, NA), font.size = ifelse(type == "Protein", 50, 10), font.color = ifelse(type == "Protein", "black", "transparent"), 
             shape = ifelse(type == "Protein", "diamond", "dot"), size = ifelse(type == "Protein", 40, 20), color = ifelse(type == "Protein", "black", global_lineage_colours[shape_map]), 
             original_color = color, original_font_color = font.color, title = paste0("<p><b>", id, "</b><br>Type: ", type, "<br>Module: ", module, "<br>Lineage: ", shape_map, "</p>"))
    
    
    vis_edges <- igraph::as_data_frame(g, what = "edges") %>%
      rename(from = from, to = to) %>%
      mutate(id = row_number(), value = abs(rhoE), 
             dashes = rhoE < 0, color = "#7d7d7d", original_color = "#7d7d7d", 
             title = paste0("<p><b>rhoE:</b> ", round(rhoE, 3), "<br><b>Direction:</b> ", edge_sign, "</p>"))
    list(nodes = vis_nodes, edges = vis_edges)
  })
  
  #Layout 
  output$network_graph_env <- renderVisNetwork({
    data <- network_data_env(); req(data)
    visNetwork(data$nodes, data$edges) %>%
      visPhysics(solver = "forceAtlas2Based", 
                 forceAtlas2Based = list(gravitationalConstant = -200, springLength = 150), 
                 stabilization = list(enabled = TRUE, iterations = 1000)) %>%
      visEvents(stabilizationIterationsDone = "function () {this.setOptions({physics: false});}") %>%
      visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE), 
                 nodesIdSelection = list(enabled = TRUE, main = "Search by Name", 
                                         values = data$nodes %>% arrange(desc(type == "Protein"), id) %>% 
                                           pull(id))) %>%
      
      visEdges(smooth = list(enabled = TRUE, type = "continuous")) %>% 
      visExport(type = "png", name = "IP_Protein_Env") %>%
      visInteraction(navigationButtons = TRUE, tooltipDelay = 0)
  })
  
  #Reactions based on inputs
  observe({
    data <- network_data_env(); req(data)
    mod <- input$env_mod; lin <- input$env_lin; prot <- input$env_prot
    
    valid_mod <- data$nodes %>% filter(type == "IP"); if(lin != "All") valid_mod <- valid_mod %>% filter(shape_map == lin); if(prot != "All") valid_mod <- valid_mod %>% filter(id %in% c((data$edges %>% filter(from == prot | to == prot))$from, (data$edges %>% filter(from == prot | to == prot))$to))
    valid_lin <- data$nodes %>% filter(type == "IP"); if(mod != "All") valid_lin <- valid_lin %>% filter(module == mod); if(prot != "All") valid_lin <- valid_lin %>% filter(id %in% c((data$edges %>% filter(from == prot | to == prot))$from, (data$edges %>% filter(from == prot | to == prot))$to))
    valid_prot <- data$nodes %>% filter(type == "IP"); if(mod != "All") valid_prot <- valid_prot %>% filter(module == mod); if(lin != "All") valid_prot <- valid_prot %>% filter(shape_map == lin)
    eprot <- data$edges %>% filter(from %in% valid_prot$id | to %in% valid_prot$id)
    
    updateSelectInput(session, "env_mod", choices = c("All", as.character(sort(unique(valid_mod$module[!is.na(valid_mod$module)])))), selected = ifelse(mod %in% valid_mod$module, mod, "All"))
    updateSelectInput(session, "env_lin", choices = c("All", as.character(sort(unique(valid_lin$shape_map)))), selected = ifelse(lin %in% valid_lin$shape_map, lin, "All"))
    updateSelectInput(session, "env_prot", choices = c("All", as.character(sort(unique(data$nodes$id[data$nodes$type == "Protein" & data$nodes$id %in% c(eprot$from, eprot$to)])))), selected = ifelse(prot %in% c(eprot$from, eprot$to), prot, "All"))
  })
  
  observe({
    data <- network_data_env(); req(data, input$env_mod, input$env_lin, input$env_prot)
    a_ips <- data$nodes %>% filter(type == "IP"); if(input$env_mod != "All") a_ips <- a_ips %>% filter(module == input$env_mod); if(input$env_lin != "All") a_ips <- a_ips %>% filter(shape_map == input$env_lin)
    a_prots <- data$nodes %>% filter(type == "Protein"); if(input$env_prot != "All") a_prots <- a_prots %>% filter(id == input$env_prot)
    a_edges <- data$edges %>% filter((from %in% a_ips$id & to %in% a_prots$id) | (to %in% a_ips$id & from %in% a_prots$id))
    
    t_ids <- if(input$env_mod == "All" && input$env_lin == "All" && input$env_prot == "All") data$nodes$id else unique(c(a_edges$from, a_edges$to))
    t_e_ids <- if(input$env_mod == "All" && input$env_lin == "All" && input$env_prot == "All") data$edges$id else a_edges$id
    
    visNetworkProxy("network_graph_env") %>%
      visUpdateNodes(nodes = data$nodes %>% select(id, original_color, original_font_color) %>% mutate(color = ifelse(id %in% t_ids, original_color, "rgba(200,200,200,0.1)"), font.color = ifelse(id %in% t_ids, original_font_color, "rgba(200,200,200,0)"))) %>%
      visUpdateEdges(edges = data$edges %>% select(id, original_color) %>% mutate(color = ifelse(id %in% t_e_ids, original_color, "rgba(200,200,200,0.05)")))
  })
  
  #Module summary
  output$env_module_summary <- renderUI({
    data <- network_data_env(); req(data)
    mod <- input$env_mod
    
    if(mod == "All") {
      return(HTML("<p style='color: #555; font-style: italic;'>Select a Module for information.</p>"))
    }
    
    mod_ips <- data$nodes %>% filter(module == mod, type == "IP")
    if(nrow(mod_ips) == 0) return(HTML("<p>No IPs in this module.</p>"))
    
    overall_lin_counts <- mod_ips %>% count(shape_map) %>% arrange(desc(n))
    overall_lin_text <- paste0(overall_lin_counts$shape_map, " (", overall_lin_counts$n, ")", collapse = ", ")
    
    mod_edges <- data$edges %>% filter(from %in% mod_ips$id | to %in% mod_ips$id)
    
    edge_details <- mod_edges %>%
      left_join(data$nodes %>% select(id, type, shape_map), by = c("from" = "id")) %>% rename(from_type = type, from_lineage = shape_map) %>%
      left_join(data$nodes %>% select(id, type, shape_map), by = c("to" = "id")) %>% rename(to_type = type, to_lineage = shape_map)
    
    connections <- edge_details %>%
      mutate(
        Protein = ifelse(from_type == "Protein", from, to),
        IP = ifelse(from_type == "IP", from, to),
        Lineage = ifelse(from_type == "IP", from_lineage, to_lineage)
      ) %>%
      filter(IP %in% mod_ips$id)
    
    if(nrow(connections) > 0) {
      prot_summary <- connections %>%
        group_by(Protein, Lineage) %>% summarise(n = n(), .groups = "drop") %>%
        group_by(Protein) %>% summarise(details = paste0(Lineage, " (", n, ")", collapse = ", "), .groups = "drop") %>%
        mutate(html_string = paste0("<li><b>", Protein, ":</b> ", details, "</li>")) %>%
        pull(html_string) %>% paste(collapse = "")
      
      prot_text <- paste0("<ul style='margin-top: 5px; margin-bottom: 0; padding-left: 20px;'>", prot_summary, "</ul>")
      unique_prots_count <- length(unique(connections$Protein))
    } else {
      prot_text <- "None"
      unique_prots_count <- 0
    }
    
    HTML(paste0(
      "<div style='font-size: 15px;'>",
      "<span style='color: #1e88e5;'><b>Total IPs in Module ", mod, " (", nrow(mod_ips), "):</b></span><br>", overall_lin_text, "<br><br>",
      "<span style='color: #1e88e5;'><b>Total Proteins (", unique_prots_count, ")</b></span>", prot_text,
      "</div>"
    ))
  })
  
  #Shared stats observed
  observe({
    data <- network_data_gen(); req(data)
    mod <- input$gen_mod; lin <- input$gen_lin; prot <- input$gen_prot
    valid_mod <- data$nodes %>% filter(type == "IP"); if(lin != "All") valid_mod <- valid_mod %>% filter(shape_map == lin); if(prot != "All") valid_mod <- valid_mod %>% filter(id %in% c((data$edges %>% filter(from == prot | to == prot))$from, (data$edges %>% filter(from == prot | to == prot))$to))
    valid_lin <- data$nodes %>% filter(type == "IP"); if(mod != "All") valid_lin <- valid_lin %>% filter(module == mod); if(prot != "All") valid_lin <- valid_lin %>% filter(id %in% c((data$edges %>% filter(from == prot | to == prot))$from, (data$edges %>% filter(from == prot | to == prot))$to))
    valid_prot <- data$nodes %>% filter(type == "IP"); if(mod != "All") valid_prot <- valid_prot %>% filter(module == mod); if(lin != "All") valid_prot <- valid_prot %>% filter(shape_map == lin)
    eprot <- data$edges %>% filter(from %in% valid_prot$id | to %in% valid_prot$id)
    updateSelectInput(session, "gen_mod", choices = c("All", as.character(sort(unique(valid_mod$module[!is.na(valid_mod$module)])))), selected = ifelse(mod %in% valid_mod$module, mod, "All"))
    updateSelectInput(session, "gen_lin", choices = c("All", as.character(sort(unique(valid_lin$shape_map)))), selected = ifelse(lin %in% valid_lin$shape_map, lin, "All"))
    updateSelectInput(session, "gen_prot", choices = c("All", as.character(sort(unique(data$nodes$id[data$nodes$type == "Protein" & data$nodes$id %in% c(eprot$from, eprot$to)])))), selected = ifelse(prot %in% c(eprot$from, eprot$to), prot, "All"))
  })
  
  observe({
    data <- network_data_env(); req(data)
    mod <- input$env_mod; lin <- input$env_lin; prot <- input$env_prot
    valid_mod <- data$nodes %>% filter(type == "IP"); if(lin != "All") valid_mod <- valid_mod %>% filter(shape_map == lin); if(prot != "All") valid_mod <- valid_mod %>% filter(id %in% c((data$edges %>% filter(from == prot | to == prot))$from, (data$edges %>% filter(from == prot | to == prot))$to))
    valid_lin <- data$nodes %>% filter(type == "IP"); if(mod != "All") valid_lin <- valid_lin %>% filter(module == mod); if(prot != "All") valid_lin <- valid_lin %>% filter(id %in% c((data$edges %>% filter(from == prot | to == prot))$from, (data$edges %>% filter(from == prot | to == prot))$to))
    valid_prot <- data$nodes %>% filter(type == "IP"); if(mod != "All") valid_prot <- valid_prot %>% filter(module == mod); if(lin != "All") valid_prot <- valid_prot %>% filter(shape_map == lin)
    eprot <- data$edges %>% filter(from %in% valid_prot$id | to %in% valid_prot$id)
    updateSelectInput(session, "env_mod", choices = c("All", as.character(sort(unique(valid_mod$module[!is.na(valid_mod$module)])))), selected = ifelse(mod %in% valid_mod$module, mod, "All"))
    updateSelectInput(session, "env_lin", choices = c("All", as.character(sort(unique(valid_lin$shape_map)))), selected = ifelse(lin %in% valid_lin$shape_map, lin, "All"))
    updateSelectInput(session, "env_prot", choices = c("All", as.character(sort(unique(data$nodes$id[data$nodes$type == "Protein" & data$nodes$id %in% c(eprot$from, eprot$to)])))), selected = ifelse(prot %in% c(eprot$from, eprot$to), prot, "All"))
  })
  
  observe({
    data <- network_data_gen(); req(data, input$gen_mod, input$gen_lin, input$gen_prot)
    a_ips <- data$nodes %>% filter(type == "IP"); if(input$gen_mod != "All") a_ips <- a_ips %>% filter(module == input$gen_mod); if(input$gen_lin != "All") a_ips <- a_ips %>% filter(shape_map == input$gen_lin)
    a_prots <- data$nodes %>% filter(type == "Protein"); if(input$gen_prot != "All") a_prots <- a_prots %>% filter(id == input$gen_prot)
    a_edges <- data$edges %>% filter((from %in% a_ips$id & to %in% a_prots$id) | (to %in% a_ips$id & from %in% a_prots$id))
    t_ids <- if(input$gen_mod == "All" && input$gen_lin == "All" && input$gen_prot == "All") data$nodes$id else unique(c(a_edges$from, a_edges$to))
    t_e_ids <- if(input$gen_mod == "All" && input$gen_lin == "All" && input$gen_prot == "All") data$edges$id else a_edges$id
    visNetworkProxy("network_graph_gen") %>% visUpdateNodes(nodes = data$nodes %>% select(id, original_color, original_font_color) %>% mutate(color = ifelse(id %in% t_ids, original_color, "rgba(200,200,200,0.1)"), font.color = ifelse(id %in% t_ids, original_font_color, "rgba(200,200,200,0)"))) %>% visUpdateEdges(edges = data$edges %>% select(id, original_color) %>% mutate(color = ifelse(id %in% t_e_ids, original_color, "rgba(200,200,200,0.05)")))
  })
  
  observe({
    data <- network_data_env(); req(data, input$env_mod, input$env_lin, input$env_prot)
    a_ips <- data$nodes %>% filter(type == "IP"); if(input$env_mod != "All") a_ips <- a_ips %>% filter(module == input$env_mod); if(input$env_lin != "All") a_ips <- a_ips %>% filter(shape_map == input$env_lin)
    a_prots <- data$nodes %>% filter(type == "Protein"); if(input$env_prot != "All") a_prots <- a_prots %>% filter(id == input$env_prot)
    a_edges <- data$edges %>% filter((from %in% a_ips$id & to %in% a_prots$id) | (to %in% a_ips$id & from %in% a_prots$id))
    t_ids <- if(input$env_mod == "All" && input$env_lin == "All" && input$env_prot == "All") data$nodes$id else unique(c(a_edges$from, a_edges$to))
    t_e_ids <- if(input$env_mod == "All" && input$env_lin == "All" && input$env_prot == "All") data$edges$id else a_edges$id
    visNetworkProxy("network_graph_env") %>% visUpdateNodes(nodes = data$nodes %>% select(id, original_color, original_font_color) %>% mutate(color = ifelse(id %in% t_ids, original_color, "rgba(200,200,200,0.1)"), font.color = ifelse(id %in% t_ids, original_font_color, "rgba(200,200,200,0)"))) %>% visUpdateEdges(edges = data$edges %>% select(id, original_color) %>% mutate(color = ifelse(id %in% t_e_ids, original_color, "rgba(200,200,200,0.05)")))
  })
}







#Run App
shinyApp(ui = ui, server = server)