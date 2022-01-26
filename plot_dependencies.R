# Plot the dependency graph of analysis scripts (R and Rmd files).
# The graph is saved in 'plots/dependency_graph.pdf'.

library(igraph)
library(tidyverse)
library(fs)



# Prepare edge and vertex lists -------------------------------------------

get_dependencies <- function(file) {
  contents <- read_lines(file)
  if (path_ext(file) == "R") {
    regex_in <- "@DEPI\\s+(.+)$"
    regex_out <- "@DEPO\\s+(.+)$"
  } else if (path_ext(file) == "Rmd") {
    regex_in <- "@DEPI\\s+(.+)\\s+-->$"
    regex_out <- "@DEPO\\s+(.+)\\s+-->$"
  } else {
    stop("Unsupported file type.")
  }
  
  
  in_files <- 
    str_match(contents, regex_in) %>% 
    magrittr::extract(, 2) %>% 
    discard(is.na)
  
  out_files <- 
    str_match(contents, regex_out) %>% 
    magrittr::extract(, 2) %>% 
    discard(is.na)
  
  bind_rows(
    tibble(
      source = in_files,
      sink = as.character(file)
    ),
    tibble(
      source = as.character(file),
      sink = out_files
    )
  )
}

script_files <- dir_ls(regex = "\\.R(md)?$")

edge_list <- map_dfr(script_files, get_dependencies)

generated_files <- 
  c(edge_list$source, edge_list$sink) %>% 
  unique() %>% 
  setdiff(script_files)

vertex_list <-
  bind_rows(
    tibble(
      name = as.character(script_files),
      color = "#80b3ff",
      shape = "box"
    ),
    tibble(
      name = generated_files,
      color = "#ffe680",
      shape = "ellipse"
    )
  ) %>% 
  mutate(
    label = name,
    style = "filled",
    fontname = "Helvetica"
  ) %>% 
  filter(!name %>% str_starts("wip_"))
  


# Make graph --------------------------------------------------------------

G <- graph_from_data_frame(edge_list, vertices = vertex_list)
vertex_attr(G, "color", index = V(G)[degree(G) == 0]) <- "#cccccc"

dotfile <- tempfile()
write_graph(G, dotfile, format = "dot")
system(str_glue("dot -Tpng {dotfile} > plots/dependency_graph.png"))
file.remove(dotfile)
