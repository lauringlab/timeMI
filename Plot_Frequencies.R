#Plotting functions for amino acid frequencies

get_index <- function(pos){
  return(mi_network_info[mi_network_info$index == pos,]$id %>%
           str_replace("_", " "))
}
get_index_v <- Vectorize(get_index)

get_position_freq <- function(position, seqs){
  #takes a position and gets the frequency of residues
  temp <- mutate(seqs, `Amino Acid` = str_sub(x, start = position, end = position)) %>%
    count(`Amino Acid`, date, .drop = FALSE) %>%
    group_by(date) %>%
    mutate(freq = n/sum(n)) %>%
    ungroup() %>%
    complete(date, `Amino Acid`, fill = list(n = 0, freq = 0))
  temp
}

get_subgraph_freqs <- function(elements){
  names(elements) <- elements
  m <- lapply(elements, get_position_freq, seq_times) %>%
    bind_rows(.id = "position") %>%
    mutate(position = get_index_v(position))
  return(m)
}

pal2 <- c("A" = "#8CFF8C", "G" = "#FFFFFF", "L" = "#455E45", "S" = "#FF7042",
          "V" = "#D695C2", "T" = "#B84C00", "K" = "#4747B8", "D" = "#A00042",
          "I" = "#004C00", "N" = "#FF7C70", "E" = "#660000", "P" = "#525252",
          "R" = "#00007C", "F" = "#543C42", "Q" = "#FF4C4C", "Y" = "#8C704C",
          "H" = "#7070FF", "C" = "#FFFF70", "M" = "#B8A042", "W" = "#4F4600",
          "X" = "#B8B8B8")

plot_example_freqs <- function(nodes){
  ggplot(get_subgraph_freqs(nodes), aes(x = date, y = freq, fill = `Amino Acid`,
                                        color = `Amino Acid`, group = `Amino Acid`)) +
    geom_area(alpha = 0.4) +
    facet_wrap(~position, ncol = 1) +
    scale_fill_manual(values = pal2) +
    scale_color_manual(values = pal2) +
    scale_y_continuous(breaks = c(0, 0.5, 1.0)) +
    theme_minimal() +
    labs(title = "Amino Acid Frequencies",
         x = "Date", y = "Frequency") +
    theme(text = element_text(size = 12))
}
