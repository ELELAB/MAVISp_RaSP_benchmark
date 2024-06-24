#!/usr/bin/env Rscript

# rasp_ros_bench.r - script for benchmarking RAsP data
# Copyright (C) 2024 Katrine Meldgard, Danish Cancer Society
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# User must supply the path to the mavisp data as command line argument
# (not with simple mode or ensemble mode specified)


#------------------------------------------
#--------------Initilisation---------------
#------------------------------------------
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library('tidymodels'))
suppressPackageStartupMessages(library('corrr'))
library(forcats)
library(ggnewscale)

args <- commandArgs()
summary_content <- c()

if (length(args) == 0) {
    stop('Please provide the path to the mavisp data')
}

mavisp_index_path <- str_c(args,'simple_mode')
index_file <- list.files(mavisp_index_path,
                         pattern = 'index.csv',
                         full.names = TRUE)[1]
mavisp_path <- str_c(args,'simple_mode/dataset_tables/')

mavisp_files <- list.files(mavisp_path,
                            pattern = '.csv',
                            full.names = TRUE)
                            

#----------------------------------
#-------------Load data------------
#----------------------------------
# Load mavisp data into one tibble
mavisp_data <- mavisp_files |>
      map(\(mavisp_files) read_csv(file = mavisp_files,
                                    col_select = c(1,
                                                    matches('(Stability \\(Rosetta)'),
                                                    matches('(Stability \\(RaSP)'),
                                                    matches('(Stability classification, [A-Za-z0-9]+, \\(Rosetta, FoldX\\))'),
                                                    matches('(Stability classification, [A-Za-z0-9]+, \\(RaSP, FoldX\\))'),
                                                    matches('(Solvent Accessibility)')),
                                    col_types = 'c',
                                    id = 'protein') |>
                            rename('mutation' = 2,
                                    'rosetta_value' = matches('(Stability \\(Rosetta)'),
                                    'rasp_value' = matches('(Stability \\(RaSP)'),
                                    'rosetta' = matches('(Stability classification, [A-Za-z0-9]+, \\(Rosetta, FoldX\\))'),
                                    'rasp' = matches('(Stability classification, [A-Za-z0-9]+, \\(RaSP, FoldX\\))'),
                                    'sasa' = matches('(Solvent Accessibility)')) |>
                            mutate('protein' = str_remove(basename(protein),'-simple_mode.csv'),
                                    across(matches('revel'), as.double))) |>
                                
      list_rbind() 

index <- read_csv(file = index_file,
                  col_select = c('Protein','Uniprot AC'))

print("INDEX LOADED")

#----------------------------------
# ---------Augment data------------
#----------------------------------

augmented_mavisp_data <- mavisp_data |>
    mutate(aa_from = str_split_i(mutation, '\\d+', 1),
            aa_to = str_split_i(mutation, '\\d+', 2),
            aa_change = str_c(aa_from,aa_to),
            rosetta = factor(rosetta,
                             levels = c('Destabilizing', 'Neutral', 'Stabilizing', 'Uncertain')),
            rasp = factor(rasp,
                          levels = c('Destabilizing', 'Neutral', 'Stabilizing', 'Uncertain')),
            consensus = case_when(rosetta == rasp ~ 'TRUE',
                                    rosetta != rasp ~ 'FALSE')) |>
    drop_na(consensus)
    
#-------------------------------------------
#------- Overall evaluation metrics---------
#-------------------------------------------

# Number of proteins
n_proteins <- augmented_mavisp_data |>
    count(protein) |>
    nrow()

# Number of observations
n_observations <- augmented_mavisp_data |>
    nrow()

# Pearson correlation
pearson_cor <- augmented_mavisp_data |>
    select(rosetta_value, rasp_value) |>
    correlate() |> 
    # Output is correlation matrix (2x2). The last value 
    # of rosetta column corresponds to rasp
    pull(rosetta_value) |>
    last()

# MAE
mae <- augmented_mavisp_data |>
    mutate(diff = rosetta_value - rasp_value) |>
    pull(diff) |>
    mean()

# Max Rosetta value
max_ros_value <- augmented_mavisp_data |>
    select(rosetta_value) |>
    max()

# Max RasP value
max_rasp_value <- augmented_mavisp_data |>
    select(rasp_value) |>
    max()

# Number of outliers
nr_outliers <- augmented_mavisp_data |>
    filter(rosetta_value >= 15) |>
    nrow()

# Pearson correlation without outliers
pearson_cor_no_outlier <- augmented_mavisp_data |>
    select(rosetta_value, rasp_value) |>
    filter(rosetta_value < 15) |>
    correlate() |> 
    # Output is correlation matrix (2x2). The last value 
    # of rosetta column corresponds to rasp
    pull(rosetta_value) |>
    last()

# MAE without outliers
mae_no_outliers <- augmented_mavisp_data |>
    filter(rosetta_value < 15) |>
    mutate(diff = rosetta_value - rasp_value) |>
    pull(diff) |>
    mean()

# Accuracy
n_correct <- augmented_mavisp_data |>
    filter(consensus == 'TRUE') |>
    count(consensus) |>
    pull(n)

# Precision
precision <- augmented_mavisp_data |> 
    precision(truth = rosetta,
              estimate = rasp,
              estimator = 'macro')

precision_3 <- precision$.estimate

precision_4 <- precision_3 * 3/4

#Recall
recall <- augmented_mavisp_data |> 
    recall(truth = rosetta,
              estimate = rasp,
              estimator = 'macro')

recall_4 <- recall$.estimate

recall_3 <- recall_4 * 4/3
 
all_accuracy <- n_correct/n_observations
summary_content <- append(summary_content, c('Number of proteins',n_proteins))
summary_content <- append(summary_content, c('Number of observations',n_observations))
summary_content <- append(summary_content, c('Pearson correlation',pearson_cor))
summary_content <- append(summary_content, c('MAE',mae))
summary_content <- append(summary_content, c('Maximum Rosetta value',max_ros_value))
summary_content <- append(summary_content, c('Maximum RaSP',max_rasp_value))
summary_content <- append(summary_content, c('Number of outliers (rosetta values >= 15)',nr_outliers))
summary_content <- append(summary_content, c('Pearson correlation without rosetta values >= 15',pearson_cor_no_outlier))
summary_content <- append(summary_content, c('MAE without rosetta values >= 15',mae_no_outliers))
summary_content <- append(summary_content, c('Accuracy of RaSP with Rosetta as truth',all_accuracy))
summary_content <- append(summary_content, c('Precision of RaSP with Rosetta as truth (3 classes)',precision_3))
summary_content <- append(summary_content, c('Recall of RaSP with Rosetta as truth (3 classes)',recall_3))
summary_content <- append(summary_content, c('Precision of RaSP with Rosetta as truth (4 classes)',precision_4))
summary_content <- append(summary_content, c('Recall of RaSP with Rosetta as truth (4 classes)',recall_4))

# Make confusion matrix
confusion_matrix <- augmented_mavisp_data |>
    conf_mat(truth = rosetta,
             estimate = rasp)

tidy_confusion_matrix <- confusion_matrix |>
    tidy() |>
    mutate(rasp = str_split_i(name, pattern = '_', i = 2),
            rosetta = str_split_i(name, pattern = '_', i = 3)) |>
    mutate(rosetta = case_when(rosetta == 1 ~ 'Destabilizing',
                                rosetta == 2 ~ 'Neutral',
                                rosetta == 3 ~ 'Stabilizing',
                                rosetta == 4 ~ 'Uncertain'),
            rasp = case_when(rasp == 1 ~ 'Destabilizing',
                            rasp == 2 ~ 'Neutral',
                            rasp == 3 ~ 'Stabilizing',
                            rasp == 4 ~ 'Uncertain'))


confusion_plot <- tidy_confusion_matrix |>
    ggplot(mapping = aes(x = rosetta,
                         y = rasp,
                         fill = value)) +
    geom_tile() +
    geom_text(aes(label=value)) +
    theme_light() +
    scale_fill_gradient2(low = 'white', mid = 'tomato3', high = 'cornflowerblue', midpoint = 22361) +
    labs(x = 'Rosetta',
         y = 'RaSP',
         fill = 'Value') +
    theme(axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10)) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0))   

ggsave(plot = confusion_plot,
        file = 'polished_plots/confusion_matrix.png',
        width = 6,
        height = 4.5)

#-----------------------------------------------------------
#--------- Bar plot of the counts in each protein ----------
#-----------------------------------------------------------

protein_counts <- augmented_mavisp_data |>
    group_by(protein) |>
    count(protein) |>
    arrange(desc(n)) |>
    ungroup() |>
    mutate(accumulated = cumsum(n),
            prot_num = row_number())

protein_counts |> 
    left_join(index, 
              by = join_by(protein == Protein)) |>
    select(-prot_num) |>
    rename(nr_mutations = n,
           accumulated_count = accumulated) |>
    write_csv(file = 'protein_counts.csv')

max_below_half <- protein_counts |>
    filter(accumulated < n_observations/2) |>
    last() |>
    pull(prot_num)

summary_content <- append(summary_content, '\n')
summary_content <- append(summary_content, 'Number of proteins constituting half of the observations.')
summary_content <- append(summary_content, max_below_half)

accumulated_counts <- protein_counts |>
    mutate(accumulated = cumsum(n),
            prot_num = row_number()) |>
    ggplot(mapping = aes(x = reorder(prot_num, -n),
                         y = accumulated)) +
    geom_point() +
    theme_light() +
    geom_hline(yintercept = c(n_observations, n_observations*0.5)) +
    geom_text(aes(0,n_observations,label = '100%', vjust = -1, hjust = -0.25)) +
    geom_text(aes(0,n_observations*0.5,label = '50%', vjust = -1, hjust = -0.5)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          plot.title = element_text(hjust = 0.5, size = 20)) +
    labs(x = 'Proteins in dataset sorted from highest to lowest number of observations',
         y = 'Accumulated number of observations',
         title = 'Accumulated number of observations') +
    scale_x_discrete(breaks = seq(10,n_observations, by = 10))

ggsave(plot = accumulated_counts,
        file = 'polished_plots/accumulated_counts.png',
        width = 7,
        height = 8)



# Removing the proteins with most observations
filtered_data <- augmented_mavisp_data |>
    anti_join(protein_counts |>
              filter(!(prot_num %in% c(1:max_below_half))))


# Pearson correlation
pearson_cor_filt <- filtered_data |>
    select(rosetta_value, rasp_value) |>
    correlate() |> 
    # Output is correlation matrix (2x2). The last value 
    # of rosetta column corresponds to rasp
    pull(rosetta_value) |>
    last()

# MAE
mae_filt <- filtered_data |>
    mutate(diff = rosetta_value - rasp_value) |>
    pull(diff) |>
    mean()  

summary_content <- append(summary_content, 'Pearson correlation without those 23 proteins')
summary_content <- append(summary_content, pearson_cor_filt)
summary_content <- append(summary_content, 'MAE without those 23 proteins')
summary_content <- append(summary_content, mae_filt)

#------------------------------------------------------------------
#--------- Bar plot of the counts in each classification ----------
#------------------------------------------------------------------

class_diff_bar <- augmented_mavisp_data |>
    select(rosetta, rasp) |>
    rename('Rosetta' = rosetta,
            'RaSP' = rasp) |>
    pivot_longer(cols = c(Rosetta, RaSP),
                 names_to = 'method',
                 values_to = 'class') |>
    count(method, class) |>
    complete(method = unique(method),
             class = unique(class),
             fill = list(n = 0)) |>
    ggplot(mapping = aes(x = class,
                         y = n,
                         fill = method)) +
    geom_col(position='dodge') +
    theme_light() +
    geom_text(aes(label = n),  
                vjust = -0.5, 
                size = 5,
                colour = "black", 
                position=position_dodge(width=0.9)) +
    labs(x = 'Class',
            y = 'Count',
            fill = 'Method') +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=18),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          plot.margin = margin(t = 1, r = 1, b = 0.5, l = 1, unit = "cm")) +
    scale_fill_manual(labels = c('RaSP', 'Rosetta'),
                      values = c('skyblue','orange')) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 51000))

ggsave(plot = class_diff_bar,
        file = 'polished_plots/class_diff_bar.png',
        width = 10)



# ------------------------------------------------------------------
# -------------- Scatter plots for rosetta vs rasp values ----------
# ------------------------------------------------------------------

# Scatter plot for all
all_scatter <- augmented_mavisp_data |>
    ggplot(mapping = aes(x = rosetta_value,
                         y = rasp_value)) +
    geom_point(mapping = aes(color = consensus),
                alpha = 0.5) +
    theme_light() +
    coord_fixed() +
    labs(x = "Rosetta \u0394\u0394G (kcal/mol)",
         y = "RaSP \u0394\u0394G (kcal/mol)", 
         color = 'Consensus') +
    scale_color_discrete(labels = c('No','Yes')) +
    new_scale_color() +
    geom_smooth(mapping = aes(color = 'Linear regression'),
            method = "lm",
            se = FALSE,
            linewidth = 0.5) +
    scale_color_manual(' ', values = 'black') +
    geom_vline(xintercept = c(-3, -2, 2, 3),
                linetype = 'dotted') +
    geom_hline(yintercept = c(-3, -2, 2, 3),
                linetype = 'dotted') +
    geom_hline(mapping = aes(yintercept = -3, 
                             linetype = 'MAVISp decision boundaries')) +
    scale_linetype_manual(' ', values = 'dotted') +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          plot.margin = margin(t = 1, r = 1, b = 0.5, l = 1, unit = "cm")) +
    annotate("text", x=20, y=22, size = 6, label= str_c("Pearson correlation: ",round(pearson_cor,digits = 2)))

ggsave(plot = all_scatter,
        file = 'polished_plots/rasp_rosetta_scatter.png', 
        width = 14)


#-------------------------------------------------
#-------------- Density plot for sasa ------------
#-------------------------------------------------

sasa_density <- augmented_mavisp_data |>
    group_by(consensus) |>
    ggplot(mapping = aes(x = sasa,
                         color = consensus)) +
    geom_density(size = 2) +
    theme_light() +
    labs(x = "Relative side-chain solvent accessible area (%)",
         y = "Density", 
         color = 'Consensus') +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          plot.margin = margin(t = 1, r = 1, b = 0.5, l = 1, unit = "cm"))

ggsave(plot = sasa_density,
        file = 'polished_plots/sasa_density.png')



#------------------------------------------------------------------
#------------------- Bar plot of amino acid accuracy --------------
#------------------------------------------------------------------

aa_accuracy_plot <- augmented_mavisp_data |>
    pivot_longer(cols = c(aa_from, aa_to),
                 names_to = 'sub_class',
                 values_to = 'aa') |>
    group_by(aa, sub_class) |>
    count(sub_class, aa, consensus) |>
    pivot_wider(id_cols = c(sub_class, aa), 
                names_from = consensus,
                values_from = n) |>
    rename('false' = `FALSE`,
            'true' = `TRUE`) |>
    mutate(total = false+true,
            accuracy = true/total) |>
    ggplot(mapping = aes(x = aa, 
                         y = accuracy,
                         fill = sub_class)) +
    geom_col(position = 'dodge') +
    theme_light() +
    labs(x = 'Amino Acid Type',
         y = 'Accuracy',
         fill = 'Substitution type') +
    scale_fill_manual(labels = c('From', 'To'),
                      values = c('skyblue','orange')) +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          plot.margin = margin(t = 1, r = 1, b = 0.5, l = 1, unit = "cm")) +
    scale_y_continuous(expand = c(0,0), 
                       limits = c(0, 1))

ggsave(plot = aa_accuracy_plot,
       file = 'polished_plots/aa_accuracy.png',
       width = 14)


aa_distribution_plot <- augmented_mavisp_data |>
    pivot_longer(cols = c(aa_from, aa_to),
                 names_to = 'sub_class',
                 values_to = 'aa') |>
    group_by(aa, sub_class) |>
    count(sub_class, aa, consensus) |>
    pivot_wider(id_cols = c(sub_class, aa), 
                names_from = consensus,
                values_from = n) |>
    rename('false' = `FALSE`,
            'true' = `TRUE`) |>
    mutate(total = false+true,
            accuracy = true/total) |>
    ggplot(mapping = aes(x = aa, 
                         y = total,
                         fill = sub_class)) +
    geom_col(position = 'dodge') +
    theme_light() +
    labs(x = 'Amino Acid Type',
         y = 'Count',
         fill = 'Substitution type') +
    scale_fill_manual(labels = c('From', 'To'),
                      values = c('skyblue','orange')) +
    theme(axis.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm")) +
    scale_y_continuous(expand = c(0,0), 
                        limits = c(0, 7500),
                        breaks = seq(0,7500,by = 500)) 

ggsave(plot = aa_distribution_plot,
       file = 'polished_plots/aa_distribution.png',
       width = 14)


# Write summary file
filecon <- file('summary.txt')
writeLines(summary_content, filecon)
close(filecon)




