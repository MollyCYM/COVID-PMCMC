library(fs)
library(purrr)
library(rlang)
library(whisker)

# Read in template script
template <- readLines("para_template.R") #prepare base file
#output <- whisker.render(template, template_variables) #output store file? what is templete_variables

# Parameter options //list all the tunning parameters
param_options <- list(
    initial_z = c(0, 1, 2),
    dxdt = c("sigma*e", "-sigma*e")
)

# Produce combinations
combinations_df <- expand.grid(param_options)

#' Function to generate a simulation file from template  // not understand yet
generate <- function(param_set) {
    # Put parameters into template
    output <- whisker.render(template, param_set)

    path_reduce <- function(path, part) {
        file.path(path, path_sanitize(part))
    }

    # Create output directory  //for example /media/yc22/SEAGATE_1/MollyCOVID/data  para3_alpha28.csv
    path <- reduce(param_set, path_reduce, .init = ".")
    filename <- file.path(path, "para.R")
    dir.create(path, recursive = TRUE)

    # Save output R script
    file_conn <- file(filename)
    writeLines(output, file_conn)
    close(file_conn)
}

# Apply above function to dataframe of parameter combinations
apply(combinations_df, 1, generate)

#How to run them on Create
