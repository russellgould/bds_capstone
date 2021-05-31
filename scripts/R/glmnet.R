#!/usr/bin/env Rscript

library(dplyr)
library(randomForest)
library(caret)
library(glmnet)

args <- commandArgs(trailingOnly = TRUE)

# input phyloseq object
ps_object_path <- file.path(args[1])

# where to save glmnet stuff
output_path <- file.path(args[2])
dir.create(output_path)

objects_output_path <- file.path(output_path, "glmnet_output")
dir.create(objects_output_path)

data_test_path <- file.path(objects_output_path, "data_test.Rdata")
fit_object_path <- file.path(objects_output_path, "fit.Rdata")

# load ps object
load(ps_object_path)

# format dataset
depth <- sample_data(ps)$Depth
lat <- sample_data(ps)$Latitude
hab <- sample_data(ps)$Habitat
data <- otu_table(ps)
data <- data.frame(depth, data)
# bin depth to create categorical variables
data <- data %>% mutate(Bins = cut(depth, breaks = c(-Inf, 400, 700, 1000, Inf)))
# data$Bins <- recode_factor(as.factor(data$Bins) ('(-Inf,400]' = "shallow",
#                                                  '(400,700]' = "shal_middle",
#                                                  '(700,1e+03]' = "middle_deep",
#                                                  '(1e+03, Inf]' = "deep"))
# rename bins
levels(data$Bins) <- c("shallow", "shal_middle", "middle_deep", "deep")

# partition data into train/test sets
inTraining <- createDataPartition(data$depth, p = .75, list = FALSE)
data_train <- data[inTraining, ]
data_test <- data[-inTraining, ]

save(data_test, data_test_path)

# create x and y for glm
x_train <- data_train[, -1]
y_train <- data_train[, 1]
fit <- glmnet(x_train, y_train) # can add more parameters to make more specific to this dataset

save(fit, fit_object_path)
