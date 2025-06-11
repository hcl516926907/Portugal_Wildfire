dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Porgual_Wildfire"
dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Portugal_Wildfire"

library(ggplot2)

library(dplyr)
library(RColorBrewer)

library(xgboost)
library(rBayesianOptimization)
library(caret) 
library(fastDummies)
library(zoo)
# bru_safe_sp(force = TRUE)

load(file.path(dir.data, "Data_For_Fitting_Council.RData"))
data.fit <- data.fit.council
data.fit$ba <- exp(data.fit$log_ba)
data.fit <- data.fit %>% replace(is.na(.), 0)

write.csv(data.fit, "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Porgual_Wildfire/data_fit.csv", row.names = FALSE, fileEncoding = "UTF-8")

data.fit$xgb_cnt <- data.fit$y
data.fit$xgb_ba <- sqrt(data.fit$ba)



w <- 9
h <- 3


df <- data.fit %>%
  arrange(NAME_2, time.idx) %>%   # Ensure data is ordered by location and time
  group_by(NAME_2)                # Apply the lag per location

# Loop over the required lag steps and add columns dynamically
for(i in 0:(w  - 1)) {
  df <- df %>%
    mutate(!!paste0("ba_lag", i) := lag(xgb_ba, n = i),
           !!paste0("cnt_lag", i) := lag(xgb_cnt, n = i))
}

for(i in 1:h) {
  df <- df %>%
    mutate(!!paste0("xgb_ba_h", i) := lead(xgb_ba, n = i),
           !!paste0("xgb_cnt_h", i) := lead(xgb_cnt, n = i))
}

df <- ungroup(df)



dist.stat <- data.fit %>% group_by(NAME_1, time.idx) %>% summarise(dist_xgb_ba = sqrt(mean(ba)), 
                                                                   dist_xgb_cnt = mean(y)) 

dist.stat %>% 
  arrange(NAME_1, time.idx) %>%   # Ensure data is ordered by location and time
  group_by(NAME_1)  

for(i in 0:(w  - 1)) {
  dist.stat <- dist.stat %>%
    mutate(!!paste0("dist_ba_lag", i) := lag(dist_xgb_ba, n = i),
           !!paste0("dist_cnt_lag", i) := lag(dist_xgb_cnt, n = i))
}

dist.stat <- ungroup(dist.stat)



######################################## month ######################
df <- df %>%
  mutate(
    month_sin = sin(2 * pi * month / 12),
    month_cos = cos(2 * pi * month / 12)
  )

##################################### council level moving average #########################

df <- df %>% arrange(NAME_2, time.idx) %>%
  group_by(NAME_2) %>%
  mutate(
    ba_ma_3 = rollmean(xgb_ba, k = 3, fill = NA, align = "right"),
    ba_ma_6 = rollmean(xgb_ba, k = 6, fill = NA, align = "right"),
    ba_ma_9 = rollmean(xgb_ba, k = 9, fill = NA, align = "right"),
    ba_ma_12 = rollmean(xgb_ba, k = 12, fill = NA, align = "right"),
    ba_ma_24 = rollmean(xgb_ba, k = 24, fill = NA, align = "right"),
    ba_ma_36 = rollmean(xgb_ba, k = 36, fill = NA, align = "right"),
  ) %>%
  ungroup()

# df[,c('NAME_2','time.idx', 'xgb_ba','ba_ma_3','ba_ma_6')]

df <- df %>% arrange(NAME_2, time.idx) %>%
  group_by(NAME_2) %>%
  mutate(
    cnt_ma_3 = rollmean(xgb_cnt, k = 3, fill = NA, align = "right"),
    cnt_ma_6 = rollmean(xgb_cnt, k = 6, fill = NA, align = "right"),
    cnt_ma_9 = rollmean(xgb_cnt, k = 9, fill = NA, align = "right"),
    cnt_ma_12 = rollmean(xgb_cnt, k = 12, fill = NA, align = "right"),
    cnt_ma_24 = rollmean(xgb_cnt, k = 24, fill = NA, align = "right"),
    cnt_ma_36 = rollmean(xgb_cnt, k = 36, fill = NA, align = "right"),
  ) %>%
  ungroup()


# df[,c('NAME_2','time.idx', 'xgb_cnt','cnt_ma_3','cnt_ma_6')]

##################################### district level moving average #########################
dist.stat <- dist.stat %>% arrange(NAME_1, time.idx) %>%
  group_by(NAME_1) %>%
  mutate(
    dist_ba_ma_3 = rollmean(dist_xgb_ba, k = 3, fill = NA, align = "right"),
    dist_ba_ma_6 = rollmean(dist_xgb_ba, k = 6, fill = NA, align = "right"),
    dist_ba_ma_9 = rollmean(dist_xgb_ba, k = 9, fill = NA, align = "right"),
    dist_ba_ma_12 = rollmean(dist_xgb_ba, k = 12, fill = NA, align = "right"),
    dist_ba_ma_24 = rollmean(dist_xgb_ba, k = 24, fill = NA, align = "right"),
    dist_ba_ma_36 = rollmean(dist_xgb_ba, k = 36, fill = NA, align = "right"),
  ) %>%
  ungroup()

# dist.stat[,c('NAME_1','time.idx', 'dist_xgb_ba','dist_ba_lag1','dist_ba_ma_3','dist_ba_ma_6')]

dist.stat <- dist.stat %>% arrange(NAME_1, time.idx) %>%
  group_by(NAME_1) %>%
  mutate(
    dist_cnt_ma_3 = rollmean(dist_xgb_cnt, k = 3, fill = NA, align = "right"),
    dist_cnt_ma_6 = rollmean(dist_xgb_cnt, k = 6, fill = NA, align = "right"),
    dist_cnt_ma_9 = rollmean(dist_xgb_cnt, k = 9, fill = NA, align = "right"),
    dist_cnt_ma_12 = rollmean(dist_xgb_cnt, k = 12, fill = NA, align = "right"),
    dist_cnt_ma_24 = rollmean(dist_xgb_cnt, k = 24, fill = NA, align = "right"),
    dist_cnt_ma_36 = rollmean(dist_xgb_cnt, k = 36, fill = NA, align = "right"),
  ) %>%
  ungroup()
# dist.stat[,c('NAME_1','time.idx', 'dist_xgb_cnt','dist_cnt_lag1','dist_cnt_ma_3','dist_cnt_ma_6')]

##################################### council level historical information 1 #########################
# Only values with a positive index are used.
df <- df %>% arrange(NAME_2, time.idx) %>%
  group_by(NAME_2) %>%
  mutate(hist_ba_h1_1m = sapply(time.idx, function(t) {
    kmax <- floor((t - 11) / 12)  # maximum k that gives a valid index
    if(kmax < 0) return(NA_real_)
    # Generate indices for all k
    inds <- sapply(0:kmax, function(k) t - 11 - 12 * k)
    # Extract the corresponding values from xgb_ba (relying on order)
    valid_vals <- xgb_ba[inds]
    if(length(valid_vals) == 0) return(NA_real_)
    mean(valid_vals, na.rm = TRUE)
  })) %>%
  ungroup()

df <- df %>% arrange(NAME_2, time.idx) %>%
  group_by(NAME_2) %>%
  mutate(hist_cnt_h1_1m = sapply(time.idx, function(t) {
    kmax <- floor((t - 11) / 12)  # maximum k that gives a valid index
    if(kmax < 0) return(NA_real_)
    # Generate indices for all k
    inds <- sapply(0:kmax, function(k) t - 11 - 12 * k)
    # Extract the corresponding values from xgb_ba (relying on order)
    valid_vals <- xgb_cnt[inds]
    if(length(valid_vals) == 0) return(NA_real_)
    mean(valid_vals, na.rm = TRUE)
  })) %>%
  ungroup()

# as.data.frame(df[,c('NAME_2','time.idx','year','month',  'xgb_cnt', "hist_cnt_h1_1m")])

##################################### council level historical information 1 #########################
dist.stat <- dist.stat %>% arrange(NAME_1, time.idx) %>%
  group_by(NAME_1) %>%
  mutate(dist_hist_ba_h1_1m = sapply(time.idx, function(t) {
    kmax <- floor((t - 11) / 12)  # maximum k that gives a valid index
    if(kmax < 0) return(NA_real_)
    # Generate indices for all k
    inds <- sapply(0:kmax, function(k) t - 11 - 12 * k)
    # Extract the corresponding values from xgb_ba (relying on order)
    valid_vals <- dist_xgb_ba[inds]
    if(length(valid_vals) == 0) return(NA_real_)
    mean(valid_vals, na.rm = TRUE)
  })) %>%
  ungroup()

dist.stat <- dist.stat %>% arrange(NAME_1, time.idx) %>%
  group_by(NAME_1) %>%
  mutate(dist_hist_cnt_h1_1m = sapply(time.idx, function(t) {
    kmax <- floor((t - 11) / 12)  # maximum k that gives a valid index
    if(kmax < 0) return(NA_real_)
    # Generate indices for all k
    inds <- sapply(0:kmax, function(k) t - 11 - 12 * k)
    # Extract the corresponding values from xgb_ba (relying on order)
    valid_vals <- dist_xgb_cnt[inds]
    if(length(valid_vals) == 0) return(NA_real_)
    mean(valid_vals, na.rm = TRUE)
  })) %>%
  ungroup()

# as.data.frame(dist.stat[,c('NAME_1','time.idx', 'dist_xgb_cnt', "dist_hist_ba_h1_1m")])



##################################### council level historical information 2 #########################


### 4. Create Historical Average 2
# For each time.idx, compute the mean over indices: 
# time.idx - 11 - 12*k, time.idx - 10 - 12*k, and time.idx - 12 - 12*k for k = 0, 1, 2, ...
df <- df %>% arrange(NAME_2, time.idx) %>%
  group_by(NAME_2) %>%
  mutate(hist_ba_h1_3m = sapply(time.idx, function(t) {
    # Use the smallest offset (10) to determine maximum valid k.
    kmax <- floor((t - 10) / 12)
    if(kmax < 0) return(NA_real_)
    # For each k, generate the three indices and combine them
    inds <- unlist(lapply(0:kmax, function(k) {
      t - c(11, 10, 12) - 12 * k
    }))
    # Keep only positive indices
    inds <- inds[inds > 0]
    valid_vals <- xgb_ba[inds]
    if(length(valid_vals) == 0) return(NA_real_)
    mean(valid_vals, na.rm = TRUE)
  })) %>%
  ungroup()

df <- df %>% arrange(NAME_2, time.idx) %>%
  group_by(NAME_2) %>%
  mutate(hist_cnt_h1_3m = sapply(time.idx, function(t) {
    # Use the smallest offset (10) to determine maximum valid k.
    kmax <- floor((t - 10) / 12)
    if(kmax < 0) return(NA_real_)
    # For each k, generate the three indices and combine them
    inds <- unlist(lapply(0:kmax, function(k) {
      t - c(11, 10, 12) - 12 * k
    }))
    # Keep only positive indices
    inds <- inds[inds > 0]
    valid_vals <- xgb_cnt[inds]
    if(length(valid_vals) == 0) return(NA_real_)
    mean(valid_vals, na.rm = TRUE)
  })) %>%
  ungroup()
# as.data.frame(df[,c('NAME_2','time.idx','year','month',  'xgb_cnt', "hist_cnt_h1_3m")])




##################################### district level historical information 2 #########################

dist.stat <- dist.stat %>% arrange(NAME_1, time.idx) %>%
  group_by(NAME_1) %>%
  mutate(dist_hist_ba_h1_3m = sapply(time.idx, function(t) {
    # Use the smallest offset (10) to determine maximum valid k.
    kmax <- floor((t - 10) / 12)
    if(kmax < 0) return(NA_real_)
    # For each k, generate the three indices and combine them
    inds <- unlist(lapply(0:kmax, function(k) {
      t - c(11, 10, 12) - 12 * k
    }))
    # Keep only positive indices
    inds <- inds[inds > 0]
    valid_vals <- dist_xgb_ba[inds]
    if(length(valid_vals) == 0) return(NA_real_)
    mean(valid_vals, na.rm = TRUE)
  })) %>%
  ungroup()

dist.stat <- dist.stat %>% arrange(NAME_1, time.idx) %>%
  group_by(NAME_1) %>%
  mutate(dist_hist_cnt_h1_3m = sapply(time.idx, function(t) {
    # Use the smallest offset (10) to determine maximum valid k.
    kmax <- floor((t - 10) / 12)
    if(kmax < 0) return(NA_real_)
    # For each k, generate the three indices and combine them
    inds <- unlist(lapply(0:kmax, function(k) {
      t - c(11, 10, 12) - 12 * k
    }))
    # Keep only positive indices
    inds <- inds[inds > 0]
    valid_vals <- dist_xgb_cnt[inds]
    if(length(valid_vals) == 0) return(NA_real_)
    mean(valid_vals, na.rm = TRUE)
  })) %>%
  ungroup()
# as.data.frame(dist.stat[,c('NAME_1','time.idx', 'dist_xgb_cnt', "dist_hist_cnt_h1_3m")])



##################################### council level historical information 3 #########################


### 5. Create Historical Average 3
# For each time.idx, compute the mean over indices:
# time.idx - 11 - 12*k, time.idx - 10 - 12*k, time.idx - 12 - 12*k, time.idx - 9 - 12*k, and time.idx - 13 - 12*k for k = 0, 1, 2, ...
df <- df %>% arrange(NAME_2, time.idx) %>%
  group_by(NAME_2) %>%
  mutate(hist_ba_h1_5m = sapply(time.idx, function(t) {
    # Use the smallest offset (9) to determine the maximum valid k.
    kmax <- floor((t - 9) / 12)
    if(kmax < 0) return(NA_real_)
    # For each k, generate the five indices and combine them
    inds <- unlist(lapply(0:kmax, function(k) {
      t - c(11, 10, 12, 9, 13) - 12 * k
    }))
    # Only consider indices that are positive
    inds <- inds[inds > 0]
    valid_vals <- xgb_ba[inds]
    if(length(valid_vals) == 0) return(NA_real_)
    mean(valid_vals, na.rm = TRUE)
  })) %>%
  ungroup()


df <- df %>% arrange(NAME_2, time.idx) %>%
  group_by(NAME_2) %>%
  mutate(hist_cnt_h1_5m = sapply(time.idx, function(t) {
    # Use the smallest offset (9) to determine the maximum valid k.
    kmax <- floor((t - 9) / 12)
    if(kmax < 0) return(NA_real_)
    # For each k, generate the five indices and combine them
    inds <- unlist(lapply(0:kmax, function(k) {
      t - c(11, 10, 12, 9, 13) - 12 * k
    }))
    # Only consider indices that are positive
    inds <- inds[inds > 0]
    valid_vals <- xgb_cnt[inds]
    if(length(valid_vals) == 0) return(NA_real_)
    mean(valid_vals, na.rm = TRUE)
  })) %>%
  ungroup()

##################################### district level historical information 3 #########################

dist.stat <- dist.stat %>% arrange(NAME_1, time.idx) %>%
  group_by(NAME_1) %>%
  mutate(dist_hist_ba_h1_5m = sapply(time.idx, function(t) {
    # Use the smallest offset (9) to determine the maximum valid k.
    kmax <- floor((t - 9) / 12)
    if(kmax < 0) return(NA_real_)
    # For each k, generate the five indices and combine them
    inds <- unlist(lapply(0:kmax, function(k) {
      t - c(11, 10, 12, 9, 13) - 12 * k
    }))
    # Only consider indices that are positive
    inds <- inds[inds > 0]
    valid_vals <- dist_xgb_ba[inds]
    if(length(valid_vals) == 0) return(NA_real_)
    mean(valid_vals, na.rm = TRUE)
  })) %>%
  ungroup()


dist.stat <- dist.stat %>% arrange(NAME_1, time.idx) %>%
  group_by(NAME_1) %>%
  mutate(dist_hist_cnt_h1_5m = sapply(time.idx, function(t) {
    # Use the smallest offset (9) to determine the maximum valid k.
    kmax <- floor((t - 9) / 12)
    if(kmax < 0) return(NA_real_)
    # For each k, generate the five indices and combine them
    inds <- unlist(lapply(0:kmax, function(k) {
      t - c(11, 10, 12, 9, 13) - 12 * k
    }))
    # Only consider indices that are positive
    inds <- inds[inds > 0]
    valid_vals <- dist_xgb_cnt[inds]
    if(length(valid_vals) == 0) return(NA_real_)
    mean(valid_vals, na.rm = TRUE)
  })) %>%
  ungroup()



dist.stat[,c('dist_ba_lag1')]

df <- merge(df, dist.stat, by=c('NAME_1','time.idx'))


df <- df %>% arrange(NAME_2, time.idx)

# as.data.frame(df[,c( 'NAME_1','NAME_2','time.idx','year','month',  'xgb_cnt', "hist_cnt_h1_5m",'dist_cnt_lag1')])
# df[df$NAME_1=='Santarém' & df$time.idx==1,c( 'NAME_1','NAME_2','time.idx','year','month',  'xgb_cnt', "hist_cnt_h1_5m",'dist_cnt_lag1')]


df$NAME_1 <- as.character(df$NAME_1)

df$HVegTyp <- as.factor(df$HVegTyp)
# levels(data.fit$HVegTyp) <-  c('NoCover','Evergreen_broadleaf_trees','Mixed_forest','Interrupted_forest')
df$LVegTyp <- as.factor(df$LVegTyp)
# levels(data.fit$LVegTyp) <- c('NoCover','Crops','Tall_grass','Semidesert','Evergreen_shrubs')
df <- dummy_cols(df,   select_columns = c('HVegTyp','LVegTyp'),remove_first_dummy=TRUE)


################################### define training set and test set #################

df <- as.data.frame(df)
ba.covar.names.h1 <- c( paste('ba_lag',0:(w-1),sep=''),
                  paste('dist_ba_lag',0:(w-1),sep=''),
                  'ba_ma_3', 'ba_ma_6','ba_ma_9','ba_ma_12','ba_ma_24','ba_ma_36',
                  'dist_ba_ma_3', 'dist_ba_ma_6','dist_ba_ma_9','dist_ba_ma_12','dist_ba_ma_24','dist_ba_ma_36',
                  'hist_ba_h1_1m','hist_ba_h1_3m','hist_ba_h1_5m',
                  'dist_hist_ba_h1_1m','dist_hist_ba_h1_3m','dist_hist_ba_h1_5m',

                 'lon.grid','lat.grid','year',
                 'month_sin','month_cos',
                 'FWI','HVegCov','HVegLAI','LVegCov','LVegLAI',
                 'Pricp','RHumi','Temp','UComp','VComp'
)
# df[,ba.covar.names.h1]

cnt.covar.names.h1 <- c( paste('cnt_lag',0:(w-1),sep=''),
                     paste('dist_cnt_lag',0:(w-1),sep=''),
                     'cnt_ma_3', 'cnt_ma_6','cnt_ma_9','cnt_ma_12','cnt_ma_24','cnt_ma_36',
                     'dist_cnt_ma_3', 'dist_cnt_ma_6','dist_cnt_ma_9','dist_cnt_ma_12','dist_cnt_ma_24','dist_cnt_ma_36',
                     'hist_cnt_h1_1m','hist_cnt_h1_3m','hist_cnt_h1_5m',
                     'dist_hist_cnt_h1_1m','dist_hist_cnt_h1_3m','dist_hist_cnt_h1_5m',
                     
                     'lon.grid','lat.grid','year',
                     'month_sin','month_cos',
                     'FWI','HVegCov','HVegLAI','LVegCov','LVegLAI',
                     'Pricp','RHumi','Temp','UComp','VComp'
)
# df[,cnt.covar.names.h1]

df = na.omit(df)

data.fit.train <- df[(df$year <= 2022) ,]

data.fit.test <- df[df$year > 2022,]

acf(data.fit.train$xgb_ba,lag.max=100)
acf(data.fit.train$xgb_cnt,lag.max=100)


overall_ts <- data.fit.train %>%
  group_by(time.idx) %>%
  summarize(mean_xgb_ba = mean(xgb_ba, na.rm = TRUE)) %>%
  arrange(time.idx)

# 2. Compute the ACF (does NOT plot)
acf_res <- acf(overall_ts$mean_xgb_ba, plot = FALSE, lag.max=70 )

# 3. Turn it into a data frame
acf_df <- with(acf_res, 
               data.frame(
                 lag = as.numeric(lag),
                 acf = as.numeric(acf)
               )) %>%
  filter(lag >= 0)    # keep non‐negative lags only

# 4. Plot with ggplot2
p <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_col(width = 0.1) +
  labs(x     = "Lag",
       y     = "ACF of Burn Area") +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15),
    axis.title.x=element_text(size=15),
    axis.title.y=element_text(size=15)) 

png(filename = file.path(dir.out, "ACF_BA.pdf"), width =4000 , height = 2000, res=300)
print(p)
dev.off()



overall_ts <- data.fit.train %>%
  group_by(time.idx) %>%
  summarize(mean_xgb_cnt = mean(xgb_cnt, na.rm = TRUE)) %>%
  arrange(time.idx)

# 2. Compute the ACF (does NOT plot)
acf_res <- acf(overall_ts$mean_xgb_cnt, plot = FALSE, lag.max=70 )

# 3. Turn it into a data frame
acf_df <- with(acf_res, 
               data.frame(
                 lag = as.numeric(lag),
                 acf = as.numeric(acf)
               )) %>%
  filter(lag >= 0)    # keep non‐negative lags only

# 4. Plot with ggplot2
ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_col(width = 0.1) +
  labs(x     = "Lag",
       y     = "ACF of Fire Count") +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15),
    axis.title.x=element_text(size=15),
    axis.title.y=element_text(size=15)) 
# Assuming your dataset 'df.train' is already available in the environment

png(filename = file.path(dir.out, "ACF_Count.pdf"), width =4000 , height = 2000, res=300)
print(p)
dev.off()


XGBoostAR <- function(data, covar.names, target.name, eta){
  
  
  # Prepare the dataset: separate predictors and response
  label <- data[,target.name]
  features <- data[, covar.names]
  dtrain <- xgb.DMatrix(data = as.matrix(features), label = label)
  # Define the objective function for Bayesian optimization
  # This function performs 6-fold cross-validation with early stopping and returns the negative RMSE.
  # (BayesianOptimization maximizes the objective, so we use -RMSE.)
  xgb_cv_bayes <- function(max_depth, subsample, colsample_bytree) {
    
    params <- list(
      booster = "gbtree",
      eta = eta,  # Learning rate; adjust if necessary.
      max_depth = max_depth,
      subsample = subsample,
      colsample_bytree = colsample_bytree,
      objective = "reg:squarederror",  # For a continuous regression target.
      eval_metric = "rmse"
    )
    
    cv_res <- xgb.cv(
      params = params,
      data = dtrain,
      nrounds = 1000,
      nfold = 6,
      early_stopping_rounds = 10,
      verbose = 0,
      maximize = FALSE # Lower RMSE is better.
    )
    
    best_rmse <- min(cv_res$evaluation_log$test_rmse_mean)
    best_iter <- cv_res$best_iteration
    
    # Return a list with Score (to be maximized) and the best number of rounds.
    list(Score = -best_rmse, Pred = 0)
  }
  
  
  # Run Bayesian optimization to tune max_depth, subsample, and colsample_bytree
  opt_results <- BayesianOptimization(
    FUN = xgb_cv_bayes,
    bounds = list(
      max_depth = c(7L, 13L),           # Tune over reasonable integer values (converted inside the function).
      subsample = c(0.3, 0.9),        # Proportion of observations to sample.
      colsample_bytree = c(0.3, 0.9)   # Proportion of features to sample.
    ),
    init_points = 5,  # Initial random explorations.
    n_iter = 25,      # Number of iterations for Bayesian optimization.
    acq = "ucb",      # Acquisition function (upper confidence bound).
    kappa = 2.576,
    eps = 0.0,
    verbose = TRUE
  )
  
  # Extract the best hyperparameters and the corresponding best nrounds from the optimization history.
  best_params <- opt_results$Best_Par
  
  # Here we extract the nrounds corresponding to the best score observed in the optimization.
  final_params <- list(
    booster = "gbtree",
    eta = eta,
    max_depth = best_params['max_depth'],
    subsample = best_params['subsample'],
    colsample_bytree = best_params['colsample_bytree'],
    objective = "reg:squarederror",
    eval_metric = "rmse"
  )
  
  
  cv_final <- xgb.cv(
    params = final_params,
    data = dtrain,
    nrounds = 1000,
    nfold = 6,
    early_stopping_rounds = 10,
    verbose = 0,
    maximize = FALSE
  )
  
  best_nrounds <- cv_final$best_iteration
  cat("Best number of rounds: ", best_nrounds, "\n")
  
  
  # If not already present, create a column to store predictions on the validation set
  data$score <- NA
  
  # Create 6 folds using the caret package (returns indices for the validation set in each fold)
  folds <- createFolds(data[,target.name], k = 6, list = TRUE, returnTrain = FALSE)
  
  # Loop through each fold to train the model and generate predictions on the held-out set
  for (i in seq_along(folds)) {
    
    # Get validation indices for the current fold
    val_idx <- folds[[i]]
    
    # The training indices are all rows not in the validation fold
    train_idx <- setdiff(seq_len(nrow(data)), val_idx)
    
    # Prepare training data: remove the response and prediction column from predictors
    predictors <- covar.names
    
    dtrain <- xgb.DMatrix(data = as.matrix(data[train_idx, predictors]), 
                          label = data[, target.name][train_idx])
    dval <- xgb.DMatrix(data = as.matrix(data[val_idx, predictors]), 
                        label = data[, target.name][val_idx])
    
  
    model <- xgb.train(
      params = final_params,
      data = dtrain,
      nrounds = best_nrounds,
      watchlist = list(val = dval),
      early_stopping_rounds = 10,
      verbose = 0   # change to 1 if you want progress output
    )
    
    # Predict on the validation set
    preds <- predict(model, newdata = dval)
    
    # Save the predictions back into the original data as the out-of-fold prediction
    data$score[val_idx] <- preds
  }

  return(list('best_params'=final_params, 'best_nrounds'=best_nrounds, 'score'=data$score))
}

set.seed(123)
res_cnt_h1 = XGBoostAR(as.data.frame(data.fit.train), cnt.covar.names.h1, 'xgb_cnt_h1', eta=0.01)
save(res_cnt_h1, 
     file=file.path(dir.out,'AutoRegressive_XGBoost_Hyperparameters_CNT.RData'))

res_ba_h1 = XGBoostAR(as.data.frame(data.fit.train), ba.covar.names.h1, 'xgb_ba_h1', eta=0.005 )
save(res_ba_h1, 
     file=file.path(dir.out,'AutoRegressive_XGBoost_Hyperparameters_BA.RData'))

data.fit.train$pred_ba_h1 <- res_ba_h1$score
data.fit.train$pred_cnt_h1 <- res_cnt_h1$score

acf(data.fit.train$xgb_ba)
# Check the first few predictions
head(data.fit.train$pred_ba_h1)

data <- data.fit.train
council = 'Odemira'
plot(1:(dim(data[data$NAME_2==council,])[1]), data[data$NAME_2==council, ]$xgb_ba,type='l')
lines(data[data$NAME_2==council, ]$pred_ba_h1,col='red')

# Train the final XGBoost model using the tuned hyperparameters and best number of rounds.


final_model <- function(data, covar.names, target.name, tuning_res){

  label <- data[,target.name]
  features <- data[, covar.names]
  dtrain <- xgb.DMatrix(data = as.matrix(features), label = label)
  
  model <- xgb.train(
    params = tuning_res$best_params,
    data = dtrain,
    nrounds = tuning_res$best_nrounds ,
    verbose = 0  )
  return(model)
}

# Optionally, inspect your final model or its performance.

model_cnt_h1 <- final_model(as.data.frame(data.fit.train), cnt.covar.names.h1, 'xgb_cnt_h1',res_cnt_h1)
model_ba_h1 <- final_model(as.data.frame(data.fit.train), ba.covar.names.h1, 'xgb_ba_h1', res_ba_h1)



print(xgb.importance(model = model_cnt_h1))

print(xgb.importance(model = model_ba_h1))


data.fit.test$pred_ba_h1 <- predict(model_ba_h1, newdata = as.matrix(data.fit.test[,ba.covar.names.h1]))

data.fit.test$pred_cnt_h1 <- predict(model_cnt_h1, newdata = as.matrix(data.fit.test[,cnt.covar.names.h1]))


plot(1:(dim(data.fit.test[data.fit.test$NAME_2==council,])[1]), data.fit.test[data.fit.test$NAME_2==council, ]$xgb_cnt,type='l')
lines(data.fit.test[data.fit.test$NAME_2==council, ]$pred_cnt_h1,col='red')



data.fit <- rbind(data.fit.train, data.fit.test)

plot(1:(dim(data.fit[data.fit$NAME_2==council,])[1]), data.fit[data.fit$NAME_2==council, ]$xgb_cnt,type='l')
lines(data.fit[data.fit$NAME_2==council, ]$pred_cnt_h1,col='red')

plot(1:(dim(data.fit[data.fit$NAME_2==council,])[1]), data.fit[data.fit$NAME_2==council, ]$xgb_ba,type='l')
lines(data.fit[data.fit$NAME_2==council, ]$pred_ba_h1,col='red')


save(data.fit, 
     file=file.path(dir.out,'AutoRegressive_XGBoost_Predictions.RData'))