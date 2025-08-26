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

library(lubridate)
library(pROC)

# bru_safe_sp(force = TRUE)

load(file.path(dir.data, "Data_For_Fitting_Council.RData"))
data.fit <- data.fit.council
data.fit$ba <- exp(data.fit$log_ba)
data.fit <- data.fit %>% replace(is.na(.), 0)

write.csv(data.fit, "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Porgual_Wildfire/data_fit.csv", row.names = FALSE, fileEncoding = "UTF-8")

data.fit$xgb_cnt <- data.fit$y
data.fit$xgb_ba <- sqrt(data.fit$ba)
data.fit$xgb_z <- as.vector(data.fit$y>0) + 0

####################################################################
#Create autoregressive covariates 
####################################################################
w <- 9
h <- 1


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
           !!paste0("xgb_cnt_h", i) := lead(xgb_cnt, n = i),
           !!paste0("xgb_z_h", i) := lead(xgb_z, n = i))
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



#################### angular representation of month ######################
df <- df %>%
  mutate(
    month_sin = sin(2 * pi * month / 12),
    month_cos = cos(2 * pi * month / 12)
  )

#################### council level moving average #########################

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

###################### district level moving average #########################
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

#################### council level historical information 1 #########################
# Only values with a positive index are used.
df <- df %>% arrange(NAME_2, time.idx) %>%
  group_by(NAME_2) %>%
  mutate(hist_ba_h1_1m = sapply(time.idx, function(t) {
    kmax <- floor((t - 11) / 12)  # maximum k that gives a valid index
    if(kmax < 0) return(NA_real_)
    # Generate indices for all k
    kmax <- min(kmax,2)
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
    kmax <- min(kmax,2)
    inds <- sapply(0:kmax, function(k) t - 11 - 12 * k)
    # Extract the corresponding values from xgb_ba (relying on order)
    valid_vals <- xgb_cnt[inds]
    if(length(valid_vals) == 0) return(NA_real_)
    mean(valid_vals, na.rm = TRUE)
  })) %>%
  ungroup()


###################### district level historical information 1 ####################
dist.stat <- dist.stat %>% arrange(NAME_1, time.idx) %>%
  group_by(NAME_1) %>%
  mutate(dist_hist_ba_h1_1m = sapply(time.idx, function(t) {
    kmax <- floor((t - 11) / 12)  # maximum k that gives a valid index
    if(kmax < 0) return(NA_real_)
    # Generate indices for all k
    kmax <- min(kmax,2)
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
    kmax <- min(kmax,2)
    inds <- sapply(0:kmax, function(k) t - 11 - 12 * k)
    # Extract the corresponding values from xgb_ba (relying on order)
    valid_vals <- dist_xgb_cnt[inds]
    if(length(valid_vals) == 0) return(NA_real_)
    mean(valid_vals, na.rm = TRUE)
  })) %>%
  ungroup()

#################### council level historical information 2 #########################

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
    kmax <- min(kmax,2)
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
    kmax <- min(kmax,2)
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


##################################### district level historical information 2 #########################

dist.stat <- dist.stat %>% arrange(NAME_1, time.idx) %>%
  group_by(NAME_1) %>%
  mutate(dist_hist_ba_h1_3m = sapply(time.idx, function(t) {
    # Use the smallest offset (10) to determine maximum valid k.
    kmax <- floor((t - 10) / 12)
    if(kmax < 0) return(NA_real_)
    # For each k, generate the three indices and combine them
    kmax <- min(kmax,2)
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
    kmax <- min(kmax,2)
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



###################### council level historical information 3 #########################


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
    kmax <- min(kmax,2)
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
    kmax <- min(kmax,2)
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

###################### district level historical information 3 #########################

dist.stat <- dist.stat %>% arrange(NAME_1, time.idx) %>%
  group_by(NAME_1) %>%
  mutate(dist_hist_ba_h1_5m = sapply(time.idx, function(t) {
    # Use the smallest offset (9) to determine the maximum valid k.
    kmax <- floor((t - 9) / 12)
    if(kmax < 0) return(NA_real_)
    # For each k, generate the five indices and combine them
    kmax <- min(kmax,2)
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
    kmax <- min(kmax,2)
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


df$NAME_1 <- as.character(df$NAME_1)

df$HVegTyp <- as.factor(df$HVegTyp)
# levels(data.fit$HVegTyp) <-  c('NoCover','Evergreen_broadleaf_trees','Mixed_forest','Interrupted_forest')
df$LVegTyp <- as.factor(df$LVegTyp)
# levels(data.fit$LVegTyp) <- c('NoCover','Crops','Tall_grass','Semidesert','Evergreen_shrubs')
df <- dummy_cols(df,   select_columns = c('HVegTyp','LVegTyp','NAME_1'),remove_first_dummy=TRUE)


################################### define training set and test set #################

df <- as.data.frame(df)
ba.covar.names.h1 <- c(
                  paste('ba_lag',0:(w-1),sep=''),
                  paste('dist_ba_lag',0:(w-1),sep=''),
                  'ba_ma_3', 'ba_ma_6','ba_ma_9',
                  'ba_ma_12', 'ba_ma_24','ba_ma_36',
                  'dist_ba_ma_3', 'dist_ba_ma_6','dist_ba_ma_9',
                  'dist_ba_ma_12','dist_ba_ma_24','dist_ba_ma_36',
                  'hist_ba_h1_1m','hist_ba_h1_3m','hist_ba_h1_5m',
                  'dist_hist_ba_h1_1m','dist_hist_ba_h1_3m','dist_hist_ba_h1_5m',

                  'lon.grid','lat.grid','year',
                  'FWI','HVegCov','HVegLAI','LVegCov','LVegLAI',
                  'Pricp','RHumi','Temp','UComp','VComp','DewPoint',
                  'HVegTyp_6','HVegTyp_18','HVegTyp_19',
                  'LVegTyp_1','LVegTyp_7','LVegTyp_11','LVegTyp_16',
                  'month_sin','month_cos'
)

cnt.covar.names.h1 <- c( 
                     paste('cnt_lag',0:(w-1),sep=''),
                     paste('dist_cnt_lag',0:(w-1),sep=''),
                     'cnt_ma_3', 'cnt_ma_6','cnt_ma_9',
                     'cnt_ma_12', 'cnt_ma_24','cnt_ma_36',
                     'dist_cnt_ma_3', 'dist_cnt_ma_6','dist_cnt_ma_9',
                     'dist_cnt_ma_12','dist_cnt_ma_24','dist_cnt_ma_36',
                     'hist_cnt_h1_1m','hist_cnt_h1_3m','hist_cnt_h1_5m',
                     'dist_hist_cnt_h1_1m','dist_hist_cnt_h1_3m','dist_hist_cnt_h1_5m',

                     'lon.grid','lat.grid','year',
                     'FWI','HVegCov','HVegLAI','LVegCov','LVegLAI',
                     'Pricp','RHumi','Temp','UComp','VComp','DewPoint',
                     'HVegTyp_6','HVegTyp_18','HVegTyp_19',
                     'LVegTyp_1','LVegTyp_7','LVegTyp_11','LVegTyp_16',
                     'month_sin','month_cos'
)


combined.covar.names.h1 <- unique(c(cnt.covar.names.h1, ba.covar.names.h1))

# time <= 2022-12
data.fit.train <- df[(df$time.idx <= 144) & (df$year >= 2014),]

data.fit.test <- df[df$year > 2022,]


####################################################################
# ACF plots
####################################################################

acf(data.fit.train$xgb_ba,lag.max=120)
acf(data.fit.train$xgb_cnt,lag.max=120)

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





####################################################################
# Tune hyperparameters of XGBoost using Bayesian optimization
####################################################################

# make_purged_year_folds <- function(data, horizon_months = 1, embargo_months = 0,
#                                    warmup_years = 2, min_train_rows = 1L) {
#   dvec   <- as_date(paste(data[['year']],data[['month']],01, sep='-'))
#   yrs    <- sort(unique(year(dvec)))
#   folds  <- list()
#   
#   # label end = date + horizon months
#   label_end <- dvec %m+% months(horizon_months)
#   
#   for (y in yrs) {
#     # validation block = whole calendar year y
#     v_mask   <- year(dvec) == y
#     if (!any(v_mask)) next
#     v_start  <- ymd(paste0(y, "-01-01"))
#     
#     # respect warm-up: need at least 'warmup_years' fully past years
#     if (y <= min(yrs) + warmup_years - 1) next
#     
#     # purge/embargo: remove from training any rows whose label window touches the val start
#     # i.e., keep rows with label_end < (v_start - embargo)
#     cutoff   <- v_start %m-% months(embargo_months)
#     t_mask   <- (label_end < cutoff)
#     
#     tr_idx   <- which(t_mask)
#     va_idx   <- which(v_mask)
#     
#     if (length(tr_idx) >= min_train_rows && length(va_idx) > 0) {
#       folds[[length(folds) + 1]] <- list(train_idx = tr_idx,
#                                          val_idx   = va_idx,
#                                          val_year  = y)
#     }
#   }
#   return(folds)
# }

make_year_folds <- function(data, horizon_months = 1, 
                                   warmup_years = 2, min_train_rows = 1L) {
  dvec   <- as_date(paste(data[['year']],data[['month']],01, sep='-'))
  yrs    <- sort(unique(year(dvec)))
  folds  <- list()
  
  # label end = date + horizon months
  label_end <- dvec %m+% months(horizon_months)
  
  for (y in yrs) {
    # validation block = whole calendar year y
    v_mask   <- year(dvec) == y
    if (!any(v_mask)) next
    v_start  <- ymd(paste0(y, "-01-01"))
    
    # respect warm-up: need at least 'warmup_years' fully past years
    if (y <= min(yrs) + warmup_years - 1) next
    
    # purge/embargo: remove from training any rows whose label window touches the val start
    # i.e., keep rows with label_end < (v_start - embargo)
    cutoff   <- v_start
    t_mask   <- (label_end <= cutoff)
    
    tr_idx   <- which(t_mask)
    va_idx   <- which(v_mask)
    
    if (length(tr_idx) >= min_train_rows && length(va_idx) > 0) {
      folds[[length(folds) + 1]] <- list(train_idx = tr_idx,
                                         val_idx   = va_idx,
                                         val_year  = y)
    }
  }
  return(folds)
}



val_metric <- function(y_true,
                       y_pred, 
                       metric = c("auc", "rmse", "poisson-nloglik", "tweedie-nloglik@1.5"),
                       positive_label = 1,
                       eps = 1e-12) {
  na.ind = is.na(y_pred) # due to warming up
  y_true = y_true[!na.ind]
  y_pred = y_pred[!na.ind] 
  metric <- if (length(metric) == 1) metric else match.arg(metric)
  
  if (metric == "auc") {
    # y_true should be 0/1
    return(as.numeric(pROC::auc(y_true == positive_label, y_pred)))
  }
  
  if (metric == "rmse") {
    return(sqrt(mean((y_true - y_pred)^2)))
  }
  
  if (metric == "poisson-nloglik") {
    # Poisson negative log-likelihood (up to an additive constant):
    # mean( mu - y * log(mu) ), where mu > 0
    mu <- pmax(y_pred, eps)
    if (any(y_true < 0, na.rm = TRUE)) {
      stop("poisson-nloglik requires non-negative targets.")
    }
    return(mean(mu - y_true * log(mu)))
  }
  
  if (grepl("^tweedie-nloglik@", metric)) {
    # Tweedie negative log-likelihood (power rho in (1,2) usually):
    # mean( mu^(2-rho)/(2-rho) - y * mu^(1-rho)/(1-rho) ), ignoring constants in y and phi
    rho_str <- sub("^tweedie-nloglik@", "", metric)
    rho <- as.numeric(rho_str)
    if (!is.finite(rho)) stop("Could not parse Tweedie power from metric string.")
    if (rho <= 1 || rho >= 2) warning("rho is typically in (1,2) for compound Poisson-gamma Tweedie.")
    
    mu <- pmax(y_pred, eps)
    if (any(y_true < 0, na.rm = TRUE)) {
      stop("tweedie-nloglik requires non-negative targets.")
    }
    term <- (mu^(2 - rho)) / (2 - rho) - (y_true * mu^(1 - rho)) / (1 - rho)
    # Guard against numerical issues
    term[!is.finite(term)] <- NA_real_
    return(mean(term, na.rm = TRUE))
  }
  
  stop(sprintf("Unknown metric: %s", metric))
}



XGBoostAR <- function(data,
                      covar.names, 
                      target.name, 
                      eta,
                      max_depth, 
                      objective, 
                      evaluation_metric,
                      horizon_months = 1,
                      embargo_months = 0,
                      warmup_years = 2,
                      init_points = 10,
                      n_iter = 20,
                      early_stopping_rounds = 50,
                      max_rounds = 1000){
  
  set.seed(42)
  # Prepare the dataset: separate predictors and response
  
  folds <- make_year_folds(
    data             = data,
    horizon_months   = horizon_months,
    warmup_years     = warmup_years
  )
  
  if (length(folds) == 0) stop("No valid purged folds: increase warmup_years or check date column.")
  
  label <- data[,target.name]
  features <- data[, covar.names]
  maximize_flag <- identical(evaluation_metric, "auc")
  # Define the objective function for Bayesian optimization
  # This function performs 6-fold cross-validation with early stopping and returns the negative RMSE.
  # (BayesianOptimization maximizes the objective, so we use -RMSE.)
  xgb_cv_bayes <- function(max_depth, subsample, colsample_bytree, min_child_weight, gamma) {
    
    params <- list(
      booster = "gbtree",
      eta = eta,  # Learning rate; adjust if necessary.
      max_depth = max_depth,
      subsample = subsample,
      colsample_bytree = colsample_bytree,
      objective = objective,  # For a continuous regression target.
      eval_metric = evaluation_metric,
      # scale_pos_weight = if (evaluation_metric == 'auc') 8.34 else 1,
      min_child_weight = min_child_weight,
      gamma = gamma
    )
    
    oof <- rep(NA_real_, nrow(data))
    
    for (i in seq_along(folds)) {
      tr <- folds[[i]]$train_idx
      va <- folds[[i]]$val_idx
      
      dtr <- xgb.DMatrix(data = as.matrix(features[tr, ]), label = label[tr])
      dva <- xgb.DMatrix(data = as.matrix(features[va, ]), label = label[va])
      
      mdl <- xgb.train(params = params,
                       data   = dtr,
                       nrounds = max_rounds,
                       watchlist = list(train = dtr, eval = dva),
                       early_stopping_rounds = early_stopping_rounds, verbose = 0)
      
      oof[va] <- predict(mdl, dva, ntreelimit = mdl$best_ntreelimit)
    }
    
    best_metric <- val_metric(label, oof, evaluation_metric)
    
    
    list(Score = if (maximize_flag) best_metric else -best_metric, Pred = 0)
  }
  
  
  # Run Bayesian optimization to tune max_depth, subsample, and colsample_bytree
  opt_results <- BayesianOptimization(
    FUN = xgb_cv_bayes,
    bounds = list(
      max_depth = max_depth,           # Tune over reasonable integer values (converted inside the function).
      subsample = c(0.6, 1),        # Proportion of observations to sample.
      colsample_bytree = c(0.6, 1),   # Proportion of features to sample.
      min_child_weight = c(1,20),
      gamma = c(0,5)
    ),
    init_points = init_points,  # Initial random explorations.
    n_iter = n_iter,      # Number of iterations for Bayesian optimization.
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
    objective = objective,
    eval_metric = evaluation_metric,
    # scale_pos_weight = if (evaluation_metric == 'auc') 8.34 else 1,
    min_child_weight = best_params['min_child_weight'],
    gamma = best_params['gamma']
  )

  return( list(
                best_params   = final_params,
                folds_info    = folds,         # returned for transparency/debugging
                horizon_months = horizon_months,
                embargo_months = embargo_months,
                warmup_years   = warmup_years
  ))
}

xgb_super_learner_predict  <- function(data_train,
                                       data_test,
                           tuning_res,
                           target.name,
                           covar.names,
                           max_rounds=1000,
                           early_stopping_rounds=50) {
  set.seed(42)
  horizon_months  <- tuning_res$horizon_months
  embargo_months  <- tuning_res$embargo_months
  warmup_years    <- tuning_res$warmup_years
  
  folds <- make_year_folds(
    data            = data_train,
    horizon_months  = horizon_months,
    warmup_years    = warmup_years
  )

  
  params <- tuning_res$best_params
  
  if (length(folds) == 0) stop("No valid purged folds on training set.")
  
  y_tr   <- data_train[[target.name]]
  X_tr   <- as.matrix(data_train[, covar.names, drop = FALSE])
  X_te   <- as.matrix(data_test[,  covar.names, drop = FALSE])
  
  oof <- rep(NA_real_, nrow(data_train))
  best_iter <- integer(length(folds))
  
  for (i in seq_along(folds)) {
    tr <- folds[[i]]$train_idx
    va <- folds[[i]]$val_idx
    
    dtr <- xgb.DMatrix(data = as.matrix(X_tr[tr, ]), label = y_tr[tr])
    dva <- xgb.DMatrix(data = as.matrix(X_tr[va, ]), label = y_tr[va])
    
    mdl <- xgb.train(params = params,
                     data   = dtr,
                     nrounds = max_rounds,
                     watchlist = list(train = dtr, eval = dva),
                     early_stopping_rounds = early_stopping_rounds, verbose = 0)
    
    oof[va] <- predict(mdl, dva, nrounds = mdl$best_iteration)
    best_iter[i] <- mdl$best_iteration
  }
  
  dtr_full <- xgb.DMatrix(X_tr, label = y_tr)
  mdl_full <- xgb.train(
    params  = params,
    data    = dtr_full,
    nrounds = as.integer(median(best_iter, na.rm = TRUE)),
    verbose = 0
  )
  
  dte <- xgb.DMatrix(X_te)
  test_pred <- predict(mdl_full, newdata = dte)
  
  return(list(
    oof_pred   = oof,        # for meta-learner training
    test_pred  = test_pred,  # final forecast on test set
    final_model = mdl_full
  ))
}



######################Fit the model############################

res_ba_h1 = XGBoostAR(as.data.frame(data.fit.train),
                      ba.covar.names.h1 , 
                      'xgb_ba_h1',
                      eta=0.1, 
                      max_depth=c(2L,5L),
                      objective='reg:tweedie', 
                      evaluation_metric='tweedie-nloglik@1.5',
                      init_points = 10,
                      n_iter = 30)

save(res_ba_h1, file=file.path(dir.out,'AutoRegressive_XGBoost_Hyperparameters_BA.RData'))


res_cnt_h1 =  XGBoostAR(as.data.frame(data.fit.train), 
                        cnt.covar.names.h1 , 
                        'xgb_cnt_h1',
                        eta=0.1, 
                        max_depth=c(2L,5L), 
                        objective="count:poisson",
                        evaluation_metric="poisson-nloglik",
                        init_points = 10,
                        n_iter = 30)

save(res_cnt_h1, file=file.path(dir.out,'AutoRegressive_XGBoost_Hyperparameters_CNT.RData'))



sl_preds_ba_h1 <- xgb_super_learner_predict(
  data_train  = data.fit.train,     # past-only
  data_test   = data.fit.test,      # future period
  tuning_res  = res_ba_h1,
  covar.names = ba.covar.names.h1,
  target.name = "xgb_ba_h1"
)

sl_preds_cnt_h1 <- xgb_super_learner_predict(
  data_train  = data.fit.train,     # past-only
  data_test   = data.fit.test,      # future period
  tuning_res  = res_cnt_h1,
  covar.names = cnt.covar.names.h1,
  target.name = "xgb_cnt_h1"
)


data.fit.train['pred_ba_h1'] = sl_preds_ba_h1$oof_pred
data.fit.test['pred_ba_h1'] = sl_preds_ba_h1$test_pred

data.fit.train['pred_cnt_h1'] = sl_preds_cnt_h1$oof_pred
data.fit.test['pred_cnt_h1'] = sl_preds_cnt_h1$test_pred




data <- data.fit.train
council = 'Alcanena'
plot(1:(dim(data[data$NAME_2==council,])[1]), data[data$NAME_2==council, ]$xgb_ba_h1,type='l')
points(data[data$NAME_2==council, ]$pred_ba_h1,col='red')

plot(1:(dim(data[data$NAME_2==council,])[1]), data[data$NAME_2==council, ]$xgb_cnt_h1,type='l')
points(data[data$NAME_2==council, ]$pred_cnt_h1,col='red')






####################################################################
# Prepare the dataset for INLA modelling
####################################################################



data.fit <- rbind(data.fit.train, data.fit.test)

plot(1:(dim(data.fit[data.fit$NAME_2==council,])[1]), data.fit[data.fit$NAME_2==council, ]$xgb_cnt_h1,type='l')
lines(data.fit[data.fit$NAME_2==council, ]$pred_cnt_h1,col='red')

plot(1:(dim(data.fit[data.fit$NAME_2==council,])[1]), data.fit[data.fit$NAME_2==council, ]$xgb_ba_h1,type='l')
lines(data.fit[data.fit$NAME_2==council, ]$pred_ba_h1,col='red')


data.fit <- data.fit %>%
  arrange(NAME_2, time.idx) %>%   # Ensure data is ordered by location and time
  group_by(NAME_2) %>%               # Apply the lag per location
    mutate(!!paste0("score_ba_h", 1) := lag(pred_ba_h1, n = 1),
           !!paste0("score_cnt_h",1) := lag(pred_cnt_h1, n = 1)) %>% ungroup() |> as.data.frame()



data.fit[data.fit$NAME_2==council,c('time.idx','NAME_2','xgb_ba','xgb_ba_h1','pred_ba_h1','score_ba_h1','year','month')]
data.fit[data.fit$NAME_2=='Águeda',c('time.idx','NAME_2','xgb_cnt','xgb_cnt_h1','pred_cnt_h1','score_cnt_h1','year','month')]

first_month = min(data.fit$time.idx)
data.fit <- data.fit[!is.na(data.fit$score_ba_h1),]


mean((data.fit[data.fit$year<=2022,'score_cnt_h1']-data.fit[data.fit$year<=2022,'xgb_cnt'])^2)
mean((data.fit[data.fit$year>2022,'score_cnt_h1']-data.fit[data.fit$year>2022,'xgb_cnt'])^2)

mean((data.fit[data.fit$year<=2022,'score_ba_h1']-data.fit[data.fit$year<=2022,'xgb_ba'])^2)
mean((data.fit[data.fit$year>2022,'score_ba_h1']-data.fit[data.fit$year>2022,'xgb_ba'])^2)


save(data.fit, 
     file=file.path(dir.out,'AutoRegressive_XGBoost_Predictions.RData'))