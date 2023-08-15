##############################################
#        Analysis master dissertation        #
##############################################


# Set up ####
# clear console
rm(list=ls()); cat("\014")


# load packages
library(stats)
library(ggplot2)
library(dplyr)
library(yarrr)
library(ggplot2)
library(lme4)
library(car)
library(rstatix)


#####################################
### Load and preprocess the data ####
#####################################
# NB. S001 and S002 no valid data, so are skipped here

setwd("C:/Users/louis/OneDrive - UGent/Documents/Psychologie/THESIS/Datacollectie/Data/Computertaakdata")
S003  <- read.csv(file = 'SST_data_S003.csv', sep=',', dec= '.')
S004  <- read.csv(file = 'SST_data_S004.csv', sep=';', dec= '.')
S005  <- read.csv(file = 'SST_data_S005.csv', sep=',', dec= '.', skip = 1)
S006  <- read.csv(file = 'SST_data_S006.csv', sep=';', dec= '.')
S007  <- read.csv(file = 'SST_data_S007.csv', sep=',', dec= '.')
S008  <- read.csv(file = 'SST_data_S008.csv', sep=',', dec= '.')
S009  <- read.csv(file = 'SST_data_S009.csv', sep=',', dec= '.')
S010  <- read.csv(file = 'SST_data_S010.csv', sep=',', dec= '.')
S011  <- read.csv(file = 'SST_data_S011.csv', sep=',', dec= '.')
S012  <- read.csv(file = 'SST_data_S012.csv', sep=',', dec= '.')
S013  <- read.csv(file = 'SST_data_S013.csv', sep=',', dec= '.')
S014  <- read.csv(file = 'SST_data_S014.csv', sep=',', dec= '.')
S015  <- read.csv(file = 'SST_data_S015.csv', sep=',', dec= '.')
S016  <- read.csv(file = 'SST_data_S016.csv', sep=',', dec= '.')
S017  <- read.csv(file = 'SST_data_S017.csv', sep=';', dec= '.')
S018  <- read.csv(file = 'SST_data_S018.csv', sep=';', dec= '.')
S019  <- read.csv(file = 'SST_data_S019.csv', sep=',', dec= '.')
S020  <- read.csv(file = 'SST_data_S020.csv', sep=',', dec= '.')

# error with file H002 so skip this
# H007: undefined responses, so also removed
H001 <- read.csv(file = 'SST_data_H001.csv', sep=';', dec= '.')
H003 <- read.csv(file = 'SST_data_H003.csv', sep=';', dec= '.')
H004 <- read.csv(file = 'SST_data_H004.csv', sep=';', dec= '.')
H005 <- read.csv(file = 'SST_data_H005.csv', sep=';', dec= '.')
H006 <- read.csv(file = 'SST_data_H006.csv', sep=';', dec= '.', nrows = 1288)
H008 <- read.csv(file = 'SST_data_H008.csv', sep=';', dec= '.')
H009 <- read.csv(file = 'SST_data_H009.csv', sep=';', dec= '.')
H010 <- read.csv(file = 'SST_data_H010.csv', sep=';', dec= '.')
H011 <- read.csv(file = 'SST_data_H011.csv', sep=';', dec= '.')
H012 <- read.csv(file = 'SST_data_H012.csv', sep=';', dec= '.')
H013 <- read.csv(file = 'SST_data_H013.csv', sep=';', dec= '.')
H014 <- read.csv(file = 'SST_data_H014.csv', sep=';', dec= '.')
H015 <- read.csv(file = 'SST_data_H015.csv', sep=';', dec= '.', nrows = 1288)
H016 <- read.csv(file = 'SST_data_H016.csv', sep=';', dec= '.')
H017 <- read.csv(file = 'SST_data_H017.csv', sep=';', dec= '.')
H018 <- read.csv(file = 'SST_data_H018.csv', sep=';', dec= '.', nrows = 1288)
H019 <- read.csv(file = 'SST_data_H019.csv', sep=';', dec= '.')
H020 <- read.csv(file = 'SST_data_H020.csv', sep=';', dec= '.')
H021 <- read.csv(file = 'SST_data_H021.csv', sep=';', dec= '.')
H022 <- read.csv(file = 'SST_data_H022.csv', sep=';', dec= '.')
H023 <- read.csv(file = 'SST_data_H023.csv', sep=';', dec= '.')
H024 <- read.csv(file = 'SST_data_H024.csv', sep=';', dec= '.')
H025 <- read.csv(file = 'SST_data_H025.csv', sep=';', dec= '.')
H026 <- read.csv(file = 'SST_data_H026.csv', sep=';', dec= '.')
H027 <- read.csv(file = 'SST_data_H027.csv', sep=';', dec= '.')
H028 <- read.csv(file = 'SST_data_H028.csv', sep=';', dec= '.')
H029 <- read.csv(file = 'SST_data_H029.csv', sep=';', dec= '.')
H030 <- read.csv(file = 'SST_data_H030.csv', sep=';', dec= '.')
H031 <- read.csv(file = 'SST_data_H031.csv', sep=';', dec= '.', nrows = 1288)
H032 <- read.csv(file = 'SST_data_H032.csv', sep=';', dec= '.')
H033 <- read.csv(file = 'SST_data_H033.csv', sep=';', dec= '.')
H034 <- read.csv(file = 'SST_data_H034.csv', sep=';', dec= '.')
H035 <- read.csv(file = 'SST_data_H035.csv', sep=';', dec= '.')
H036 <- read.csv(file = 'SST_data_H036.csv', sep=';', dec= '.')
H037 <- read.csv(file = 'SST_data_H037.csv', sep=';', dec= '.')
H038 <- read.csv(file = 'SST_data_H038.csv', sep=';', dec= '.')
H039 <- read.csv(file = 'SST_data_H039.csv', sep=';', dec= '.')
H040 <- read.csv(file = 'SST_data_H040.csv', sep=';', dec= '.')
H041 <- read.csv(file = 'SST_data_H041.csv', sep=';', dec= '.')
H042 <- read.csv(file = 'SST_data_H042.csv', sep=';', dec= '.')
H043 <- read.csv(file = 'SST_data_H043.csv', sep=';', dec= '.', skip = 1)
H044 <- read.csv(file = 'SST_data_H044.csv', sep=';', dec= '.')
H045 <- read.csv(file = 'SST_data_H045.csv', sep=';', dec= '.')

Scores <- read.csv(file='Scores.csv', sep = ';')

# Create full dataset
Fulldata <- rbind(S003,S004,S005,S006,S007,S008,S009,S010,
                  S011,S012,S013,S014,S015,S016,S017,S018,S019,S020,
                  H001, H003, H004, H005, H006, H008, H009, H010, 
                  H011, H012, H013, H014, H015, H016, H017, H018, H019, H020, 
                  H021, H022, H023, H024, H025, H026, H027, H028, H029, H030, 
                  H031, H032, H033, H034, H035, H036, H037, H038, H039, H040,
                  H041, H042, H043, H044, H045)
Fulldata <- merge(Fulldata, Scores, by = "participantID")

# Create groups 'Healthy or Tics'
Fulldata <- Fulldata %>%
  mutate(group = ifelse(substr(participantID, 1, 1) == "S", "tics","healthy"))
Fulldata <- Fulldata[!duplicated(Fulldata), ]


# Preprocess the data (reclass, rename, and redefine some variables)
Fulldata$age          <- as.numeric(Fulldata$age.x)
Fulldata$gender       <- as.factor(Fulldata$gender.x)
Fulldata$handedness   <- as.factor(Fulldata$handedness)
Fulldata$handedness   <- relevel(factor(Fulldata$handedness),ref = "Rechts") ## reference = Right 
Fulldata$group        <- as.factor(Fulldata$group)
Fulldata$blocknumber  <- as.factor(Fulldata$block_i)
Fulldata$rt           <- as.numeric(Fulldata$rt)
Fulldata$SSD          <- as.numeric(Fulldata$SSD)
Fulldata$score        <- as.numeric(Fulldata$score)
Fulldata$correct      <- Fulldata$correct == 'true' # change to boolean
Fulldata$ticscore     <- as.numeric(Fulldata$score)

# Set ticscore for healthy participants to 0
Fulldata$ticscore <- ifelse(is.na(Fulldata$score), 0, Fulldata$score)

# Finally, select only the relevant columns (first) and rows (second)
# to obtain the finale data set
data <- subset(Fulldata, select = c(participantID, age, gender, 
                                    handedness, ticscore, group, block_i, 
                                    trial_i, stim, signal, SSD, response, 
                                    rt, correct))
data <- subset(data, block_i!= 0) # remove practice block
# clear up the console a bit; keep only the final data set
rm(Fulldata, list = grep("H0|S0|Scores",ls(), value=T) )


# R: this is for me only! rename some variables so it complies with the code below
#colnames(data) <- gsub("YTSS", "ticscore",colnames(data))
#data$group <- ifelse(data$group == "tic", "tics", "healthy")


################
# demographics #
################

datademo <- subset(data, select = c(participantID, age, gender, handedness, ticscore, group))
datademo <- datademo[!duplicated(datademo), ]
# add rtmean per participant in the demographic dataframe
datademo <- merge(datademo, aggregate(rt ~ participantID, data = data, mean), by = "participantID", suffixes = c("", "_mean"))
datademo$rt <- as.numeric(datademo$rt)

# check whether there is a significant difference for rt between healthy and tic pax
rthealthy <- subset(datademo, group == "healthy")
rttics    <- subset(datademo, group == "tics")

t.test(rthealthy$rt, rttics$rt)
## p = 0.0578 --> no significant difference, healthy and tic participants don't have a difference in rt BUT a trend is visible


# Calculate proportions
# general gender
table(datademo$gender, datademo$group) 

# handedness
table(datademo$handedness, datademo$group)

# Plot a distribution of the Tic Severity Scores in participants with tics
hist(datademo2$ticscore)

#Plot a distribution of the reaction time over Group
pirateplot(rt ~ group, data = datademo, main = 'Distribution of RT over Group', ylab = 'Reaction time (in ms)', xlab = 'Group', plot = TRUE)



######################################
#### Analysis 1: SSRT differences ####
######################################

###############
# Preparation #
############### 
# Subset go trials and remove missing reaction time (rt) values
# The signal column does not use 'go' and 'stop' but 'no' and 'yes' respectively;

#create an empty list to store SSRT values for each participant
ssrt_list <-list()

# Iterate over unique participant IDs
for (participantID in unique(data$participantID)) {
  # subset data for current participant
  data1 <- data[data$participantID == participantID, ]
  
  # subset go trials for the current participant
  go_trials <- data1[data1$signal == "no" & !is.na(data1$rt), ]
  
  #subset stop trials for the current participant
  # R: leave incorrect stop trials (i.e., those in which RT = NA) in the data,
  # this is important for the SSRT calculation
  stop_trials <- data1[data1$signal == "yes" & !is.na(data1$SSD), ]
  
  #STEP 1: Sort go trials in ascending order of reaction times
  sorted_go_trials <- go_trials[order(go_trials$rt), ]
  
  # STEP 2: Calculate the middle 50% of go reaction times
  # n is an index value based on the probability of being incorrect on a stop trial
  # so if we have 100 go trials, and were inaccurate on 50% of the stop trials, n = 
  # 100 * 0.50 = 50. The "nth Go RT" is then the RT on trials n, thus on the 50th one
  # n <- nrow(sorted_go_trials)
  probab_respond_to_signal <- 1-mean(stop_trials$correct)
  n <- round(nrow(sorted_go_trials) * probab_respond_to_signal)
  nth_gotrial <- sorted_go_trials[n, ]
  
  # Step 3: Calculate the SSRT
  # middle_rt <- median(middle_go_trials$rt)
  nth_rt <- nth_gotrial$rt
  mean_ssd <- round(mean(stop_trials$SSD))
  ssrt <- nth_rt - mean_ssd # R: there you go!
  
  # Store the SSRT value for the current participant
  ssrt_list[[as.character(participantID)]] <- ssrt
  
  # Print the calculated SSRT
  cat("Participant:", participantID, "SSRT:", ssrt, "ms\n") 
  
}


# Recommendation 7: Checking violation of independence assumption
# Check the independence assumption and compare mean RTs
# Note. stop trials are trials in which signal == "yes"!

stop_trials <- data[data$signal == "yes" & !is.na(data$rt), ]
if (nrow(stop_trials) > 0) {
  mean_rt_stop <- aggregate(rt ~ participantID, data = stop_trials, mean)
  mean_rt_go <- aggregate(rt ~ participantID, data = data[data$signal == "yes" & data$response != "" & !is.na(data$rt), ], mean)
  violation <- subset(merge(mean_rt_stop, mean_rt_go, by = "participantID"), rt.x > rt.y)
} else {
  violation <- data.frame()  # No violation as there are no valid stop trials
}


# Recommendation 8: Estimating SSRT using the integration method:++
calculate_SSRT <- function(data) {
  go_trials <- subset(data, signal == "yes" & !is.na(rt))
  stop_trials <- subset(data, signal == "no" & !is.na(rt) & !is.na(SSD))
  
  # Calculate the mean and standard deviation of go trial reaction times
  go_rt_mean <- mean(go_trials$rt)
  go_rt_sd <- sd(go_trials$rt)
  
  # Calculate the mean stop signal delay (SSD)
  ssd_mean <- mean(stop_trials$SSD)
  
  # Calculate the estimated SSRT using the integration method
  ssrt <- go_rt_mean + (go_rt_sd * qnorm((sum(go_trials$rt > ssd_mean) / length(go_trials$rt))))
  
  return(ssrt)
}

# Call the calculate_SSRT function on your data
ssrt_estimate <- calculate_SSRT(data)
ssrt_estimate

# Recommendation 9: Refraining from estimating SSRT under certain conditions
# Filter out SSRT estimates based on recommendation 9
filtered_SSRT_estimates <- ssrt_list[!(probab_respond_to_signal < 0.25 | probab_respond_to_signal > 0.75)]

# Print the filtered SSRT estimates that are nonzero
filtered_SSRT_estimates <- filtered_SSRT_estimates[sapply(filtered_SSRT_estimates,length)>0]

# Create a new dataframe dataFull by adding a column with filtered SSRT estimates to use in the further analyses
datademo$filtered_SSRT <- as.numeric(filtered_SSRT_estimates)

length(datademo$filtered_SSRT)

# Subset the SSRT values for the tics and healthy groups
tics_ssrt <- datademo$filtered_SSRT[datademo$group == "tics"]
healthy_ssrt <- datademo$filtered_SSRT[datademo$group == "healthy"]


##############
## Analysis ##
##############
# H0: SSRT healthy = SSRT tic participants
# H1: SSRT heatlhy != SSRT tic participants

# Perform independent t-test to test whether there is a difference in ssrt for healthy and tic pax
result <- t.test(tics_ssrt, healthy_ssrt)

# Print the results
cat("Independent t-test results:\n")
cat("Tics Group Mean SSRT:", mean(tics_ssrt), "\n")
cat("Healthy Group Mean SSRT:", mean(healthy_ssrt), "\n")
cat("t-value:", result$statistic, "\n")
cat("p-value:", result$p.value, "\n")
cat("Conclusion:\n")
if (result$p.value < 0.05) {
  cat("There is a significant difference in SSRT between tics and healthy participants.\n")
} else {
  cat("There is no significant difference in SSRT between tics and healthy participants.\n")
}


# H0: SSRT low tic severity = SSRT high tic severity
# H1: SSRT low tic severity != SSRT high tic severity

# linear model to test whether people with higher tic score have a different SSRT from people with a lower tic score.

# create a subdataset for participant with tics only (see above)
datademo2 <- subset(datademo, group == "tics")

myLM1 <- lm(formula = ticscore ~ filtered_SSRT, data = datademo2)
summary(myLM1)



######################################
### Analysis 2: Sequential effects ###
######################################

###############
# Preparation #
############### 

# Do not REMOVE any trials from the data, rather label them as invalid for this
# analysis in particular; removing them will erroneously change the order/sequence 
# of the trials left over

# Determine whether the previous go stimulus was an alteration (alt) or 
# repetition (rep) with respect to the current go stimulus.
dataSE <- data
dataSE$type <- c(NA, as.character(data$stim[1:(length(data$stim)-1)]))
dataSE$type <- as.factor(ifelse(dataSE$stim == dataSE$type, "rep", "alt"))

###############
# SE analysis #
###############
# Create the 'previoustrialtype' column
dataSE$previoustrialtype <- c(NA, ifelse(
  dataSE$signal[-1] == "no" & dataSE$correct[-1], 
  "SGO", ifelse(
    dataSE$signal[-1] == "yes" & dataSE$correct[-1], 
    "SS", "US"))) 

# Add columns needed to determine sequential effects
dataSE <- data
dataSE$prevstim <- as.factor(c(NA, as.character(dataSE$stim[-nrow(dataSE)])))
dataSE$repORalt <- ifelse(dataSE$stim == dataSE$prevstim, "rep","alt")

dataSE$prevsignal <- as.factor(c(NA, as.character(dataSE$signal[-nrow(dataSE)])))
dataSE$preccorrect <- c(NA, dataSE$correct[-nrow(dataSE)])
dataSE$previoustrialtype <- paste0(ifelse(dataSE$preccorrect,'S','U'),
                                   ifelse(dataSE$prevsignal == 'no', 'GO', 'S'))

# Add a column indicating whether the data needs to be included in the analysis
# I.e., allowed current trials: correct go
# Allowed previous trials: correct go (SGo), correct stop (SS), incorrect stop (US)
dataSE <- dataSE[which(dataSE$signal == 'no' & dataSE$correct 
                       & dataSE$previoustrialtype %in% c("SGO","SS","US")),]


# Subset the data for alteration and repetition trials
alteration_data <- data.frame(dataSE[dataSE$repORalt == "alt", ])
repetition_data <- data.frame(dataSE[dataSE$repORalt == "rep", ])

# Calculate means for alteration and repetition trials per participant and previoustrialtype
alteration_mean <- aggregate(rt ~ participantID + previoustrialtype, data = alteration_data, mean)
repetition_mean <- aggregate(rt ~ participantID + previoustrialtype, data = repetition_data, mean)

# Merge alteration and repetition means by participant and previoustrialtype
dataSE2 <- left_join(alteration_mean, repetition_mean, by = c("participantID", "previoustrialtype"), suffix = c("_alteration", "_repetition"))

# Calculate the difference scores between alteration and repetition means per participant and previoustrialtype
dataSE2$difference_score <- dataSE2$rt_alteration - dataSE2$rt_repetition

# Create the final dataframe with ParticipantID, previoustrialtype, and DifferenceScore columns
dataSE3 <- dataSE2[, c("participantID", "previoustrialtype", "difference_score")]

# Merge with demographic data
dataSE3 <- merge(dataSE2, datademo, by = "participantID")

# Print the final dataframe
print(dataSE3) # Note. the previoustrialtype now also considers previous successful stops (SS)



# Plot a distribution of the difference scores
hist(dataSE2$difference_score, main = "Distribution of the difference scores", xlab = "Difference scores")

# Plot a distribution of the difference score by condition and group
pirateplot(difference_score ~ previoustrialtype + group, data = dataSE3,main = "Difference score by previous trial type and group", ylab = "Difference score",xlab = "Previous trial type", plot = TRUE)


# --> difference score (between repetition and alternation trials) per participant per previous trial type (SGO, SS and US)
# H0: difference score tic participants = healthy participants
# H1: difference score tic participants > healthy participants

# With this data, we are trying to predict a numerical value (RT, or RT difference) 
# according to categorical variables (i.e., previoustrialtype and healthy/tic group);
# therefore, we need to do an ANOVA instead of a linear mixed model:
aovdata <- by(dataSE, dataSE[c('repORalt', 'previoustrialtype','participantID')], function(x){
  x$rt <- round(mean(x$rt)); 
  x[1, c('participantID','group','ticscore','previoustrialtype','repORalt', 'rt')]
}); aovdata <- do.call("rbind", aovdata) 
# have a look at the data we are working with:
# head(aovdata, 20); 
View(aovdata)

analysis2 <- anova_test(data = aovdata, 
                        dv = rt, 
                        wid = participantID,
                        between = group,
                        within = c(repORalt, previoustrialtype)
                        #covariate = # add age/gender/handedness here if you like
)
get_anova_table(analysis2) # look at the results table
# result: there is no significant interaction effect between group and any of the
# sequential effects (repORalt, prevtrialtype, nor the three-way interaction)
# conclusion: whichever way you look at sequential effects, the tic versus control
# groups do not perform differently
# note that the above function is nearly the same as the following formula:
# (with the exception of the error definition)
summary( aov(
  rt ~ group*repORalt*previoustrialtype + Error(participantID/1), 
  data= aovdata
) )
# and that we can make the model simpler by omitting the main effects of prevtrialtype and repORalt
summary( aov(
  rt ~ group*previoustrialtype:repORalt + Error(participantID/1), 
  data= aovdata
) ) # the same conclusions hold
# Also note that we find the main effect of seq. effects that we established above as well.



# H0: difference score higher tic severity  = lower tic severity
# H1: difference score higher tic severity  > lower tic severity

# Create dataset with only people with tics
dataSE_tics <- subset(dataSE3, group == "tics")

# Linear model to check whether ticscore has an effect on the difference score when looking at different previous trial types (only for people with tics)
model2 <- lm(formula = ticscore ~ difference_score:previoustrialtype, data = dataSE_tics)
summary(model2) 


