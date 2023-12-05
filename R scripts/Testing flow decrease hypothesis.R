# Load packages
library(dplyr)
library(ggplot2)
library(lubridate)

# Set seed for reproducibility
set.seed(123)

# Define time period
start_date <- as.Date("2009-01-01")
end_date <- as.Date("2020-12-31")
dates <- seq(start_date, end_date, by = "months")

# Simulate monthly flow data with seasonal variations and random noise
monthly_flows <- 100 + 50 * sin(2 * pi * (seq_along(dates) %% 12) / 12) + rnorm(length(dates), mean = 0, sd = 10)

# Create a data frame
flow_data <- data.frame(Date = dates, Flow = monthly_flows)

# Plot the flow data
ggplot(flow_data, aes(x = Date, y = Flow)) +
  geom_line() +
  geom_point() +
  labs(title = "Annual Mean Flow (August to August)",
       x = "Year",
       y = "Flow (cubic meters per second)")

# Calculate annual means
# Randomly select the starting month for each year
start_months <- sample(c("Sep", "Oct", "Nov"), length(unique(year(dates))), replace = TRUE)
# Ensure only one instance of November
start_months[which(start_months == "Nov")[2:length(start_months)]] <- "Oct"
flow_data$Start_Month <- rep(start_months, each = 12)[seq_along(dates)]

# Calculate annual means based on the randomly selected starting month
annual_means <- flow_data %>%
  mutate(Year = lubridate::year(Date),
         Month_Num = lubridate::month(Date),
         Year_Start_Month = match(Start_Month, month.abb)) %>%
  filter(between(Month_Num, Year_Start_Month, Year_Start_Month + 11)) %>%
  group_by(Year) %>%
  summarise(Annual_Mean_Flow = mean(Flow)) 

# Remove the first row (year 2009) from the annual means dataset
annual_means <- annual_means[-1, ]

# Plot the annual means
ggplot(annual_means, aes(x = Year, y = Annual_Mean_Flow)) +
  geom_line() +
  geom_point() +
  labs(title = "Annual Mean Flow (August to August)",
       x = "Year",
       y = "Flow (cubic meters per second)")

# generate shifted data
# Set seed for reproducibility
set.seed(123)

# Define time period
start_date <- as.Date("2009-01-01")
end_date <- as.Date("2020-12-31")
dates <- seq(start_date, end_date, by = "months")

# Simulate monthly flow data with seasonal variations and random noise
monthly_flows_shifted <- 100 + 50 * sin(2 * pi * ((seq_along(dates) + 1) %% 12) / 12) + rnorm(length(dates), mean = 0, sd = 10)

# Create a data frame
flow_data_shifted <- data.frame(Date = dates, Flow = monthly_flows_shifted)

# Plot the shifted flow data
ggplot(flow_data_shifted, aes(x = Date, y = Flow)) +
  geom_line() +
  geom_point() +
  labs(title = "Seasonally Shifted Flow Data (One Month Later)",
       x = "Year",
       y = "Flow (cubic meters per second)")

# Calculate annual means for shifted data
# Randomly select the starting month for each year
start_months_shifted <- sample(c("Sep", "Oct", "Nov"), length(unique(year(dates))), replace = TRUE)
# Ensure only one instance of November
start_months_shifted[which(start_months_shifted == "Nov")[2:length(start_months_shifted)]] <- "Oct"
flow_data_shifted$Start_Month <- rep(start_months_shifted, each = 12)[seq_along(dates)]

# Calculate annual means based on the randomly selected starting month for shifted data
annual_means_shifted <- flow_data_shifted %>%
  mutate(Year = lubridate::year(Date),
         Month_Num = lubridate::month(Date),
         Year_Start_Month = match(Start_Month, month.abb)) %>%
  filter(between(Month_Num, Year_Start_Month, Year_Start_Month + 11)) %>%
  group_by(Year) %>%
  summarise(Annual_Mean_Flow_Shifted = mean(Flow)) 

# Remove the first row (year 2009) from the annual means dataset for shifted data
annual_means_shifted <- annual_means_shifted[-1, ]

# Plot the annual means for shifted data
ggplot(annual_means_shifted, aes(x = Year, y = Annual_Mean_Flow_Shifted)) +
  geom_line() +
  geom_point() +
  labs(title = "Annual Mean Flow (One-Month Seasonally Shifted)",
       x = "Year",
       y = "Flow (cubic meters per second)")

# Linear regression on the unshifted annual means
lm_unshifted <- lm(Annual_Mean_Flow ~ Year, data = annual_means)

# Linear regression on the shifted annual means
lm_shifted <- lm(Annual_Mean_Flow_Shifted ~ Year, data = annual_means_shifted)

# Display the regression coefficients (slopes)
cat("Slope for Unshifted Annual Means:", coef(lm_unshifted)[2], "\n")
cat("Slope for Shifted Annual Means:", coef(lm_shifted)[2], "\n")
