library(readODS)

# Read data
df = read_ods("applications/oecd/data/OECDData.ods")

# Save countries for later
countries = df$...3

# Remove irrelevant row and columns
df$Country = NULL
df$...2 = NULL
df$...3 = NULL
df = df[-1, ]

# Transform df into numeric df
for (colname in colnames(df)) {
  df[[colname]] = as.numeric(df[[colname]])
}

# Add countries back to the dataframe as index
df = data.frame(df)
rownames(df) = countries[-1]

# Change column names
df$Housing = df$Housing...12
df$Community = df$Housing...13
df$Life.satisfaction = df$Housing...14
df$Housing...12 = NULL
df$Housing...13 = NULL
df$Housing...14 = NULL

# Remove NAs
df = na.omit(df)

# Clusters found in Cavicchia et al. (2022)
c1 = c("AUS", "AUT", "BEL", "CAN", "DNK", "FIN", "FRA", "DEU", "ISL", "IRL",
       "ITA", "JPN", "LUX", "NLD", "NZL", "NOR", "ESP", "SWE", "CHE", "GBR",
       "USA")
c2 = c("CHL", "CZE", "EST", "GRC", "HUN", "ISR", "KOR", "LVA", "MEX", "POL",
       "PRT", "SVK", "SVN", "TUR")

df1 = df[c1, ]
df2 = df[c2, ]

# Save df
save(df1, df2, file = "applications/oecd/data/OECDData.RData")
