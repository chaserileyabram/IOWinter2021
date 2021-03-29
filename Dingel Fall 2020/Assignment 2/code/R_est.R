# Chase Abram
# R_est
# Assignment 2, International Macro and Trade
# Fall 2020


# packages <- c("foreign", "fixest", "dplyr", "xtable")
# lapply(packages, library, character.only = TRUE)

# df <- read.dta("col_regfile09.dta")
# df <- df %>% filter(flow > 0)

# Make logflow
lf = log(df$flow)
df <- cbind(df, logflow = lf)

# Make logdistance
ld <- log(df$distw)
df <- cbind(df, logdist = ld)

# Time and run feols with exp-year and imp-year FE
start_time <- Sys.time()
my_reg <- feols(fml = logflow ~ logdist + contig 
                + comlang_off|year^iso_o+year^iso_d, 
                data = df)
end_time <- Sys.time()

# Output time taken
print(end_time - start_time)

# Save table
summary(my_reg, se = "white")
etable(my_reg, se = "white", tex = TRUE, file = "Table_3_R.tex")






