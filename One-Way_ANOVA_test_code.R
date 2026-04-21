#===============================================================================
## One-way ANOVA test to compare Boyce Index scores based on species status, ##
##      widespread species and endemic species to Madagascar.                ##
#===============================================================================
# The subset of 14 species from the Archbold-Vernay Expedition (AVE) is split into
# two groups: Wide-spread species and Endemic species using the Boyce-Index figures to
# carry out the comparison. 
#===================================
# Part 1: One-way ANOVA test using Random Sampling Boyce Index values
# 1. Load file holding Boyce Index Vales data
sdm_data_random <- file.choose()

# 2. Run the One-way ANOVA test using the 'aov()'function 
# The formula: Dependent_Variable ~ Independent_Variable, data = file name
random_sampling_model <- aov(boyce_index ~ species_type, data = sdm_data_random)

# 3. View the results using 'summary()' function
summary(random_sampling_model)

# 4. Check for Homogeneity of Variance using Levene's Test, this Requires 'car' package
install.packages("car")
library(car)
leveneTest(boyce_index ~ species_type, data = sdm_data_random)

# 4. Visualization (Boxplot)
boxplot(boyce_index ~ species_type, data = sdm_data_random,
        main = "Boyce Index by Species Type (Random Sampling)",
        xlab = "Species Type",
        ylab = "Boyce Index Score",
        col = c("lightblue", "lightgreen"))

#===================================
# Part 2: One-way ANOVA test using Sampling Bias Boyce Index values
# 1. Load file holding Boyce Index Vales data
sdm_data_bias <- file.choose()

# 2. Run the One-way ANOVA test using the 'aov()'function 
# The formula: Dependent_Variable ~ Independent_Variable, data = file name
bias_sampling_model <- aov(boyce_index ~ species_type, data = sdm_data_bias)

# 3. View the results using 'summary()' function
summary(bias_sampling_model)

# 4. Check for Homogeneity of Variance using Levene's Test, this Requires 'car' package
install.packages("car")
library(car)
leveneTest(boyce_index ~ species_type, data = sdm_data_bias)

# 4. Visualization (Boxplot)
boxplot(boyce_index ~ species_type, data = sdm_data_bias,
        main = "Boyce Index by Species Type (Sampling Bias)",
        xlab = "Species Type",
        ylab = "Boyce Index Score",
        col = c("lightblue", "lightgreen"))
