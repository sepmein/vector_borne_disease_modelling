# ===============================================================================
# Package Installation Script for SEIR-SEI MCMC Analysis
# ===============================================================================

cat("Installing required packages for SEIR-SEI MCMC analysis...\n")

# List of required packages
required_packages <- c(
    "rstan", # Stan interface for R
    "tidybayes", # Bayesian analysis tools
    "tidyverse", # Data manipulation suite
    "ggplot2", # Plotting
    "posterior", # Posterior analysis tools
    "bayesplot", # Bayesian plotting
    "lubridate", # Date handling
    "dplyr", # Data manipulation
    "tidyr", # Data tidying
    "stringr" # String manipulation
)

# Function to install packages safely
install_if_missing <- function(package_name) {
    if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
        cat(sprintf("Installing %s...\n", package_name))
        install.packages(package_name, repos = "https://cran.r-project.org")

        # Verify installation
        if (require(package_name, character.only = TRUE, quietly = TRUE)) {
            cat(sprintf("✅ %s installed successfully\n", package_name))
        } else {
            cat(sprintf("❌ Failed to install %s\n", package_name))
        }
    } else {
        cat(sprintf("✅ %s already installed\n", package_name))
    }
}

# Install Stan first (special handling)
cat("\n=== Installing RStan (may take several minutes) ===\n")
if (!require("rstan", quietly = TRUE)) {
    # Configure for faster compilation
    options(repos = c(CRAN = "https://cran.r-project.org"))

    # Install dependencies first
    install.packages(c("V8", "StanHeaders"), dependencies = TRUE)

    # Install rstan
    install.packages("rstan", dependencies = TRUE)

    # Configure rstan
    library(rstan)
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
} else {
    cat("✅ rstan already installed\n")
}

# Install other packages
cat("\n=== Installing other required packages ===\n")
for (pkg in required_packages[-1]) { # Skip rstan as it's already handled
    install_if_missing(pkg)
}

# Final verification
cat("\n=== Final Package Verification ===\n")
missing_packages <- c()
for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        missing_packages <- c(missing_packages, pkg)
    }
}

if (length(missing_packages) == 0) {
    cat("\n🎉 All packages installed successfully!\n")
    cat("You can now run the MCMC analysis.\n")
} else {
    cat("\n⚠️  The following packages failed to install:\n")
    for (pkg in missing_packages) {
        cat(sprintf("   - %s\n", pkg))
    }
    cat("\nPlease install these manually or check your R installation.\n")
}

cat("\n===============================================================================\n")
cat("Installation complete! You can now run: results <- main_analysis()\n")
cat("===============================================================================\n")
