# Base R Shiny image
FROM rocker/shiny

# Make a directory in the container
RUN mkdir /home/shiny-app

# Install R dependencies
RUN R -e "install.packages(c('plotly', 'dplyr', 'shinyWidgets'))"

# Copy the Shiny app code
COPY Rshiny_app_data-driven.R /home/shiny-app/app.R
ADD images_rshiny /home/shiny-app/images_rshiny

# Expose the application port
EXPOSE 8180

# Run the R Shiny app
CMD Rscript /home/shiny-app/app.R