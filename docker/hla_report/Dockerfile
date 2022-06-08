FROM rocker/r-rmd

USER root

RUN apt-get update --allow-releaseinfo-change && apt-get install -y \
    libssl-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    librsvg2-dev \
    libavfilter-dev \
    libwebp-dev \
    cargo \
    libpoppler-cpp-dev \
    libmagick++-dev

RUN Rscript -e 'install.packages(c("rmarkdown", "rjson", "jsonlite", "dplyr", "knitr", "kableExtra"), repos = "http://cran.rstudio.com/");'

# WORKDIR /app
# RUN chmod 777 -R /app

# ADD ./clinical_hla_report_generator.Rmd clinical_hla_report_generator.Rmd
# ADD ./research_hla_report_generator.Rmd research_hla_report_generator.Rmd
# ADD ./run.R run.R

ADD ./hla_report_generator.Rmd hla_report_generator.Rmd

CMD /bin/bash