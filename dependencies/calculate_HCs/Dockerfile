FROM rocker/r-ver:4.5.2

RUN Rscript -e 'install.packages( \
        c( \
            "data.table", \
            "Matrix", \
            "sparsesvd", \
            "optparse", \
            "R.utils" \
        ), \
        repos="https://www.stats.bris.ac.uk/R/")'

COPY calculate_haplotype_components.R /usr/local/bin/calculate_haplotype_components

RUN chmod +xr /usr/local/bin/calculate_haplotype_components

CMD ["calculate_haplotype_components", "-h"]
