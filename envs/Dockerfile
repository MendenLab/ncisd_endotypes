# Use Miniconda as base
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /app

# Copy your environment YAML file
COPY python38_conda_env.yml .

# Create conda environment
RUN conda install -n base -c conda-forge mamba
RUN mamba env create -f python38_conda_env.yml
RUN conda clean -afy

# Port für Jupyter freigeben
EXPOSE 8888 8787


# After conda env is created and activated...
RUN conda run -n r_py_env R -e "\
  if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes'); \
  remotes::install_github('jmw86069/multienrichjam', dependencies=TRUE); \
  remotes::install_github('davidsjoberg/ggsankey', dependencies=TRUE)"


# Copy your code into the container
COPY . /app

# Set environment path (so conda activate works)
ENV PATH=/opt/conda/envs/analysis_env/bin:$PATH

# Start both JupyterLab and RStudio (supervisord alternative)
CMD ["conda", "run", "-n", "r_py_env", "bash", "-c", "\
    jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root & \
    rstudio-server start && \
    tail -f /dev/null"]
