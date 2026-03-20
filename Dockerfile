# CRISPRArchitect — Docker Container
# ====================================
# Packages the full toolkit + Streamlit web app into a
# reproducible container that runs anywhere.
#
# BUILD:
#   docker build -t crisprarchitect .
#
# RUN:
#   docker run -p 8501:8501 crisprarchitect
#
# Then open http://localhost:8501 in your browser.

FROM python:3.11-slim

# Metadata
LABEL maintainer="CRISPRArchitect Team"
LABEL description="Multi-site genome editing strategy optimizer"
LABEL version="0.1.0"

# Avoid Python buffering issues in Docker
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# Set working directory
WORKDIR /app

# Install system dependencies (minimal)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    curl && \
    rm -rf /var/lib/apt/lists/*

# Copy requirements first (cache layer)
COPY requirements.txt .
COPY pyproject.toml .

# Install Python dependencies
RUN pip install --no-cache-dir \
    numpy>=1.24.0 \
    scipy>=1.10.0 \
    matplotlib>=3.7.0 \
    pandas>=2.0.0 \
    streamlit>=1.30.0 \
    plotly>=5.18.0 \
    requests>=2.28.0

# Copy the entire project
COPY . .

# Install the package itself
RUN pip install --no-cache-dir -e .

# Expose Streamlit port
EXPOSE 8501

# Health check
HEALTHCHECK --interval=30s --timeout=10s --retries=3 \
    CMD curl -f http://localhost:8501/_stcore/health || exit 1

# Streamlit configuration for Docker
RUN mkdir -p /app/.streamlit
RUN echo '[server]\n\
port = 8501\n\
headless = true\n\
enableCORS = false\n\
enableXsrfProtection = false\n\
address = "0.0.0.0"\n\
\n\
[theme]\n\
primaryColor = "#2196F3"\n\
' > /app/.streamlit/config.toml

# Run the app
ENTRYPOINT ["streamlit", "run", "webapp/app.py", "--server.port=8501", "--server.address=0.0.0.0"]
