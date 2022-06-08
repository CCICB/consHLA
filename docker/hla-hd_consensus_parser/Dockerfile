# Use a Python runtime as parent image
FROM python:3.8.6-slim

RUN pip install pandas numpy

# Copy Python script into container
ADD ./hlahd_consensus_parser_v2.py /app/hlahd_consensus_parser_v2.py
WORKDIR /app

CMD ["/bin/bash"]