import logging

# Create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Set the logging output to the desired format.
# https://docs.python.org/3.7/library/logging.html#logrecord-attributes
formatter = logging.Formatter(
    "{asctime} [{levelname}]: {message}", datefmt="%Y-%m-%d %H:%M:%S", style="{"
)

# Log to a file. Log only warning to the log file.
file_handler = logging.FileHandler("karma.log")
file_handler.setFormatter(formatter)
file_handler.setLevel(logging.WARNING)
logger.addHandler(file_handler)

# Log to the stdout.
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)
