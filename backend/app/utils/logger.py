import logging
from pathlib import Path


def setup_logging(log_dir: str = "logs"):
    Path(log_dir).mkdir(exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.handlers.RotatingFileHandler(f"{log_dir}/app.log", maxBytes=1024 * 1024, backupCount=3),  # 1MB
            logging.StreamHandler(),
        ],
    )

    return logging.getLogger(__name__)
