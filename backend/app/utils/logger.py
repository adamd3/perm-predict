import logging
from pathlib import Path

LOG_DIR = "logs"
Path(LOG_DIR).mkdir(exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.handlers.RotatingFileHandler(f"{LOG_DIR}/app.log", maxBytes=1024 * 1024, backupCount=3),  # 1MB
        logging.StreamHandler(),
    ],
)

logger = logging.getLogger(__name__)
