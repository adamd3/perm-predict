from pydantic_settings import BaseSettings
from pathlib import Path


class Settings(BaseSettings):
    MODEL_PATH: str = "app/model/permeability_model.onnx"
    FEATURE_COUNT: int = 4200
    BATCH_SIZE: int = 1000
    FINGERPRINT_RADIUS: int = 2
    LOG_DIR: str = "logs"

    class Config:
        env_file = ".env"


settings = Settings()
