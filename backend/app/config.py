from pydantic_settings import BaseSettings
from pathlib import Path


class Settings(BaseSettings):
    # Celery Configuration
    CELERY_BROKER_URL: str = "redis://localhost:6379/0"
    CELERY_RESULT_BACKEND: str = "redis://localhost:6379/0"
    
    # ML Model Paths
    MODEL_CLASSIFIER_PATH: str = "app/ml_models/classifier_model.pkl"
    MODEL_REGRESSOR_PATH: str = "app/ml_models/regressor_model.pkl"
    
    # Feature Extraction Settings
    FEATURE_COUNT: int = 4200
    FINGERPRINT_RADIUS: int = 2
    
    # Processing Settings
    BATCH_SIZE: int = 1000
    
    # Logging
    LOG_DIR: str = "logs"
    LOG_LEVEL: str = "INFO"
    
    # API Settings
    API_TITLE: str = "Perm-Predict GraphQL API"
    API_VERSION: str = "2.0.0"
    API_DESCRIPTION: str = "Machine learning-based prediction of chemical accumulation in bacteria"

    class Config:
        env_file = ".env"


settings = Settings()
