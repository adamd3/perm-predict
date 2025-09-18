from pydantic_settings import BaseSettings
from pathlib import Path


class Settings(BaseSettings):
    # Celery Configuration
    CELERY_BROKER_URL: str
    CELERY_RESULT_BACKEND: str
    
    # ML Model Paths
    MODEL_CLASSIFIER_PATH: str = "app/ml_models/classifier_model.pkl"
    
    # Regression/Ensemble Model Paths (commented out for classification-only mode)
    # MODEL_REGRESSOR_PATH: str = "app/ml_models/regressor_model.pkl"  # Legacy - kept for compatibility
    # MODEL_XGBOOST_REGRESSOR_PATH: str = "app/ml_models/xgboost_regressor.pkl"
    # MODEL_ATTENTIVEFP_PATH: str = "app/ml_models/attentivefp_regressor.pt"
    # MODEL_DIMENET_PATH: str = "app/ml_models/dimenet_regressor.pt"
    # MODEL_BLENDER_PATH: str = "app/ml_models/blender_model.pkl"
    
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
