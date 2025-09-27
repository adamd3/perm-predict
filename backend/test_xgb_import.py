import xgboost as xgb
import sys

try:
    print(f"XGBoost import successful. Version: {xgb.__version__}")
except Exception as e:
    print(f"Error importing XGBoost: {e}", file=sys.stderr)
    sys.exit(1)
