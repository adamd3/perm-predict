import joblib
from onnxmltools import convert_xgboost
from onnxconverter_common.data_types import FloatTensorType


def convert_to_onnx(xgboost_model_path: str, onnx_model_path: str, feature_count: int = 4200):
    model = joblib.load(xgboost_model_path)
    booster = model.get_booster()

    # Fix feature names
    original_feature_names = booster.feature_names
    if original_feature_names is not None:
        booster.feature_names = [f"f{i}" for i in range(len(original_feature_names))]

    initial_type = [("float_input", FloatTensorType([None, feature_count]))]
    onnx_model = convert_xgboost(booster, initial_types=initial_type, target_opset=15)

    with open(onnx_model_path, "wb") as f:
        f.write(onnx_model.SerializeToString())


if __name__ == "__main__":
    convert_to_onnx(
        xgboost_model_path="app/model/xgboost_model.pkl", onnx_model_path="app/model/permeability_model.onnx"
    )
