import os
import torch

# Configuration
USE_CPU = os.getenv("USE_CPU", "true").lower() == "true"
DEVICE = torch.device("cpu") if USE_CPU else torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Model paths
BASE_PATH = os.path.dirname(os.path.dirname(__file__))
MODELS_PATH = os.path.join(BASE_PATH, "models")
DTI_MODEL_PATH = os.path.join(MODELS_PATH, "best_dti_model.pth")
GENERATOR_MODEL_PATH = os.path.join(MODELS_PATH, "generator_model.pth")

# Feature flags
LOAD_MODELS_ON_STARTUP = os.getenv("LOAD_MODELS_ON_STARTUP", "false").lower() == "true"