//export const API_URL = "http://localhost:8000";

// For production (update this after deploying backend)
export const API_URL = import.meta.env.VITE_API_URL || "http://localhost:8000";

export const PROJECT_TYPES = {
  DTI: "dti",
  DRUG_DISCOVERY: "drug-discovery",
  BOTH: "both"
};

export const PROJECT_TYPE_LABELS = {
  [PROJECT_TYPES.DTI]: "DTI Prediction",
  [PROJECT_TYPES.DRUG_DISCOVERY]: "Drug Discovery",
  [PROJECT_TYPES.BOTH]: "Combined Research"
};

export const ROUTES = {
  HOME: "/",
  LOGIN: "/login",
  DTI: "/tools/dti-prediction",
  DRUG_DISCOVERY: "/tools/drug-discovery",
  PROJECTS: "/projects",
  PROJECT_DETAIL: "/projects/:id"
};
