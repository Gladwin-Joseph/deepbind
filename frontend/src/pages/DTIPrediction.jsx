import { useState, useEffect } from "react";
import { useLocation } from "react-router-dom";
import axios from "axios";
import {
  Send,
  Loader2,
  AlertCircle,
  CheckCircle2,
  Activity
} from "lucide-react";
import useStore from "../store/useStore";
import toast from "react-hot-toast";
import "./DTIPrediction.css";

const API_URL = "http://localhost:8000";

function DTIPrediction() {
  const location = useLocation();
  const {
    currentProjectId,
    addResultToProject,
    setCurrentProject,
    getProjectViewState,
    saveProjectViewState
  } = useStore();
  const [smiles, setSmiles] = useState("");
  const [proteinSeq, setProteinSeq] = useState("");
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState(null);
  const [error, setError] = useState(null);
  const [loadedFromState, setLoadedFromState] = useState(false);

  const exampleData = {
    smiles: "COC1=C(C=C2C(=C1)CCN=C2C3=CC(=C(C=C3)Cl)Cl)Cl",
    protein: "MTVKTEAAKGTLTYSRMRGMVAILIAFMKQRRMGLNDFIQKIANNSYISFRFKSELVSSGC"
  };

  // Load previous state when project is selected
  useEffect(() => {
    if (location.state?.projectId && !loadedFromState) {
      const projectId = location.state.projectId;
      setCurrentProject(projectId);

      // Restore previous view state
      const viewState = getProjectViewState(projectId);
      console.log("Loading view state for DTI project:", projectId, viewState);

      if (viewState && viewState.dti) {
        setSmiles(viewState.dti.smiles || "");
        setProteinSeq(viewState.dti.proteinSeq || "");
        setResult(viewState.dti.result || null);
        setLoadedFromState(true);
        toast.success("Previous session restored!");
      }
    }
  }, [location.state, setCurrentProject, getProjectViewState, loadedFromState]);

  // Save view state whenever inputs or results change
  useEffect(() => {
    if (currentProjectId && loadedFromState) {
      console.log("Saving view state for DTI project:", currentProjectId);
      saveProjectViewState(currentProjectId, {
        dti: {
          smiles,
          proteinSeq,
          result
        }
      });
    }
  }, [
    currentProjectId,
    smiles,
    proteinSeq,
    result,
    saveProjectViewState,
    loadedFromState
  ]);

  const handleSubmit = async (e) => {
    e.preventDefault();
    setLoading(true);
    setError(null);
    setResult(null);
    setLoadedFromState(true); // Enable auto-save after first submission

    try {
      const response = await axios.post(`${API_URL}/api/predict-dti`, {
        smiles: smiles.trim(),
        protein_sequence: proteinSeq.trim()
      });

      setResult(response.data);

      // Save to current project if exists
      if (currentProjectId) {
        addResultToProject(currentProjectId, {
          type: "dti",
          smiles: smiles.trim(),
          proteinSequence: proteinSeq.trim(),
          bindingAffinity: response.data.binding_affinity,
          proteinLength: response.data.protein_length
        });
        toast.success("Result saved to project!");
      }
    } catch (err) {
      setError(
        err.response?.data?.detail ||
          "Prediction failed. Please check your inputs."
      );
      toast.error("Prediction failed!");
    } finally {
      setLoading(false);
    }
  };

  const loadExample = () => {
    setSmiles(exampleData.smiles);
    setProteinSeq(exampleData.protein);
    setLoadedFromState(true);
  };

  return (
    <div className="content-grid">
      <div className="input-section">
        <form onSubmit={handleSubmit} className="prediction-form">
          {currentProjectId && (
            <div className="project-indicator">
              <Activity size={16} />
              <span>Working on project</span>
            </div>
          )}

          <div className="form-group">
            <label htmlFor="smiles">
              Drug SMILES String
              <span className="required">*</span>
            </label>
            <input
              id="smiles"
              type="text"
              value={smiles}
              onChange={(e) => setSmiles(e.target.value)}
              placeholder="e.g., COC1=C(C=C2C(=C1)CCN=C2C3=CC(=C(C=C3)Cl)Cl)Cl"
              required
              className="input-field"
            />
            <span className="input-hint">
              Enter molecular structure in SMILES notation
            </span>
          </div>

          <div className="form-group">
            <label htmlFor="protein">
              Protein Sequence (FASTA)
              <span className="required">*</span>
            </label>
            <textarea
              id="protein"
              value={proteinSeq}
              onChange={(e) => setProteinSeq(e.target.value)}
              placeholder="e.g., MTVKTEAAKGTLTYSRMRGMVAILIAFMKQRRMGLNDFIQKIANNS..."
              required
              rows={6}
              className="input-field"
            />
            <span className="input-hint">Enter amino acid sequence</span>
          </div>

          <div className="button-group">
            <button
              type="button"
              onClick={loadExample}
              className="btn btn-secondary"
              disabled={loading}
            >
              Load Example
            </button>

            <button
              type="submit"
              className="btn btn-primary"
              disabled={loading || !smiles || !proteinSeq}
            >
              {loading ? (
                <>
                  <Loader2 size={20} className="spin" />
                  Predicting...
                </>
              ) : (
                <>
                  <Send size={20} />
                  Predict Affinity
                </>
              )}
            </button>
          </div>
        </form>
      </div>

      <div className="result-section">
        {error && (
          <div className="alert alert-error">
            <AlertCircle size={20} />
            <div>
              <strong>Error</strong>
              <p>{error}</p>
            </div>
          </div>
        )}

        {result && (
          <div className="result-card">
            <div className="result-header">
              <CheckCircle2 size={24} className="success-icon" />
              <h2>Prediction Result</h2>
            </div>

            <div className="affinity-display">
              <span className="affinity-label">Binding Affinity Score</span>
              <span className="affinity-value">
                {result.binding_affinity.toFixed(4)}
              </span>
            </div>

            <div className="result-details">
              <div className="detail-item">
                <span className="detail-label">SMILES</span>
                <span className="detail-value mono">{result.smiles}</span>
              </div>

              <div className="detail-item">
                <span className="detail-label">Protein Length</span>
                <span className="detail-value">
                  {result.protein_length} amino acids
                </span>
              </div>
            </div>

            <div className="interpretation">
              <h3>Interpretation</h3>
              <p>
                {result.binding_affinity > 10
                  ? "Strong binding affinity - High potential for interaction"
                  : result.binding_affinity > 8
                  ? "Moderate binding affinity - Potential for interaction"
                  : "Weak binding affinity - Low potential for interaction"}
              </p>
            </div>
          </div>
        )}

        {!error && !result && (
          <div className="placeholder">
            <Activity size={48} className="placeholder-icon" />
            <p>
              Enter drug and protein information to predict binding affinity
            </p>
          </div>
        )}
      </div>
    </div>
  );
}

export default DTIPrediction;
