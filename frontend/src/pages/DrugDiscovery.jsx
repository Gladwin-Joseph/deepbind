import { useState, useEffect } from "react";
import { useLocation } from "react-router-dom";
import axios from "axios";
import {
  Sparkles,
  Loader2,
  AlertCircle,
  TrendingUp,
  Beaker,
  Download
} from "lucide-react";
import MoleculeViewer from "../components/MoleculeViewer";
import useStore from "../store/useStore";
import toast from "react-hot-toast";
import jsPDF from "jspdf";
import autoTable from "jspdf-autotable";
import "./DrugDiscovery.css";

const API_URL = "http://localhost:8000";

function DrugDiscovery() {
  const location = useLocation();
  const {
    currentProjectId,
    addResultToProject,
    setCurrentProject,
    getProjectViewState,
    saveProjectViewState
  } = useStore();
  const [proteinSeq, setProteinSeq] = useState("");
  const [numCandidates, setNumCandidates] = useState(10);
  const [loading, setLoading] = useState(false);
  const [results, setResults] = useState(null);
  const [error, setError] = useState(null);
  const [loadedFromState, setLoadedFromState] = useState(false);

  const exampleProtein =
    "MKTLLILALLTLVSLSGCGGVFQESPPFLSRRAKLGSGQLEVLQREQKGDYTAESLWLPEGADNLRDLLHNILLRLLNEKAASPLEEPAQLRRIYNQVDLQRLISQVMGFGKVYFSGLSDRKQAVGGVLLPQIQLHCVGDGVEYLFLPKGFGGHMRERFNRLLEERFDKFAERMERGVDLVSTTGFEGQVARQLGRLFKAKYPFLASITEKVRKICQENGLSVAIGVLVSFAGVGASAQESPLKEKLGTLQFTYDHSLEAQIIPEPEGPLEEVAVLGHPIEQKEGPMCDLLWSDPDKDSGLVPSEEVLAKRFQGEFGNMYHPHERPTLLQRFLEWLKDKKLLTPTGCRMNIETLVNSELDGEEQLAREIKILAKLMVAHTGKTVGEEAHYIEGGMVNLKKQKAVYPQDLVTGVEVLRHQGHHTFIEYKQVDVEMSIFKFPPDKIVAEERQIRKQNELEEELKAEREALNLEREALNHTEAKQKKLELFPVLEQLHKAYPLRV";

  useEffect(() => {
    if (location.state?.projectId && !loadedFromState) {
      const projectId = location.state.projectId;
      setCurrentProject(projectId);

      const viewState = getProjectViewState(projectId);
      console.log("Loading view state for project:", projectId, viewState);

      if (viewState && viewState.drugDiscovery) {
        setProteinSeq(viewState.drugDiscovery.proteinSeq || "");
        setNumCandidates(viewState.drugDiscovery.numCandidates || 10);
        setResults(viewState.drugDiscovery.results || null);
        setLoadedFromState(true);
        toast.success("Previous session restored!");
      }
    }
  }, [location.state, setCurrentProject, getProjectViewState, loadedFromState]);

  useEffect(() => {
    if (currentProjectId && loadedFromState) {
      console.log("Saving view state for project:", currentProjectId);
      saveProjectViewState(currentProjectId, {
        drugDiscovery: {
          proteinSeq,
          numCandidates,
          results
        }
      });
    }
  }, [
    currentProjectId,
    proteinSeq,
    numCandidates,
    results,
    saveProjectViewState,
    loadedFromState
  ]);

  const handleSubmit = async (e) => {
    e.preventDefault();
    setLoading(true);
    setError(null);
    setResults(null);
    setLoadedFromState(true);

    try {
      const response = await axios.post(`${API_URL}/api/discover-drugs`, {
        protein_sequence: proteinSeq.trim(),
        num_candidates: parseInt(numCandidates)
      });

      setResults(response.data);

      if (currentProjectId) {
        addResultToProject(currentProjectId, {
          type: "drug-discovery",
          proteinSequence: proteinSeq.trim(),
          candidatesCount: response.data.valid_molecules,
          topCandidates: response.data.candidates.slice(0, 5)
        });
        toast.success("Results saved to project!");
      }
    } catch (err) {
      setError(
        err.response?.data?.detail ||
          "Drug generation failed. Please try again."
      );
      toast.error("Drug generation failed!");
    } finally {
      setLoading(false);
    }
  };

  const downloadPDF = () => {
    try {
      if (!results) {
        toast.error("No results to download!");
        return;
      }

      // Create new PDF document
      const doc = new jsPDF();
      const pageWidth = doc.internal.pageSize.getWidth();
      let yPosition = 20;

      // Add Title
      doc.setFontSize(22);
      doc.setTextColor(99, 102, 241);
      doc.text("Drug Discovery Report", pageWidth / 2, yPosition, {
        align: "center"
      });

      yPosition += 8;
      doc.setFontSize(12);
      doc.setTextColor(100, 100, 100);
      doc.text("AI-Generated Drug Candidates", pageWidth / 2, yPosition, {
        align: "center"
      });

      yPosition += 10;
      doc.setFontSize(10);
      doc.setTextColor(80, 80, 80);
      const currentDate = new Date().toLocaleString();
      doc.text(`Generated: ${currentDate}`, pageWidth / 2, yPosition, {
        align: "center"
      });

      // Summary Section
      yPosition += 15;
      doc.setFontSize(16);
      doc.setTextColor(0, 0, 0);
      doc.text("Summary", 14, yPosition);

      yPosition += 8;
      doc.setFontSize(11);
      doc.setTextColor(60, 60, 60);

      doc.text(
        `Number of Candidates Requested: ${numCandidates}`,
        14,
        yPosition
      );
      yPosition += 6;

      doc.text(`Total Generated: ${results.total_generated}`, 14, yPosition);
      yPosition += 6;

      doc.text(`Valid Molecules: ${results.valid_molecules}`, 14, yPosition);
      yPosition += 6;

      doc.text(
        `Protein Length: ${results.protein_length} amino acids`,
        14,
        yPosition
      );

      // Protein Sequence Section
      yPosition += 12;
      doc.setFontSize(14);
      doc.setTextColor(0, 0, 0);
      doc.text("Target Protein Sequence", 14, yPosition);

      yPosition += 8;
      doc.setFontSize(9);
      doc.setFont(undefined, "normal");
      doc.setTextColor(60, 60, 60);

      // Split protein sequence into multiple lines
      const maxCharsPerLine = 80;
      const proteinLines = [];
      for (let i = 0; i < proteinSeq.length; i += maxCharsPerLine) {
        proteinLines.push(proteinSeq.substring(i, i + maxCharsPerLine));
      }

      // Show first 5 lines (400 chars) with ellipsis if longer
      const linesToShow = proteinLines.slice(0, 5);
      linesToShow.forEach((line, index) => {
        doc.text(line, 14, yPosition + index * 5);
      });

      if (proteinLines.length > 5) {
        yPosition += linesToShow.length * 5 + 5;
        doc.setFontSize(8);
        doc.setFont(undefined, "italic");
        doc.text(
          `... (showing first ${linesToShow.length * maxCharsPerLine} of ${
            proteinSeq.length
          } characters)`,
          14,
          yPosition
        );
      } else {
        yPosition += linesToShow.length * 5;
      }

      // Top Drug Candidates Table
      yPosition += 12;
      doc.setFontSize(16);
      doc.setFont(undefined, "bold");
      doc.setTextColor(0, 0, 0);
      doc.text("Top Drug Candidates", 14, yPosition);

      // Prepare table data
      const tableData = results.candidates.map((candidate) => [
        `#${candidate.rank}`,
        candidate.smiles.substring(0, 50) +
          (candidate.smiles.length > 50 ? "..." : ""),
        candidate.affinity_score.toFixed(4),
        candidate.molecular_weight
          ? candidate.molecular_weight.toFixed(2)
          : "N/A",
        candidate.logp ? candidate.logp.toFixed(2) : "N/A"
      ]);

      // Add table using autoTable
      autoTable(doc, {
        startY: yPosition + 5,
        head: [["Rank", "SMILES", "Affinity Score", "MW (Da)", "LogP"]],
        body: tableData,
        theme: "grid",
        headStyles: {
          fillColor: [99, 102, 241],
          textColor: [255, 255, 255],
          fontSize: 10,
          fontStyle: "bold",
          halign: "center"
        },
        bodyStyles: {
          fontSize: 9,
          cellPadding: 3
        },
        columnStyles: {
          0: { cellWidth: 15, halign: "center" },
          1: { cellWidth: 80, fontSize: 8 },
          2: { cellWidth: 25, halign: "center" },
          3: { cellWidth: 20, halign: "center" },
          4: { cellWidth: 20, halign: "center" }
        },
        alternateRowStyles: {
          fillColor: [245, 247, 250]
        },
        margin: { left: 14, right: 14 }
      });

      // Footer
      const finalY = doc.lastAutoTable.finalY + 15;
      doc.setFontSize(8);
      doc.setFont(undefined, "italic");
      doc.setTextColor(120, 120, 120);
      doc.text("DeepBind Platform", pageWidth / 2, finalY, {
        align: "center"
      });
      doc.text(
        "For research purposes only. Validate all candidates experimentally.",
        pageWidth / 2,
        finalY + 4,
        { align: "center" }
      );

      // Save the PDF
      const fileName = `drug-discovery-${
        new Date().toISOString().split("T")[0]
      }.pdf`;
      doc.save(fileName);

      toast.success("PDF downloaded successfully!");
    } catch (error) {
      console.error("PDF generation error:", error);
      toast.error("Failed to generate PDF. Please try again.");
    }
  };

  const loadExample = () => {
    setProteinSeq(exampleProtein);
    setLoadedFromState(true);
  };

  const getScoreColor = (score) => {
    if (score > 10) return "#10b981";
    if (score > 9) return "#3b82f6";
    if (score > 8) return "#f59e0b";
    return "#ef4444";
  };

  return (
    <div className="discovery-layout">
      <div className="input-panel">
        <form onSubmit={handleSubmit} className="discovery-form">
          {currentProjectId && (
            <div className="project-indicator">
              <Sparkles size={16} />
              <span>Working on project</span>
            </div>
          )}

          <div className="form-group">
            <label htmlFor="protein-seq">
              Target Protein Sequence
              <span className="required">*</span>
            </label>
            <textarea
              id="protein-seq"
              value={proteinSeq}
              onChange={(e) => setProteinSeq(e.target.value)}
              placeholder="Enter protein amino acid sequence..."
              required
              rows={8}
              className="input-field"
            />
            <span className="input-hint">
              {proteinSeq.length > 0
                ? `${proteinSeq.length} amino acids`
                : "Enter FASTA sequence"}
            </span>
          </div>

          <div className="form-group">
            <label htmlFor="num-candidates">Number of Candidates</label>
            <input
              id="num-candidates"
              type="number"
              value={numCandidates}
              onChange={(e) => setNumCandidates(e.target.value)}
              min="1"
              max="100"
              className="input-field"
            />
            <span className="input-hint">Generate 1-100 drug candidates</span>
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
              disabled={loading || !proteinSeq}
            >
              {loading ? (
                <>
                  <Loader2 size={20} className="spin" />
                  Generating...
                </>
              ) : (
                <>
                  <Sparkles size={20} />
                  Generate Drugs
                </>
              )}
            </button>
          </div>
        </form>
      </div>

      <div className="results-panel">
        {error && (
          <div className="alert alert-error">
            <AlertCircle size={20} />
            <div>
              <strong>Error</strong>
              <p>{error}</p>
            </div>
          </div>
        )}

        {results && (
          <div className="results-container">
            <div className="results-header">
              <div className="results-title-section">
                <h2>Generated Drug Candidates</h2>
                <button
                  className="btn btn-primary download-pdf-btn"
                  onClick={downloadPDF}
                >
                  <Download size={20} />
                  Download PDF
                </button>
              </div>
              <div className="stats">
                <div className="stat-item">
                  <span className="stat-label">Valid Molecules</span>
                  <span className="stat-value">
                    {results.valid_molecules}/{results.total_generated}
                  </span>
                </div>
                <div className="stat-item">
                  <span className="stat-label">Protein Length</span>
                  <span className="stat-value">
                    {results.protein_length} aa
                  </span>
                </div>
              </div>
            </div>

            <div className="candidates-list">
              {results.candidates.map((candidate) => (
                <div key={candidate.rank} className="candidate-card">
                  <div className="candidate-visual">
                    <MoleculeViewer smiles={candidate.smiles} />
                    <div
                      className="rank-badge"
                      style={{
                        backgroundColor: getScoreColor(candidate.affinity_score)
                      }}
                    >
                      #{candidate.rank}
                    </div>
                  </div>

                  <div className="candidate-body">
                    <div className="candidate-header">
                      <div className="affinity-score">
                        <TrendingUp size={16} />
                        <span>{candidate.affinity_score.toFixed(4)}</span>
                      </div>
                    </div>

                    <div className="smiles-display">
                      <Beaker size={16} />
                      <code>{candidate.smiles}</code>
                    </div>

                    <div className="properties">
                      <div className="property">
                        <span className="property-label">MW</span>
                        <span className="property-value">
                          {candidate.molecular_weight?.toFixed(2) || "N/A"}
                        </span>
                      </div>
                      <div className="property">
                        <span className="property-label">LogP</span>
                        <span className="property-value">
                          {candidate.logp?.toFixed(2) || "N/A"}
                        </span>
                      </div>
                    </div>
                  </div>
                </div>
              ))}
            </div>
          </div>
        )}

        {!error && !results && (
          <div className="placeholder">
            <Sparkles size={48} className="placeholder-icon" />
            <p>Enter target protein sequence to generate drug candidates</p>
            <p className="placeholder-hint">
              Our AI will generate and rank novel molecules based on predicted
              binding affinity
            </p>
          </div>
        )}
      </div>
    </div>
  );
}

export default DrugDiscovery;
