import { useState } from "react";
import {
  Book,
  Activity,
  Sparkles,
  Folder,
  Code,
  Database,
  Zap
} from "lucide-react";
import "./Documentation.css";

function Documentation() {
  const [activeSection, setActiveSection] = useState("getting-started");

  const sections = [
    {
      id: "getting-started",
      title: "Getting Started",
      icon: Zap,
      content: (
        <>
          <h2>Getting Started with DeepBind</h2>
          <p>
            Welcome to DeepBind! This guide will help you get started with our
            platform.
          </p>

          <h3>Quick Start</h3>
          <ol>
            <li>
              <strong>Create an Account:</strong> Sign up using your Google
              account or email.
            </li>
            <li>
              <strong>Create a Project:</strong> Click "New Project" in the
              sidebar and select your research type.
            </li>
            <li>
              <strong>Choose Your Tool:</strong> Select DTI Prediction or Drug
              Discovery based on your needs.
            </li>
            <li>
              <strong>Input Your Data:</strong> Enter molecular structures
              (SMILES) or protein sequences (FASTA).
            </li>
            <li>
              <strong>Analyze Results:</strong> View predictions, save results,
              and iterate on your research.
            </li>
          </ol>

          <h3>System Requirements</h3>
          <ul>
            <li>Modern web browser (Chrome, Firefox, Safari, or Edge)</li>
            <li>Stable internet connection</li>
            <li>JavaScript enabled</li>
            <li>Minimum 4GB RAM recommended for optimal performance</li>
          </ul>
        </>
      )
    },
    {
      id: "dti-prediction",
      title: "DTI Prediction",
      icon: Activity,
      content: (
        <>
          <h2>Drug-Target Interaction Prediction</h2>
          <p>
            Learn how to predict binding affinity between drugs and target
            proteins.
          </p>

          <h3>Input Format</h3>
          <h4>SMILES String</h4>
          <p>
            SMILES (Simplified Molecular Input Line Entry System) represents
            molecular structures:
          </p>
          <div className="code-block">
            <code>Example: CCO (Ethanol)</code>
            <br />
            <code>Example: COC1=C(C=C2C(=C1)CCN=C2C3=CC(=C(C=C3)Cl)Cl)Cl</code>
          </div>

          <h4>Protein Sequence (FASTA)</h4>
          <p>Enter the amino acid sequence of your target protein:</p>
          <div className="code-block">
            <code>
              Example: MTVKTEAAKGTLTYSRMRGMVAILIAFMKQRRMGLNDFIQKIANNS...
            </code>
          </div>

          <h3>Understanding Results</h3>
          <ul>
            <li>
              <strong>Binding Affinity Score:</strong> Higher values indicate
              stronger binding
            </li>
            <li>
              <strong>Score {">"} 10:</strong> Strong binding - High therapeutic
              potential
            </li>
            <li>
              <strong>Score 8-10:</strong> Moderate binding - Worth
              investigating
            </li>
            <li>
              <strong>Score {"<"} 8:</strong> Weak binding - Low priority
            </li>
          </ul>

          <h3>Best Practices</h3>
          <ul>
            <li>
              Validate SMILES strings using RDKit or similar tools before input
            </li>
            <li>Ensure protein sequences are in standard amino acid format</li>
            <li>
              Run multiple predictions with similar molecules for comparison
            </li>
            <li>
              Save promising candidates to your project for further analysis
            </li>
          </ul>
        </>
      )
    },
    {
      id: "drug-discovery",
      title: "Drug Discovery",
      icon: Sparkles,
      content: (
        <>
          <h2>AI-Powered Drug Discovery</h2>
          <p>Generate novel drug candidates using our Generative AI models.</p>

          <h3>How It Works</h3>
          <ol>
            <li>
              <strong>Input Target Protein:</strong> Provide the protein
              sequence you want to target
            </li>
            <li>
              <strong>Set Parameters:</strong> Choose how many drug candidates
              to generate (1-100)
            </li>
            <li>
              <strong>Generation:</strong> Our GAN model creates novel molecular
              structures
            </li>
            <li>
              <strong>Ranking:</strong> Candidates are automatically ranked by
              predicted binding affinity
            </li>
          </ol>

          <h3>Understanding Molecular Properties</h3>
          <h4>Molecular Weight (MW)</h4>
          <p>
            The total mass of the molecule. Drug-like molecules typically have
            MW between 150-500 Da.
          </p>

          <h4>LogP (Partition Coefficient)</h4>
          <p>
            Measures lipophilicity (fat-solubility). Optimal drug-like LogP is
            typically between 0-5.
          </p>

          <h3>Filtering Results</h3>
          <p>Top candidates are selected based on:</p>
          <ul>
            <li>Predicted binding affinity</li>
            <li>Drug-likeness properties (Lipinski's Rule of Five)</li>
            <li>Synthetic accessibility</li>
            <li>Structural novelty</li>
          </ul>

          <h3>Next Steps</h3>
          <p>After generation:</p>
          <ol>
            <li>Review top candidates and their properties</li>
            <li>Export SMILES strings for further analysis</li>
            <li>Use DTI Prediction to validate binding affinity</li>
            <li>Perform additional computational analysis (ADME, toxicity)</li>
          </ol>
        </>
      )
    },
    {
      id: "projects",
      title: "Project Management",
      icon: Folder,
      content: (
        <>
          <h2>Managing Your Research Projects</h2>
          <p>Organize your drug discovery research efficiently.</p>

          <h3>Creating Projects</h3>
          <ol>
            <li>Click "New Project" in the sidebar</li>
            <li>Enter project name and description</li>
            <li>
              Select project type:
              <ul>
                <li>
                  <strong>DTI Prediction:</strong> For binding affinity studies
                </li>
                <li>
                  <strong>Drug Discovery:</strong> For molecule generation
                </li>
                <li>
                  <strong>Both:</strong> For comprehensive research
                </li>
              </ul>
            </li>
            <li>Click "Create Project"</li>
          </ol>

          <h3>Working with Projects</h3>
          <ul>
            <li>
              <strong>Auto-save:</strong> All results are automatically saved to
              active projects
            </li>
            <li>
              <strong>Resume Work:</strong> Click any project to restore your
              last session
            </li>
            <li>
              <strong>View Results:</strong> Access all saved predictions and
              generated molecules
            </li>
            <li>
              <strong>Delete Projects:</strong> Remove projects you no longer
              need
            </li>
          </ul>

          <h3>Data Storage</h3>
          <p>
            Your data is stored locally in your browser using secure
            localStorage:
          </p>
          <ul>
            <li>Data persists across sessions</li>
            <li>No data sent to external servers (except for predictions)</li>
            <li>Clear browser data to reset all projects</li>
          </ul>
        </>
      )
    },
    {
      id: "api",
      title: "API Reference",
      icon: Code,
      content: (
        <>
          <h2>API Documentation</h2>
          <p>Technical reference for developers.</p>

          <h3>DTI Prediction Endpoint</h3>
          <div className="code-block">
            <code>POST /api/predict-dti</code>
          </div>
          <p>
            <strong>Request Body:</strong>
          </p>
          <div className="code-block">
            <code>{`{
  "smiles": "string",
  "protein_sequence": "string"
}`}</code>
          </div>
          <p>
            <strong>Response:</strong>
          </p>
          <div className="code-block">
            <code>{`{
  "binding_affinity": float,
  "smiles": "string",
  "protein_length": int
}`}</code>
          </div>

          <h3>Drug Discovery Endpoint</h3>
          <div className="code-block">
            <code>POST /api/discover-drugs</code>
          </div>
          <p>
            <strong>Request Body:</strong>
          </p>
          <div className="code-block">
            <code>{`{
  "protein_sequence": "string",
  "num_candidates": int (1-100)
}`}</code>
          </div>
          <p>
            <strong>Response:</strong>
          </p>
          <div className="code-block">
            <code>{`{
  "protein_length": int,
  "total_generated": int,
  "valid_molecules": int,
  "candidates": [
    {
      "rank": int,
      "smiles": "string",
      "affinity_score": float,
      "molecular_weight": float,
      "logp": float
    }
  ]
}`}</code>
          </div>

          <h3>Rate Limits</h3>
          <ul>
            <li>DTI Predictions: 100 requests per hour</li>
            <li>Drug Discovery: 20 requests per hour</li>
            <li>Contact us for higher limits</li>
          </ul>
        </>
      )
    },
    {
      id: "models",
      title: "Model Information",
      icon: Database,
      content: (
        <>
          <h2>AI Models & Architecture</h2>
          <p>Technical details about our machine learning models.</p>

          <h3>DTI Prediction Model</h3>
          <ul>
            <li>
              <strong>Architecture:</strong> Graph Convolutional Network (GCN)
            </li>
            <li>
              <strong>Protein Encoding:</strong> ProtBERT (1024-dim embeddings)
            </li>
            <li>
              <strong>Training Data:</strong> KIBA dataset (15,000+
              interactions)
            </li>
            <li>
              <strong>Performance:</strong> MSE {"<"} 0.7 on test set
            </li>
          </ul>

          <h3>Drug Generation Model</h3>
          <ul>
            <li>
              <strong>Architecture:</strong> Generative Adversarial Network
              (GAN)
            </li>
            <li>
              <strong>Latent Dimension:</strong> 56
            </li>
            <li>
              <strong>Output Dimension:</strong> 128 (molecular embeddings)
            </li>
            <li>
              <strong>Training Epochs:</strong> 1000
            </li>
          </ul>

          <h3>Technical Stack</h3>
          <ul>
            <li>
              <strong>Backend:</strong> FastAPI, PyTorch, RDKit
            </li>
            <li>
              <strong>Frontend:</strong> React, Zustand, TailwindCSS
            </li>
            <li>
              <strong>ML Libraries:</strong> PyTorch Geometric, Transformers
            </li>
            <li>
              <strong>Deployment:</strong> CUDA-enabled GPU servers
            </li>
          </ul>

          <h3>Limitations</h3>
          <ul>
            <li>
              Predictions are computational estimates, not experimental results
            </li>
            <li>Generated molecules require synthesis validation</li>
            <li>Toxicity and ADME properties not evaluated by default</li>
            <li>Model trained primarily on kinase inhibitors (KIBA dataset)</li>
          </ul>
        </>
      )
    }
  ];

  const activeContent = sections.find((s) => s.id === activeSection);

  return (
    <div className="documentation-container">
      <div className="documentation-sidebar">
        <div className="docs-header">
          <Book size={32} />
          <h2>Documentation</h2>
        </div>

        <nav className="docs-nav">
          {sections.map((section) => {
            const Icon = section.icon;
            return (
              <button
                key={section.id}
                className={`docs-nav-item ${
                  activeSection === section.id ? "active" : ""
                }`}
                onClick={() => setActiveSection(section.id)}
              >
                <Icon size={20} />
                <span>{section.title}</span>
              </button>
            );
          })}
        </nav>
      </div>

      <div className="documentation-content">
        {activeContent && activeContent.content}
      </div>
    </div>
  );
}

export default Documentation;
