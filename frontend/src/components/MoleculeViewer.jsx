import { useEffect, useState } from "react";
import axios from "axios";

const API_URL = "http://localhost:8000";

function MoleculeViewer({ smiles }) {
  const [imageData, setImageData] = useState(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(false);

  useEffect(() => {
    if (!smiles) {
      setLoading(false);
      return;
    }

    const fetchMoleculeImage = async () => {
      setLoading(true);
      setError(false);

      try {
        const response = await axios.post(
          `${API_URL}/api/render-molecule-base64`,
          {
            smiles: smiles,
            width: 200,
            height: 200
          }
        );

        setImageData(response.data.image);
      } catch (err) {
        console.error("Error rendering molecule:", err);
        setError(true);
      } finally {
        setLoading(false);
      }
    };

    fetchMoleculeImage();
  }, [smiles]);

  if (loading) {
    return (
      <div
        style={{
          width: "200px",
          height: "200px",
          display: "flex",
          alignItems: "center",
          justifyContent: "center",
          background: "#0f172a",
          border: "1px solid #334155",
          borderRadius: "8px",
          color: "#6366f1"
        }}
      >
        Loading...
      </div>
    );
  }

  if (error || !imageData) {
    return (
      <div
        style={{
          width: "200px",
          height: "200px",
          display: "flex",
          alignItems: "center",
          justifyContent: "center",
          background: "#0f172a",
          border: "1px solid #334155",
          borderRadius: "8px",
          color: "#64748b",
          fontSize: "0.875rem",
          textAlign: "center",
          padding: "1rem"
        }}
      >
        Unable to render
      </div>
    );
  }

  return (
    <img
      src={imageData}
      alt="Molecule Structure"
      style={{
        width: "200px",
        height: "200px",
        border: "1px solid #334155",
        borderRadius: "8px",
        background: "#0f172a"
      }}
    />
  );
}

export default MoleculeViewer;
