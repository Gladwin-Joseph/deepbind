import { useState } from "react";
import { useNavigate } from "react-router-dom";
import { X, Activity, Sparkles, Layers } from "lucide-react";
import useStore from "../store/useStore";
import toast from "react-hot-toast";
import "./ProjectModal.css";

function ProjectModal({ isOpen, onClose }) {
  const navigate = useNavigate();
  const createProject = useStore((state) => state.createProject);

  const [formData, setFormData] = useState({
    name: "",
    description: "",
    type: "dti"
  });

  const projectTypes = [
    {
      value: "dti",
      label: "DTI Prediction",
      icon: Activity,
      description: "Drug-target interaction analysis"
    },
    {
      value: "drug-discovery",
      label: "Drug Discovery",
      icon: Sparkles,
      description: "AI-powered molecule generation"
    },
    {
      value: "both",
      label: "Combined Research",
      icon: Layers,
      description: "Both DTI and drug generation"
    }
  ];

  const handleSubmit = (e) => {
    e.preventDefault();

    if (!formData.name.trim()) {
      toast.error("Project name is required");
      return;
    }

    const projectId = createProject(formData);
    toast.success("Project created successfully!");

    // Reset form
    setFormData({ name: "", description: "", type: "dti" });
    onClose();

    // Navigate based on project type
    if (formData.type === "dti") {
      navigate("/tools/dti-prediction", { state: { projectId } });
    } else if (formData.type === "drug-discovery") {
      navigate("/tools/drug-discovery", { state: { projectId } });
    } else {
      navigate("/projects");
    }
  };

  if (!isOpen) return null;

  return (
    <div className="modal-overlay" onClick={onClose}>
      <div className="modal-content" onClick={(e) => e.stopPropagation()}>
        <div className="modal-header">
          <h2>Create New Project</h2>
          <button className="modal-close" onClick={onClose}>
            <X size={24} />
          </button>
        </div>

        <form onSubmit={handleSubmit} className="modal-form">
          <div className="form-group">
            <label htmlFor="project-name">
              Project Name <span className="required">*</span>
            </label>
            <input
              id="project-name"
              type="text"
              value={formData.name}
              onChange={(e) =>
                setFormData({ ...formData, name: e.target.value })
              }
              placeholder="e.g., Cancer Drug Discovery"
              className="modal-input"
              required
            />
          </div>

          <div className="form-group">
            <label htmlFor="project-description">Description</label>
            <textarea
              id="project-description"
              value={formData.description}
              onChange={(e) =>
                setFormData({ ...formData, description: e.target.value })
              }
              placeholder="Brief description of your research project..."
              className="modal-textarea"
              rows={3}
            />
          </div>

          <div className="form-group">
            <label>Project Type</label>
            <div className="project-types">
              {projectTypes.map((type) => {
                const Icon = type.icon;
                return (
                  <div
                    key={type.value}
                    className={`project-type-card ${
                      formData.type === type.value ? "selected" : ""
                    }`}
                    onClick={() =>
                      setFormData({ ...formData, type: type.value })
                    }
                  >
                    <Icon size={24} />
                    <div className="project-type-info">
                      <span className="project-type-label">{type.label}</span>
                      <span className="project-type-description">
                        {type.description}
                      </span>
                    </div>
                  </div>
                );
              })}
            </div>
          </div>

          <div className="modal-actions">
            <button type="button" className="btn-secondary" onClick={onClose}>
              Cancel
            </button>
            <button type="submit" className="btn-primary">
              Create Project
            </button>
          </div>
        </form>
      </div>
    </div>
  );
}

export default ProjectModal;
