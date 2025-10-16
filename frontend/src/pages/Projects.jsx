import { useState } from "react";
import { useNavigate } from "react-router-dom";
import {
  Activity,
  Sparkles,
  Folder,
  Trash2,
  Calendar,
  FileText
} from "lucide-react";
import useStore from "../store/useStore";
import { format } from "date-fns";
import toast from "react-hot-toast";
import "./Projects.css";

function Projects() {
  const navigate = useNavigate();
  const { projects, deleteProject, setCurrentProject } = useStore();
  const [filter, setFilter] = useState("all");

  const filteredProjects = projects.filter((p) => {
    if (filter === "all") return true;
    return p.type === filter;
  });

  const handleDelete = (projectId, projectName) => {
    if (window.confirm(`Delete project "${projectName}"?`)) {
      deleteProject(projectId);
      toast.success("Project deleted successfully");
    }
  };

  const handleProjectClick = (project) => {
    setCurrentProject(project.id);

    // Navigate to the appropriate tool based on project type
    if (project.type === "dti") {
      navigate("/tools/dti-prediction", { state: { projectId: project.id } });
    } else if (project.type === "drug-discovery") {
      navigate("/tools/drug-discovery", { state: { projectId: project.id } });
    } else {
      // For 'both' type, stay on projects page or go to a specific view
      toast.info(
        "This project supports both tools. Select a tool from the sidebar."
      );
    }
  };

  return (
    <div className="projects-container">
      <div className="projects-header">
        <div>
          <h1>My Projects</h1>
          <p>Manage and organize your research projects</p>
        </div>
      </div>

      <div className="projects-filters">
        <button
          className={`filter-btn ${filter === "all" ? "active" : ""}`}
          onClick={() => setFilter("all")}
        >
          All Projects ({projects.length})
        </button>
        <button
          className={`filter-btn ${filter === "dti" ? "active" : ""}`}
          onClick={() => setFilter("dti")}
        >
          <Activity size={16} />
          DTI ({projects.filter((p) => p.type === "dti").length})
        </button>
        <button
          className={`filter-btn ${
            filter === "drug-discovery" ? "active" : ""
          }`}
          onClick={() => setFilter("drug-discovery")}
        >
          <Sparkles size={16} />
          Drug Discovery (
          {projects.filter((p) => p.type === "drug-discovery").length})
        </button>
      </div>

      {filteredProjects.length > 0 ? (
        <div className="projects-grid">
          {filteredProjects.map((project) => (
            <div key={project.id} className="project-card">
              <div className="project-card-header">
                <div className="project-type-badge">
                  {project.type === "dti" ? (
                    <Activity size={16} />
                  ) : (
                    <Sparkles size={16} />
                  )}
                  <span>
                    {project.type === "dti" ? "DTI" : "Drug Discovery"}
                  </span>
                </div>
                <button
                  className="delete-btn"
                  onClick={(e) => {
                    e.stopPropagation();
                    handleDelete(project.id, project.name);
                  }}
                >
                  <Trash2 size={16} />
                </button>
              </div>

              <div
                className="project-card-body"
                onClick={() => handleProjectClick(project)}
              >
                <h3>{project.name}</h3>
                <p>{project.description || "No description provided"}</p>

                <div className="project-meta">
                  <div className="meta-item">
                    <Calendar size={14} />
                    <span>
                      {format(new Date(project.createdAt), "MMM dd, yyyy")}
                    </span>
                  </div>
                  <div className="meta-item">
                    <FileText size={14} />
                    <span>{project.results?.length || 0} results</span>
                  </div>
                </div>
              </div>
            </div>
          ))}
        </div>
      ) : (
        <div className="empty-state">
          <Folder size={64} />
          <h2>No Projects Found</h2>
          <p>
            {filter === "all"
              ? "Create your first project to get started"
              : `No ${filter} projects yet`}
          </p>
        </div>
      )}
    </div>
  );
}

export default Projects;
