import { useState } from "react";
import { useNavigate, useLocation } from "react-router-dom";
import {
  Plus,
  Activity,
  Sparkles,
  FolderOpen,
  ChevronRight,
  Folder
} from "lucide-react";
import useStore from "../store/useStore";
import ProjectModal from "./ProjectModal";
import { format } from "date-fns";
import toast from "react-hot-toast";
import "./Sidebar.css";

function Sidebar() {
  const navigate = useNavigate();
  const location = useLocation();
  const [isModalOpen, setIsModalOpen] = useState(false);
  const { projects, setCurrentProject } = useStore();

  const recentProjects = projects.slice(0, 5);

  const handleProjectClick = (project) => {
    setCurrentProject(project.id);

    // Navigate based on project type
    if (project.type === "dti") {
      navigate("/tools/dti-prediction", { state: { projectId: project.id } });
    } else if (project.type === "drug-discovery") {
      navigate("/tools/drug-discovery", { state: { projectId: project.id } });
    } else if (project.type === "both") {
      // For combined projects, show options
      toast(
        (t) => (
          <div
            style={{ display: "flex", flexDirection: "column", gap: "0.75rem" }}
          >
            <div style={{ fontWeight: "600", fontSize: "1rem" }}>
              Choose Research Tool
            </div>
            <div style={{ display: "flex", gap: "0.5rem" }}>
              <button
                onClick={() => {
                  toast.dismiss(t.id);
                  navigate("/tools/dti-prediction", {
                    state: { projectId: project.id }
                  });
                }}
                style={{
                  padding: "0.5rem 1rem",
                  background: "#6366f1",
                  color: "white",
                  border: "none",
                  borderRadius: "6px",
                  cursor: "pointer",
                  fontSize: "0.875rem",
                  fontWeight: "500"
                }}
              >
                DTI Prediction
              </button>
              <button
                onClick={() => {
                  toast.dismiss(t.id);
                  navigate("/tools/drug-discovery", {
                    state: { projectId: project.id }
                  });
                }}
                style={{
                  padding: "0.5rem 1rem",
                  background: "#a855f7",
                  color: "white",
                  border: "none",
                  borderRadius: "6px",
                  cursor: "pointer",
                  fontSize: "0.875rem",
                  fontWeight: "500"
                }}
              >
                Drug Discovery
              </button>
            </div>
          </div>
        ),
        {
          duration: 6000,
          style: {
            background: "#1e293b",
            color: "#fff",
            border: "1px solid #334155",
            minWidth: "300px"
          }
        }
      );
    }
  };

  const tools = [
    {
      id: "dti",
      name: "DTI Prediction",
      icon: Activity,
      path: "/tools/dti-prediction"
    },
    {
      id: "drug",
      name: "Drug Discovery",
      icon: Sparkles,
      path: "/tools/drug-discovery"
    }
  ];

  const getProjectIcon = (type) => {
    if (type === "dti") return <Activity size={16} />;
    if (type === "drug-discovery") return <Sparkles size={16} />;
    return <Folder size={16} />;
  };

  return (
    <>
      <aside className="sidebar">
        {/* Create Project Button */}
        <button
          className="create-project-btn"
          onClick={() => setIsModalOpen(true)}
        >
          <Plus size={20} />
          New Project
        </button>

        {/* Recent Projects */}
        <div className="sidebar-section">
          <h3 className="sidebar-title">
            <FolderOpen size={18} />
            Recent Projects
          </h3>
          <div className="sidebar-list">
            {recentProjects.length > 0 ? (
              recentProjects.map((project) => (
                <div
                  key={project.id}
                  className="sidebar-item"
                  onClick={() => handleProjectClick(project)}
                >
                  <div className="sidebar-item-icon">
                    {getProjectIcon(project.type)}
                  </div>
                  <div className="sidebar-item-content">
                    <span className="sidebar-item-name">{project.name}</span>
                    <span className="sidebar-item-date">
                      {format(new Date(project.lastModified), "MMM dd")}
                    </span>
                  </div>
                  <ChevronRight size={16} className="sidebar-item-arrow" />
                </div>
              ))
            ) : (
              <p className="sidebar-empty">No projects yet</p>
            )}
          </div>
          {projects.length > 5 && (
            <button
              className="sidebar-view-all"
              onClick={() => navigate("/projects")}
            >
              View All Projects
            </button>
          )}
        </div>

        {/* Research Tools */}
        <div className="sidebar-section">
          <h3 className="sidebar-title">Research Tools</h3>
          <div className="sidebar-list">
            {tools.map((tool) => {
              const Icon = tool.icon;
              const isActive = location.pathname === tool.path;
              return (
                <div
                  key={tool.id}
                  className={`sidebar-tool ${isActive ? "active" : ""}`}
                  onClick={() => navigate(tool.path)}
                >
                  <Icon size={20} />
                  <span>{tool.name}</span>
                </div>
              );
            })}
          </div>
        </div>
      </aside>

      <ProjectModal
        isOpen={isModalOpen}
        onClose={() => setIsModalOpen(false)}
      />
    </>
  );
}

export default Sidebar;
