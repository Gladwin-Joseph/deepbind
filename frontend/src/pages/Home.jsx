import { useUser } from "@clerk/clerk-react";
import { useNavigate } from "react-router-dom";
import {
  Activity,
  Sparkles,
  Folder,
  TrendingUp,
  Plus,
  ArrowRight
} from "lucide-react";
import useStore from "../store/useStore";
import { format } from "date-fns";
import toast from "react-hot-toast";
import "./Home.css";

function Home() {
  const { user } = useUser();
  const navigate = useNavigate();
  const { stats, projects, setCurrentProject } = useStore();

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

  const statCards = [
    {
      icon: Activity,
      label: "DTI Predictions",
      value: stats.totalDTIPredictions,
      color: "#6366f1",
      bgColor: "rgba(99, 102, 241, 0.1)"
    },
    {
      icon: Sparkles,
      label: "Molecules Generated",
      value: stats.totalDrugsGenerated,
      color: "#a855f7",
      bgColor: "rgba(168, 85, 247, 0.1)"
    },
    {
      icon: Folder,
      label: "Active Projects",
      value: stats.activeProjects,
      color: "#10b981",
      bgColor: "rgba(16, 185, 129, 0.1)"
    },
    {
      icon: TrendingUp,
      label: "Success Rate",
      value: "94%",
      color: "#f59e0b",
      bgColor: "rgba(245, 158, 11, 0.1)"
    }
  ];

  return (
    <div className="home-container">
      {/* Hero Section */}
      <div className="home-hero">
        <h1>Hello, {user?.firstName || "Researcher"}!</h1>
        <p>
          Your molecular research command center is ready. Explore advanced
          AI-driven drug discovery tools and accelerate breakthrough
          discoveries.
        </p>
      </div>

      {/* Stats Grid */}
      <div className="stats-grid">
        {statCards.map((stat, index) => {
          const Icon = stat.icon;
          return (
            <div
              key={index}
              className="stat-card"
              style={{ borderLeftColor: stat.color }}
            >
              <div
                className="stat-icon"
                style={{ backgroundColor: stat.bgColor }}
              >
                <Icon size={24} style={{ color: stat.color }} />
              </div>
              <div className="stat-content">
                <span className="stat-value">{stat.value}</span>
                <span className="stat-label">{stat.label}</span>
              </div>
            </div>
          );
        })}
      </div>

      {/* Quick Actions */}
      <div className="quick-actions-section">
        <h2>Quick Actions</h2>
        <div className="action-cards">
          <div
            className="action-card"
            onClick={() => navigate("/tools/dti-prediction")}
          >
            <div
              className="action-icon"
              style={{
                background: "linear-gradient(135deg, #6366f1 0%, #4f46e5 100%)"
              }}
            >
              <Activity size={32} />
            </div>
            <h3>DTI Prediction</h3>
            <p>
              Predict drug-target binding affinity using advanced graph neural
              networks
            </p>
            <button className="action-button">
              Start Analysis
              <ArrowRight size={16} />
            </button>
          </div>

          <div
            className="action-card"
            onClick={() => navigate("/tools/drug-discovery")}
          >
            <div
              className="action-icon"
              style={{
                background: "linear-gradient(135deg, #a855f7 0%, #9333ea 100%)"
              }}
            >
              <Sparkles size={32} />
            </div>
            <h3>Drug Discovery</h3>
            <p>Generate novel drug candidates using generative AI models</p>
            <button className="action-button">
              Generate Molecules
              <ArrowRight size={16} />
            </button>
          </div>
        </div>
      </div>

      {/* Recent Projects */}
      <div className="recent-projects-section">
        <div className="section-header">
          <h2>Recent Projects</h2>
          <button
            className="view-all-btn"
            onClick={() => navigate("/projects")}
          >
            View All
            <ArrowRight size={16} />
          </button>
        </div>

        {recentProjects.length > 0 ? (
          <div className="projects-list">
            {recentProjects.map((project) => (
              <div
                key={project.id}
                className="project-item"
                onClick={() => handleProjectClick(project)}
              >
                <div className="project-icon">
                  {project.type === "dti" ? (
                    <Activity size={20} />
                  ) : project.type === "drug-discovery" ? (
                    <Sparkles size={20} />
                  ) : (
                    <Folder size={20} />
                  )}
                </div>
                <div className="project-info">
                  <h4>{project.name}</h4>
                  <p>{project.description || "No description"}</p>
                  <span className="project-date">
                    Last modified:{" "}
                    {format(new Date(project.lastModified), "MMM dd, yyyy")}
                  </span>
                </div>
                <div className="project-stats">
                  <span className="result-count">
                    {project.results?.length || 0} results
                  </span>
                </div>
              </div>
            ))}
          </div>
        ) : (
          <div className="empty-state">
            <Folder size={48} />
            <p>No projects yet. Create your first project to get started!</p>
          </div>
        )}
      </div>
    </div>
  );
}

export default Home;
