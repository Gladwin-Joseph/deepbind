import { create } from "zustand";
import { persist } from "zustand/middleware";
import { v4 as uuidv4 } from "uuid";

const useStore = create(
  persist(
    (set, get) => ({
      // Projects
      projects: [],
      currentProjectId: null,

      // Last view state for each project (stores form inputs and results)
      projectViewStates: {}, // { projectId: { dti: {...}, drugDiscovery: {...} } }

      // Statistics
      stats: {
        totalDTIPredictions: 0,
        totalDrugsGenerated: 0,
        activeProjects: 0
      },

      // Set current project
      setCurrentProject: (projectId) => {
        set({ currentProjectId: projectId });
      },

      // Save project view state (called whenever inputs/results change)
      saveProjectViewState: (projectId, viewState) => {
        set((state) => ({
          projectViewStates: {
            ...state.projectViewStates,
            [projectId]: {
              ...state.projectViewStates[projectId],
              ...viewState,
              lastViewed: new Date().toISOString()
            }
          }
        }));
      },

      // Get project view state
      getProjectViewState: (projectId) => {
        return get().projectViewStates[projectId] || null;
      },

      // Create new project
      createProject: (projectData) => {
        const newProject = {
          id: uuidv4(),
          ...projectData,
          createdAt: new Date().toISOString(),
          lastModified: new Date().toISOString(),
          results: []
        };

        set((state) => ({
          projects: [newProject, ...state.projects],
          stats: {
            ...state.stats,
            activeProjects: state.projects.length + 1
          }
        }));

        return newProject.id;
      },

      // Add result to project
      addResultToProject: (projectId, result) => {
        set((state) => {
          const updatedProjects = state.projects.map((project) => {
            if (project.id === projectId) {
              return {
                ...project,
                results: [
                  ...project.results,
                  { ...result, timestamp: new Date().toISOString() }
                ],
                lastModified: new Date().toISOString()
              };
            }
            return project;
          });

          // Update stats
          const newStats = { ...state.stats };
          if (result.type === "dti") {
            newStats.totalDTIPredictions += 1;
          } else if (result.type === "drug-discovery") {
            newStats.totalDrugsGenerated += result.candidatesCount || 0;
          }

          return {
            projects: updatedProjects,
            stats: newStats
          };
        });
      },

      // Get project by ID
      getProject: (projectId) => {
        return get().projects.find((p) => p.id === projectId);
      },

      // Get projects by type
      getProjectsByType: (type) => {
        return get().projects.filter((p) => p.type === type);
      },

      // Delete project
      deleteProject: (projectId) => {
        set((state) => {
          // Remove project view state
          const newViewStates = { ...state.projectViewStates };
          delete newViewStates[projectId];

          return {
            projects: state.projects.filter((p) => p.id !== projectId),
            projectViewStates: newViewStates,
            currentProjectId:
              state.currentProjectId === projectId
                ? null
                : state.currentProjectId,
            stats: {
              ...state.stats,
              activeProjects: Math.max(0, state.stats.activeProjects - 1)
            }
          };
        });
      },

      // Update project
      updateProject: (projectId, updates) => {
        set((state) => ({
          projects: state.projects.map((p) =>
            p.id === projectId
              ? { ...p, ...updates, lastModified: new Date().toISOString() }
              : p
          )
        }));
      },

      // Clear all data (for testing)
      clearAll: () => {
        set({
          projects: [],
          currentProjectId: null,
          projectViewStates: {},
          stats: {
            totalDTIPredictions: 0,
            totalDrugsGenerated: 0,
            activeProjects: 0
          }
        });
      }
    }),
    {
      name: "drug-discovery-storage"
    }
  )
);

export default useStore;
