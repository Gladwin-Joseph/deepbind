import {
  BrowserRouter as Router,
  Routes,
  Route,
  Navigate
} from "react-router-dom";
import { useAuth } from "@clerk/clerk-react";
import Navbar from "./components/Navbar";
import Sidebar from "./components/Sidebar";
import AIChatbot from "./components/AiChatbot";
import ProtectedRoute from "./components/ProtectedRoute";
import Login from "./pages/Login";
import Home from "./pages/Home";
import DTIPrediction from "./pages/DTIPrediction";
import DrugDiscovery from "./pages/DrugDiscovery";
import Projects from "./pages/Projects";
import FAQ from "./pages/FAQ";
import Contact from "./pages/Contact";
import Documentation from "./pages/Documentation";
import "./App.css";

function AppLayout({ children }) {
  return (
    <div className="app-layout">
      <Navbar />
      <div className="app-main">
        <Sidebar />
        <div className="app-content">{children}</div>
      </div>
      {/* AI Chatbot - Always visible when logged in */}
      <AIChatbot />
    </div>
  );
}

function App() {
  const { isSignedIn, isLoaded } = useAuth();

  if (!isLoaded) {
    return (
      <div className="loading-screen">
        <div className="spinner"></div>
      </div>
    );
  }

  return (
    <Router>
      <Routes>
        {/* Public Route */}
        <Route path="/login" element={<Login />} />

        {/* Protected Routes */}
        <Route
          path="/"
          element={
            <ProtectedRoute>
              <AppLayout>
                <Home />
              </AppLayout>
            </ProtectedRoute>
          }
        />

        <Route
          path="/tools/dti-prediction"
          element={
            <ProtectedRoute>
              <AppLayout>
                <DTIPrediction />
              </AppLayout>
            </ProtectedRoute>
          }
        />

        <Route
          path="/tools/drug-discovery"
          element={
            <ProtectedRoute>
              <AppLayout>
                <DrugDiscovery />
              </AppLayout>
            </ProtectedRoute>
          }
        />

        <Route
          path="/projects"
          element={
            <ProtectedRoute>
              <AppLayout>
                <Projects />
              </AppLayout>
            </ProtectedRoute>
          }
        />

        <Route
          path="/documentation"
          element={
            <ProtectedRoute>
              <AppLayout>
                <Documentation />
              </AppLayout>
            </ProtectedRoute>
          }
        />

        <Route
          path="/faq"
          element={
            <ProtectedRoute>
              <AppLayout>
                <FAQ />
              </AppLayout>
            </ProtectedRoute>
          }
        />

        <Route
          path="/contact"
          element={
            <ProtectedRoute>
              <AppLayout>
                <Contact />
              </AppLayout>
            </ProtectedRoute>
          }
        />

        {/* Redirect */}
        <Route path="*" element={<Navigate to="/" replace />} />
      </Routes>
    </Router>
  );
}

export default App;
