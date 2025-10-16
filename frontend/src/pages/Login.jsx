import { SignIn, SignUp, useAuth } from "@clerk/clerk-react";
import { useState } from "react";
import { Navigate } from "react-router-dom";
import { Atom, Sparkles, Activity, Beaker } from "lucide-react";
import "./Login.css";

function Login() {
  const { isSignedIn } = useAuth();
  const [mode, setMode] = useState("sign-in"); // 'sign-in' or 'sign-up'

  if (isSignedIn) {
    return <Navigate to="/" replace />;
  }

  return (
    <div className="login-container">
      <div className="login-left">
        <div className="login-brand">
          <Atom size={48} />
          <h1>DeepBind</h1>
          <p>
            Advanced molecular research platform powered by artificial
            intelligence
          </p>
        </div>

        <div className="login-features">
          <div className="feature-item">
            <Activity size={24} />
            <div>
              <h3>DTI Prediction</h3>
              <p>
                Predict drug-target interactions with state-of-the-art GNN
                models
              </p>
            </div>
          </div>

          <div className="feature-item">
            <Sparkles size={24} />
            <div>
              <h3>AI Drug Generation</h3>
              <p>Generate novel drug candidates using generative AI</p>
            </div>
          </div>

          <div className="feature-item">
            <Beaker size={24} />
            <div>
              <h3>Project Management</h3>
              <p>Organize and track your research projects efficiently</p>
            </div>
          </div>
        </div>

        <div className="login-footer">
          <p>© 2025 DeepBind. All rights reserved.</p>
        </div>
      </div>

      <div className="login-right">
        <div className="clerk-wrapper">
          {/* Toggle buttons */}
          <div className="auth-toggle">
            <button
              className={`toggle-btn ${mode === "sign-in" ? "active" : ""}`}
              onClick={() => setMode("sign-in")}
            >
              Sign In
            </button>
            <button
              className={`toggle-btn ${mode === "sign-up" ? "active" : ""}`}
              onClick={() => setMode("sign-up")}
            >
              Sign Up
            </button>
          </div>

          {mode === "sign-in" ? (
            <SignIn
              appearance={{
                elements: {
                  rootBox: {
                    width: "100%"
                  },
                  card: {
                    backgroundColor: "#1e293b",
                    border: "1px solid #334155",
                    boxShadow: "0 25px 50px -12px rgba(0, 0, 0, 0.5)",
                    borderRadius: "16px",
                    padding: "2rem"
                  },
                  headerTitle: {
                    color: "#ffffff",
                    fontSize: "1.5rem",
                    fontWeight: "700"
                  },
                  headerSubtitle: {
                    color: "#94a3b8"
                  },
                  socialButtonsBlockButton: {
                    backgroundColor: "#334155",
                    border: "1px solid #475569",
                    color: "#ffffff",
                    borderRadius: "8px",
                    padding: "0.75rem"
                  },
                  dividerLine: {
                    backgroundColor: "#334155"
                  },
                  dividerText: {
                    color: "#94a3b8"
                  },
                  formFieldLabel: {
                    color: "#e2e8f0"
                  },
                  formFieldInput: {
                    backgroundColor: "#0f172a",
                    border: "1px solid #334155",
                    color: "#ffffff"
                  },
                  formButtonPrimary: {
                    backgroundColor: "#6366f1",
                    color: "#ffffff",
                    "&:hover": {
                      backgroundColor: "#4f46e5"
                    }
                  },
                  footerActionLink: {
                    color: "#6366f1"
                  }
                }
              }}
              redirectUrl="/"
              afterSignInUrl="/"
              signUpUrl="#"
              // 👇 this line makes Clerk show no redirect link
              routing="virtual"
              // 👇 we use this callback to handle link clicks manually
              onSignUpClick={(e) => {
                e.preventDefault();
                setMode("sign-up");
              }}
            />
          ) : (
            <SignUp
              appearance={{
                elements: {
                  rootBox: {
                    width: "100%",
                    display: "flex",
                    justifyContent: "center"
                  },
                  card: {
                    backgroundColor: "#1e293b",
                    border: "1px solid #334155",
                    boxShadow: "0 25px 50px -12px rgba(0, 0, 0, 0.5)",
                    borderRadius: "16px",
                    padding: "2rem"
                  },
                  headerTitle: {
                    color: "#ffffff",
                    fontSize: "1.5rem",
                    fontWeight: "700"
                  },
                  headerSubtitle: {
                    color: "#94a3b8"
                  },
                  socialButtonsBlockButton: {
                    backgroundColor: "#334155",
                    border: "1px solid #475569",
                    color: "#ffffff",
                    borderRadius: "8px",
                    padding: "0.75rem"
                  },
                  dividerLine: {
                    backgroundColor: "#334155"
                  },
                  dividerText: {
                    color: "#94a3b8"
                  },
                  formFieldLabel: {
                    color: "#e2e8f0"
                  },
                  formFieldInput: {
                    backgroundColor: "#0f172a",
                    border: "1px solid #334155",
                    color: "#ffffff"
                  },
                  formButtonPrimary: {
                    backgroundColor: "#6366f1",
                    color: "#ffffff",
                    "&:hover": {
                      backgroundColor: "#4f46e5"
                    }
                  },
                  footerActionLink: {
                    color: "#6366f1"
                  }
                }
              }}
              redirectUrl="/"
              afterSignUpUrl="/"
              routing="virtual"
              // Custom handler for "Already have an account? Sign in"
              onSignInClick={(e) => {
                e.preventDefault();
                setMode("sign-in");
              }}
            />
          )}
        </div>
      </div>
    </div>
  );
}

export default Login;
