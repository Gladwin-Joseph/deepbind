import React from "react";
import ReactDOM from "react-dom/client";
import { ClerkProvider } from "@clerk/clerk-react";
import { Toaster } from "react-hot-toast";
import App from "./App.jsx";
import "./index.css";

const PUBLISHABLE_KEY = import.meta.env.VITE_CLERK_PUBLISHABLE_KEY;

if (!PUBLISHABLE_KEY) {
  throw new Error("Missing Clerk Publishable Key");
}

// Dark theme customization for Clerk
const clerkAppearance = {
  baseTheme: undefined,
  variables: {
    colorPrimary: "#6366f1",
    colorBackground: "#0f172a",
    colorInputBackground: "#1e293b",
    colorInputText: "#ffffff",
    colorText: "#ffffff",
    colorTextSecondary: "#94a3b8",
    colorDanger: "#ef4444",
    colorSuccess: "#10b981",
    colorWarning: "#f59e0b",
    borderRadius: "0.5rem",
    fontFamily:
      '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif'
  },
  elements: {
    // Root container
    rootBox: {
      backgroundColor: "#0f172a"
    },
    // Card containers
    card: {
      backgroundColor: "#1e293b",
      border: "1px solid #334155",
      boxShadow: "0 25px 50px -12px rgba(0, 0, 0, 0.5)"
    },
    // Header
    headerTitle: {
      color: "#ffffff"
    },
    headerSubtitle: {
      color: "#94a3b8"
    },
    // Navigation
    navbar: {
      backgroundColor: "#1e293b",
      borderBottom: "1px solid #334155"
    },
    navbarButton: {
      color: "#94a3b8"
    },
    navbarButtonActive: {
      color: "#6366f1"
    },
    // Page wrapper
    page: {
      backgroundColor: "#0f172a"
    },
    pageScrollBox: {
      backgroundColor: "#0f172a"
    },
    // Form elements
    formFieldLabel: {
      color: "#e2e8f0"
    },
    formFieldInput: {
      backgroundColor: "#0f172a",
      border: "1px solid #334155",
      color: "#ffffff"
    },
    formFieldInput__focused: {
      borderColor: "#6366f1"
    },
    formButtonPrimary: {
      backgroundColor: "#6366f1",
      color: "#ffffff"
    },
    formButtonPrimary__hover: {
      backgroundColor: "#4f46e5"
    },
    // Profile section
    profileSection: {
      backgroundColor: "#1e293b",
      border: "1px solid #334155",
      borderRadius: "0.75rem"
    },
    profileSectionTitle: {
      color: "#ffffff"
    },
    profileSectionContent: {
      color: "#94a3b8"
    },
    // Badges and tags
    badge: {
      backgroundColor: "rgba(99, 102, 241, 0.1)",
      color: "#6366f1",
      border: "1px solid rgba(99, 102, 241, 0.3)"
    },
    // Buttons
    button: {
      color: "#94a3b8",
      backgroundColor: "transparent",
      border: "1px solid #334155"
    },
    buttonPrimary: {
      backgroundColor: "#6366f1",
      color: "#ffffff"
    },
    buttonPrimary__hover: {
      backgroundColor: "#4f46e5"
    },
    buttonSecondary: {
      backgroundColor: "transparent",
      color: "#94a3b8",
      border: "1px solid #334155"
    },
    buttonSecondary__hover: {
      backgroundColor: "#334155",
      color: "#ffffff"
    },
    // Avatar
    avatarBox: {
      border: "2px solid #6366f1"
    },
    // Social buttons
    socialButtonsBlockButton: {
      backgroundColor: "#334155",
      border: "1px solid #475569",
      color: "#ffffff"
    },
    socialButtonsBlockButton__hover: {
      backgroundColor: "#475569"
    },
    // Dividers
    dividerLine: {
      backgroundColor: "#334155"
    },
    dividerText: {
      color: "#94a3b8"
    },
    // Footer
    footerActionLink: {
      color: "#6366f1"
    },
    footerActionLink__hover: {
      color: "#4f46e5"
    },
    // Alert/Error messages
    alert: {
      backgroundColor: "rgba(239, 68, 68, 0.1)",
      border: "1px solid rgba(239, 68, 68, 0.3)",
      color: "#fca5a5"
    },
    alertText: {
      color: "#fca5a5"
    },
    // Table elements (for connected accounts)
    table: {
      backgroundColor: "#1e293b",
      border: "1px solid #334155"
    },
    tableHead: {
      backgroundColor: "#0f172a",
      color: "#94a3b8"
    },
    tableBody: {
      color: "#e2e8f0"
    },
    tableRow: {
      borderBottom: "1px solid #334155"
    },
    // Menu items
    menuItem: {
      color: "#94a3b8"
    },
    menuItem__hover: {
      backgroundColor: "#334155",
      color: "#ffffff"
    },
    menuItemActive: {
      backgroundColor: "#6366f1",
      color: "#ffffff"
    },
    // Identity preview
    identityPreviewText: {
      color: "#94a3b8"
    },
    identityPreviewEditButton: {
      color: "#6366f1"
    },
    identityPreviewEditButton__hover: {
      color: "#4f46e5"
    },
    // OTP input
    otpCodeFieldInput: {
      backgroundColor: "#0f172a",
      border: "1px solid #334155",
      color: "#ffffff"
    },
    // Modal
    modalBackdrop: {
      backgroundColor: "rgba(0, 0, 0, 0.75)"
    },
    modalContent: {
      backgroundColor: "#1e293b",
      border: "1px solid #334155"
    }
  }
};

ReactDOM.createRoot(document.getElementById("root")).render(
  <React.StrictMode>
    <ClerkProvider
      publishableKey={PUBLISHABLE_KEY}
      appearance={clerkAppearance}
    >
      <App />
      <Toaster
        position="top-right"
        toastOptions={{
          duration: 3000,
          style: {
            background: "#1e293b",
            color: "#fff",
            border: "1px solid #334155"
          },
          success: {
            iconTheme: {
              primary: "#10b981",
              secondary: "#fff"
            }
          },
          error: {
            iconTheme: {
              primary: "#ef4444",
              secondary: "#fff"
            }
          }
        }}
      />
    </ClerkProvider>
  </React.StrictMode>
);
