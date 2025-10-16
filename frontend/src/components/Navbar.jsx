import { Link, useLocation } from "react-router-dom";
import { UserButton, useUser } from "@clerk/clerk-react";
import {
  Atom,
  Home,
  Folder,
  Bell,
  HelpCircle,
  FileText,
  Mail
} from "lucide-react";
import "./Navbar.css";

function Navbar() {
  const location = useLocation();
  const { user } = useUser();

  return (
    <nav className="navbar-dark">
      <div className="nav-container">
        <div className="nav-brand">
          <Link to="/" className="brand-link">
            <Atom size={32} />
            <span>DeepBind</span>
          </Link>
        </div>

        <div className="nav-center">
          <Link
            to="/"
            className={`nav-link ${location.pathname === "/" ? "active" : ""}`}
          >
            <Home size={20} />
            <span>Dashboard</span>
          </Link>

          <Link
            to="/projects"
            className={`nav-link ${
              location.pathname.startsWith("/projects") ? "active" : ""
            }`}
          >
            <Folder size={20} />
            <span>Projects</span>
          </Link>

          <Link
            to="/documentation"
            className={`nav-link ${
              location.pathname === "/documentation" ? "active" : ""
            }`}
          >
            <FileText size={20} />
            <span>Docs</span>
          </Link>

          <Link
            to="/faq"
            className={`nav-link ${
              location.pathname === "/faq" ? "active" : ""
            }`}
          >
            <HelpCircle size={20} />
            <span>FAQ</span>
          </Link>

          <Link
            to="/contact"
            className={`nav-link ${
              location.pathname === "/contact" ? "active" : ""
            }`}
          >
            <Mail size={20} />
            <span>Contact</span>
          </Link>
        </div>

        <div className="nav-right">
          <div className="user-section">
            <div className="user-info">
              <span className="user-name">
                {user?.firstName} {user?.lastName}
              </span>
              <span className="user-email">
                {user?.primaryEmailAddress?.emailAddress}
              </span>
            </div>
            <UserButton
              appearance={{
                elements: {
                  userButtonAvatarBox: "user-avatar-box",
                  userButtonPopoverCard: "user-popover-card",
                  userButtonPopoverActionButton: "user-popover-button",
                  userButtonPopoverActionButtonText: "user-popover-text",
                  userButtonPopoverActionButtonIcon: "user-popover-icon",
                  userButtonPopoverFooter: "user-popover-footer"
                }
              }}
            />
          </div>
        </div>
      </div>
    </nav>
  );
}

export default Navbar;
