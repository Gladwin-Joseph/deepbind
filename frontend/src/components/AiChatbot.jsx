import { GoogleGenerativeAI } from "@google/generative-ai";
import { useState, useRef, useEffect } from "react";
import { MessageCircle, X, Send, Loader2, Sparkles } from "lucide-react";
import useStore from "../store/useStore";
import "./AIChatbot.css";

const apiKey = import.meta.env.VITE_GEMINI_API_KEY;

// ✅ Pass the API key when initializing
const genAI = new GoogleGenerativeAI(apiKey);

function AIChatbot() {
  const [isOpen, setIsOpen] = useState(false);
  const [messages, setMessages] = useState([
    {
      role: "assistant",
      content:
        "Hello! I'm your AI research assistant. I can help you with drug discovery, molecular biology, and answer questions about your projects. How can I assist you today?"
    }
  ]);
  const [input, setInput] = useState("");
  const [loading, setLoading] = useState(false);
  const messagesEndRef = useRef(null);
  const { currentProjectId, getProject } = useStore();

  const scrollToBottom = () => {
    messagesEndRef.current?.scrollIntoView({ behavior: "smooth" });
  };

  useEffect(() => {
    scrollToBottom();
  }, [messages]);

  const getContextPrompt = () => {
    let context =
      "You are an expert AI assistant specialized in drug discovery, molecular biology, bioinformatics, and pharmaceutical research. ";

    if (currentProjectId) {
      const project = getProject(currentProjectId);
      if (project) {
        context += `The user is currently working on a project named "${project.name}" (${project.type}). `;
        if (project.description) {
          context += `Project description: ${project.description}. `;
        }
        if (project.results && project.results.length > 0) {
          context += `The project has ${project.results.length} results. `;
        }
      }
    }

    context +=
      "Provide helpful, accurate, and concise answers about drug discovery, protein-drug interactions, molecular structures, SMILES notation, and related topics.";
    return context;
  };

  const handleSend = async () => {
    if (!input.trim() || loading) return;

    const userMessage = { role: "user", content: input };
    setMessages((prev) => [...prev, userMessage]);
    setInput("");
    setLoading(true);

    try {
      // ✅ Use the updated model name
      const model = genAI.getGenerativeModel({
        model: "gemini-2.5-flash" // or "gemini-1.5-pro" for more advanced responses
      });

      const context = getContextPrompt();
      const prompt = `${context}\n\nUser question: ${input}`;

      const result = await model.generateContent(prompt);
      const response = await result.response;
      const text = response.text();

      setMessages((prev) => [...prev, { role: "assistant", content: text }]);
    } catch (error) {
      console.error("Gemini API Error:", error);
      setMessages((prev) => [
        ...prev,
        {
          role: "assistant",
          content:
            "Sorry, I encountered an error. Please make sure your Gemini API key is configured correctly in the .env file."
        }
      ]);
    } finally {
      setLoading(false);
    }
  };

  const handleKeyPress = (e) => {
    if (e.key === "Enter" && !e.shiftKey) {
      e.preventDefault();
      handleSend();
    }
  };

  return (
    <>
      {/* Chatbot Toggle Button */}
      <button
        className="chatbot-toggle"
        onClick={() => setIsOpen(!isOpen)}
        aria-label="Toggle AI Assistant"
      >
        {isOpen ? <X size={24} /> : <MessageCircle size={24} />}
      </button>

      {/* Chatbot Window */}
      {isOpen && (
        <div className="chatbot-window">
          <div className="chatbot-header">
            <div className="chatbot-header-content">
              <Sparkles size={20} />
              <div>
                <h3>AI Research Assistant</h3>
              </div>
            </div>
            <button onClick={() => setIsOpen(false)} className="chatbot-close">
              <X size={20} />
            </button>
          </div>

          <div className="chatbot-messages">
            {messages.map((msg, idx) => (
              <div key={idx} className={`message ${msg.role}`}>
                <div className="message-content">{msg.content}</div>
              </div>
            ))}
            {loading && (
              <div className="message assistant">
                <div className="message-content">
                  <Loader2 size={16} className="spin" />
                  <span>Thinking...</span>
                </div>
              </div>
            )}
            <div ref={messagesEndRef} />
          </div>

          <div className="chatbot-input">
            <textarea
              value={input}
              onChange={(e) => setInput(e.target.value)}
              onKeyPress={handleKeyPress}
              placeholder="Ask me anything about drug discovery..."
              rows={1}
              disabled={loading}
            />
            <button
              onClick={handleSend}
              disabled={!input.trim() || loading}
              className="send-button"
            >
              <Send size={20} />
            </button>
          </div>
        </div>
      )}
    </>
  );
}

export default AIChatbot;
