import { useState } from "react";
import { ChevronDown, HelpCircle } from "lucide-react";
import "./FAQ.css";

function FAQ() {
  const [openIndex, setOpenIndex] = useState(null);

  const faqs = [
    {
      category: "General",
      questions: [
        {
          q: "What is DeepBind?",
          a: "DeepBind is an advanced molecular research platform powered by artificial intelligence. It combines Graph Neural Networks (GNN) for drug-target interaction predictions and Generative AI for novel drug candidate generation, accelerating the drug discovery process."
        },
        {
          q: "Who can use this platform?",
          a: "This platform is designed for researchers, pharmaceutical scientists, bioinformaticians, and students working in drug discovery, molecular biology, and computational chemistry fields."
        },
        {
          q: "Do I need programming knowledge?",
          a: "No programming knowledge is required! Our platform provides an intuitive interface where you simply input molecular data (SMILES strings and protein sequences) and get results instantly."
        }
      ]
    },
    {
      category: "DTI Prediction",
      questions: [
        {
          q: "What is DTI Prediction?",
          a: "Drug-Target Interaction (DTI) prediction analyzes the binding affinity between a drug molecule and a target protein using advanced Graph Neural Networks. It helps determine how strongly a drug candidate will interact with its biological target."
        },
        {
          q: "What is a SMILES string?",
          a: 'SMILES (Simplified Molecular Input Line Entry System) is a notation that describes the structure of chemical molecules using short ASCII strings. For example, "CCO" represents ethanol.'
        },
        {
          q: "How accurate are the predictions?",
          a: "Our DTI model is trained on the KIBA dataset with thousands of validated drug-protein interactions. The model achieves strong performance with mean squared error below 0.7 on test data. However, predictions should be validated with experimental methods."
        },
        {
          q: "What do the binding affinity scores mean?",
          a: "Binding affinity scores indicate the strength of interaction: Scores > 10 suggest strong binding (high potential), 8-10 indicate moderate binding, and < 8 suggest weak binding. Higher scores generally indicate better drug candidates."
        }
      ]
    },
    {
      category: "Drug Discovery",
      questions: [
        {
          q: "How does DeepBind work?",
          a: "Our system uses Generative Adversarial Networks (GANs) to create novel drug-like molecules. The generator learns patterns from known drugs and creates new molecular structures, while the discriminator ensures they have drug-like properties."
        },
        {
          q: "Are generated molecules safe?",
          a: "Generated molecules are computational predictions based on chemical similarity to known drugs. They require extensive validation, including toxicity testing, ADME studies, and clinical trials before any real-world application."
        },
        {
          q: "Can I download the generated molecules?",
          a: "Yes! All generated SMILES strings can be copied from the results. You can use them in other molecular modeling software or for further analysis."
        },
        {
          q: "What are MW and LogP values?",
          a: "MW (Molecular Weight) indicates the size of the molecule. LogP (Partition Coefficient) measures how lipophilic (fat-loving) a molecule is, affecting its ability to cross cell membranes and its bioavailability."
        }
      ]
    },
    {
      category: "Projects",
      questions: [
        {
          q: "How do I create a project?",
          a: 'Click the "New Project" button in the sidebar, enter a project name and description, select the project type (DTI Prediction, Drug Discovery, or Both), and click "Create Project". Your project will be saved automatically.'
        },
        {
          q: "Are my projects saved?",
          a: "Yes! All projects are saved locally in your browser using persistent storage. Your data remains private and is not sent to any server except for predictions through our API."
        },
        {
          q: "Can I export my project data?",
          a: "Currently, project data is stored locally. We recommend copying important SMILES strings and results for your records. Export functionality is coming in future updates."
        }
      ]
    },
    {
      category: "Technical",
      questions: [
        {
          q: "What models power this platform?",
          a: "We use Graph Convolutional Networks (GCN) for DTI prediction, ProtBERT for protein embeddings, and Generative Adversarial Networks (GAN) for drug generation. The DTI model is trained on KIBA dataset with 15,000+ interactions."
        },
        {
          q: "How long do predictions take?",
          a: "DTI predictions typically take 2-5 seconds. Drug generation takes 1-2 minutes for 10 candidates, depending on your hardware and protein sequence length."
        },
        {
          q: "Is my data secure?",
          a: "Yes! Your molecular data is processed securely. Project data is stored locally in your browser. Only the molecular structures needed for prediction are sent to our secure backend API."
        },
        {
          q: "Can I use this for commercial purposes?",
          a: "Please review our terms of service. For commercial applications, contact us for licensing information. Academic and research use is generally permitted with proper attribution."
        }
      ]
    }
  ];

  const toggleFAQ = (index) => {
    setOpenIndex(openIndex === index ? null : index);
  };

  let globalIndex = 0;

  return (
    <div className="faq-container">
      <div className="faq-header">
        <HelpCircle size={48} />
        <h1>Frequently Asked Questions</h1>
        <p>Find answers to common questions about DeepBind</p>
      </div>

      <div className="faq-content">
        {faqs.map((category, catIdx) => (
          <div key={catIdx} className="faq-category">
            <h2 className="category-title">{category.category}</h2>
            <div className="faq-list">
              {category.questions.map((faq) => {
                const currentIndex = globalIndex++;
                const isOpen = openIndex === currentIndex;

                return (
                  <div
                    key={currentIndex}
                    className={`faq-item ${isOpen ? "open" : ""}`}
                  >
                    <button
                      className="faq-question"
                      onClick={() => toggleFAQ(currentIndex)}
                    >
                      <span>{faq.q}</span>
                      <ChevronDown
                        size={20}
                        className={`chevron ${isOpen ? "rotate" : ""}`}
                      />
                    </button>
                    {isOpen && (
                      <div className="faq-answer">
                        <p>{faq.a}</p>
                      </div>
                    )}
                  </div>
                );
              })}
            </div>
          </div>
        ))}
      </div>

      <div className="faq-footer">
        <h3>Still have questions?</h3>
        <p>Can't find what you're looking for? Contact our support team.</p>
        <button
          className="btn btn-primary"
          onClick={() => (window.location.href = "/contact")}
        >
          Contact Support
        </button>
      </div>
    </div>
  );
}

export default FAQ;
