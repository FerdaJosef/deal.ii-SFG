

<html lang="en">
<head>
<meta charset="UTF-8">
<style>
    @page {
        size: A4;
        margin: 20mm;
        background-color: #fcfcfc;
    }
    body {
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        color: #333;
        line-height: 1.6;
        font-size: 11pt;
        background-color: #fcfcfc;
    }
    .header-banner {
        background-color: #2c3e50;
        color: white;
        padding: 20px;
        margin: -20mm -20mm 20mm -20mm;
        text-align: center;
    }
    h1 { margin: 0; font-size: 22pt; }
    h2 { 
        color: #2980b9; 
        border-left: 5px solid #2980b9; 
        padding-left: 10px; 
        margin-top: 30px;
        font-size: 16pt;
    }
    h3 { color: #34495e; font-size: 13pt; margin-top: 20px; }
    code {
        background-color: #f0f0f0;
        padding: 2px 4px;
        border-radius: 4px;
        font-family: 'Courier New', Courier, monospace;
        font-size: 10pt;
    }
    pre {
        background-color: #2d2d2d;
        color: #ccc;
        padding: 15px;
        border-radius: 6px;
        overflow-x: auto;
        font-size: 9.5pt;
        line-height: 1.4;
    }
    .step-box {
        background-color: #edf2f7;
        border: 1px solid #cbd5e0;
        padding: 15px;
        margin: 15px 0;
        border-radius: 8px;
    }
    .tip {
        border-left: 4px solid #38a169;
        background-color: #f0fff4;
        padding: 10px 15px;
        margin: 20px 0;
    }
    .math {
        font-family: 'Times New Roman', serif;
        font-style: italic;
        font-weight: bold;
        color: #2c3e50;
    }
</style>
</head>

    <div class="header-banner">
        <h1>deal.II-SFG Simulation Framework</h1>
        <p>Technical Documentation & Implementation Guide</p>
    </div>

    <h2>1. Introduction</h2>
    <p>This document provides the standard operating procedure for setting up, configuring, and running stochastic phase-field simulations using the deal.II-SFG framework.</p>

    <h2>2. Setup & Execution Pipeline</h2>

    <div class="step-box">
        <h3>Step I: Repository Initialization</h3>
        <p>Clone the repository and enter the project root.</p>
        <pre>git clone &lt;repository-url&gt;
cd deal.ii-SFG</pre>
    </div>

    <div class="step-box">
        <h3>Step II: Model Customization</h3>
        <p>The model logic is decentralized into specific source files to optimize compilation times.</p>
        <ul>
            <li><strong>Initial Conditions:</strong> Modify <code>source/Model/InitialValues.cc</code>.</li>
            <li><strong>Right-Hand Side:</strong> Define source terms and noise in <code>source/Model/RandomField.cc</code>.</li>
            <li><strong>Templates:</strong> You MUST append the explicit template instantiation at the end of each <code>.cc</code> file:
                <pre>template class RandomField&lt;1, 2&gt;;</pre>
            </li>
        </ul>
    </div>

    <div class="step-box">
        <h3>Step III: Symbolic Differentiation (AceGen)</h3>
        <p>The framework utilizes AceGen for automated residual and Jacobian generation.</p>
        <ol>
            <li>Place AceGen C-output in <code>source/Model/AceGen/</code>.</li>
            <li>Run the bridge script:
                <pre>python3 script.py source/Model/AceGen/output.c</pre>
            </li>
        </ol>
    </div>

    <div class="step-box">
        <h3>Step IV: Configuration</h3>
        <p>Ensure <code>main.cc</code> includes the correct model header and <code>CMakeLists.txt</code> points to the relevant source files.</p>
    </div>

    <div class="step-box">
        <h3>Step V: Build & Run</h3>
        <pre>mkdir build && cd build
cmake ..
make -j$(nproc)
make run</pre>
    </div>

    <div class="tip">
        <strong>Technical Note:</strong> The separation of <code>RandomField.cc</code> ensures that changes to the core FEM algorithm in <code>model.cc</code> do not require re-generating the stochastic noise logic, significantly reducing iterative development time.
    </div>



