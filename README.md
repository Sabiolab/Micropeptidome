text
           _____           _   _     _                      
 _   _   |  __ \         | | (_)   | |                     
| | | |  | |__) |__ _ __ | |_ _  __| | ___  _ __ ___   ___ 
| | | |  |  ___/ _ \ '_ \| __| |/ _` |/ _ \| '_ ` _ \ / _ \
| |_| |  | |  |  __/ |_) | |_| | (_| | (_) | | | | | |  __/
| ___/   |_|   \___| .__/ \__|_|\__,_|\___/|_| |_| |_|\___|
| |                | |                                       
|_|                |_|                  
"""


## What is Micropeptidome?

**Micropeptidome** is a framework for identifying smORF-encoded microproteins (<150aa) from both proteomic ans transcriptomic experiments. It distinguishes between:

- **SAMs** (Swiss-Prot Analog Microproteins): smORFs that resemble microproteins in Swiss-Prot, the reviewed  section of the UniProt Knowledgebase.
- **PRISMs** (Physicochemically Resembling In Silico Microproteins): smORFs that resemble synthetic sequences.

ShortStop also supports the generation of matched negative control microprotein sequences (i.e., PRISMS) for downstream applications.

---

## Why Use ShortStop?

Thousands of smORFs are actively translated in the human transcriptome, yet the functions of their encoded microproteins remain unknown. 

Protein function is conventionally inferred from evidence of evolutionary selection. However, most encoded microproteins are evolutionarily young and lack detectable sequence conservation across species. 

ShortStop provides a homology-independent approach to identify microproteins that share key physicochemical properties with known proteins, even in the absence of sequence similarity. 

ShortStop classifications supports hypothesis-driven studies to "deorphanize" the function of translating smORFs. It also provides a foundation for future benchmarking of new algorithms and methods as data emerges from ongoing studies.

---

## Key Features

- Classifies smORFs as SAMs or PRISMs using a pre-trained ML model
- Generates composition-matched negative control microprotein sequences
- Extracts physicochemical features from input smORFs
- Supports model training, prediction, and feature extraction as CLI modules

---

## Requirements

You’ll need:

1. A GTF file of smORFs that must contain CDS and transcripts features
2. A matched reference genome (e.g., hg38, which automatically downloads upon initiating demo mode).

---

## Installation

> ✅ We recommend the creation of a conda environment:
> ```bash
> conda create -n micropeptidome python=3.9
> conda activate micropeptidome
> ```

### Option 1 – Direct from GitHub (recommended)
```bash
pip install git+https://github.com/Sabiolab/Micropeptidome.git
```

### Option 2 – Clone and Install Locally
```bash
git clone https://github.com/Sabiolab/Micropeptidome.git
cd Micropeptidome
pip install .
```

### ⚠️ If you get a C compilation error during install...
Install a C compiler for your system:

- **Ubuntu/Debian**
  ```bash
  sudo apt-get install build-essential
  ```

- **Fedora/CentOS**
  ```bash
  sudo dnf install gcc
  ```

- **Arch Linux**
  ```bash
  sudo pacman -S base-devel
  ```

- **Windows**  
  Download and install: [Microsoft C++ Build Tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/)



## License and Contributions

This project is licensed for **non-commercial academic research use only**.  
See [LICENSE.md](./LICENSE.md) for full terms.

By contributing to this repository, you agree to the [Contributor License Agreement (CLA)](./CLA.md).  

Please read our [CONTRIBUTING.md](./CONTRIBUTING.md) before submitting code or issues.

By downloading or using this tool, you agree to the terms in LICENSE.md and CLA.md.
