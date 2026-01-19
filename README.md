```text
          ____            _   _     _                      
 _   _   |  __ \         | | (_)   | |                     
| | | |  | |__) |__ _ __ | |_ _  __| | ___  _ __ ___   ___ 
| | | |  |  ___/ _ \ '_ \| __| |/ _` |/ _ \| '_ ` _ \ / _ \
| |_| |  | |  |  __/ |_) | |_| | (_| | (_) | | | | | |  __/
| ___/   |_|   \___| .__/ \__|_|\__,_|\___/|_| |_| |_|\___|
| |                | |                                       
|_|                |_|                  

```

## What is Micropeptidome?

**Micropeptidome** is a framework for identifying microproteins (<150 aa) from both proteomic and transcriptomic experiments. It inludes several tools:

- **getefear**: transform your list (.csv) of microproteins in a .gtf doc which can be used to classify later with ShortStop.
- **ShortStop**: Classifies smORFs as SAMs or PRISMs using a pre-trained ML model ([click for detailed documentation](https://github.com/brendan-miller-salk/ShortStop/blob/master/README.md)). 

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
pip install git+https://github.com/Sabiolab/Micropeptidome/ShortStop.git
```

### Option 2 – Clone and Install Locally
```bash
git clone https://github.com/Sabiolab/Micropeptidome/ShortStop.git
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

---


## Folder structure

```text
Micropeptidome
├── README.md
│  
├── scripts/
│   ├── fastear.py
│   └── getefear.py
│   ├── cuentaSAM.py
│   ├── heterogenicitySAMs.py
│   ├── probabilidaes_SAMs.R
│   ├── Venn_SAMs.R
│   └── getefear.py
├── RNAseq/
│   ├── README_RNAseq.md
│   ├── De_novo_transcripts.py
│   ├── filetr_smorf_pep.py
│   └── smorfs_transcript_to_genome_gtf.py
└── ShortStop/  --> Clone from 'brendan-miller-salk/ShortStop'
    ├── README.md
    └── src/shortstop/
        └── shortStop.py
```


---

## License and Contributions

This project is licensed for **non-commercial academic research use only**.  
See [LICENSE.md](./LICENSE.md) for full terms.

By contributing to this repository, you agree to the [Contributor License Agreement (CLA)](./CLA.md).  


By downloading or using this tool, you agree to the terms in LICENSE.md and CLA.md.
