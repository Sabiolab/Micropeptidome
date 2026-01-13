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

## Do you have data from transcriptomic experiments?

This workflow processes **RNA-seq** data to identify small open reading frames (smORFs) encoding putative microproteins (10-150 amino acids).
The pipeline performs transcriptome assembly, ORF prediction, size-based filtering, and coordinate transformation to generate genome-based annotations compatible with **ShortStop**. 

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



## License and Contributions

This project is licensed for **non-commercial academic research use only**.  
See [LICENSE.md](./LICENSE.md) for full terms.

By contributing to this repository, you agree to the [Contributor License Agreement (CLA)](./CLA.md).  


By downloading or using this tool, you agree to the terms in LICENSE.md and CLA.md.
