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

1. BAM and FASTQ files. Preferrably paired end.
2. A matched reference genome (e.g., hg38, which automatically downloads upon initiating ShortStop demo mode)
3. 

---
