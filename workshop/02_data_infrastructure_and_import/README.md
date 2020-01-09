- The components of a single cell experiment
  - counts
  - sample metadata
  - feature metadata
- Constructing and importing elements of a single-cell experiment (counts, sample metadata, feature metadata)
  - Simplified example 
    - Before workshop: Take an subsamples SCE and write each element to disk. Provide students with files.
    - During workshop: Read files back in to construct the SCE.
  - Realistic example
    - During workshop: A more realistic example, e.g., using **DropletUtils** to load a 10X dataset.
  - Working with SCORE, you'll receive an SCE.
- The *SingleCellExperiment* class
  - `show`
    - Explain what is displayed when you print an SCE
  - The starting point
    - `assays`
      - `counts()`, `logcounts()`
    - `colData`
      - Subsetting
    - `rowData`
      - Subsetting
  - The add ons
    - `reducedDims()`
    - `spikeNames()`(?)
    - `altExps()`
    - `sizeFactors()`(?)
    
# TODOs

- [ ] Should I show an example SCE before constructing one?
