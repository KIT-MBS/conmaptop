# Plots and Data for Contact Map Analysis
Git repository with Data and Analysis for the publication *Influence of Contact Map Topology on RNA Structure Prediction*.

## Requirements
To execute all the scripts you need:
- pyDCA [^1]
- self written Helper library see [here](https://gitlab.jsc.fz-juelich.de/faber1/biohelpers)
- for effective sequence number calculation `sequeff`, see [here](https://gitlab.jsc.fz-juelich.de/faber1/sequeff)
- Testset from Pucci et al.[^2]

[^1]: Zerihun,M.B., Pucci,F., Peter,E.K. and Schug,A. (2020) pydca v1.0: a comprehensive software for direct coupling analysis of RNA and protein sequences. Bioinformatics, 36, 2264–2265.
[^2]: Pucci,F., Zerihun,M.B., Peter,E.K. and Schug,A. (2020) Evaluating DCA-based method performances for RNA contact prediction by a well-curated data set.


## Structure
### `FinalConfig`
Simulation data w/o any restraints

**Columns `results.csv`:**
- **RNA**: ID of the representative from https://www.rcsb.org/
- **Clustering File**: Clustered file, where `all` represents a REMC with 10 Replicas and `single` only one Replica, postfix `A` for an energy threshold of $0.01$ and `B` for $0.005$, i.e. only the $1\%$ or $0.05\%$ frames with the lowest energy are considered for clustering (for more details see xxx)
- **Threshold**: Distance threshold for the clustering in $\mathring{A}$
- **Cluster**: Cluster number
- **Configuration**: $1$ for Standard Config
- **size**: Number of residues

### `contact_distribution`
Comparison of different contact map topologies.

### `false_contacts`
Influence of false positives on the structure prediction

### `applications`
Application of the different scores to an additional validation set (included in the folder).

### `Plots`
Figures for the publication.

