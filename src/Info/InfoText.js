const infoText = {
  UMAP: {
    title: "Cdr3 Amino Acid Sequence UMAP",
    text:
      "UMAP embedding generated from <br>the subset of CD8/CD4 T cells<br> identified with CellAssign. <br>The top 10 TRB sequences by <br>frequency are mapped to cells <br>in the embedding. The radius slider <br>allows the user to highlight areas <br>of high density for a given sequence <br>that can be referenced against<br> the subtype UMAP."
  },
  SUBTYPEUMAP: {
    title: "Subtype UMAP",
    text:
      "UMAP embedding generated <br>from the subset of CD8/CD4 T-cells<br> identified with CellAssign. Subtypes are assigned from <br>differential gene expression <br>derived from Leiden clustering. <br>The list of current subtypes and known gene markers is<br> provided below."
  },
  HEATMAP: {
    title: "Subtype to Sequence Heatmap",
    text:
      "Frequency of the top 10 TRB <br>sequences is shown for each of the <br>identified subtypes. <br>The heatmap intensity captures <br>the relative proportion of the <br>given sequence with respect to <br>total number of TCR+ sequences within each subtype <br>(Frequency of Sequence in Subtype / Total <br>TCR+ Cells in Subtype)."
  },
  BARPLOT: {
    title: "Generation Probabilities",
    text:
      "Distribution of log10 generation probabilities generated from OLGA for all TRB sequences. Description of the OLGA method can be found here: https://academic.oup.com/bioinformatics/article/35/17/2974/5292315/. The log10 probability of each top 10 TRB sequence can be referenced with respect to the entire distribution of the sample. Sequences with a lower probability of generation will be found in the left tail of the distribution."
  },
  HISTOGRAM: {
    title: "Clonotype Expansion",
    text:
      "Visualization of clonal expansion in each subtype by plotting the fraction of cells that belong to sequences expanded to N cells. The fraction of cells with sequences that belong to a single cell (labeled 1 in legend) describes the relative proportion of unexpanded clonotypes. The fraction of cells with sequences that belong to >10 cells describes the relative proportion of highly expanded clonotypes within the subtype."
  }
};
export default infoText;
