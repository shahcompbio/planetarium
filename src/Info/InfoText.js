const infoText = {
  UMAP: {
    title: "Cdr3 Amino Acid Sequence UMAP",
    text:
      "UMAP embedding generated from <br>the subset of CD8/CD4 T cells<br> identified with CellAssign. <br>The top 10 TRB sequences by <br>frequency are mapped to cells <br>in the embedding. The radius slider <br>allows the user to highlight areas <br>of high density for a given sequence <br>that can be referenced against<br> the subtype UMAP."
  },
  SUBTYPEUMAP: {
    title: "Subtype UMAP",
    text:
      "UMAP embedding generated from the subset <br>of CD8/CD4 T-cells identified with CellAssign. <br>Subtypes are assigned from differential gene <br>expression derived from Leiden clustering. <br>The list of current subtypes and known gene<br>markers is provided below."
  },
  HEATMAP: {
    title: "Subtype to Sequence Heatmap",
    text:
      "Frequency of the top 10 TRB sequences is shown <br>for each of the identified subtypes. The heatmap <br>intensity captures the relative proportion of the <br>given sequence with respect to total number of<br> TCR+ sequences within each subtype (Frequency<br> of Sequence in Subtype / Total TCR+ Cells in Subtype)."
  },
  HISTOGRAM: {
    title: "Generation Probabilities",
    text:
      "Distribution of log10 generation probabilities <br>generated from OLGA for all TRB sequences.<br>Description of the OLGA method can be found here: <br>https://academic.oup.com/bioinformatics/article/35/17/2974/5292315/. <br>The log10 probability of each top 10 TRB sequence <br>can be referenced with respect to the entire distribution<br> of the sample. Sequences with a lower probability of <br>generation will be found in the left tail of the distribution."
  },
  BARPLOT: {
    title: "Clonotype Expansion",
    text:
      "Visualization of clonal expansion in each subtype<br> by plotting the fraction of cells that belong to sequences <br>expanded to N cells. The fraction of cells with sequences<br> that belong to a single cell (labeled 1 in legend) describes<br> the relative proportion of unexpanded clonotypes. <br>The fraction of cells with sequences that belong <br>to >10 cells describes the relative proportion of highly<br> expanded clonotypes within the subtype."
  },
  TABLE: {
    title: "Differentially Expressed Genes",
    text:
      "Differentially expressed genes for each subtype (1 vs. All)<br> using the Wilcoxon test. The top 50 statistically significant <br>(p < 0.001) with a minimum log fold change of 0.25 and <br>expression found in >= 50% of the subtype cells are included<br> for each subtype. It is possible to sort DEGs by both<br> p-value and log fold change."
  }
};
export default infoText;
