export const isHighlighted = (
  highlightedColumn,
  highlightedRow,
  column,
  row
) => {
  if (highlightedColumn !== null || highlightedRow !== null) {
    if (highlightedColumn === column || highlightedRow === row) {
      return true;
    } else {
      return false;
    }
  } else {
    return true;
  }
};

export const isValueHighlighted = (value, highlighted) =>
  highlighted === null || value === highlighted;
