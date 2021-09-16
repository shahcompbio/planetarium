export default (value, highlighted, compareFn = (a, b) => a === b) =>
  highlighted === null || compareFn(value, highlighted);
