const initialState = {
  xParam: "UMAP_1",
  yParam: "UMAP_2",
  clonotypeParam: "cdr3s_aa",
  cellIdParam: "cell_id",
  subtypeParam: "subtype",
  logXParam: "log10_probability",
  logYParam: "logProbability"
};

const dashboardReducer = (state, action) => {
  switch (action.type) {
    case "OVERRIDE": {
      return {
        ...state,
        ...action.value
      };
    }
    default:
      return state;
  }
};
export { initialState };
export default dashboardReducer;
