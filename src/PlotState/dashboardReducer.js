const initialState = {
  xParam: "UMAP_1",
  yParam: "UMAP_2",
  clonotypeParam: "cdr3s_aa",
  cellIdParam: "cell_id"
};

const dashboardReducer = (state, action) => {
  switch (action.type) {
    case "OVERRIDE": {
      return {
        ...state
      };
    }
    default:
      return state;
  }
};
export { initialState };
export default dashboardReducer;
