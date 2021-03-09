const initialState = {
  xParam: "UMAP_1",
  yParam: "UMAP_2",
  clonotypeParam: "clonotype",
  cellIdParam: "cell_id",
  subtypeParam: "subtype"
};

const dashboardReducer = (state, action) => {
  console.log(action);
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
