const { Parser, Drawer } = require('smiles-drawer');

function clean2d(smiles) {
    let parsed = Parser.parse(smiles)
    let drawer = new Drawer({});
    drawer.initDraw(parsed, 'light', true);
    drawer.processGraph();

    let vertices = drawer.graph.vertices;
    let xy = Array();
    for (let i = 0; i < vertices.length; i++) {
      let position = vertices[i].position;
      xy.push([position.x, position.y]);
    }
    return xy;
}

export { clean2d };
