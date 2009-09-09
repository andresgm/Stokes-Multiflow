% calcula los vectores unitarios de un arreglo de vectores fila vector
function uvect = univect(vector)

uvect = vector./repmat(normesp(vector),[1 3]);
