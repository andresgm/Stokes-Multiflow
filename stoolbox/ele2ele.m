% determina los elementos vecinos a cada elemento
function ele2elecell = ele2ele(elements)

numelements = size(elements,1);
ele2elecell = cell(numelements,1);
elelist = repmat([1:1:numelements]',[1 3]);

for i=1:numelements
    nodevertex = elements(i,:);
      
    indl1 = elelist(elements(:) == nodevertex(1));
    indl2 = elelist(elements(:) == nodevertex(2));
    indl3 = elelist(elements(:) == nodevertex(3));

    ind = [indl1;indl2;indl3];
    
    % elimine el elemento i
    ind(ind == i) = [];
    ind = sort(unique(ind));
    % elimine el 
    ele2elecell{i} = ind;
end
