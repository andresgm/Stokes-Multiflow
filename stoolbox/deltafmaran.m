% calcula los esfuerzos de marangoni generales (tangenciales) y debido a
% tension supeficial variable
function [rdeltafmaran,rdeltafcurv,rsigmavar] = deltafmaran(struct,gamma,maran,rkcurv)


[rmaranstress,rsigmavar] = marangoni(struct,gamma,maran);
rdeltafmaran = rmaranstress.*maran.rkmaran;
rdeltafcurv = deltafcurv(struct.curv,rsigmavar).*rkcurv;
%%