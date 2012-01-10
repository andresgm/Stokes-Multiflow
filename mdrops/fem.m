% FEM algorithm for the calculating the tensions on a lipid membrane. 
% This tensions arise in order to mantain constant surface area.

function [ftensionn,ftensiont] = fem(geom,reference)

% The function fem has the following input arguments:
% geom: is a structure that contains the coordinates and coordination of the
% elements, along with other properties of the membrane


