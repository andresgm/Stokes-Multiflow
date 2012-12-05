function [xg,wg] = quadrature
% Returns 2 lists, the first are the function evaluation points for the
% domain [0,1] in 2D. The second list are the accompanying weights.

A0 = 1/3;
A1 = 0.059715871789770;
B1 = 0.470142064105115;
A2 = 0.797426985353087;
B2 = 0.101286507323456;
W0 = 0.1125;
W1 = 0.066197076394253;
W2 = 0.062969590272413;

xg = [  A0, A0;
        A1, B1;
        B1, A1;
        B1, B1;
        B2, A2;
        B2, B2;
        A2, B2];
wg = [  W0, W1, W1, W1, W2, W2, W2];