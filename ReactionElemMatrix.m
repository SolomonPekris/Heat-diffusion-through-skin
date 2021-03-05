function [elemat] = ReactionElemMatrix(lamda, eID, msh)
%%REACTIONELEMMMATRIX calculates the local 2-by-2 element matrix for the linear
%%reaction operator, for any element in the finite element mesh for a given
%%reaction coefficient lamda, for an element eID, in a linear 1D msh using.
%%a gauss scheme. For lamda = 1, this is also equal to the local element
%%Mass matrix

J = msh.elem(eID).J; %set Jacobian
elemat = 0; %initialise local element matrix
gq = CreateGQScheme(2); %Creates Gauss scheme

for i=1:2 % loop and calculate element matrix using gauss scheme
    xi = gq.xipts(i); % set xi to gauss points
    psim=[(1-xi)/2; (xi+1)/2]; % first weighting function
    psin=[(1-xi)/2  (xi+1)/2]; % second weighting function
    elemat = elemat + gq.gsw(i) * psim * psin; % sum gauss points
end

elemat = elemat * J * lamda; % multiply entire matrix by J and lamda
