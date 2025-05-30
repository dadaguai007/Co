function output = mzm(Ai,V, Vdc, Vpi)
% MZM pull-pull
% scalar 
%Ai equal optical signal
% optical filed E is output
if isscalar(Ai)
    Ai = Ai * ones(size(V));
end
    output =Ai.* cos(pi/2 * (V + Vdc) / Vpi);
end


