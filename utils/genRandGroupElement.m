function R = genRandGroupElement(ROffType,d)
%GENRANDGROUPELEMENT Summary of this function goes here
%   Detailed explanation goes here

switch ROffType
    case 'O'
        R = orth(rand(d)); %%% "orth" returns an orthonormal basis (det(R)=1)
        if rand(1) < 0.5
            R(:,1) = -R(:,1);
        end
    case 'SO'
        R = orth(rand(d)); %%% "orth" returns an orthonormal basis (det(R)=1)
    case 'Perm'
        RVec = assignmentoptimal(rand(d));
        R = full(sparse(1:length(RVec), RVec, ones(1,length(RVec))));
    otherwise
        error(['unknown ROffType: ' ROffType])
end


end

