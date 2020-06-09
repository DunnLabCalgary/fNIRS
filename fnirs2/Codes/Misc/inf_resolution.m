function [dod] = inf_resolution (d,dod)

%%%%%%%%%%%remove infs%%%%%%%%%%%
[a,b]=find(isinf(dod));
if isempty(a)==0
    for i=1:length(a)
        d(a(i),b(i))=nan;
        %d(:,b(i))=inpaint_nans(d(:,b(i)));
        if a(i)~=1
            dod(a(i),b(i))=dod(a(i)-1,b(i));
        else
            dod(a(i),b(i))=dod(a(i)+1,b(i));
        end
    end
end
dod=inpaint_nans(dod);
%ending inf resolution

end