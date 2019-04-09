% Finds the flipped, conjugated, circ-shifted, and phase-rotated version 
% of each column of "in" that best matches each column of "ref"
% 
% [out,fxn]=disambig1Dfft(in,ref)

function [out,fxn]=disambig1Dfft(in,ref)

 % extract input size
 [N,L] = size(in);

 % disambiguate the intput signals
 out = nan(size(in));
 kk_best = nan(1,L);
 flip_best = nan(1,L);
 con_best = nan(1,L);
 angl_best = nan(1,L);
 for ll=1:L
   minErr = inf;
   for flip=0:1
     if flip,
       in_flip = flipud(in(:,ll));
     else
       in_flip = in(:,ll);
     end
     for con=0:1
       if con,
         in_con = conj(in_flip);
       else
         in_con = in_flip;
       end
       for kk=0:N-1
         in_shift = circshift(in_con,kk);
         angl = sign(in_shift'*ref(:,ll));
         in_rot = in_shift*angl;
         err = norm(in_rot-ref(:,ll));
         if err<minErr
           minErr = err;
           out(:,ll) = in_rot;
           kk_best(ll) = kk;
           flip_best(ll) = flip;
           con_best(ll) = con;
           angl_best(ll) = angl;
         end
       end
     end
   end
 end

 J = speye(L); 
 % this fxn applies flipud to selected columns
 flipud1 = @(inmat,cols) [flipud(inmat(:,cols)),inmat(:,setdiff(1:L,cols))]...
                *J([cols,setdiff(1:L,cols)],:); 
 % this fxn applies conj to selected columns
 conj1 = @(inmat,cols) [conj(inmat(:,cols)),inmat(:,setdiff(1:L,cols))]...
                *J([cols,setdiff(1:L,cols)],:); 
 fxn = @(invec) bsxfun(@times,...
                   circshift(...
                     conj1(...
                       flipud1(invec,find(flip_best)),...
                     find(con_best)),...
                   kk_best),...
                 angl_best);

%if nargout>1
%  if flip_best
%    if con_best
%      fxn = @(invec) circshift(conj(flipud(invec)),kk_best)*angl_best;
%    else
%      fxn = @(invec) circshift(flipud(invec),kk_best)*angl_best;
%    end
%  else
%    if con_best
%      fxn = @(invec) circshift(conj(invec),kk_best)*angl_best;
%    else
%      fxn = @(invec) circshift(invec,kk_best)*angl_best;
%    end
%  end
%end
