function [M, inliers] = ransac(x,modelfunct,distfunct,degenfunct,s,t,feedback,maxdatatrials,maxtrials)


%Arguments: 
 %%x is a stack of X1 and X2 are 3byn homogenous 2d points (normalized)..only for our code call by ransacfund.
%     modelfunct - Handle to a function that fits a model to s
%                 data from x.  
%                 Note it is possible that the model fitting function can return
%                 multiple models (for example up to 3 fundamental matrices
%                 can be fitted to 7 matched points).  In this case it is
%                 assumed that the fitting function returns a cell array of
%                 models.
%                 If this function cannot fit a model it should return M as
%                 an empty matrix.
%
%     distfunct    - Handle to a function that evaluates the
%                 distances from the model to data x.
%                 
%                 This function must evaluate the distances between points
%                 and the model returning the indices of elements in x that
%                 are inliers, that is, the points that are within distance
%                 't' of the model.  Additionally, if M is a cell array of
%                 possible models 'distfunct' will return the model that has the
%                 most inliers.  If there is only one model this function
%                 must still copy the model to the output.  After this call M
%                 will be a non-cell object representing only one model. 
%
%     degenfunct   - Handle to a function that determines whether a
%                 set of datapoints will produce a degenerate model.
%                 This is used to discard random samples that do not
%                 result in useful models.
%                 It is assumed that degenfn is a boolean function of
%                 the form: 
%                    r = degenfunct(x)
%                 It may be that you cannot devise a test for degeneracy in
%                 which case you should write a dummy function that always
%                 returns a value of 1 (true) and rely on 'modelfunct' to return
%                 an empty model should the data set be degenerate.
%
%     s         - The minimum number of samples from x required by
%                 fittingfn to fit a model.
%
%     t         - The distance threshold between a data point and the model
%                 used to decide whether the point is an inlier or not.
%
%     feedback  - An optional flag 0/1. If set to one the trial count and the
%                 estimated total number of trials required is printed out at
%                 each step.  Defaults to 0.
%
%     maxDataTrials - Maximum number of attempts to select a non-degenerate
%                     data set. This parameter is optional and defaults to 100.
%
%     maxTrials - Maximum number of iterations. This parameter is optional and
%                 defaults to 1000.
%
% Returns:
%     M         - The model having the greatest number of inliers.
%     inliers   - An array of indices of the elements of x that were
%                 the inliers for the best model.
%Reference: http://www.peterkovesi.com and Andrew Zisserman's Multiple View
%Geometry


%test of parameters
error(nargchk(6,9,nargin));
fprintf('ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo1');
 if nargin<9; maxtrials= 1000; end;
 if nargin<8; maxdatatrials = 100; end;
 if nargin<7; feedback=1; end;
 
 [rows npts] = size(x);
 

 
p = 0.99;

bestM = NaN; 
trialcount = 0;
bestscore = 0;
N=1;

while trialcount<100
    %N>trialcount
    degenerate =1;
    count=1;
    
    while degenerate
       
        ind = randsample(npts,s);
       
       degenerate  = feval(degenfunct, x(:,ind));
       
       if ~degenerate
          
           M = feval(modelfunct,x(:,ind));
           
           if isempty(M)
               degenerate = 1; 
           end
       end
       
       
       count = count+1;
       if count > maxdatatrials
           warning('unable to select a nondegenerate data set')
           break 
       end
       
    end
       
    [inliers] = feval(distfunct,M,x,t);
    
    
    ninliers = length(inliers);
    
    if ninliers > bestscore
        bestscore = ninliers;
        bestinliers =inliers;
        bestM = M;
        
   
     
     fracinliers = ninliers/npts;
     pNoOutliers = 1-fracinliers
     N = log(1-p)/log(pNoOutliers)
     
     
    end
     
    trialcount = trialcount+1
    
    if feedback
        fprintf('trial %d out of %d     \r',trialcount,ceil(N));
         fprintf('ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo');
    end
 
    %safegaurd against being stuck in this loop forever
    if trialcount > maxtrials
        warning(sprintf('ransac reached the maximum number of %d trials', maxtrials));
        break
    end
end
 
 if ~isnan(bestM)
     M = bestM;
     inliers = bestinliers;
     fprintf('ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo');
 else
     M = [];
     inliers = [];
     warning('ransac was unable to find a useful solution');
 end
 

