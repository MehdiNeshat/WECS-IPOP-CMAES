function [stop]=Stop_output_Func(x, optimValues, state)
   stop = false;
   if (optimValues.fval==0) %&& optimValues.iteration > 10
      stop = true; 
      % Here, I set max. iterations to 10 if still detecting Inf
      % function value; This is specific to my own research problem.
   end

end