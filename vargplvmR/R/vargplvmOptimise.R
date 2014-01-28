vargplvmOptimise <-
function (model, display, iters, ...)
{
# % VARGPLVMOPTIMISE Optimise the VARGPLVM.
# % FORMAT
# % DESC takes a given GP-LVM model structure and optimises with
# % respect to parameters and latent positions. 
# % ARG model : the model to be optimised.
# % ARG display : flag dictating whether or not to display
# % optimisation progress (set to greater than zero) (default value 1). 
# % ARG iters : number of iterations to run the optimiser
# % for (default value 2000).
# % RETURN model : the optimised model.
# %
# % SEEALSO : vargplvmCreate, vargplvmLogLikelihood,
# % vargplvmLogLikeGradients, vargplvmObjective, vargplvmGradient
# % 
# % COPYRIGHT : Michalis K. Titsias, 2009
# % 
# % COPYRIGHT : Neil D. Lawrence, 2005, 2006
# 
# % VARGPLVM

if (nargs() < 3)
{
  iters <- 2000 
  if (nargs() < 2)
    display <- 1 
}

varargin <-list(...)
options <- optOptions() 
if (length(varargin) == 2)
{
  cat("vargplvmOptimise \n")
#   to do
#     if (varargin[[1]] == "gradcheck")
#     {
#         assert(islogical(varargin{2})) 
#         options(9) <- varargin{2} 
#         if options(9)
#             [params, names] <- modelExtractParam(model) 
#             for i<-1:length(names)
#                 fprintf('%d\t%s\n', i, names{i}) 
#             end
#             feval('gradchek2', params, @vargplvmObjective, @vargplvmGradient, model, names) 
#         else
#             params <- modelExtractParam(model) 
#         end
#     }
} else {
 
    params <- modelExtractParam(model)
}

options[2] <- 0.1*options[2]  
options[3] <- 0.1*options[3] 

if (display)
{
  options[1] <- 1 
  if (length(params) <= 100)
    options[9] <- 1 
}
options[14] <- iters 

if (("optimiser" %in% names(model))) # && ~isa(model.optimiser, 'function_handle')
{    
    if ("optimiser" %in% names(model))
      optim <- model$optimiser 
    else
      optim <- "scg" 
    
    
    if (optim == "optimiMinimize")
    {
#         % Carl Rasmussen's minimize function 
        params <- do.call(optim, list("vargplvmObjectiveGradient", params, options, model)) 
    } else if (optim == "scg2") {
#         % NETLAB style optimization with a slight modification so that an
#         % objectiveGradient can be used where applicable, in order to re-use
#         % precomputed quantities.
        params <- do.call(optim, list("vargplvmObjectiveGradient", params,  options,  "vargplvmGradient", model)) 
    } else {
#         % NETLAB style optimization.
        params <- do.call(optim, list("vargplvmObjective", params,  options,  "vargplvmGradient", model)) 
        }
    
} 
# else if (("optimiser" %in% names(model)) && isa(model.optimiser, 'function_handle')
# {    
#     f <- fcnchk(model.optimiser) 
#     params <- f(model) 
#     
# } else {
#     error('vargplvmOptimise: Invalid optimiser setting.') 
# }

# %model <- vargplvmExpandParam(model, params) 
model <- modelExpandParam(model, params) 

return (model)
}
