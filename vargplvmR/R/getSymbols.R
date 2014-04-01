getSymbols <-
function (number)
{
# % GETSYMBOLS Get a cell array of different plot symbols.
# % FORMAT
# % DESC returns a cell array of different plot symbols. A maximum of
# % 66 distinct symbols will be created.
# % ARG number : the number of different plot symbols required.
# % RETURN symbol : cell array of the different symbols.
# %
# % SEEALSO : plot
# %
# % COPYRIGHT : Neil D. Lawrence, 2005
  # 
# % NDLUTIL
  
#   symbolColour <- c("r", "g", "b", "c", "m")  #  %, "y"
  symbolColour <- c(2:6)
#   symbolShape <- c("x", "o", "+", "*", "s", "d", "v", "^", "<", ">", "p") 
  symbolShape <-c(4,1,3,8,0,5,6,2,15,16,17)
  counter <- 0 
  symbol <- list()
  while (counter < number)
  {
    symbol[[counter+1]] <- c(symbolColour[(counter %% length(symbolColour))+1],
                             symbolShape[(counter %% length(symbolShape))+1]) 
    counter <- counter +1 
  }
  return (symbol)
}
