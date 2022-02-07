setGeneric('%Q%', function(x, ...) standardGeneric('%Q%'))
setMethod("%Q%", signature(x = "GRanges"), function(x, y) {
  condition_call  = substitute(y)
  ## serious R voodoo gymnastics .. but I think finally hacked it to remove ghosts
  ## create environment that combines the calling env with the granges env
  env = as(c(as.list(parent.frame(2)), as.list(as.data.frame(x))), 'environment')
  parent.env(env) = parent.frame()
  ix = tryCatch(eval(condition_call, env), error = function(e) NULL)
  if (is.null(ix))
  {
    condition_call  = substitute(y)
    ix = eval(condition_call, GenomicRanges::as.data.frame(x))
  }
  return(x[ix])
})