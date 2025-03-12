#' @name alpha
#' @title alpha
#' @description
#' Give transparency value to colors
#'
#' Takes provided colors and gives them the specified alpha (ie transparency) value
#'
#' @author Marcin Imielinski
#' @param col RGB color
#' @keywords internal
#' @export
alpha = function(col, alpha)
{
  col.rgb = col2rgb(col)
  out = rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
  names(out) = names(col)
  return(out)
}

#' @name dodo.call
#' @title dodo.call
#' @description
#' do.call implemented using eval parse for those pesky (e.g. S4) case when do.call does not work
dodo.call = function(FUN, args)
{
    if (!is.character(FUN))
        FUN = substitute(FUN)
    cmd = paste(FUN, '(', paste('args[[', 1:length(args), ']]', collapse = ','), ')', sep = '')
    return(eval(parse(text = cmd)))
}


#' @name dedup
#' @title dedup
#'
#' @description
#' relabels duplicates in a character vector with .1, .2, .3
#' (where "." can be replaced by any user specified suffix)
#'
#' @param x input vector to dedup
#' @param suffix suffix separator to use before adding integer for dups in x
#' @return length(x) vector of input + suffix separator + integer for dups and no suffix for "originals"
#' @author Marcin Imielinski
dedup = function(x, suffix = '.')
{
  dup = duplicated(x);
  udup = setdiff(unique(x[dup]), NA)
  udup.ix = lapply(udup, function(y) which(x==y))
  udup.suffices = lapply(udup.ix, function(y) c('', paste(suffix, 2:length(y), sep = '')))
  out = x;
  out[unlist(udup.ix)] = paste(out[unlist(udup.ix)], unlist(udup.suffices), sep = '');
  return(out)
}

#' @name deanchor
#' @title deanchor
#' @author Kevin Hadi
deanchor = function(jj, gm, flip = FALSE) {
  if (flip)
    bps = gr2dt(c(gr.flipstrand(jj$left), jj$right))[, .(sn = seqnames, st = start, en = end, str = strand, key = c('Left', 'Right'))]
  else
    bps = gr2dt(c(jj$left, jj$right))[, .(sn = seqnames, st = start, en = end, str = strand, key = c('Left', 'Right'))]
  setkey(bps, key)
  newgr = dt2gr(gr2dt(gm$gr)[, .(seqnames = bps[.(seqnames), sn],
                                 start = ifelse(bps[.(seqnames), str] == '+',
                                                start + bps[.(seqnames), st],
                                                bps[.(seqnames), st]-end),
                                 end = ifelse(bps[.(seqnames), str] == '+',
                                              end + bps[.(seqnames), st],
                                              bps[.(seqnames), st]-start))])
  new.gm = gM(newgr, gm$dat)
  return(new.gm)
}
